/**
 * GNNVertexProducerAlpaka - Full PyTorchAlpaka-based GNN vertex producer
 *
 * Uses AlpakaModel with:
 * - TensorCollection for inputs (track features)
 * - Raw tensor outputs for 2D assignment matrix A
 * - Manual copy to SoA for outputs
 *
 * Model outputs: A[N,K], z_hat[K], t_hat[K], p[K], pi[K,4]
 */

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "PhysicsTools/PyTorchAlpaka/interface/TensorCollection.h"
#include "PhysicsTools/PyTorchAlpaka/interface/alpaka/AlpakaModel.h"
#include "PhysicsTools/PyTorch/interface/TorchInterface.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexGNNSoA.h"
#include "RecoVertex/PrimaryVertexProducer/interface/alpaka/VertexGNNDeviceCollection.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::vertexgnn {

  class GNNVertexProducerAlpaka : public stream::EDProducer<> {
  public:
    explicit GNNVertexProducerAlpaka(const edm::ParameterSet& params)
        : EDProducer<>(params),
          trackFeaturesToken_(consumes(params.getParameter<edm::InputTag>("trackFeatures"))),
          slotPredictionsToken_{produces()},
          assignmentToken_{produces()},
          model_(params.getParameter<edm::FileInPath>("model").fullPath()),
          numSlots_(params.getParameter<int>("numSlots")),
          verbose_(params.getUntrackedParameter<bool>("verbose", false)) {
      edm::LogInfo("GNNVertexProducerAlpaka") << "Loaded TorchScript model";
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::FileInPath>("model", 
          edm::FileInPath("RecoVertex/PrimaryVertexProducer/data/dummy_vertex_slot.pt"));
      desc.add<edm::InputTag>("trackFeatures", edm::InputTag("trackFeatureSource"));
      desc.add<int>("numSlots", 200);
      desc.addUntracked<bool>("verbose", false);
      descriptions.addWithDefaultLabel(desc);
    }

    void produce(device::Event& event, const device::EventSetup& eventSetup) override {
      // Get input
      const auto& trackFeatures = event.get(trackFeaturesToken_);
      const int N = trackFeatures.const_view().metadata().size();
      const int K = numSlots_;

      if (verbose_) {
        edm::LogInfo("GNNVertexProducerAlpaka") << "N=" << N << " tracks, K=" << K << " slots";
      }

      // Allocate outputs
      auto slotPredictions = SlotPredictionsDeviceCollection(K, event.queue());
      auto assignments = AssignmentDeviceCollection(N * K, event.queue());

      // Get SoA views for writing outputs
      auto slotView = slotPredictions.view();
      auto assignView = assignments.view();

      // =========================================================================
      // BUILD INPUT TENSOR FROM SOA
      // =========================================================================
      // Extract features from SoA to flat vector
      auto inputView = trackFeatures.const_view();
      std::vector<float> features(N * 13);
      
      for (int i = 0; i < N; ++i) {
        features[i * 13 + 0] = inputView[i].vz();
        features[i * 13 + 1] = inputView[i].dz();
        features[i * 13 + 2] = inputView[i].pt();
        features[i * 13 + 3] = inputView[i].eta();
        features[i * 13 + 4] = inputView[i].mva();
        features[i * 13 + 5] = inputView[i].pl();
        features[i * 13 + 6] = inputView[i].t_pi();
        features[i * 13 + 7] = inputView[i].t_k();
        features[i * 13 + 8] = inputView[i].t_p();
        features[i * 13 + 9] = inputView[i].s_pi();
        features[i * 13 + 10] = inputView[i].s_k();
        features[i * 13 + 11] = inputView[i].s_p();
        features[i * 13 + 12] = inputView[i].has_time();
      }

      // Create input tensor [1, N, 13]
      auto options = ::torch::TensorOptions().dtype(::torch::kFloat32);
      auto inputTensor = ::torch::from_blob(features.data(), {1, N, 13}, options);

      // =========================================================================
      // RUN INFERENCE (raw forward, not via TensorCollection)
      // =========================================================================
      std::vector<::torch::jit::IValue> inputs;
      inputs.push_back(inputTensor);

      if (verbose_) {
        edm::LogInfo("GNNVertexProducerAlpaka") << "Running inference...";
      }

      // Get base model and run forward
      auto output = model_.forward(inputs);

      // Parse tuple outputs: (A, z_hat, t_hat, p, pi)
      auto tuple = output.toTuple();
      auto A_tensor = tuple->elements()[0].toTensor();       // [N, K]
      auto z_hat_tensor = tuple->elements()[1].toTensor();   // [K]
      auto t_hat_tensor = tuple->elements()[2].toTensor();   // [K]
      auto p_tensor = tuple->elements()[3].toTensor();       // [K]
      auto pi_tensor = tuple->elements()[4].toTensor();      // [K, 4]

      // Get accessors
      auto A_acc = A_tensor.accessor<float, 2>();
      auto z_acc = z_hat_tensor.accessor<float, 1>();
      auto t_acc = t_hat_tensor.accessor<float, 1>();
      auto p_acc = p_tensor.accessor<float, 1>();
      auto pi_acc = pi_tensor.accessor<float, 2>();

      if (verbose_) {
        edm::LogInfo("GNNVertexProducerAlpaka") << "Output shapes: A[" 
            << A_tensor.size(0) << "," << A_tensor.size(1) << "], z_hat["
            << z_hat_tensor.size(0) << "], p[" << p_tensor.size(0) << "]";
        
        // Show first few vertex slots with high p
        edm::LogInfo("GNNVertexProducerAlpaka") << "--- Top vertex predictions ---";
        for (int k = 0; k < std::min(K, 10); ++k) {
          if (p_acc[k] > 0.1f) {
            edm::LogInfo("GNNVertexProducerAlpaka") 
                << "  Slot " << k << ": z=" << z_acc[k] << " t=" << t_acc[k] 
                << " p=" << p_acc[k];
          }
        }
        
        // Show first track assignments
        edm::LogInfo("GNNVertexProducerAlpaka") << "--- Track assignments (first 5 tracks) ---";
        for (int i = 0; i < std::min(N, 5); ++i) {
          float max_p = 0;
          int best_k = -1;
          for (int k = 0; k < K; ++k) {
            if (A_acc[i][k] > max_p) {
              max_p = A_acc[i][k];
              best_k = k;
            }
          }
          edm::LogInfo("GNNVertexProducerAlpaka") 
              << "  Track " << i << " -> Slot " << best_k << " (prob=" << max_p << ")";
        }
      }

      // =========================================================================
      // COPY OUTPUTS TO SOA
      // =========================================================================
      
      // Slot predictions: z_hat, t_hat, p, pi
      for (int k = 0; k < K; ++k) {
        slotView[k].z_hat() = z_acc[k];
        slotView[k].t_hat() = t_acc[k];
        slotView[k].p() = p_acc[k];
        slotView[k].pi_0() = pi_acc[k][0];
        slotView[k].pi_1() = pi_acc[k][1];
        slotView[k].pi_2() = pi_acc[k][2];
        slotView[k].pi_3() = pi_acc[k][3];
      }

      // Assignment matrix A: flatten [N, K] -> [N*K]
      for (int i = 0; i < N; ++i) {
        for (int k = 0; k < K; ++k) {
          assignView[i * K + k].prob() = A_acc[i][k];
        }
      }

      if (verbose_) {
        edm::LogInfo("GNNVertexProducerAlpaka") << "Copied outputs to SoA";
      }

      // =========================================================================
      // PUT OUTPUTS
      // =========================================================================
      event.emplace(slotPredictionsToken_, std::move(slotPredictions));
      event.emplace(assignmentToken_, std::move(assignments));
    }

  private:
    const device::EDGetToken<TrackFeaturesDeviceCollection> trackFeaturesToken_;
    const device::EDPutToken<SlotPredictionsDeviceCollection> slotPredictionsToken_;
    const device::EDPutToken<AssignmentDeviceCollection> assignmentToken_;
    torch::AlpakaModel model_;
    const int numSlots_;
    const bool verbose_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vertexgnn

DEFINE_FWK_ALPAKA_MODULE(vertexgnn::GNNVertexProducerAlpaka);
