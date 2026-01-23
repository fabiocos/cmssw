/**
 * GNNVertexProducerAlpaka - Proper Alpaka-based GNN vertex producer
 *
 * Uses TensorCollection for BOTH input and output, following SimpleNet.cc pattern:
 *   - view.records().column() accessor for TensorCollection::add()
 *   - model_.forward(event.queue(), inputs, outputs)
 *
 * Input: TrackFeaturesDeviceCollection [N, 13]
 * Output: GNNOutputDeviceCollection [N, K] with Eigen columns
 */

#include <Eigen/Core>
#include <Eigen/Dense>

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
#include "RecoVertex/PrimaryVertexProducer/interface/VertexGNNSoA.h"
#include "RecoVertex/PrimaryVertexProducer/interface/alpaka/VertexGNNDeviceCollection.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::vertexgnn {

  class GNNVertexProducerAlpaka : public stream::EDProducer<> {
  public:
    explicit GNNVertexProducerAlpaka(const edm::ParameterSet& params)
        : EDProducer<>(params),
          trackFeaturesToken_(consumes(params.getParameter<edm::InputTag>("trackFeatures"))),
          gnnOutputToken_{produces()},
          model_(params.getParameter<edm::FileInPath>("model").fullPath()),
          verbose_(params.getUntrackedParameter<bool>("verbose", false)) {
      edm::LogInfo("GNNVertexProducerAlpaka") << "Loaded TorchScript model";
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::FileInPath>("model", 
          edm::FileInPath("RecoVertex/PrimaryVertexProducer/data/dummy_vertex_slot.pt"));
      desc.add<edm::InputTag>("trackFeatures", edm::InputTag("trackFeatureSource"));
      desc.addUntracked<bool>("verbose", false);
      descriptions.addWithDefaultLabel(desc);
    }

    void produce(device::Event& event, const device::EventSetup& eventSetup) override {
      // Get input collection
      const auto& trackFeatures = event.get(trackFeaturesToken_);
      const auto N = trackFeatures.const_view().metadata().size();

      if (verbose_) {
        edm::LogInfo("GNNVertexProducerAlpaka") << "N=" << N << " tracks, K=" << ::vertexgnn::kNumSlots << " slots";
      }

      // Allocate output collection on device (batch = N)
      auto gnnOutput = GNNOutputDeviceCollection(N, event.queue());

      // Get SoA records (returns tuple-wrapped accessors for TensorCollection)
      auto inputRecords = trackFeatures.const_view().records();
      auto outputRecords = gnnOutput.view().records();

      // =========================================================================
      // INPUT: TensorCollection from SoA records → [N, 13]
      // =========================================================================
      cms::torch::alpakatools::TensorCollection<Queue> inputs(N);
      inputs.add<::vertexgnn::TrackFeaturesSoA>("features",
          inputRecords.vz(), inputRecords.dz(), inputRecords.pt(), inputRecords.eta(),
          inputRecords.mva(), inputRecords.pl(), inputRecords.t_pi(), inputRecords.t_k(),
          inputRecords.t_p(), inputRecords.s_pi(), inputRecords.s_k(), inputRecords.s_p(),
          inputRecords.has_time());

      // =========================================================================
      // OUTPUT: TensorCollection to SoA Eigen columns
      // =========================================================================
      // Model outputs 5 tensors: (A[N,K], z_hat[N,K], t_hat[N,K], p[N,K], pi[N,3])
      // Order must match model output tuple order!
      cms::torch::alpakatools::TensorCollection<Queue> outputs(N);
      outputs.add<::vertexgnn::GNNOutputSoA>("A", outputRecords.A());
      outputs.add<::vertexgnn::GNNOutputSoA>("z_hat", outputRecords.z_hat());
      outputs.add<::vertexgnn::GNNOutputSoA>("t_hat", outputRecords.t_hat());
      outputs.add<::vertexgnn::GNNOutputSoA>("p", outputRecords.p());
      outputs.add<::vertexgnn::GNNOutputSoA>("pi", outputRecords.pi());

      // =========================================================================
      // INFERENCE: Proper Alpaka forward
      // =========================================================================
      if (verbose_) {
        edm::LogInfo("GNNVertexProducerAlpaka") << "Running model_.forward(queue, inputs, outputs)...";
      }

      model_.forward(event.queue(), inputs, outputs);

      if (verbose_) {
        edm::LogInfo("GNNVertexProducerAlpaka") << "Inference complete";
      }

      // =========================================================================
      // PUT OUTPUT
      // =========================================================================
      event.emplace(gnnOutputToken_, std::move(gnnOutput));
    }

  private:
    const device::EDGetToken<TrackFeaturesDeviceCollection> trackFeaturesToken_;
    const device::EDPutToken<GNNOutputDeviceCollection> gnnOutputToken_;
    torch::AlpakaModel model_;
    const bool verbose_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vertexgnn

DEFINE_FWK_ALPAKA_MODULE(vertexgnn::GNNVertexProducerAlpaka);
