/**
 * GNNVertexProducerAlpaka - Alpaka-based GNN vertex producer with GPU inference
 *
 * Consumes TrackFeaturesHostCollection from standard EDProducer, copies to device,
 * runs inference on GPU, and produces GNNOutputDeviceCollection.
 *
 * The framework handles automatic D2H transfer when downstream consumers need HostCollection.
 *
 * Input: TrackFeaturesHostCollection [N, 13] (from TrackFeatureProducer)
 * Output: GNNOutputDeviceCollection [N, K] + [N, 3]
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
#include "HeterogeneousCore/AlpakaInterface/interface/host.h"
#include "PhysicsTools/PyTorchAlpaka/interface/TensorCollection.h"
#include "PhysicsTools/PyTorchAlpaka/interface/alpaka/AlpakaModel.h"
#include "DataFormats/VertexGNNReco/interface/VertexGNNSoA.h"
#include "DataFormats/VertexGNNReco/interface/VertexGNNHostCollection.h"
#include "DataFormats/VertexGNNReco/interface/alpaka/VertexGNNDeviceCollection.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::vertexgnn {

  class GNNVertexProducerAlpaka : public stream::EDProducer<> {
  public:
    explicit GNNVertexProducerAlpaka(const edm::ParameterSet& params)
        : EDProducer<>(params),
          // Consume HostCollection from TrackFeatureProducer (standard EDProducer)
          trackFeaturesToken_(consumes<::vertexgnn::TrackFeaturesHostCollection>(
              params.getParameter<edm::InputTag>("trackFeatures"))),
          gnnOutputToken_{produces()},
          model_(params.getParameter<edm::FileInPath>("model").fullPath()),
          verbose_(params.getUntrackedParameter<bool>("verbose", false)) {
      edm::LogInfo("GNNVertexProducerAlpaka") << "Loaded TorchScript model (GPU inference enabled)";
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
      // Get input HOST collection from standard EDProducer
      const auto& hostInput = event.get(trackFeaturesToken_);
      const auto N = hostInput.const_view().metadata().size();

      if (verbose_) {
        edm::LogInfo("GNNVertexProducerAlpaka") << "N=" << N << " tracks (from HostCollection), K=" << ::vertexgnn::kNumSlots << " slots";
      }

      // Copy input to device for GPU inference
      TrackFeaturesDeviceCollection deviceInput(N, event.queue());
      alpaka::memcpy(event.queue(), deviceInput.buffer(), hostInput.buffer());

      // Allocate output on device
      GNNOutputDeviceCollection deviceOutput(N, event.queue());

      // Get SoA records for TensorCollection
      auto inputRecords = deviceInput.const_view().records();
      auto outputRecords = deviceOutput.view().records();

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
      cms::torch::alpakatools::TensorCollection<Queue> outputs(N);
      outputs.add<::vertexgnn::GNNOutputSoA>("A", outputRecords.A());
      outputs.add<::vertexgnn::GNNOutputSoA>("z_hat", outputRecords.z_hat());
      outputs.add<::vertexgnn::GNNOutputSoA>("t_hat", outputRecords.t_hat());
      outputs.add<::vertexgnn::GNNOutputSoA>("p", outputRecords.p());
      outputs.add<::vertexgnn::GNNOutputSoA>("pi", outputRecords.pi());

      // =========================================================================
      // INFERENCE: Run on device (GPU if available)
      // =========================================================================
      if (verbose_) {
        edm::LogInfo("GNNVertexProducerAlpaka") << "Running model_.forward(queue, inputs, outputs) on device...";
      }

      model_.forward(event.queue(), inputs, outputs);

      if (verbose_) {
        edm::LogInfo("GNNVertexProducerAlpaka") << "Inference complete";
      }

      // =========================================================================
      // PUT OUTPUT: DeviceCollection - framework handles D2H for consumers
      // =========================================================================
      event.emplace(gnnOutputToken_, std::move(deviceOutput));
    }

  private:
    const edm::EDGetTokenT<::vertexgnn::TrackFeaturesHostCollection> trackFeaturesToken_;
    const device::EDPutToken<GNNOutputDeviceCollection> gnnOutputToken_;
    torch::AlpakaModel model_;
    const bool verbose_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vertexgnn

DEFINE_FWK_ALPAKA_MODULE(vertexgnn::GNNVertexProducerAlpaka);
