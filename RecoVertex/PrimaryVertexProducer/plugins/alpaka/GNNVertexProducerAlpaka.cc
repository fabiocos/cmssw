/**
 * GNNVertexProducerAlpaka - PyTorchAlpaka-based GNN vertex clustering
 *
 * Uses TorchScript model for direct inference via PhysicsTools/PyTorchAlpaka.
 * Produces vertex clusters from track features using the VertexSlotModel.
 */

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
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
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include <cmath>
#include <vector>

namespace ALPAKA_ACCELERATOR_NAMESPACE::vertexgnn {

  class GNNVertexProducerAlpaka : public stream::EDProducer<> {
  public:
    explicit GNNVertexProducerAlpaka(const edm::ParameterSet& params);

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    void produce(device::Event& event, const device::EventSetup& eventSetup) override;

  private:
    // Configuration
    const int numSlots_;
    const float existenceThreshold_;
    const int verbosity_;

    // Tokens
    const device::EDGetToken<TrackFeaturesDeviceCollection> trackFeaturesToken_;
    const device::EDPutToken<SlotPredictionsDeviceCollection> slotPredictionsToken_;
    const device::EDPutToken<AssignmentDeviceCollection> assignmentToken_;

    // PyTorch model
    torch::AlpakaModel model_;
  };

  GNNVertexProducerAlpaka::GNNVertexProducerAlpaka(const edm::ParameterSet& params)
      : EDProducer<>(params),
        numSlots_(params.getParameter<int>("numSlots")),
        existenceThreshold_(params.getParameter<double>("existenceThreshold")),
        verbosity_(params.getUntrackedParameter<int>("verbosity", 0)),
        trackFeaturesToken_(consumes(params.getParameter<edm::InputTag>("trackFeatures"))),
        slotPredictionsToken_{produces()},
        assignmentToken_{produces()},
        model_(params.getParameter<edm::FileInPath>("model").fullPath()) {}

  void GNNVertexProducerAlpaka::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::FileInPath>("model",
                              edm::FileInPath("RecoVertex/PrimaryVertexProducer/data/best_vertex_slot.pt"));
    desc.add<edm::InputTag>("trackFeatures", edm::InputTag("trackFeatureProducer"));
    desc.add<int>("numSlots", 200);
    desc.add<double>("existenceThreshold", 0.5);
    desc.addUntracked<int>("verbosity", 0);
    descriptions.addWithDefaultLabel(desc);
  }

  void GNNVertexProducerAlpaka::produce(device::Event& event, const device::EventSetup& eventSetup) {
    // Get input track features
    const auto& trackFeatures = event.get(trackFeaturesToken_);
    const auto N = trackFeatures.const_view().metadata().size();
    const auto K = numSlots_;

    if (verbosity_ > 0) {
      edm::LogInfo("GNNVertexProducerAlpaka") << "Processing " << N << " tracks with " << K << " slots";
    }

    // Allocate output collections
    auto slotPredictions = SlotPredictionsDeviceCollection(K, event.queue());
    auto assignments = AssignmentDeviceCollection(N * K, event.queue());

    // Get records for tensor binding
    auto inputRecords = trackFeatures.const_view().records();
    auto slotRecords = slotPredictions.view().records();
    auto assignRecords = assignments.view().records();

    // Build input tensor collection (13 features)
    cms::torch::alpakatools::TensorCollection<Queue> inputs(N);
    inputs.add<::vertexgnn::TrackFeaturesSoA>("x",
                                               inputRecords.vz(),
                                               inputRecords.dz(),
                                               inputRecords.pt(),
                                               inputRecords.eta(),
                                               inputRecords.mva(),
                                               inputRecords.pl(),
                                               inputRecords.t_pi(),
                                               inputRecords.t_k(),
                                               inputRecords.t_p(),
                                               inputRecords.s_pi(),
                                               inputRecords.s_k(),
                                               inputRecords.s_p(),
                                               inputRecords.has_time());

    // Build output tensor collection (5 outputs)
    // Note: The model outputs are:
    //   A: [N, K] assignment probabilities
    //   z_hat: [K] vertex z positions
    //   t_hat: [K] vertex t positions
    //   p: [K] existence probabilities
    //   pi: [K, 4] PID weights
    cms::torch::alpakatools::TensorCollection<Queue> outputs(K);
    outputs.add<::vertexgnn::AssignmentSoA>("A", assignRecords.prob());
    outputs.add<::vertexgnn::SlotPredictionsSoA>("z_hat", slotRecords.z_hat());
    outputs.add<::vertexgnn::SlotPredictionsSoA>("t_hat", slotRecords.t_hat());
    outputs.add<::vertexgnn::SlotPredictionsSoA>("p", slotRecords.p());
    outputs.add<::vertexgnn::SlotPredictionsSoA>("pi", 
                                                  slotRecords.pi_0(),
                                                  slotRecords.pi_1(),
                                                  slotRecords.pi_2(),
                                                  slotRecords.pi_3());

    // Run inference
    model_.forward(event.queue(), inputs, outputs);

    // Put outputs into event
    event.emplace(slotPredictionsToken_, std::move(slotPredictions));
    event.emplace(assignmentToken_, std::move(assignments));
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vertexgnn

DEFINE_FWK_ALPAKA_MODULE(vertexgnn::GNNVertexProducerAlpaka);
