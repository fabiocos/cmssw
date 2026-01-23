/**
 * GNNVertexBuilderFromAlpaka - Consumes Alpaka SoA outputs, produces reco::VertexCollection
 *
 * This is a STANDARD EDProducer that:
 *   1. Consumes TrackCollection (standard EDM)
 *   2. Consumes SlotPredictionsDeviceCollection from GNNVertexProducerAlpaka (Alpaka SoA)
 *   3. Consumes AssignmentDeviceCollection from GNNVertexProducerAlpaka (Alpaka SoA)
 *   4. Builds reco::VertexCollection using GNNClusterizerFromAlpaka
 *
 * The INFERENCE happens in GNNVertexProducerAlpaka (Alpaka producer) upstream.
 * This module only CONSUMES the SoA results and builds vertices.
 *
 * Pipeline:
 *   TrackFeatureProducerAlpaka (Alpaka) → TrackFeaturesDeviceCollection
 *         ↓
 *   GNNVertexProducerAlpaka (Alpaka) → SlotPredictions + Assignments
 *         ↓
 *   GNNVertexBuilderFromAlpaka (this, standard) → reco::VertexCollection
 */

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "RecoVertex/PrimaryVertexProducer/interface/GNNClusterizerFromAlpaka.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexGNNSoA.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexGNNHostCollection.h"

class GNNVertexBuilderFromAlpaka : public edm::stream::EDProducer<> {
public:
  explicit GNNVertexBuilderFromAlpaka(const edm::ParameterSet& params);
  ~GNNVertexBuilderFromAlpaka() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void produce(edm::Event& event, const edm::EventSetup& eventSetup) override;

private:
  // Input tokens
  const edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  const edm::EDGetTokenT<vertexgnn::SlotPredictionsHostCollection> slotPredictionsToken_;
  const edm::EDGetTokenT<vertexgnn::AssignmentHostCollection> assignmentsToken_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbToken_;
  
  // Clusterizer
  std::unique_ptr<vertexgnn::GNNClusterizerFromAlpaka> clusterizer_;
  
  // Configuration
  const int numSlots_;
  const bool verbose_;
};

GNNVertexBuilderFromAlpaka::GNNVertexBuilderFromAlpaka(const edm::ParameterSet& params)
    : trackToken_(consumes<reco::TrackCollection>(params.getParameter<edm::InputTag>("tracks"))),
      slotPredictionsToken_(consumes<vertexgnn::SlotPredictionsHostCollection>(
          params.getParameter<edm::InputTag>("slotPredictions"))),
      assignmentsToken_(consumes<vertexgnn::AssignmentHostCollection>(
          params.getParameter<edm::InputTag>("assignments"))),
      ttbToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      numSlots_(params.getParameter<int>("numSlots")),
      verbose_(params.getUntrackedParameter<bool>("verbose", false)) {
  
  produces<reco::VertexCollection>();
  
  // Create clusterizer
  clusterizer_ = std::make_unique<vertexgnn::GNNClusterizerFromAlpaka>(params);
  
  edm::LogInfo("GNNVertexBuilderFromAlpaka") << "Initialized - consumes Alpaka SoA outputs";
}

void GNNVertexBuilderFromAlpaka::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("slotPredictions", edm::InputTag("gnnVertexProducer"));
  desc.add<edm::InputTag>("assignments", edm::InputTag("gnnVertexProducer"));
  desc.add<int>("numSlots", 200);
  desc.add<double>("existenceThreshold", 0.5);
  desc.add<double>("trackAssignmentThreshold", 0.0);
  desc.addUntracked<bool>("verbose", false);
  
  descriptions.addWithDefaultLabel(desc);
}

void GNNVertexBuilderFromAlpaka::produce(edm::Event& event, const edm::EventSetup& eventSetup) {
  // Get inputs
  auto tracks = event.getHandle(trackToken_);
  auto slotPredictions = event.getHandle(slotPredictionsToken_);
  auto assignments = event.getHandle(assignmentsToken_);
  auto ttBuilder = eventSetup.getHandle(ttbToken_);
  
  auto output = std::make_unique<reco::VertexCollection>();
  
  if (!tracks.isValid() || tracks->empty() || 
      !slotPredictions.isValid() || !assignments.isValid()) {
    event.put(std::move(output));
    return;
  }
  
  // Build TransientTracks
  std::vector<reco::TransientTrack> ttracks;
  ttracks.reserve(tracks->size());
  for (const auto& track : *tracks) {
    ttracks.push_back(ttBuilder->build(track));
  }
  
  const int K = slotPredictions->const_view().metadata().size();
  // N must come from assignment collection, not tracks - to match upstream producer
  const int assignmentSize = assignments->const_view().metadata().size();
  const int N = (K > 0) ? (assignmentSize / K) : 0;
  
  // Warn if track count doesn't match - this indicates TrackFeatureSource isn't consuming real tracks
  if (N != static_cast<int>(ttracks.size())) {
    edm::LogWarning("GNNVertexBuilderFromAlpaka") 
        << "Track count mismatch! GNN processed " << N << " tracks but generalTracks has " 
        << ttracks.size() << ". Using N=" << N << " from assignment matrix.";
  }
  
  if (verbose_) {
    edm::LogInfo("GNNVertexBuilderFromAlpaka") 
        << "Consuming Alpaka outputs: N=" << N << " (GNN), K=" << K << " slots";
  }
  
  // Get SoA views - these are the ALREADY COMPUTED Alpaka outputs
  auto slotView = slotPredictions->const_view();
  auto assignView = assignments->const_view();
  
  // Extract pointers to data from SoA
  // Need to build flat arrays for the clusterizer interface
  std::vector<float> A_flat(N * K);
  std::vector<float> z_hat(K), t_hat(K), p(K);
  std::vector<float> pi_0(K), pi_1(K), pi_2(K), pi_3(K);
  
  for (int k = 0; k < K; ++k) {
    z_hat[k] = slotView[k].z_hat();
    t_hat[k] = slotView[k].t_hat();
    p[k] = slotView[k].p();
    pi_0[k] = slotView[k].pi_0();
    pi_1[k] = slotView[k].pi_1();
    pi_2[k] = slotView[k].pi_2();
    pi_3[k] = slotView[k].pi_3();
  }
  
  for (int i = 0; i < N; ++i) {
    for (int k = 0; k < K; ++k) {
      A_flat[i * K + k] = assignView[i * K + k].prob();
    }
  }
  
  // Build vertices using clusterizer
  // If N doesn't match ttracks.size(), we can only use first N tracks
  std::vector<reco::TransientTrack> tracksToUse;
  if (N <= static_cast<int>(ttracks.size())) {
    tracksToUse.assign(ttracks.begin(), ttracks.begin() + N);
  } else {
    tracksToUse = ttracks;
    edm::LogWarning("GNNVertexBuilderFromAlpaka") 
        << "GNN processed more tracks than available: " << N << " vs " << ttracks.size();
  }
  
  std::vector<TransientVertex> transVertices = clusterizer_->vertices(
      tracksToUse,
      A_flat.data(),
      z_hat.data(),
      t_hat.data(),
      p.data(),
      pi_0.data(),
      pi_1.data(),
      pi_2.data(),
      pi_3.data(),
      N, K);
  
  // Convert TransientVertex to reco::Vertex
  for (const auto& tv : transVertices) {
    // Convert GlobalPoint to reco::Vertex::Point
    reco::Vertex::Point pos(tv.position().x(), tv.position().y(), tv.position().z());
    
    // Convert GlobalError to reco::Vertex::Error (6-element symmetric matrix)
    reco::Vertex::Error err;
    err(0, 0) = tv.positionError().cxx();
    err(1, 0) = tv.positionError().cyx();
    err(1, 1) = tv.positionError().cyy();
    err(2, 0) = tv.positionError().czx();
    err(2, 1) = tv.positionError().czy();
    err(2, 2) = tv.positionError().czz();
    
    reco::Vertex rv(pos, err, tv.totalChiSquared(), tv.degreesOfFreedom(), 
                    tv.originalTracks().size());
    output->push_back(rv);
  }
  
  if (verbose_) {
    edm::LogInfo("GNNVertexBuilderFromAlpaka") 
        << "Built " << output->size() << " vertices from Alpaka SoA";
  }
  
  event.put(std::move(output));
}

DEFINE_FWK_MODULE(GNNVertexBuilderFromAlpaka);
