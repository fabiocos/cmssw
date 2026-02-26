/**
 * GNNVertexBuilderFromAlpaka - Consumes Alpaka SoA outputs, produces reco::VertexCollection
 *
 * This is a STANDARD EDProducer that:
 *   1. Consumes TrackCollection (standard EDM)
 *   2. Consumes GNNOutputHostCollection from GNNVertexProducerAlpaka (unified Alpaka SoA)
 *   3. Builds reco::VertexCollection using GNNClusterizerFromAlpaka
 *
 * The GNNOutputHostCollection contains all outputs with batch=N:
 *   - A[N, K]: Assignment probabilities (Eigen::Vector per track)
 *   - z_hat[N, K], t_hat[N, K], p[N, K]: Slot predictions (replicated per track)
 *   - pi[N, 3]: PID weights per track (Eigen::Vector)
 *
 * Since slot predictions are replicated across all tracks, we only need
 * to read from the first track to get the slot values.
 */

#include <Eigen/Core>
#include <Eigen/Dense>

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
  const edm::EDGetTokenT<vertexgnn::GNNOutputHostCollection> gnnOutputToken_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbToken_;
  
  // Clusterizer
  std::unique_ptr<vertexgnn::GNNClusterizerFromAlpaka> clusterizer_;
  
  // Configuration
  const bool verbose_;
};

GNNVertexBuilderFromAlpaka::GNNVertexBuilderFromAlpaka(const edm::ParameterSet& params)
    : trackToken_(consumes<reco::TrackCollection>(params.getParameter<edm::InputTag>("tracks"))),
      gnnOutputToken_(consumes<vertexgnn::GNNOutputHostCollection>(
          params.getParameter<edm::InputTag>("gnnOutput"))),
      ttbToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      verbose_(params.getUntrackedParameter<bool>("verbose", false)) {
  
  produces<reco::VertexCollection>();
  
  // Create clusterizer
  clusterizer_ = std::make_unique<vertexgnn::GNNClusterizerFromAlpaka>(params);
  
  edm::LogInfo("GNNVertexBuilderFromAlpaka") << "Initialized - consumes unified GNNOutputHostCollection";
}

void GNNVertexBuilderFromAlpaka::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("gnnOutput", edm::InputTag("gnnVertexProducer"));
  desc.add<double>("existenceThreshold", 0.5);
  desc.add<double>("trackAssignmentThreshold", 0.0);
  desc.addUntracked<bool>("verbose", false);
  
  descriptions.addWithDefaultLabel(desc);
}

void GNNVertexBuilderFromAlpaka::produce(edm::Event& event, const edm::EventSetup& eventSetup) {
  // Get inputs
  auto tracks = event.getHandle(trackToken_);
  auto gnnOutput = event.getHandle(gnnOutputToken_);
  auto ttBuilder = eventSetup.getHandle(ttbToken_);
  
  auto output = std::make_unique<reco::VertexCollection>();
  
  if (!tracks.isValid() || tracks->empty() || !gnnOutput.isValid()) {
    event.put(std::move(output));
    return;
  }
  
  // Build TransientTracks
  std::vector<reco::TransientTrack> ttracks;
  ttracks.reserve(tracks->size());
  for (const auto& track : *tracks) {
    ttracks.push_back(ttBuilder->build(track));
  }
  
  // Get SoA view
  auto gnnView = gnnOutput->const_view();
  const int N = gnnView.metadata().size();
  constexpr int K = vertexgnn::kNumSlots;  // 180 (matches v17p1 model)
  
  if (verbose_) {
    edm::LogInfo("GNNVertexBuilderFromAlpaka") 
        << "Consuming unified GNN output: N=" << N << " tracks, K=" << K << " slots";
  }
  
  // Extract data from Eigen columns
  // Slot predictions are replicated per track, so read from track 0
  std::vector<float> A_flat(N * K);
  std::vector<float> z_hat(K), t_hat(K), p(K);
  std::vector<float> pi_0(N), pi_1(N), pi_2(N);  // Per track, not per slot!
  
  // Read slot predictions from first track (they're replicated)
  // NOTE: Copy element data immediately to avoid dangling references
  if (N > 0) {
    auto elem = gnnView[0];  // Copy the element to avoid temporary issues
    for (int k = 0; k < K; ++k) {
      z_hat[k] = elem.z_hat()[k];
      t_hat[k] = elem.t_hat()[k];
      p[k] = elem.p()[k];
    }
  }
  
  // Read assignment matrix A[N, K] and per-track PID weights
  for (int i = 0; i < N; ++i) {
    auto elem = gnnView[i];  // Copy element to avoid dangling reference
    for (int k = 0; k < K; ++k) {
      A_flat[i * K + k] = elem.A()[k];
    }
    // Per-track PID weights from Eigen column [3]
    pi_0[i] = elem.pi()[0];
    pi_1[i] = elem.pi()[1];
    pi_2[i] = elem.pi()[2];
  }
  
  // Warn if track count doesn't match
  if (N != static_cast<int>(ttracks.size())) {
    edm::LogWarning("GNNVertexBuilderFromAlpaka") 
        << "Track count mismatch! GNN processed " << N << " tracks but generalTracks has " 
        << ttracks.size() << ". Using N=" << N << " from GNN output.";
  }
  
  // Build vertices using clusterizer
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
      N, K);
  
  // Convert TransientVertex to reco::Vertex
  for (const auto& tv : transVertices) {
    reco::Vertex::Point pos(tv.position().x(), tv.position().y(), tv.position().z());
    
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
        << "Built " << output->size() << " vertices from unified GNN output";
  }
  
  event.put(std::move(output));
}

DEFINE_FWK_MODULE(GNNVertexBuilderFromAlpaka);
