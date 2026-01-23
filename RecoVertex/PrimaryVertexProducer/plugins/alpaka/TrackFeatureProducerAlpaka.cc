/**
 * TrackFeatureProducerAlpaka - Prepares track features for GNN inference
 *
 * Converts TransientTracks to SoA format with 13 features per track.
 */

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
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
#include "RecoVertex/PrimaryVertexProducer/interface/VertexGNNSoA.h"
#include "RecoVertex/PrimaryVertexProducer/interface/alpaka/VertexGNNDeviceCollection.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include <cmath>
#include <vector>

namespace ALPAKA_ACCELERATOR_NAMESPACE::vertexgnn {

  class TrackFeatureProducerAlpaka : public stream::EDProducer<> {
  public:
    explicit TrackFeatureProducerAlpaka(const edm::ParameterSet& params);

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    void produce(device::Event& event, const device::EventSetup& eventSetup) override;

  private:
    // Helper to combine sigma values
    static float combineSigma(float s1, float s2) {
      return std::sqrt(s1 * s1 + s2 * s2);
    }

    // Configuration
    const int verbosity_;

    // Standard EDM tokens (non-Alpaka inputs)
    const edm::EDGetTokenT<std::vector<reco::TransientTrack>> tracksToken_;
    
    // Alpaka output token
    const device::EDPutToken<TrackFeaturesDeviceCollection> featuresToken_;
  };

  TrackFeatureProducerAlpaka::TrackFeatureProducerAlpaka(const edm::ParameterSet& params)
      : EDProducer<>(params),
        verbosity_(params.getUntrackedParameter<int>("verbosity", 0)),
        tracksToken_(consumes<std::vector<reco::TransientTrack>>(params.getParameter<edm::InputTag>("tracks"))),
        featuresToken_{produces()} {}

  void TrackFeatureProducerAlpaka::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("tracks", edm::InputTag("unsortedOfflinePrimaryVertices", ""));
    desc.addUntracked<int>("verbosity", 0);
    descriptions.addWithDefaultLabel(desc);
  }

  void TrackFeatureProducerAlpaka::produce(device::Event& event, const device::EventSetup& eventSetup) {
    // Get tracks from standard EDM event
    edm::Handle<std::vector<reco::TransientTrack>> tracksHandle;
    event.get(tracksToken_, tracksHandle);
    const auto& tracks = *tracksHandle;
    const auto N = tracks.size();

    if (verbosity_ > 0) {
      edm::LogInfo("TrackFeatureProducerAlpaka") << "Processing " << N << " tracks";
    }

    // Allocate output collection
    auto features = TrackFeaturesDeviceCollection(N, event.queue());
    auto view = features.view();

    // Fill features from tracks (on host, then transferred)
    for (size_t i = 0; i < N; ++i) {
      const auto& tt = tracks[i];
      const auto& track = tt.track();

      // Basic track properties
      float track_vz = track.vz();
      float track_dz = track.dzError();
      float track_pt = track.pt();
      float track_eta = track.eta();
      
      // Quality indicators
      float track_mva = 1.0f;  // Placeholder - would come from MVA selector
      float track_pl = std::abs(track.dz() / track.dzError());

      // Timing features (placeholders - require MTD info)
      // In full implementation, these would come from MTD track extras
      float t_pi = 0.0f, t_k = 0.0f, t_p = 0.0f;
      float s_pi = 0.2f, s_k = 0.2f, s_p = 0.2f;
      float has_time = 0.0f;

      // Store features
      view[i].vz() = std::isfinite(track_vz) ? track_vz : 1e9f;
      view[i].dz() = std::isfinite(track_dz) ? track_dz : 1e9f;
      view[i].pt() = std::isfinite(track_pt) ? track_pt : 1e9f;
      view[i].eta() = std::isfinite(track_eta) ? track_eta : 1e9f;
      view[i].mva() = track_mva;
      view[i].pl() = std::isfinite(track_pl) ? track_pl : 1e9f;
      view[i].t_pi() = t_pi;
      view[i].t_k() = t_k;
      view[i].t_p() = t_p;
      view[i].s_pi() = s_pi;
      view[i].s_k() = s_k;
      view[i].s_p() = s_p;
      view[i].has_time() = has_time;
    }

    // Put output into event
    event.emplace(featuresToken_, std::move(features));
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::vertexgnn

DEFINE_FWK_ALPAKA_MODULE(vertexgnn::TrackFeatureProducerAlpaka);
