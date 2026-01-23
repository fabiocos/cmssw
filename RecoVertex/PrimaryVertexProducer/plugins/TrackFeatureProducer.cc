/**
 * TrackFeatureProducer - Standard EDProducer for track features
 *
 * Consumes real generalTracks with MTD timing and produces TrackFeaturesHostCollection
 * with 13 features per track (identical to GNNClusterizer.cc feature extraction).
 *
 * Features: vz, dz, pt, eta, mva, pl, t_pi, t_k, t_p, s_pi, s_k, s_p, has_time
 *
 * Note: This is a standard EDProducer (not Alpaka) because:
 *   1. Feature extraction is CPU-bound, no GPU needed
 *   2. Standard EDM products (TrackCollection, ValueMaps) can't be consumed by Alpaka device::Event
 *
 * The output HostCollection can be consumed by downstream Alpaka producers.
 */

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "HeterogeneousCore/AlpakaInterface/interface/host.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexGNNSoA.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexGNNHostCollection.h"

#include <cmath>

namespace vertexgnn {

  /**
   * Combine timing uncertainties in quadrature (same as GNNClusterizer.cc)
   */
  inline float combine_sigma(float sig_h, float sig_tmtd, float eps = 1e-12f) {
    float h = std::max(sig_h, 0.0f);
    float t = std::max(sig_tmtd, 0.0f);
    float out = std::sqrt(h * h + t * t + eps);
    if (!std::isfinite(out)) {
      out = 1e9f;
    }
    return out;
  }

  class TrackFeatureProducer : public edm::stream::EDProducer<> {
  public:
    explicit TrackFeatureProducer(const edm::ParameterSet& params)
        : trackToken_(consumes<reco::TrackCollection>(params.getParameter<edm::InputTag>("tracks"))),
          beamSpotToken_(consumes<reco::BeamSpot>(params.getParameter<edm::InputTag>("beamSpot"))),
          trkTimesToken_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("trackTimesLabel"))),
          trkTimeResosToken_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("trackTimeResosLabel"))),
          trkMTDAssocToken_(consumes<edm::ValueMap<int>>(params.getParameter<edm::InputTag>("trackAssocSrc"))),
          MTDtimeToken_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("tmtdSrc"))),
          sigmaMTDtimeToken_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("sigmatmtdSrc"))),
          pathLengthToken_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("pathmtd"))),
          btlMatchChi2Token_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("btlMatchChi2Src"))),
          btlMatchTime_Chi2Token_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("btlMatchTimeChi2Src"))),
          etlMatchChi2Token_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("etlMatchChi2Src"))),
          etlMatchTime_Chi2Token_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("etlMatchTimeChi2Src"))),
          trkTimePiToken_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("tofPi"))),
          trkTimeKToken_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("tofK"))),
          trkTimePToken_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("tofP"))),
          sigmaTrkTimePiToken_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("sigmatofpiSrc"))),
          sigmaTrkTimeKToken_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("sigmatofkSrc"))),
          sigmaTrkTimePToken_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("sigmatofpSrc"))),
          npixBarrelToken_(consumes<edm::ValueMap<int>>(params.getParameter<edm::InputTag>("npixBarrelSrc"))),
          npixEndcapToken_(consumes<edm::ValueMap<int>>(params.getParameter<edm::InputTag>("npixEndcapSrc"))),
          trackMTDTimeQualityToken_(consumes<edm::ValueMap<float>>(params.getParameter<edm::InputTag>("trackMTDTimeQualityVMapTag"))),
          ttbToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
          minTrackTimeQuality_(params.getParameter<double>("minTrackTimeQuality")),
          useMVASelection_(params.getParameter<bool>("useMVACut")),
          verbose_(params.getUntrackedParameter<bool>("verbose", false)) {
      produces<TrackFeaturesHostCollection>();
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
      desc.add<edm::InputTag>("beamSpot", edm::InputTag("offlineBeamSpot"));
      desc.add<edm::InputTag>("trackTimesLabel", edm::InputTag("tofPID4DnoPID:t0safe"));
      desc.add<edm::InputTag>("trackTimeResosLabel", edm::InputTag("tofPID4DnoPID:sigmat0safe"));
      desc.add<edm::InputTag>("trackMTDTimeQualityVMapTag", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"));
      desc.add<edm::InputTag>("trackAssocSrc", edm::InputTag("trackExtenderWithMTD:generalTrackassoc"));
      desc.add<edm::InputTag>("tmtdSrc", edm::InputTag("trackExtenderWithMTD:generalTracktmtd"));
      desc.add<edm::InputTag>("sigmatmtdSrc", edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd"));
      desc.add<edm::InputTag>("pathmtd", edm::InputTag("trackExtenderWithMTD:generalTrackPathLength"));
      desc.add<edm::InputTag>("btlMatchChi2Src", edm::InputTag("trackExtenderWithMTD", "btlMatchChi2"));
      desc.add<edm::InputTag>("btlMatchTimeChi2Src", edm::InputTag("trackExtenderWithMTD", "btlMatchTimeChi2"));
      desc.add<edm::InputTag>("etlMatchChi2Src", edm::InputTag("trackExtenderWithMTD", "etlMatchChi2"));
      desc.add<edm::InputTag>("etlMatchTimeChi2Src", edm::InputTag("trackExtenderWithMTD", "etlMatchTimeChi2"));
      desc.add<edm::InputTag>("tofPi", edm::InputTag("trackExtenderWithMTD:generalTrackTofPi"));
      desc.add<edm::InputTag>("tofK", edm::InputTag("trackExtenderWithMTD:generalTrackTofK"));
      desc.add<edm::InputTag>("tofP", edm::InputTag("trackExtenderWithMTD:generalTrackTofP"));
      desc.add<edm::InputTag>("sigmatofpiSrc", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofPi"));
      desc.add<edm::InputTag>("sigmatofkSrc", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofK"));
      desc.add<edm::InputTag>("sigmatofpSrc", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofP"));
      desc.add<edm::InputTag>("npixBarrelSrc", edm::InputTag("trackExtenderWithMTD", "npixBarrel"));
      desc.add<edm::InputTag>("npixEndcapSrc", edm::InputTag("trackExtenderWithMTD", "npixEndcap"));
      desc.add<double>("minTrackTimeQuality", 0.8);
      desc.add<bool>("useMVACut", false);
      desc.addUntracked<bool>("verbose", false);
      descriptions.addWithDefaultLabel(desc);
    }

    void produce(edm::Event& event, const edm::EventSetup& eventSetup) override {
      // Get track handle
      auto trackHandle = event.getHandle(trackToken_);
      auto const& beamSpot = event.get(beamSpotToken_);
      
      // Get timing ValueMaps
      auto const& trackTimes = event.get(trkTimesToken_);
      auto const& trackTimeResos = event.get(trkTimeResosToken_);
      auto const& trackMTDAssoc = event.get(trkMTDAssocToken_);
      auto const& trackMTDTimes = event.get(MTDtimeToken_);
      auto const& trackMTDTimesRes = event.get(sigmaMTDtimeToken_);
      auto const& pathlength = event.get(pathLengthToken_);
      auto const& btlMatchChi2 = event.get(btlMatchChi2Token_);
      auto const& btlMatchTimeChi2 = event.get(btlMatchTime_Chi2Token_);
      auto const& etlMatchChi2 = event.get(etlMatchChi2Token_);
      auto const& etlMatchTimeChi2 = event.get(etlMatchTime_Chi2Token_);
      auto const& trkPiTime = event.get(trkTimePiToken_);
      auto const& trkKTime = event.get(trkTimeKToken_);
      auto const& trkPTime = event.get(trkTimePToken_);
      auto const& sigmatrkPiTime = event.get(sigmaTrkTimePiToken_);
      auto const& sigmatrkKTime = event.get(sigmaTrkTimeKToken_);
      auto const& sigmatrkPTime = event.get(sigmaTrkTimePToken_);
      auto const& npixBarrel = event.get(npixBarrelToken_);
      auto const& npixEndcap = event.get(npixEndcapToken_);
      auto const& trackMTDTimeQualities = event.get(trackMTDTimeQualityToken_);

      // Get TransientTrackBuilder
      auto const& ttBuilder = eventSetup.getData(ttbToken_);

      // Build all TransientTracks with timing info using the batch method
      std::vector<reco::TransientTrack> t_tks = ttBuilder.build(
          trackHandle,
          beamSpot,
          trackTimes,
          trackTimeResos,
          trackMTDAssoc,
          trackMTDTimes,
          trackMTDTimesRes,
          trackMTDTimeQualities,
          pathlength,
          btlMatchChi2,
          btlMatchTimeChi2,
          etlMatchChi2,
          etlMatchTimeChi2,
          trkPiTime,
          trkKTime,
          trkPTime,
          sigmatrkPiTime,
          sigmatrkKTime,
          sigmatrkPTime,
          npixBarrel,
          npixEndcap);

      const int N = t_tks.size();

      if (verbose_) {
        edm::LogInfo("TrackFeatureProducer") << "Processing " << N << " tracks";
      }

      // Allocate output HostCollection
      auto trackFeatures = std::make_unique<TrackFeaturesHostCollection>(N, cms::alpakatools::host());
      auto view = trackFeatures->view();

      // Extract features from each TransientTrack
      for (int i = 0; i < N; ++i) {
        const reco::TransientTrack& ttrack = t_tks[i];

        // Extract features (same logic as GNNClusterizer.cc)
        float track_vz = ttrack.track().vz();
        float track_dz = ttrack.track().dzError();
        float track_pt = ttrack.track().pt();
        float track_eta = ttrack.track().eta();
        float track_mva = ttrack.MVAquality();
        float track_pl = ttrack.pathLength();

        // Ensure missing MVA/PL are -1.0f
        if (!std::isfinite(track_mva)) track_mva = -1.0f;
        if (!std::isfinite(track_pl)) track_pl = -1.0f;

        // Timing information
        float track_t_pi = ttrack.MTDtime() - ttrack.trackTime_pi();
        float track_t_k = ttrack.MTDtime() - ttrack.trackTime_k();
        float track_t_p = ttrack.MTDtime() - ttrack.trackTime_p();
        float track_s_tmtd = ttrack.MTDtimeErr();
        float track_s_pi = ttrack.sigma_time_pi();
        float track_s_k = ttrack.sigma_time_k();
        float track_s_p = ttrack.sigma_time_p();

        // Build timing masks
        bool tmtd_ok = (track_s_tmtd >= 0.0f);
        bool pi_ok = (track_s_pi >= 0.0f) && tmtd_ok;
        bool k_ok = (track_s_k >= 0.0f) && tmtd_ok;
        bool p_ok = (track_s_p >= 0.0f) && tmtd_ok;

        // Apply masks to timing features
        float feat_t_pi = 0.0f, feat_t_k = 0.0f, feat_t_p = 0.0f;
        float feat_s_pi = 0.2f, feat_s_k = 0.2f, feat_s_p = 0.2f;

        if (pi_ok) {
          feat_t_pi = track_t_pi;
          feat_s_pi = combine_sigma(track_s_pi, track_s_tmtd);
        }
        if (k_ok) {
          feat_t_k = track_t_k;
          feat_s_k = combine_sigma(track_s_k, track_s_tmtd);
        }
        if (p_ok) {
          feat_t_p = track_t_p;
          feat_s_p = combine_sigma(track_s_p, track_s_tmtd);
        }

        // has_time flag
        float has_time = tmtd_ok ? 1.0f : 0.0f;

        // Replace NaN/Inf
        auto sanitize = [](float v) { return std::isfinite(v) ? v : 1e9f; };

        // Store features in SoA
        view[i].vz() = sanitize(track_vz);
        view[i].dz() = sanitize(track_dz);
        view[i].pt() = sanitize(track_pt);
        view[i].eta() = sanitize(track_eta);
        view[i].mva() = sanitize(track_mva);
        view[i].pl() = sanitize(track_pl);
        view[i].t_pi() = sanitize(feat_t_pi);
        view[i].t_k() = sanitize(feat_t_k);
        view[i].t_p() = sanitize(feat_t_p);
        view[i].s_pi() = sanitize(feat_s_pi);
        view[i].s_k() = sanitize(feat_s_k);
        view[i].s_p() = sanitize(feat_s_p);
        view[i].has_time() = has_time;
      }

      if (verbose_) {
        edm::LogInfo("TrackFeatureProducer") << "Produced " << N << " track features";
      }

      // Put output into event
      event.put(std::move(trackFeatures));
    }

  private:
    // Input tokens
    const edm::EDGetTokenT<reco::TrackCollection> trackToken_;
    const edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
    const edm::EDGetTokenT<edm::ValueMap<float>> trkTimesToken_;
    const edm::EDGetTokenT<edm::ValueMap<float>> trkTimeResosToken_;
    const edm::EDGetTokenT<edm::ValueMap<int>> trkMTDAssocToken_;
    const edm::EDGetTokenT<edm::ValueMap<float>> MTDtimeToken_;
    const edm::EDGetTokenT<edm::ValueMap<float>> sigmaMTDtimeToken_;
    const edm::EDGetTokenT<edm::ValueMap<float>> pathLengthToken_;
    const edm::EDGetTokenT<edm::ValueMap<float>> btlMatchChi2Token_;
    const edm::EDGetTokenT<edm::ValueMap<float>> btlMatchTime_Chi2Token_;
    const edm::EDGetTokenT<edm::ValueMap<float>> etlMatchChi2Token_;
    const edm::EDGetTokenT<edm::ValueMap<float>> etlMatchTime_Chi2Token_;
    const edm::EDGetTokenT<edm::ValueMap<float>> trkTimePiToken_;
    const edm::EDGetTokenT<edm::ValueMap<float>> trkTimeKToken_;
    const edm::EDGetTokenT<edm::ValueMap<float>> trkTimePToken_;
    const edm::EDGetTokenT<edm::ValueMap<float>> sigmaTrkTimePiToken_;
    const edm::EDGetTokenT<edm::ValueMap<float>> sigmaTrkTimeKToken_;
    const edm::EDGetTokenT<edm::ValueMap<float>> sigmaTrkTimePToken_;
    const edm::EDGetTokenT<edm::ValueMap<int>> npixBarrelToken_;
    const edm::EDGetTokenT<edm::ValueMap<int>> npixEndcapToken_;
    const edm::EDGetTokenT<edm::ValueMap<float>> trackMTDTimeQualityToken_;
    const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbToken_;
    
    // Configuration
    const double minTrackTimeQuality_;
    const bool useMVASelection_;
    const bool verbose_;
  };

}  // namespace vertexgnn

DEFINE_FWK_MODULE(vertexgnn::TrackFeatureProducer);
