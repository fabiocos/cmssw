#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "DataFormats/Common/interface/Handle.h"
#include "TrackingTools/TransientTrack/interface/GsfTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTS.h"
#include "FWCore/Utilities/interface/isFinite.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

using namespace reco;
using namespace std;
using namespace edm;

TransientTrack TransientTrackBuilder::build(const Track* t) const {
  return TransientTrack(*t, theField, theTrackingGeometry);
}

TransientTrack TransientTrackBuilder::build(const Track& t) const {
  return TransientTrack(t, theField, theTrackingGeometry);
}

TransientTrack TransientTrackBuilder::build(const GsfTrack* t) const {
  return TransientTrack(new GsfTransientTrack(*t, theField, theTrackingGeometry));
}

TransientTrack TransientTrackBuilder::build(const GsfTrack& t) const {
  return TransientTrack(new GsfTransientTrack(t, theField, theTrackingGeometry));
}

TransientTrack TransientTrackBuilder::build(const CandidatePtr* t) const {
  reco::PFCandidatePtr tryPF(*t);
  edm::Ptr<pat::PackedCandidate> tryPacked(*t);
  if (tryPF.get() != nullptr && tryPF->isTimeValid()) {
    return TransientTrack(*t, tryPF->time(), tryPF->timeError(), theField, theTrackingGeometry);
  } else if (tryPacked.get() != nullptr && tryPacked->timeError() > 0.f) {
    return TransientTrack(*t, (double)tryPacked->time(), (double)tryPacked->timeError(), theField, theTrackingGeometry);
  }
  return TransientTrack(*t, theField, theTrackingGeometry);
}

TransientTrack TransientTrackBuilder::build(const CandidatePtr& t) const { return this->build(&t); }

TransientTrack TransientTrackBuilder::build(const TrackRef* t) const {
  return TransientTrack(*t, theField, theTrackingGeometry);
}

TransientTrack TransientTrackBuilder::build(const TrackRef& t) const {
  return TransientTrack(t, theField, theTrackingGeometry);
}

TransientTrack TransientTrackBuilder::build(const GsfTrackRef* t) const {
  return TransientTrack(new GsfTransientTrack(*t, theField, theTrackingGeometry));
}

TransientTrack TransientTrackBuilder::build(const GsfTrackRef& t) const {
  return TransientTrack(new GsfTransientTrack(t, theField, theTrackingGeometry));
}

vector<TransientTrack> TransientTrackBuilder::build(const edm::Handle<reco::TrackCollection>& trkColl) const {
  vector<TransientTrack> ttVect;
  ttVect.reserve((*trkColl).size());
  for (unsigned int i = 0; i < (*trkColl).size(); i++) {
    ttVect.push_back(TransientTrack(TrackRef(trkColl, i), theField, theTrackingGeometry));
  }
  return ttVect;
}

vector<TransientTrack> TransientTrackBuilder::build(const edm::Handle<reco::GsfTrackCollection>& trkColl) const {
  vector<TransientTrack> ttVect;
  ttVect.reserve((*trkColl).size());
  for (unsigned int i = 0; i < (*trkColl).size(); i++) {
    ttVect.push_back(TransientTrack(new GsfTransientTrack(GsfTrackRef(trkColl, i), theField, theTrackingGeometry)));
  }
  return ttVect;
}

vector<TransientTrack> TransientTrackBuilder::build(const edm::Handle<edm::View<Track> >& trkColl) const {
  vector<TransientTrack> ttVect;
  ttVect.reserve((*trkColl).size());
  for (unsigned int i = 0; i < (*trkColl).size(); i++) {
    const Track* trk = &(*trkColl)[i];
    const GsfTrack* gsfTrack = dynamic_cast<const GsfTrack*>(trk);
    if (gsfTrack) {
      ttVect.push_back(TransientTrack(
          new GsfTransientTrack(RefToBase<Track>(trkColl, i).castTo<GsfTrackRef>(), theField, theTrackingGeometry)));
    } else {  // gsf
      ttVect.push_back(TransientTrack(RefToBase<Track>(trkColl, i).castTo<TrackRef>(), theField, theTrackingGeometry));
    }
  }
  return ttVect;
}

vector<TransientTrack> TransientTrackBuilder::build(const edm::Handle<reco::TrackCollection>& trkColl,
                                                    const edm::ValueMap<float>& trackTimes,
                                                    const edm::ValueMap<float>& trackTimeResos) const {
  vector<TransientTrack> ttVect;
  ttVect.reserve((*trkColl).size());
  for (unsigned int i = 0; i < (*trkColl).size(); i++) {
    TrackRef ref(trkColl, i);
    double time = trackTimes[ref];
    double timeReso = trackTimeResos[ref];
    timeReso = (timeReso > 1e-6 ? timeReso
                                : defaultInvalidTrackTimeReso);  // make the error much larger than the BS time width
    if (edm::isNotFinite(time)) {
      time = 0.0;
      timeReso = defaultInvalidTrackTimeReso;
    }
    ttVect.push_back(TransientTrack(ref, time, timeReso, theField, theTrackingGeometry));
  }
  return ttVect;
}

vector<TransientTrack> TransientTrackBuilder::build(const edm::Handle<reco::GsfTrackCollection>& trkColl,
                                                    const edm::ValueMap<float>& trackTimes,
                                                    const edm::ValueMap<float>& trackTimeResos) const {
  vector<TransientTrack> ttVect;
  ttVect.reserve((*trkColl).size());
  for (unsigned int i = 0; i < (*trkColl).size(); i++) {
    GsfTrackRef ref(trkColl, i);
    double time = trackTimes[ref];
    double timeReso = trackTimeResos[ref];
    timeReso = (timeReso > 1e-6 ? timeReso
                                : defaultInvalidTrackTimeReso);  // make the error much larger than the BS time width
    if (edm::isNotFinite(time)) {
      time = 0.0;
      timeReso = defaultInvalidTrackTimeReso;
    }
    ttVect.push_back(TransientTrack(new GsfTransientTrack(ref, time, timeReso, theField, theTrackingGeometry)));
  }
  return ttVect;
}

vector<TransientTrack> TransientTrackBuilder::build(const edm::Handle<edm::View<Track> >& trkColl,
                                                    const edm::ValueMap<float>& trackTimes,
                                                    const edm::ValueMap<float>& trackTimeResos) const {
  vector<TransientTrack> ttVect;
  ttVect.reserve((*trkColl).size());
  for (unsigned int i = 0; i < (*trkColl).size(); i++) {
    const Track* trk = &(*trkColl)[i];
    const GsfTrack* gsfTrack = dynamic_cast<const GsfTrack*>(trk);
    if (gsfTrack) {
      GsfTrackRef ref = RefToBase<Track>(trkColl, i).castTo<GsfTrackRef>();
      double time = trackTimes[ref];
      double timeReso = trackTimeResos[ref];
      timeReso = (timeReso > 1e-6 ? timeReso
                                  : defaultInvalidTrackTimeReso);  // make the error much larger than the BS time width
      if (edm::isNotFinite(time)) {
        time = 0.0;
        timeReso = defaultInvalidTrackTimeReso;
      }
      ttVect.push_back(TransientTrack(new GsfTransientTrack(
          RefToBase<Track>(trkColl, i).castTo<GsfTrackRef>(), time, timeReso, theField, theTrackingGeometry)));
    } else {  // gsf
      TrackRef ref = RefToBase<Track>(trkColl, i).castTo<TrackRef>();
      double time = trackTimes[ref];
      double timeReso = trackTimeResos[ref];
      timeReso = (timeReso > 1e-6 ? timeReso
                                  : defaultInvalidTrackTimeReso);  // make the error much larger than the BS time width
      if (edm::isNotFinite(time)) {
        time = 0.0;
        timeReso = defaultInvalidTrackTimeReso;
      }
      ttVect.push_back(TransientTrack(
          RefToBase<Track>(trkColl, i).castTo<TrackRef>(), time, timeReso, theField, theTrackingGeometry));
    }
  }
  return ttVect;
}

vector<TransientTrack> TransientTrackBuilder::build(const edm::Handle<reco::TrackCollection>& trkColl,
                                                    const reco::BeamSpot& beamSpot) const {
  vector<TransientTrack> ttVect = build(trkColl);
  for (unsigned int i = 0; i < ttVect.size(); i++) {
    ttVect[i].setBeamSpot(beamSpot);
  }
  return ttVect;
}

vector<TransientTrack> TransientTrackBuilder::build(const edm::Handle<reco::GsfTrackCollection>& trkColl,
                                                    const reco::BeamSpot& beamSpot) const {
  vector<TransientTrack> ttVect = build(trkColl);
  for (unsigned int i = 0; i < ttVect.size(); i++) {
    ttVect[i].setBeamSpot(beamSpot);
  }
  return ttVect;
}

vector<TransientTrack> TransientTrackBuilder::build(const edm::Handle<edm::View<Track> >& trkColl,
                                                    const reco::BeamSpot& beamSpot) const {
  vector<TransientTrack> ttVect = build(trkColl);
  for (unsigned int i = 0; i < ttVect.size(); i++) {
    ttVect[i].setBeamSpot(beamSpot);
  }
  return ttVect;
}

vector<TransientTrack> TransientTrackBuilder::build(const edm::Handle<reco::TrackCollection>& trkColl,
                                                    const reco::BeamSpot& beamSpot,
                                                    const edm::ValueMap<float>& trackTimes,
                                                    const edm::ValueMap<float>& trackTimeResos) const {
  vector<TransientTrack> ttVect = build(trkColl, trackTimes, trackTimeResos);
  for (unsigned int i = 0; i < ttVect.size(); i++) {
    ttVect[i].setBeamSpot(beamSpot);
  }
  return ttVect;
}

vector<TransientTrack> TransientTrackBuilder::build(const edm::Handle<reco::GsfTrackCollection>& trkColl,
                                                    const reco::BeamSpot& beamSpot,
                                                    const edm::ValueMap<float>& trackTimes,
                                                    const edm::ValueMap<float>& trackTimeResos) const {
  vector<TransientTrack> ttVect = build(trkColl, trackTimes, trackTimeResos);
  for (unsigned int i = 0; i < ttVect.size(); i++) {
    ttVect[i].setBeamSpot(beamSpot);
  }
  return ttVect;
}

vector<TransientTrack> TransientTrackBuilder::build(const edm::Handle<edm::View<Track> >& trkColl,
                                                    const reco::BeamSpot& beamSpot,
                                                    const edm::ValueMap<float>& trackTimes,
                                                    const edm::ValueMap<float>& trackTimeResos) const {
  vector<TransientTrack> ttVect = build(trkColl, trackTimes, trackTimeResos);
  for (unsigned int i = 0; i < ttVect.size(); i++) {
    ttVect[i].setBeamSpot(beamSpot);
  }
  return ttVect;
}

TransientTrack TransientTrackBuilder::build(const FreeTrajectoryState& fts) const {
  return TransientTrack(new TransientTrackFromFTS(fts));
}


vector<TransientTrack> TransientTrackBuilder::build(const edm::Handle<reco::TrackCollection>& trkColl,
						    const reco::BeamSpot& beamSpot,
                                                    const edm::ValueMap<float>& trackTimes,
                                                    const edm::ValueMap<float>& trackTimeResos,
                                                    const edm::ValueMap<int>& trackMTDAssoc,
                                                    const edm::ValueMap<float>& MTDtime,
                                                    const edm::ValueMap<float>& sigmaMTDtime,
                                                    const edm::ValueMap<float>& MVAquality,
                                                    const edm::ValueMap<float>& pathLength,
                                                    const edm::ValueMap<float>& btlmatchChi2,
                                                    const edm::ValueMap<float>& btlmatchTime_Chi2,
                                                    const edm::ValueMap<float>& etlmatchChi2,
                                                    const edm::ValueMap<float>& etlmatchTime_Chi2,
                                                    const edm::ValueMap<float>& trackTime_pi,
                                                    const edm::ValueMap<float>& trackTime_k,
                                                    const edm::ValueMap<float>& trackTime_p,
                                                    const edm::ValueMap<float>& sigmatrackTime_pi,
                                                    const edm::ValueMap<float>& sigmatrackTime_k,
                                                    const edm::ValueMap<float>& sigmatrackTime_p,
                                                    const edm::ValueMap<int>& tracknpixBarrel,
                                                    const edm::ValueMap<int>& tracknpixEndcap ) const {
  vector<TransientTrack> ttVect;
  ttVect.reserve((*trkColl).size());
  for (unsigned int i = 0; i < (*trkColl).size(); i++) {
    TrackRef ref(trkColl, i);
    double time = trackTimes[ref];
    double timeReso = trackTimeResos[ref];
    int trackAsoc_mtd = trackMTDAssoc[ref];
    float MTDtime_mtd = MTDtime[ref];
    float MTDtimeErr_mtd = sigmaMTDtime[ref];
    float MVAquality_mtd = MVAquality[ref];
    float pathLength_mtd = pathLength[ref];
    float btlMatch_chi2_mtd = btlmatchChi2[ref];
    float btlMatchTime_chi2_mtd = btlmatchTime_Chi2[ref];
    float etlMatch_chi2_mtd = etlmatchChi2[ref];
    float etlMatchTime_chi2_mtd = etlmatchTime_Chi2[ref];
    float trackTime_pi_mtd = trackTime_pi[ref];
    float trackTime_k_mtd = trackTime_k[ref];
    float trackTime_p_mtd = trackTime_p[ref];
    float track_sigmaTime_pi_mtd = sigmatrackTime_pi[ref];
    float track_sigmaTime_k_mtd = sigmatrackTime_k[ref];
    float track_sigmaTime_p_mtd = sigmatrackTime_p[ref];
    int npixBarrel = tracknpixBarrel[ref];
    int npixEndcap = tracknpixEndcap[ref];


    timeReso = (timeReso > 1e-6 ? timeReso
                                : defaultInvalidTrackTimeReso);  // make the error much larger than the BS time width
    if (edm::isNotFinite(time)) {
      time = 0.0;
      timeReso = defaultInvalidTrackTimeReso;
    }
    ttVect.push_back(TransientTrack(ref, time, timeReso, theField,  theTrackingGeometry, trackAsoc_mtd, MTDtime_mtd, MTDtimeErr_mtd, MVAquality_mtd, pathLength_mtd, btlMatch_chi2_mtd, btlMatchTime_chi2_mtd, etlMatch_chi2_mtd, etlMatchTime_chi2_mtd, trackTime_pi_mtd, trackTime_k_mtd, trackTime_p_mtd, track_sigmaTime_pi_mtd, track_sigmaTime_k_mtd, track_sigmaTime_p_mtd,npixBarrel, npixEndcap));
  }
  for (unsigned int i = 0; i < ttVect.size(); i++) {
    ttVect[i].setBeamSpot(beamSpot);
  }
  return ttVect;
}
