#include "TrackingTools/TransientTrack/interface/TrackTransientTrack.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include <iostream>

/*
 * ThreadSafe statement:
 * This class is using mutable member data: initialTSOS, initialTSCP,
 * trajectoryStateClosestToBeamLine. To guarantee thread safeness we
 * rely on helper member data: m_TSOS, m_TSCP and m_SCTBL, respectively.
 * Each time we'll change mutable member data we rely on specific order of the
 * operator= and the store. It is important since C++11 will guarantee that
 * the value changed by the operator= will be seen by all threads as occuring
 * before the call to store and therefore the kSet == m_TSOS.load is always
 * guaranteed to be true if and only if the thread will see the most recent
 * value of initialTSOS
 */

using namespace reco;

TrackTransientTrack::TrackTransientTrack()
    : Track(),
      tkr_(),
      hasTime(false),
      timeExt_(0.),
      dtErrorExt_(0.),
      theField(nullptr),
      m_TSOS(kUnset),
      m_TSCP(kUnset),
      m_SCTBL(kUnset),
      trkAssoc_(0),
      mtdtime_(0.), mtdtimeErr_(0.), mva_(0.), pathlength_(0.), btlchi2_(0.), btltimechi2_(0.), etlchi2_(0.), etltimechi2_(0.), time_pi_(0.), time_k_(0.), time_p_(0.), sigma_time_pi_(0.), sigma_time_k_(0.), sigma_time_p_(0.), npixbarrel_(0), npixendcap_(0)
{}

TrackTransientTrack::TrackTransientTrack(const Track& tk, const MagneticField* field)
    : Track(tk),
      tkr_(),
      hasTime(false),
      timeExt_(0.),
      dtErrorExt_(0.),
      theField(field),
      m_TSOS(kUnset),
      m_TSCP(kUnset),
      m_SCTBL(kUnset),
      trkAssoc_(0),
      mtdtime_(0.), mtdtimeErr_(0.), mva_(0.), pathlength_(0.), btlchi2_(0.), btltimechi2_(0.), etlchi2_(0.), etltimechi2_(0.), time_pi_(0.), time_k_(0.), time_p_(0.), sigma_time_pi_(0.), sigma_time_k_(0.), sigma_time_p_(0.), npixbarrel_(0), npixendcap_(0){
  initialFTS = trajectoryStateTransform::initialFreeState(tk, field);
}

TrackTransientTrack::TrackTransientTrack(const Track& tk,
                                         const double time,
                                         const double dtime,
                                         const MagneticField* field)
    : Track(tk),
      tkr_(),
      hasTime(true),
      timeExt_(time),
      dtErrorExt_(dtime),
      theField(field),
      m_TSOS(kUnset),
      m_TSCP(kUnset),
      m_SCTBL(kUnset),
      trkAssoc_(0),
      mtdtime_(0.), mtdtimeErr_(0.), mva_(0.), pathlength_(0.), btlchi2_(0.), btltimechi2_(0.), etlchi2_(0.), etltimechi2_(0.), time_pi_(0.), time_k_(0.), time_p_(0.), sigma_time_pi_(0.), sigma_time_k_(0.), sigma_time_p_(0.), npixbarrel_(0), npixendcap_(0){
  initialFTS = trajectoryStateTransform::initialFreeState(tk, field);
}

TrackTransientTrack::TrackTransientTrack(const TrackRef& tk, const MagneticField* field)
    : Track(*tk),
      tkr_(tk),
      hasTime(false),
      timeExt_(0.),
      dtErrorExt_(0.),
      theField(field),
      m_TSOS(kUnset),
      m_TSCP(kUnset),
      m_SCTBL(kUnset),
      trkAssoc_(0),
      mtdtime_(0.), mtdtimeErr_(0.), mva_(0.), pathlength_(0.), btlchi2_(0.), btltimechi2_(0.), etlchi2_(0.), etltimechi2_(0.), time_pi_(0.), time_k_(0.), time_p_(0.), sigma_time_pi_(0.), sigma_time_k_(0.), sigma_time_p_(0.), npixbarrel_(0), npixendcap_(0){
  initialFTS = trajectoryStateTransform::initialFreeState(*tk, field);
}

TrackTransientTrack::TrackTransientTrack(const TrackRef& tk,
                                         const double time,
                                         const double dtime,
                                         const MagneticField* field)
    : Track(*tk),
      tkr_(tk),
      hasTime(true),
      timeExt_(time),
      dtErrorExt_(dtime),
      theField(field),
      m_TSOS(kUnset),
      m_TSCP(kUnset),
      m_SCTBL(kUnset),
      trkAssoc_(0),
      mtdtime_(0.), mtdtimeErr_(0.), mva_(0.), pathlength_(0.), btlchi2_(0.), btltimechi2_(0.), etlchi2_(0.), etltimechi2_(0.), time_pi_(0.), time_k_(0.), time_p_(0.), sigma_time_pi_(0.), sigma_time_k_(0.), sigma_time_p_(0.), npixbarrel_(0), npixendcap_(0){
  initialFTS = trajectoryStateTransform::initialFreeState(*tk, field);
}

TrackTransientTrack::TrackTransientTrack(const Track& tk,
                                         const MagneticField* field,
                                         const edm::ESHandle<GlobalTrackingGeometry>& tg)
    : Track(tk),
      tkr_(),
      hasTime(false),
      timeExt_(0.),
      dtErrorExt_(0.),
      theField(field),
      m_TSOS(kUnset),
      m_TSCP(kUnset),
      m_SCTBL(kUnset),
      theTrackingGeometry(tg),
      trkAssoc_(0),
      mtdtime_(0.), mtdtimeErr_(0.), mva_(0.), pathlength_(0.), btlchi2_(0.), btltimechi2_(0.), etlchi2_(0.), etltimechi2_(0.), time_pi_(0.), time_k_(0.), time_p_(0.), sigma_time_pi_(0.), sigma_time_k_(0.), sigma_time_p_(0.), npixbarrel_(0), npixendcap_(0){
 initialFTS = trajectoryStateTransform::initialFreeState(tk, field);
}

TrackTransientTrack::TrackTransientTrack(const TrackRef& tk,
                                         const MagneticField* field,
                                         const edm::ESHandle<GlobalTrackingGeometry>& tg)
    : Track(*tk),
      tkr_(tk),
      hasTime(false),
      timeExt_(0.),
      dtErrorExt_(0.),
      theField(field),
      m_TSOS(kUnset),
      m_TSCP(kUnset),
      m_SCTBL(kUnset),
      theTrackingGeometry(tg),
      trkAssoc_(0),
      mtdtime_(0.), mtdtimeErr_(0.), mva_(0.), pathlength_(0.), btlchi2_(0.), btltimechi2_(0.), etlchi2_(0.), etltimechi2_(0.), time_pi_(0.), time_k_(0.), time_p_(0.), sigma_time_pi_(0.), sigma_time_k_(0.), sigma_time_p_(0.), npixbarrel_(0), npixendcap_(0){
  initialFTS = trajectoryStateTransform::initialFreeState(*tk, field);
}

TrackTransientTrack::TrackTransientTrack(const Track& tk,
                                         const double time,
                                         const double dtime,
                                         const MagneticField* field,
                                         const edm::ESHandle<GlobalTrackingGeometry>& tg)
    : Track(tk),
      tkr_(),
      hasTime(true),
      timeExt_(time),
      dtErrorExt_(dtime),
      theField(field),
      m_TSOS(kUnset),
      m_TSCP(kUnset),
      m_SCTBL(kUnset),
      theTrackingGeometry(tg),
      trkAssoc_(0),
      mtdtime_(0.), mtdtimeErr_(0.), mva_(0.), pathlength_(0.), btlchi2_(0.), btltimechi2_(0.), etlchi2_(0.), etltimechi2_(0.), time_pi_(0.), time_k_(0.), time_p_(0.), sigma_time_pi_(0.), sigma_time_k_(0.), sigma_time_p_(0.), npixbarrel_(0), npixendcap_(0){
  initialFTS = trajectoryStateTransform::initialFreeState(tk, field);
}

TrackTransientTrack::TrackTransientTrack(const TrackRef& tk,
                                         const double time,
                                         const double dtime,
                                         const MagneticField* field,
                                         const edm::ESHandle<GlobalTrackingGeometry>& tg)
    : Track(*tk),
      tkr_(tk),
      hasTime(true),
      timeExt_(time),
      dtErrorExt_(dtime),
      theField(field),
      m_TSOS(kUnset),
      m_TSCP(kUnset),
      m_SCTBL(kUnset),
      theTrackingGeometry(tg),
      trkAssoc_(0),
      mtdtime_(0.), mtdtimeErr_(0.), mva_(0.), pathlength_(0.), btlchi2_(0.), btltimechi2_(0.), etlchi2_(0.), etltimechi2_(0.), time_pi_(0.), time_k_(0.), time_p_(0.), sigma_time_pi_(0.), sigma_time_k_(0.), sigma_time_p_(0.), npixbarrel_(0), npixendcap_(0){
  initialFTS = trajectoryStateTransform::initialFreeState(*tk, field);
}

TrackTransientTrack::TrackTransientTrack(const TrackTransientTrack& tt)
    : Track(tt),
      tkr_(tt.persistentTrackRef()),
      hasTime(tt.hasTime),
      timeExt_(tt.timeExt()),
      dtErrorExt_(tt.dtErrorExt()),
      theField(tt.field()),
      initialFTS(tt.initialFreeState()),
      m_TSOS(kUnset),
      m_TSCP(kUnset),
      trkAssoc_(0),
      mtdtime_(0), mtdtimeErr_(0), mva_(0), pathlength_(0), btlchi2_(0), btltimechi2_(0), etlchi2_(0), etltimechi2_(0), time_pi_(0), time_k_(0), time_p_(0), sigma_time_pi_(0), sigma_time_k_(0), sigma_time_p_(0), npixbarrel_(0), npixendcap_(0){ 
      //trkAssoc_(tt.trackAsocMTD()), mtdtime_(tt.MTDtime()), mtdtimeErr_(tt.MTDtimeErr()), mva_(tt.MVAquality()), pathlength_(tt.pathLength()), btlchi2_(tt.btlMatch_chi2()), btltimechi2_(tt.btlMatchTime_chi2()), etlchi2_(tt.etlMatch_chi2()), etltimechi2_(tt.etlMatchTime_chi2()), time_pi_(tt.trackTime_pi()), time_k_(tt.trackTime_k()), time_p_(tt.trackTime_p()), sigma_time_pi_(tt.sigma_time_pi()), sigma_time_k_(tt.sigma_time_k()), sigma_time_p_(tt.sigma_time_p()){
  // see ThreadSafe statement above about the order of operator= and store
  if (kSet == tt.m_TSOS.load()) {
    initialTSOS = tt.impactPointState();
    m_TSOS.store(kSet);
  }
  // see ThreadSafe statement above about the order of operator= and store
  if (kSet == tt.m_TSCP.load()) {
    initialTSCP = tt.impactPointTSCP();
    m_TSCP.store(kSet);
  }
}

void TrackTransientTrack::setTrackingGeometry(const edm::ESHandle<GlobalTrackingGeometry>& tg) {
  theTrackingGeometry = tg;
}

void TrackTransientTrack::setBeamSpot(const BeamSpot& beamSpot) {
  theBeamSpot = beamSpot;
  m_SCTBL = kUnset;
}

TrajectoryStateOnSurface TrackTransientTrack::impactPointState() const {
  // see ThreadSafe statement above about the order of operator= and store
  if (kSet == m_TSOS.load())
    return initialTSOS;
  TransverseImpactPointExtrapolator tipe(theField);
  auto tmp = tipe.extrapolate(initialFTS, initialFTS.position());
  char expected = kUnset;
  if (m_TSOS.compare_exchange_strong(expected, kSetting)) {
    initialTSOS = tmp;
    m_TSOS.store(kSet);
    return initialTSOS;
  }
  return tmp;
}

TrajectoryStateClosestToPoint TrackTransientTrack::impactPointTSCP() const {
  // see ThreadSafe statement above about the order of operator= and store
  if (kSet == m_TSCP.load())
    return initialTSCP;
  auto tmp = builder(initialFTS, initialFTS.position());
  char expected = kUnset;
  if (m_TSCP.compare_exchange_strong(expected, kSetting)) {
    initialTSCP = tmp;
    m_TSCP.store(kSet);
    return initialTSCP;
  }
  return tmp;
}

TrajectoryStateOnSurface TrackTransientTrack::outermostMeasurementState() const {
  return trajectoryStateTransform::outerStateOnSurface((*this), *theTrackingGeometry, theField);
}

TrajectoryStateOnSurface TrackTransientTrack::innermostMeasurementState() const {
  return trajectoryStateTransform::innerStateOnSurface((*this), *theTrackingGeometry, theField);
}

TrajectoryStateOnSurface TrackTransientTrack::stateOnSurface(const GlobalPoint& point) const {
  TransverseImpactPointExtrapolator tipe(theField);
  return tipe.extrapolate(initialFTS, point);
}

TrajectoryStateClosestToBeamLine TrackTransientTrack::stateAtBeamLine() const {
  // see ThreadSafe statement above about the order of operator= and store
  if (kSet == m_SCTBL.load())
    return trajectoryStateClosestToBeamLine;
  TSCBLBuilderNoMaterial blsBuilder;
  const auto tmp = blsBuilder(initialFTS, theBeamSpot);
  char expected = kUnset;
  if (m_SCTBL.compare_exchange_strong(expected, kSetting)) {
    trajectoryStateClosestToBeamLine = tmp;
    m_SCTBL.store(kSet);
    return trajectoryStateClosestToBeamLine;
  }
  return tmp;
}

TrackTransientTrack::TrackTransientTrack(const TrackRef& tk,
                                         const double time,
                                         const double dtime,
                                         const MagneticField* field,
                                         const edm::ESHandle<GlobalTrackingGeometry>& tg,
                                         const int trkAssoc,
                                         const float mtdtime, const float mtdtimeErr,
                                         const float mva, const float pathlength, const float btlchi2, const float btltimechi2, const float etlchi2, const float etltimechi2, const float time_pi, const float time_k, const float time_p, const float sigma_time_pi, const float sigma_time_k, const float sigma_time_p, const int npixbarrel, const int npixendcap)
    : Track(*tk),
      tkr_(tk),
      hasTime(true),
      timeExt_(time),
      dtErrorExt_(dtime),
      theField(field),
      m_TSOS(kUnset),
      m_TSCP(kUnset),
      m_SCTBL(kUnset),
      theTrackingGeometry(tg), trkAssoc_(trkAssoc), mtdtime_(mtdtime), mtdtimeErr_(mtdtimeErr), mva_(mva), pathlength_(pathlength), btlchi2_(btlchi2), btltimechi2_(btltimechi2), etlchi2_(etlchi2), etltimechi2_(etltimechi2), time_pi_(time_pi), time_k_(time_k), time_p_(time_p), sigma_time_pi_(sigma_time_pi), sigma_time_k_(sigma_time_k), sigma_time_p_(sigma_time_p), npixbarrel_(npixbarrel), npixendcap_(npixendcap) {
  initialFTS = trajectoryStateTransform::initialFreeState(*tk, field);
}
