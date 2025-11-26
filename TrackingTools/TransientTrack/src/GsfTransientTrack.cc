#include "TrackingTools/TransientTrack/interface/GsfTransientTrack.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "TrackingTools/GsfTools/interface/GsfPropagatorAdapter.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include <iostream>

using namespace reco;
using namespace std;

GsfTransientTrack::GsfTransientTrack()
    : GsfTrack(),
      tkr_(),
      hasTime(false),
      timeExt_(0.),
      dtErrorExt_(0.),
      theField(nullptr),
      initialTSOSAvailable(false),
      initialTSCPAvailable(false),
      blStateAvailable(false),
      theTIPExtrapolator(), 
      trkAssoc_(0.),
      mtdtime_(0.), mtdtimeErr_(0.), mva_(0.), pathlength_(0.), btlchi2_(0.), btltimechi2_(0.), etlchi2_(0.), etltimechi2_(0.), time_pi_(0.), time_k_(0.), time_p_(0.), sigma_time_pi_(0.), sigma_time_k_(0.), sigma_time_p_(0.){}

GsfTransientTrack::GsfTransientTrack(const GsfTrack& tk, const MagneticField* field)
    : GsfTrack(tk),
      tkr_(),
      hasTime(false),
      timeExt_(0.),
      dtErrorExt_(0.),
      theField(field),
      initialTSOSAvailable(false),
      initialTSCPAvailable(false),
      blStateAvailable(false), 
      trkAssoc_(0),
      mtdtime_(0.), mtdtimeErr_(0.), mva_(0.), pathlength_(0.), btlchi2_(0.), btltimechi2_(0.), etlchi2_(0.), etltimechi2_(0.), time_pi_(0.), time_k_(0.), time_p_(0.), sigma_time_pi_(0.), sigma_time_k_(0.), sigma_time_p_(0.){
  initialFTS = trajectoryStateTransform::initialFreeState(tk, field);
}

GsfTransientTrack::GsfTransientTrack(const GsfTrack& tk,
                                     const double time,
                                     const double dtime,
                                     const MagneticField* field)
    : GsfTrack(tk),
      tkr_(),
      hasTime(true),
      timeExt_(time),
      dtErrorExt_(dtime),
      theField(field),
      initialTSOSAvailable(false),
      initialTSCPAvailable(false),
      blStateAvailable(false), 
      trkAssoc_(0),
      mtdtime_(0.), mtdtimeErr_(0.), mva_(0.), pathlength_(0.), btlchi2_(0.), btltimechi2_(0.), etlchi2_(0.), etltimechi2_(0.), time_pi_(0.), time_k_(0.), time_p_(0.), sigma_time_pi_(0.), sigma_time_k_(0.), sigma_time_p_(0.){
  initialFTS = trajectoryStateTransform::initialFreeState(tk, field);
}

GsfTransientTrack::GsfTransientTrack(const GsfTrackRef& tk, const MagneticField* field)
    : GsfTrack(*tk),
      tkr_(tk),
      hasTime(false),
      timeExt_(0.),
      dtErrorExt_(0.),
      theField(field),
      initialTSOSAvailable(false),
      initialTSCPAvailable(false),
      blStateAvailable(false),
      theTIPExtrapolator(AnalyticalPropagator(field, alongMomentum)), 
      trkAssoc_(0),
      mtdtime_(0.), mtdtimeErr_(0.), mva_(0.), pathlength_(0.), btlchi2_(0.), btltimechi2_(0.), etlchi2_(0.), etltimechi2_(0.), time_pi_(0.), time_k_(0.), time_p_(0.), sigma_time_pi_(0.), sigma_time_k_(0.), sigma_time_p_(0.){
  initialFTS = trajectoryStateTransform::initialFreeState(*tk, field);
}

GsfTransientTrack::GsfTransientTrack(const GsfTrackRef& tk,
                                     const double time,
                                     const double dtime,
                                     const MagneticField* field)
    : GsfTrack(*tk),
      tkr_(tk),
      hasTime(true),
      timeExt_(time),
      dtErrorExt_(dtime),
      theField(field),
      initialTSOSAvailable(false),
      initialTSCPAvailable(false),
      blStateAvailable(false),
      theTIPExtrapolator(AnalyticalPropagator(field, alongMomentum)), 
      trkAssoc_(0),
      mtdtime_(0.), mtdtimeErr_(0.), mva_(0.), pathlength_(0.), btlchi2_(0.), btltimechi2_(0.), etlchi2_(0.), etltimechi2_(0.), time_pi_(0.), time_k_(0.), time_p_(0.), sigma_time_pi_(0.), sigma_time_k_(0.), sigma_time_p_(0.) {
  initialFTS = trajectoryStateTransform::initialFreeState(*tk, field);
}

GsfTransientTrack::GsfTransientTrack(const GsfTrack& tk,
                                     const MagneticField* field,
                                     const edm::ESHandle<GlobalTrackingGeometry>& tg)
    : GsfTrack(tk),
      tkr_(),
      hasTime(false),
      timeExt_(0.),
      dtErrorExt_(0.),
      theField(field),
      initialTSOSAvailable(false),
      initialTSCPAvailable(false),
      blStateAvailable(false),
      theTrackingGeometry(tg),
      theTIPExtrapolator(AnalyticalPropagator(field, alongMomentum)),
      trkAssoc_(0),
      mtdtime_(0.), mtdtimeErr_(0.), mva_(0.), pathlength_(0.), btlchi2_(0.), btltimechi2_(0.), etlchi2_(0.), etltimechi2_(0.), time_pi_(0.), time_k_(0.), time_p_(0.), sigma_time_pi_(0.), sigma_time_k_(0.), sigma_time_p_(0.){
  initialFTS = trajectoryStateTransform::initialFreeState(tk, field);
}

GsfTransientTrack::GsfTransientTrack(const GsfTrack& tk,
                                     const double time,
                                     const double dtime,
                                     const MagneticField* field,
                                     const edm::ESHandle<GlobalTrackingGeometry>& tg)
    : GsfTrack(tk),
      tkr_(),
      hasTime(true),
      timeExt_(time),
      dtErrorExt_(dtime),
      theField(field),
      initialTSOSAvailable(false),
      initialTSCPAvailable(false),
      blStateAvailable(false),
      theTrackingGeometry(tg),
      theTIPExtrapolator(AnalyticalPropagator(field, alongMomentum)), 
      trkAssoc_(0),
      mtdtime_(0.), mtdtimeErr_(0.), mva_(0.), pathlength_(0.), btlchi2_(0.), btltimechi2_(0.), etlchi2_(0.), etltimechi2_(0.), time_pi_(0.), time_k_(0.), time_p_(0.), sigma_time_pi_(0.), sigma_time_k_(0.), sigma_time_p_(0.){
  initialFTS = trajectoryStateTransform::initialFreeState(tk, field);
}

GsfTransientTrack::GsfTransientTrack(const GsfTrackRef& tk,
                                     const MagneticField* field,
                                     const edm::ESHandle<GlobalTrackingGeometry>& tg)
    : GsfTrack(*tk),
      tkr_(tk),
      hasTime(false),
      timeExt_(0.),
      dtErrorExt_(0.),
      theField(field),
      initialTSOSAvailable(false),
      initialTSCPAvailable(false),
      blStateAvailable(false),
      theTrackingGeometry(tg),
      theTIPExtrapolator(AnalyticalPropagator(field, alongMomentum)), 
      trkAssoc_(0),
      mtdtime_(0.), mtdtimeErr_(0.), mva_(0.), pathlength_(0.), btlchi2_(0.), btltimechi2_(0.), etlchi2_(0.), etltimechi2_(0.), time_pi_(0.), time_k_(0.), time_p_(0.), sigma_time_pi_(0.), sigma_time_k_(0.), sigma_time_p_(0.){
  initialFTS = trajectoryStateTransform::initialFreeState(*tk, field);
}

GsfTransientTrack::GsfTransientTrack(const GsfTrackRef& tk,
                                     const double time,
                                     const double dtime,
                                     const MagneticField* field,
                                     const edm::ESHandle<GlobalTrackingGeometry>& tg)
    : GsfTrack(*tk),
      tkr_(tk),
      hasTime(true),
      timeExt_(time),
      dtErrorExt_(dtime),
      theField(field),
      initialTSOSAvailable(false),
      initialTSCPAvailable(false),
      blStateAvailable(false),
      theTrackingGeometry(tg),
      theTIPExtrapolator(AnalyticalPropagator(field, alongMomentum)), 
      trkAssoc_(0),
      mtdtime_(0.), mtdtimeErr_(0.), mva_(0.), pathlength_(0.), btlchi2_(0.), btltimechi2_(0.), etlchi2_(0.), etltimechi2_(0.), time_pi_(0.), time_k_(0.), time_p_(0.), sigma_time_pi_(0.), sigma_time_k_(0.), sigma_time_p_(0.){
  initialFTS = trajectoryStateTransform::initialFreeState(*tk, field);
}

GsfTransientTrack::GsfTransientTrack(const GsfTransientTrack& tt)
    : GsfTrack(tt),
      tkr_(tt.persistentTrackRef()),
      hasTime(tt.hasTime),
      timeExt_(tt.timeExt_),
      dtErrorExt_(tt.dtErrorExt_),
      theField(tt.field()),
      initialFTS(tt.initialFreeState()),
      initialTSOSAvailable(false),
      initialTSCPAvailable(false),
      theTIPExtrapolator(AnalyticalPropagator(tt.field(), alongMomentum)),
      trkAssoc_(0),
      mtdtime_(0.), mtdtimeErr_(0.), mva_(0.), pathlength_(0.), btlchi2_(0.), btltimechi2_(0.), etlchi2_(0.), etltimechi2_(0.), time_pi_(0.), time_k_(0.), time_p_(0.), sigma_time_pi_(0.), sigma_time_k_(0.), sigma_time_p_(0.){
       if (tt.initialTSOSAvailable) {
    initialTSOS = tt.impactPointState();
    initialTSOSAvailable = true;
  }
  if (tt.initialTSCPAvailable) {
    initialTSCP = tt.impactPointTSCP();
    initialTSCPAvailable = true;
  }
}

void GsfTransientTrack::setTrackingGeometry(const edm::ESHandle<GlobalTrackingGeometry>& tg) {
  theTrackingGeometry = tg;
}

void GsfTransientTrack::setBeamSpot(const BeamSpot& beamSpot) {
  theBeamSpot = beamSpot;
  blStateAvailable = false;
}

TrajectoryStateOnSurface GsfTransientTrack::impactPointState() const {
  if (!initialTSOSAvailable)
    calculateTSOSAtVertex();
  return initialTSOS;
}

TrajectoryStateClosestToPoint GsfTransientTrack::impactPointTSCP() const {
  if (!initialTSCPAvailable) {
    initialTSCP = builder(initialFTS, initialFTS.position());
    initialTSCPAvailable = true;
  }
  return initialTSCP;
}

TrajectoryStateOnSurface GsfTransientTrack::outermostMeasurementState() const {
  return MultiTrajectoryStateTransform::outerStateOnSurface(*this, *theTrackingGeometry, theField);
}

TrajectoryStateOnSurface GsfTransientTrack::innermostMeasurementState() const {
  return MultiTrajectoryStateTransform::innerStateOnSurface(*this, *theTrackingGeometry, theField);
}

void GsfTransientTrack::calculateTSOSAtVertex() const {
  TransverseImpactPointExtrapolator tipe(theField);
  initialTSOS = tipe.extrapolate(initialFTS, initialFTS.position());
  initialTSOSAvailable = true;
}

TrajectoryStateOnSurface GsfTransientTrack::stateOnSurface(const GlobalPoint& point) const {
  return theTIPExtrapolator.extrapolate(innermostMeasurementState(), point);
}

TrajectoryStateClosestToPoint GsfTransientTrack::trajectoryStateClosestToPoint(const GlobalPoint& point) const {
  return builder(stateOnSurface(point), point);
}

TrajectoryStateClosestToBeamLine GsfTransientTrack::stateAtBeamLine() const {
  if (!blStateAvailable) {
    TSCBLBuilderNoMaterial blsBuilder;
    trajectoryStateClosestToBeamLine = blsBuilder(initialFTS, theBeamSpot);
    blStateAvailable = true;
  }
  return trajectoryStateClosestToBeamLine;
}
