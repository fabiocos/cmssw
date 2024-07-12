#ifndef TrackReco_GsfTransientTrack_h
#define TrackReco_GsfTransientTrack_h

/**
   * Concrete implementation of the TransientTrack for a multi-state reco::GsfTrack
   */

#include "TrackingTools/TransientTrack/interface/BasicTransientTrack.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"

namespace reco {

  class GsfTransientTrack : public GsfTrack, public BasicTransientTrack {
  public:
    // constructor from persistent track
    GsfTransientTrack();
    GsfTransientTrack(const GsfTrack& tk, const MagneticField* field);
    GsfTransientTrack(const GsfTrack& tk, const double time, const double dtime, const MagneticField* field);

    GsfTransientTrack(const GsfTrackRef& tk, const MagneticField* field);
    GsfTransientTrack(const GsfTrackRef& tk, const double time, const double dtime, const MagneticField* field);

    GsfTransientTrack(const GsfTrackRef& tk,
                      const MagneticField* field,
                      const edm::ESHandle<GlobalTrackingGeometry>& trackingGeometry);
    GsfTransientTrack(const GsfTrackRef& tk,
                      const double time,
                      const double dtime,
                      const MagneticField* field,
                      const edm::ESHandle<GlobalTrackingGeometry>& trackingGeometry);

    GsfTransientTrack(const GsfTrack& tk,
                      const MagneticField* field,
                      const edm::ESHandle<GlobalTrackingGeometry>& trackingGeometry);
    GsfTransientTrack(const GsfTrack& tk,
                      const double time,
                      const double dtime,
                      const MagneticField* field,
                      const edm::ESHandle<GlobalTrackingGeometry>& trackingGeometry);

    GsfTransientTrack(const GsfTransientTrack& tt);

    GsfTransientTrack& operator=(const GsfTransientTrack& tt);

    void setTrackingGeometry(const edm::ESHandle<GlobalTrackingGeometry>&) override;

    void setBeamSpot(const reco::BeamSpot& beamSpot) override;

    FreeTrajectoryState initialFreeState() const override { return initialFTS; }

    TrajectoryStateOnSurface outermostMeasurementState() const override;

    TrajectoryStateOnSurface innermostMeasurementState() const override;

    /**
    * The TrajectoryStateClosestToPoint at any point.
    * The inner-most multi-state state will be used for the propagation.
    * The TSCP will be built out of the collapsed mult-component TSOS.
    */
    TrajectoryStateClosestToPoint trajectoryStateClosestToPoint(const GlobalPoint& point) const override;

    /**
    * The TSOS at any point. The inner-most multi-state state will be used for 
    * the propagation. The TSOS will therefore have multiple components.
    */
    TrajectoryStateOnSurface stateOnSurface(const GlobalPoint& point) const override;

    /**
    * The TrajectoryStateClosestToPoint at the initial state, made from the 
    * stored singl-component state.
    */
    TrajectoryStateClosestToPoint impactPointTSCP() const override;

    /**
    * The TSOS at the initial state, made from the stored singl-component state.
    */

    TrajectoryStateOnSurface impactPointState() const override;

    bool impactPointStateAvailable() const override { return initialTSOSAvailable; }

    /**
   * access to original persistent track
   */
    GsfTrackRef persistentTrackRef() const { return tkr_; }

    TrackBaseRef trackBaseRef() const override { return TrackBaseRef(tkr_); }

    TrackCharge charge() const override { return GsfTrack::charge(); }

    const MagneticField* field() const override { return theField; }

    const Track& track() const override { return *this; }

    TrajectoryStateClosestToBeamLine stateAtBeamLine() const override;

    double timeExt() const override { return (hasTime ? timeExt_ : std::numeric_limits<double>::quiet_NaN()); }

    double dtErrorExt() const override { return (hasTime ? dtErrorExt_ : std::numeric_limits<double>::quiet_NaN()); }
    int trackAsocMTD() const override { return (hasTime ? trkAssoc_ : std::numeric_limits<double>::quiet_NaN())  ; }
    float MTDtime() const override {return (hasTime ? mtdtime_ : std::numeric_limits<double>::quiet_NaN()); }
   float MTDtimeErr() const override {return (hasTime ? mtdtimeErr_ : std::numeric_limits<double>::quiet_NaN()); }  
    float MVAquality() const override { return (hasTime ? mva_  : std::numeric_limits<double>::quiet_NaN()); }
    float pathLength() const override { return (hasTime ? pathlength_ : std::numeric_limits<double>::quiet_NaN()); }
    float btlMatch_chi2() const override { return (hasTime ? btlchi2_ : std::numeric_limits<double>::quiet_NaN()); } 
    float btlMatchTime_chi2() const override { return (hasTime ? btltimechi2_ : std::numeric_limits<double>::quiet_NaN()); }
    float etlMatch_chi2() const override { return (hasTime ? etlchi2_ : std::numeric_limits<double>::quiet_NaN()); }
    float etlMatchTime_chi2() const override { return (hasTime ? etltimechi2_ : std::numeric_limits<double>::quiet_NaN()); }
    float trackTime_pi() const override { return (hasTime ? time_pi_ : std::numeric_limits<double>::quiet_NaN()); }
    float trackTime_k() const override { return (hasTime ? time_k_ : std::numeric_limits<double>::quiet_NaN()); }
    float trackTime_p() const override { return (hasTime ? time_p_ : std::numeric_limits<double>::quiet_NaN()); }
    float sigma_time_pi() const override {return (hasTime ? sigma_time_pi_ : std::numeric_limits<double>::quiet_NaN()); }
    float sigma_time_k() const override {return (hasTime ? sigma_time_k_: std::numeric_limits<double>::quiet_NaN()); }
    float sigma_time_p() const override {return (hasTime ? sigma_time_p_ : std::numeric_limits<double>::quiet_NaN()); }
    int nPixBarrel() const override { return (hasTime ? npixbarrel_ : std::numeric_limits<double>::quiet_NaN()); }
    int nPixEndcap() const override { return (hasTime ? npixendcap_ : std::numeric_limits<double>::quiet_NaN()); }

  private:
    void init();
    void calculateTSOSAtVertex() const;

    GsfTrackRef tkr_;
    bool hasTime;
    double timeExt_, dtErrorExt_;
    const MagneticField* theField;

    FreeTrajectoryState initialFTS;
    mutable bool initialTSOSAvailable, initialTSCPAvailable, blStateAvailable;
    mutable TrajectoryStateOnSurface initialTSOS;
    mutable TrajectoryStateClosestToPoint initialTSCP;
    TSCPBuilderNoMaterial builder;
    edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
    reco::BeamSpot theBeamSpot;
      TransverseImpactPointExtrapolator theTIPExtrapolator;
    int trkAssoc_;
    float mtdtime_, mtdtimeErr_;
    float mva_;
    float pathlength_;
    float btlchi2_, btltimechi2_, etlchi2_, etltimechi2_,time_pi_,time_k_,time_p_,sigma_time_pi_, sigma_time_k_, sigma_time_p_;
    int npixbarrel_, npixendcap_ ;
    mutable TrajectoryStateClosestToBeamLine trajectoryStateClosestToBeamLine;

  };

}  // namespace reco

#endif
