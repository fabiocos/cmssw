#ifndef TrackReco_TransientTrackFromFTS_h
#define TrackReco_TransientTrackFromFTS_h

/**
   * Concrete implementation of the TransientTrack for a multi-state reco::GsfTrack
   * To be built through the factory TransientTrackFromFTSFactory or the TransientTrackBuilder
   */

#include "TrackingTools/TransientTrack/interface/BasicTransientTrack.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h"

namespace reco {

  class TransientTrackFromFTS : public BasicTransientTrack {
  public:
    TransientTrackFromFTS();

    TransientTrackFromFTS(const FreeTrajectoryState& fts);
    TransientTrackFromFTS(const FreeTrajectoryState& fts, const double time, const double dtime);

    TransientTrackFromFTS(const FreeTrajectoryState& fts,
                          const edm::ESHandle<GlobalTrackingGeometry>& trackingGeometry);
    TransientTrackFromFTS(const FreeTrajectoryState& fts,
                          const double time,
                          const double dtime,
                          const edm::ESHandle<GlobalTrackingGeometry>& trackingGeometry);

    TransientTrackFromFTS(const TransientTrackFromFTS& tt);

    TransientTrackFromFTS& operator=(const TransientTrackFromFTS& tt);

    void setTrackingGeometry(const edm::ESHandle<GlobalTrackingGeometry>&) override;

    void setBeamSpot(const reco::BeamSpot& beamSpot) override;

    FreeTrajectoryState initialFreeState() const override { return initialFTS; }

    TrajectoryStateOnSurface outermostMeasurementState() const override;

    TrajectoryStateOnSurface innermostMeasurementState() const override;

    TrajectoryStateClosestToPoint trajectoryStateClosestToPoint(const GlobalPoint& point) const override {
      return builder(initialFTS, point);
    }

    /**
    * The TSOS at any point. The initial state will be used for the propagation.
    */
    TrajectoryStateOnSurface stateOnSurface(const GlobalPoint& point) const override;

    TrajectoryStateClosestToPoint impactPointTSCP() const override;

    TrajectoryStateOnSurface impactPointState() const override;

    bool impactPointStateAvailable() const override { return initialTSOSAvailable; }

    TrackCharge charge() const override { return initialFTS.charge(); }

    const MagneticField* field() const override { return theField; }

    const Track& track() const override;

    TrackBaseRef trackBaseRef() const override { return TrackBaseRef(); }

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
    void calculateTSOSAtVertex() const;

    FreeTrajectoryState initialFTS;
    bool hasTime;
    double timeExt_;
    double dtErrorExt_;
    const MagneticField* theField;
    mutable bool initialTSOSAvailable, initialTSCPAvailable, trackAvailable, blStateAvailable;
    mutable TrajectoryStateOnSurface initialTSOS;
    mutable TrajectoryStateClosestToPoint initialTSCP;
    mutable Track theTrack;
    TSCPBuilderNoMaterial builder;
    edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
    reco::BeamSpot theBeamSpot;
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
