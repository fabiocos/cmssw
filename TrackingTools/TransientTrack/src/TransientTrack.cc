#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TrackTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/CandidatePtrTransientTrack.h"

#include <iostream>

using namespace reco;

typedef TrackTransientTrack TTT;
typedef CandidatePtrTransientTrack CTT;

TransientTrack::TransientTrack(const Track& tk, const MagneticField* field) : Base(new TTT(tk, field)) {}
TransientTrack::TransientTrack(const Track& tk, const double time, const double dtime, const MagneticField* field)
    : Base(new TTT(tk, time, dtime, field)) {}

TransientTrack::TransientTrack(const CandidatePtr& ptr, const MagneticField* field) : Base(new CTT(ptr, field)) {}

TransientTrack::TransientTrack(const TrackRef& tk, const MagneticField* field) : Base(new TTT(tk, field)) {}

TransientTrack::TransientTrack(const TrackRef& tk, const double time, const double dtime, const MagneticField* field)
    : Base(new TTT(tk, time, dtime, field)) {}

TransientTrack::TransientTrack(const Track& tk,
                               const MagneticField* field,
                               const edm::ESHandle<GlobalTrackingGeometry>& tg)
    : Base(new TTT(tk, field, tg)) {}
TransientTrack::TransientTrack(const Track& tk,
                               const double time,
                               const double dtime,
                               const MagneticField* field,
                               const edm::ESHandle<GlobalTrackingGeometry>& tg)
    : Base(new TTT(tk, time, dtime, field, tg)) {}

TransientTrack::TransientTrack(const TrackRef& tk,
                               const MagneticField* field,
                               const edm::ESHandle<GlobalTrackingGeometry>& tg)
    : Base(new TTT(tk, field, tg)) {}
TransientTrack::TransientTrack(const TrackRef& tk,
                               const double time,
                               const double dtime,
                               const MagneticField* field,
                               const edm::ESHandle<GlobalTrackingGeometry>& tg)
    : Base(new TTT(tk, time, dtime, field, tg)) {}

TransientTrack::TransientTrack(const CandidatePtr& tk,
                               const MagneticField* field,
                               const edm::ESHandle<GlobalTrackingGeometry>& tg)
    : Base(new CTT(tk, field, tg)) {}

TransientTrack::TransientTrack(const CandidatePtr& tk,
                               const double time,
                               const double dtime,
                               const MagneticField* field,
                               const edm::ESHandle<GlobalTrackingGeometry>& tg)
    : Base(new CTT(tk, time, dtime, field, tg)) {}

TransientTrack::TransientTrack(const TrackRef& tk,
                               const double time,
                               const double dtime,
                               const MagneticField* field,
                               const edm::ESHandle<GlobalTrackingGeometry>& tg,
                               const int trkAssoc,
                               const float mtdtime, const float mtdtimeErr,
                               const float mva, const float pathlength, const float btlchi2, const float btltimechi2, const float etlchi2, const float etltimechi2, const float time_pi, const float time_k, const float time_p, const float sigma_time_pi, const float sigma_time_k, const float sigma_time_p, const int npixBar, const int npixEnd)
    : Base(new TTT(tk, time, dtime, field, tg, trkAssoc, mtdtime, mtdtimeErr, mva, pathlength, btlchi2, btltimechi2, etlchi2, etltimechi2, time_pi, time_k, time_p, sigma_time_pi, sigma_time_k, sigma_time_p, npixBar, npixEnd)) {}

// TransientTrack::TransientTrack( const TransientTrack & tt ) :
//   Base( new TTT(tt)) {}

// TransientTrack& TransientTrack::operator=(const TransientTrack & tt) {
// //   std::cout << "assign op." << std::endl;
//   if (this == &tt) return *this;
//   //
//   //  std::cout << tt.tk_ << std::endl;
// //   std::cout << "assign base." << std::endl;
//   Track::operator=(tt);
// //   std::cout << "done assign base." << std::endl;
//   //  tk_ = &(tt.persistentTrack());
//   //  tk_ = tt.tk_;
// //   std::cout << "assign ref." << std::endl;
//   tkr_ = tt.persistentTrackRef();
//   initialTSOSAvailable =  tt.initialTSOSAvailable;
//   initialTSCPAvailable = tt.initialTSCPAvailable;
//   initialTSCP = tt.initialTSCP;
//   initialTSOS = tt.initialTSOS;
//   theField = tt.field();
//   initialFTS = tt.initialFreeState();
// //   std::cout << "assign op. OK" << std::endl;
//
//   return *this;
// }
