#ifndef ISO_H
#define ISO_H

#include <string>
#include <vector>

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "MagneticField/Engine/interface/MagneticField.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"

#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

class B0Isolation {
public:
  B0Isolation(reco::Vertex,
              edm::Handle<reco::TrackCollection>,
              edm::ESHandle<MagneticField>,
              ClosestApproachInRPhi,
              reco::BeamSpot,
              const pat::Muon&,
              const pat::Muon&,
              //                 reco::TrackRef , // muon -
              //                 reco::TrackRef , // muon +
              uint,  // trk minus index
              uint,  // trk plus index
              const reco::TransientTrack,
              const reco::TransientTrack,
              const reco::TransientTrack,
              const reco::TransientTrack);
  ~B0Isolation(){};

  std::vector<float> mum_isovec, mup_isovec, trkm_isovec, trkp_isovec;
  std::vector<float> mum_isopts, mup_isopts, trkm_isopts, trkp_isopts;
  std::vector<float> mum_isodr, mup_isodr, trkm_isodr, trkp_isodr;

private:
  const ParticleMass muonMass = 0.10565837;
  const ParticleMass pionMass = 0.13957018;
  const ParticleMass kaonMass = 0.493677;
  float mumasserr = 3.5e-9;

  float LHCbIsolation(
      ClosestApproachInRPhi, reco::TransientTrack, reco::TransientTrack, const ParticleMass, float, float, float);
};
#endif