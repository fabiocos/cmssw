#include "Validation/MtdValidation/interface/B0Isolation.h"
#include "TLorentzVector.h"

#define TRKMAXR 110.0 // [cm]
#define TRKMAXZ 280.0 // [cm]


B0Isolation::B0Isolation(reco::Vertex bestVtx, 
                         edm::Handle<reco::TrackCollection> tracks, 
                         edm::ESHandle<MagneticField> bFieldHandle,
                         ClosestApproachInRPhi ClosestApp,
                         reco::BeamSpot beamSpot,
                         const pat::Muon& mum,
                         const pat::Muon& mup,
//                          reco::TrackRef muTrackm, // muon -
//                          reco::TrackRef muTrackp, // muon +
                         uint itrkm, // trk minus index
                         uint itrkp,  // trk plus index
                         const reco::TransientTrack refitMumTT, 
                         const reco::TransientTrack refitMupTT,
                         const reco::TransientTrack refitTrkmTT,
                         const reco::TransientTrack refitTrkpTT
                         ){
                         
    float iso = 0;
    mum_isovec .clear();
    mup_isovec .clear();
    trkm_isovec.clear();
    trkp_isovec.clear();

//     for (std::vector<reco::Vertex>::const_iterator iVertex = vertices->begin(); iVertex != vertices->end(); iVertex++) { 
//       bestVtx = *(iVertex); if (bestVtx.isValid() == true) break; 
//     }

    for (uint itrkiso =0 ;  itrkiso < tracks->size(); itrkiso++)
    {

      reco::TrackRef tkiso(tracks,itrkiso) ;                                                
      if ( itrkiso == itrkm  || itrkiso == itrkp)           continue;
//       if (!tkiso->quality(reco::TrackBase::highPurity))     continue;
      if ( tkiso->pt() < 0.8 )                              continue;


      // check if the track is one of the two muons
      bool skip_this_track = false;              
      for (unsigned int i = 0; i < mum.numberOfSourceCandidatePtrs(); ++i) {
          const edm::Ptr<reco::Candidate> & source = mum.sourceCandidatePtr(i);
          if (! ( (mum.sourceCandidatePtr(i)).isNonnull() &&  (mum.sourceCandidatePtr(i)).isAvailable() ))   continue;
          const reco::Candidate & cand = *(source);
          if (cand.charge() == 0 || cand.bestTrack() == nullptr)      continue;
          try{ cand.bestTrack()->eta();}
          catch(...) { continue;}
          if ( deltaR(tkiso->eta(),tkiso->phi(),cand.bestTrack()->eta(),cand.bestTrack()->phi()) < 0.00001 ) {
              skip_this_track = true;
              break;
          }
      }
      if (skip_this_track) continue;
      for (unsigned int i = 0; i < mup.numberOfSourceCandidatePtrs(); ++i) {
          const edm::Ptr<reco::Candidate> & source = mup.sourceCandidatePtr(i);
          if (! ( (mup.sourceCandidatePtr(i)).isNonnull() &&  (mup.sourceCandidatePtr(i)).isAvailable() ))   continue;
          const reco::Candidate & cand = *(source);
          if (cand.charge() == 0 || cand.bestTrack() == nullptr)      continue;
          try{ cand.bestTrack()->eta();}
          catch(...) { continue;}
          if ( deltaR(tkiso->eta(),tkiso->phi(),cand.bestTrack()->eta(),cand.bestTrack()->phi()) < 0.00001 ) {
              skip_this_track = true;
              break;
          }
      }
      if (skip_this_track) continue;


      // requirement that the track is not associated to any PV
      // if (findPV(trk_index, recVtxColl) == 1 ) continue;
      
      const reco::TransientTrack TrackIsoTT((*tkiso), &(*bFieldHandle));
      // check to ensure the goodness of the track
      if (! (TrackIsoTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),
                                                                  beamSpot.position().y(),
                                                                  beamSpot.position().z())).perigeeError().transverseImpactParameterError() >0) ) continue;
 
      
      iso = LHCbIsolation(ClosestApp, refitMumTT, TrackIsoTT, muonMass, bestVtx.x(), bestVtx.y(), bestVtx.z());
      if (iso > 0 ) mum_isovec.push_back(iso);
      
      iso = LHCbIsolation(ClosestApp, refitMupTT, TrackIsoTT, muonMass, bestVtx.x(), bestVtx.y(), bestVtx.z());
      if (iso > 0 ) mup_isovec.push_back(iso);
    
      iso = LHCbIsolation(ClosestApp, refitTrkmTT, TrackIsoTT, muonMass, bestVtx.x(), bestVtx.y(), bestVtx.z());
      if (iso > 0 ) trkm_isovec.push_back(iso);
    
      iso = LHCbIsolation(ClosestApp, refitTrkpTT, TrackIsoTT, muonMass, bestVtx.x(), bestVtx.y(), bestVtx.z());
      if (iso > 0 ) trkp_isovec.push_back(iso);
    
    
      // add new iso
      ClosestApp.calculate(TrackIsoTT.initialFreeState(), refitMumTT.initialFreeState());
      if (ClosestApp.status() != false)
      {
        if ( ClosestApp.distance() < 0.1 ) {
          mum_isopts.push_back( tkiso->pt() );
//        mum_isomom.push_back( tkiso->p()  );   
          mum_isodr.push_back ( deltaR(tkiso ->eta(), tkiso->phi(), refitMumTT.track().eta(), refitMumTT.track().phi() ));   
        }                 
      }
      ClosestApp.calculate(TrackIsoTT.initialFreeState(), refitMupTT.initialFreeState());
      if (ClosestApp.status() != false)
      {
        if ( ClosestApp.distance() < 0.1 ) {
          mup_isopts.push_back( tkiso->pt() );
          mup_isodr.push_back ( deltaR(tkiso ->eta(), tkiso->phi(), refitMupTT.track().eta(), refitMupTT.track().phi() ));   
         }                 
      }
      ClosestApp.calculate(TrackIsoTT.initialFreeState(), refitTrkmTT.initialFreeState());
      if (ClosestApp.status() != false)
      {
        if ( ClosestApp.distance() < 0.1 ) {
          trkm_isopts.push_back( tkiso->pt() );
          trkm_isodr.push_back ( deltaR(tkiso ->eta(), tkiso->phi(), refitTrkmTT.track().eta(), refitTrkmTT.track().phi() ));   
        }                 
      }
      ClosestApp.calculate(TrackIsoTT.initialFreeState(), refitTrkpTT.initialFreeState());
      if (ClosestApp.status() != false)
      {
        if ( ClosestApp.distance() < 0.1 ) {
          trkp_isopts.push_back( tkiso->pt() );
          trkp_isodr.push_back ( deltaR(tkiso ->eta(), tkiso->phi(), refitTrkpTT.track().eta(), refitTrkpTT.track().phi() ));   
        }                 
      }
    
    }                         

}



float B0Isolation::LHCbIsolation( ClosestApproachInRPhi ClosestApp, 
                                  reco::TransientTrack muTrackmTT, 
                                  reco::TransientTrack TrackIsoTT,
                                  const ParticleMass muonMass,
                                  float vtx_x,
                                  float vtx_y,
                                  float vtx_z
                                )
{
  float iso  = 0;

  KinematicParticleVertexFitter KPVtxFitter; 
  KinematicParticleFactoryFromTransientTrack partFactory;

  float chi, ndf, alpha;
  TLorentzVector tmp_trk_lv, tmp_cand_lv, tmp_lv;
  ClosestApp.calculate(TrackIsoTT.initialFreeState(), muTrackmTT.initialFreeState());
  if (ClosestApp.status() == false)
  {
//      std::cout << __LINE__ << " : continue --> bad status of closest approach" << std::endl;
     return 0;
  }

  GlobalPoint XingPoint = ClosestApp.crossingPoint();
  if ((sqrt(XingPoint.x()*XingPoint.x() + XingPoint.y()*XingPoint.y()) > TRKMAXR) || (fabs(XingPoint.z()) > TRKMAXZ))
  {
//       std::cout << __LINE__ << " : continue --> closest approach crossing point outside the tracker volume" << std::endl;
      return 0;
  }

  if ( ClosestApp.distance() < 0.1 )
  {
    chi = 0.;
    ndf = 0.;
    std::vector<RefCountedKinematicParticle> tmpParticles;
    tmpParticles.push_back(partFactory.particle(muTrackmTT, muonMass,chi,ndf, mumasserr) );
    tmpParticles.push_back(partFactory.particle(TrackIsoTT, muonMass,chi,ndf, mumasserr) );


    RefCountedKinematicTree tmpVertexFitTree = KPVtxFitter.fit(tmpParticles);
    if ( ! tmpVertexFitTree->isValid()) return 0;
      tmpVertexFitTree->movePointerToTheTop();
      RefCountedKinematicVertex tmpVertex   = tmpVertexFitTree->currentDecayVertex();
      if (tmpVertexFitTree->isValid() && tmpVertex -> vertexIsValid() )
      {
      tmp_trk_lv .SetPtEtaPhiM( TrackIsoTT.track().pt(), TrackIsoTT.track().eta(), TrackIsoTT.track().phi(), kaonMass);
      tmp_cand_lv.SetPtEtaPhiM( muTrackmTT.track().pt(), muTrackmTT.track().eta(), muTrackmTT.track().phi(), kaonMass);
      tmp_lv =  tmp_trk_lv + tmp_cand_lv;
      
      alpha = tmp_lv.Angle(TVector3 ( tmpVertex->position().x() - vtx_x, 
                                      tmpVertex->position().y() - vtx_y, 
                                      tmpVertex->position().z() - vtx_z )
                                    );
      iso = (tmp_lv).P() * alpha / ( (tmp_lv).P()*alpha + tmp_cand_lv.Pt() + tmp_trk_lv.Pt() );
    }
  }

  return iso;

}


