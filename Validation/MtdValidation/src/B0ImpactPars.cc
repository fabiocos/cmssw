#include "Validation/MtdValidation/interface/B0ImpactPars.h"

B0ImpactPars::B0ImpactPars(edm::Handle<reco::VertexCollection> vertices,
                           const reco::TransientTrack refitMumTT,
                           const reco::TransientTrack refitMupTT,
                           const reco::TransientTrack refitTrkmTT,
                           const reco::TransientTrack refitTrkpTT) {
  double theMumMind0 = 100;
  double theMupMind0 = 100;
  double theTrkmMind0 = 100;
  double theTrkpMind0 = 100;
  double theMumMind0E = 100;
  double theMupMind0E = 100;
  double theTrkmMind0E = 100;
  double theTrkpMind0E = 100;
  double theMumMinIP = 100;
  double theMupMinIP = 100;
  double theTrkmMinIP = 100;
  double theTrkpMinIP = 100;
  double theMumMinIPS = 100;
  double theMupMinIPS = 100;
  double theTrkmMinIPS = 100;
  double theTrkpMinIPS = 100;
  GlobalPoint vert;

  std::pair<double, double> IPPair;
  TrajectoryStateClosestToPoint traj;

  for (std::vector<reco::Vertex>::const_iterator ipv = vertices->begin(); ipv != vertices->end(); ipv++) {
    if (!ipv->isValid())
      continue;
    vert = GlobalPoint(ipv->x(), ipv->y(), ipv->z());

    traj = refitMumTT.trajectoryStateClosestToPoint(vert);
    if (fabs(traj.perigeeParameters().transverseImpactParameter()) < theMumMind0) {
      theMumMind0 = fabs(traj.perigeeParameters().transverseImpactParameter());
      theMumMind0E = traj.perigeeError().transverseImpactParameterError();
    }

    traj = refitMupTT.trajectoryStateClosestToPoint(vert);
    if (fabs(traj.perigeeParameters().transverseImpactParameter()) < theMupMind0) {
      theMupMind0 = fabs(traj.perigeeParameters().transverseImpactParameter());
      theMupMind0E = traj.perigeeError().transverseImpactParameterError();
    }

    traj = refitTrkmTT.trajectoryStateClosestToPoint(vert);
    if (fabs(traj.perigeeParameters().transverseImpactParameter()) < theTrkmMind0) {
      theTrkmMind0 = fabs(traj.perigeeParameters().transverseImpactParameter());
      theTrkmMind0E = traj.perigeeError().transverseImpactParameterError();
    }

    traj = refitTrkpTT.trajectoryStateClosestToPoint(vert);
    if (fabs(traj.perigeeParameters().transverseImpactParameter()) < theTrkpMind0) {
      theTrkpMind0 = fabs(traj.perigeeParameters().transverseImpactParameter());
      theTrkpMind0E = traj.perigeeError().transverseImpactParameterError();
    }

    IPPair = impactParFromVtx(refitMumTT, *ipv);
    if (IPPair.first < theMumMinIP) {
      theMumMinIP = IPPair.first;
      theMumMinIPS = IPPair.second;
    }

    IPPair = impactParFromVtx(refitMupTT, *ipv);
    if (IPPair.first < theMupMinIP) {
      theMupMinIP = IPPair.first;
      theMupMinIPS = IPPair.second;
    }
    IPPair = impactParFromVtx(refitTrkmTT, *ipv);
    if (IPPair.first < theTrkmMinIP) {
      theTrkmMinIP = IPPair.first;
      theTrkmMinIPS = IPPair.second;
    }
    IPPair = impactParFromVtx(refitTrkpTT, *ipv);
    if (IPPair.first < theTrkpMinIP) {
      theTrkpMinIP = IPPair.first;
      theTrkpMinIPS = IPPair.second;
    }
  }

  mumMinIP.first = theMumMinIP;
  mumMinIP.second = theMumMinIPS;

  mupMinIP.first = theMupMinIP;
  mupMinIP.second = theMupMinIPS;

  tkmMinIP.first = theTrkmMinIP;
  tkmMinIP.second = theTrkmMinIPS;

  tkpMinIP.first = theTrkpMinIP;
  tkpMinIP.second = theTrkpMinIPS;

  mumMind0.first = theMumMind0;
  mumMind0.second = theMumMind0E;
  mupMind0.first = theMupMind0;
  mupMind0.second = theMupMind0E;
  tkmMind0.first = theTrkmMind0;
  tkmMind0.second = theTrkmMind0E;
  tkpMind0.first = theTrkpMind0;
  tkpMind0.second = theTrkpMind0E;
}

std::pair<double, double> B0ImpactPars::impactParFromVtx(reco::TransientTrack piTT, reco::Vertex myVtx) {
  std::pair<double, double> measure;
  std::pair<bool, Measurement1D> piIP_pair = IPTools::absoluteImpactParameter3D(piTT, myVtx);
  if (piIP_pair.first) {
    measure.first = piIP_pair.second.value();
    measure.second = piIP_pair.second.significance();
  } else {
    //       std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 3D" << std::endl;
    measure.first = 9999;
    measure.second = -9999.;
  }
  return measure;
}
