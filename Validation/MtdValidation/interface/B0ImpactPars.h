#ifndef IPS_H
#define IPS_H

// #include <string>
#include <vector>

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
// #include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

class B0ImpactPars {
public:
  B0ImpactPars(edm::Handle<reco::VertexCollection>,
               const reco::TransientTrack,
               const reco::TransientTrack,
               const reco::TransientTrack,
               const reco::TransientTrack);
  ~B0ImpactPars(){};

  std::pair<double, double> mumMind0, mumMinIP;
  std::pair<double, double> mupMind0, mupMinIP;
  std::pair<double, double> tkmMind0, tkmMinIP;
  std::pair<double, double> tkpMind0, tkpMinIP;

private:
  std::pair<double, double> impactParFromVtx(reco::TransientTrack, reco::Vertex);
};
#endif