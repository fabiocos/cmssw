#ifndef SequentialPrimaryVertexFitterAdapter_h
#define SequentialPrimaryVertexFitterAdapter_h

/**\class SequentialPrimaryVertexFitterAdapter
 
  Description: Adapter class for Kalman and Adaptive vertex fitters 

*/

#include <sstream>
#include <map>

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexFitterBase.h"
#include "RecoVertex/VertexPrimitives/interface/VertexFitter.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

class SequentialPrimaryVertexFitterAdapter : public PrimaryVertexFitterBase {
public:
  SequentialPrimaryVertexFitterAdapter() : fitter(nullptr), useClusterWeights_(true) {}
  SequentialPrimaryVertexFitterAdapter(const VertexFitter<5>* vertex_fitter) : fitter(vertex_fitter), useClusterWeights_(true) {}
  SequentialPrimaryVertexFitterAdapter(const VertexFitter<5>* vertex_fitter, bool useClusterWeights) 
      : fitter(vertex_fitter), useClusterWeights_(useClusterWeights) {}
  ~SequentialPrimaryVertexFitterAdapter() override = default;
  
  void setUseClusterWeights(bool use) { useClusterWeights_ = use; }

  std::vector<TransientVertex> fit(const std::vector<reco::TransientTrack>& dummy,
                                   const std::vector<TransientVertex>& clusters,
                                   const reco::BeamSpot& beamspot,
                                   const bool useBeamConstraint) override {
    std::vector<TransientVertex> pvs;
    int clusterIdx = 0;
    for (auto& cluster : clusters) {
      const std::vector<reco::TransientTrack>& tracklist = cluster.originalTracks();
      TransientVertex v;
      if (useBeamConstraint && (tracklist.size() > 1)) {
        try {
          v = fitter->vertex(tracklist, beamspot);
        } catch (VertexException& ex) {
          std::ostringstream beamspotInfo;
          beamspotInfo << "While processing SequentialPrimaryVertexFitterAdapter::fit() with BeamSpot parameters: \n"
                       << beamspot;
          ex.addContext(beamspotInfo.str());
          throw;  // rethrow the exception
        }
      } else if (!(useBeamConstraint) && (tracklist.size() > 1)) {
        v = fitter->vertex(tracklist);
      }  // else: no fit ==> v.isValid()=False

      if (v.isValid()) {
        // DEBUG: Log weight map status before copy
        edm::LogInfo("SequentialPrimaryVertexFitterAdapter") 
            << "Cluster " << clusterIdx << ": hasTrackWeight=" << cluster.hasTrackWeight()
            << ", cluster.ntracks=" << cluster.originalTracks().size()
            << ", fitter.ntracks=" << v.originalTracks().size();
        
        // Optionally preserve weight map from input cluster (e.g., GNN assignment probabilities)
        // instead of using fitter's geometric weights  
        // useClusterWeights_ = true: use GNN weights, false: use fitter geometric weights
        if (useClusterWeights_ && cluster.hasTrackWeight()) {
          const auto& clusterWM = cluster.weightMap();
          
          // DEBUG: Print ALL tracks for first cluster only (clusterIdx == 0), otherwise first 5
          bool printAll = (clusterIdx == 0);
          int maxPrint = printAll ? 9999 : 5;
          
          // DEBUG: Log cluster track pointers
          edm::LogInfo("SequentialPrimaryVertexFitterAdapter") 
              << "  [DEBUG] Cluster tracks (weight map keys): " 
              << (printAll ? "*** PRINTING ALL ***" : "first 5 only");
          int debugIdx = 0;
          for (const auto& kv : clusterWM) {
            if (debugIdx < maxPrint) {
              const reco::Track* tptr = &(kv.first.track());
              edm::LogInfo("SequentialPrimaryVertexFitterAdapter") 
                  << "    C[" << debugIdx << "]: ptr=" << tptr 
                  << ", pt=" << kv.first.track().pt()
                  << ", w=" << kv.second;
            }
            debugIdx++;
          }
          
          // DEBUG: Log fitter's original tracks
          edm::LogInfo("SequentialPrimaryVertexFitterAdapter") 
              << "  [DEBUG] Fitter tracks (v.originalTracks()): "
              << (printAll ? "*** PRINTING ALL ***" : "first 5 only");
          debugIdx = 0;
          for (const auto& tt : v.originalTracks()) {
            if (debugIdx < maxPrint) {
              const reco::Track* tptr = &(tt.track());
              edm::LogInfo("SequentialPrimaryVertexFitterAdapter") 
                  << "    F[" << debugIdx << "]: ptr=" << tptr 
                  << ", pt=" << tt.track().pt();
            }
            debugIdx++;
          }
          
          // Build a lookup from underlying track pointer to weight
          // (TransientTrack objects differ, but underlying reco::Track should be same)
          std::map<const reco::Track*, float> trackPtrToWeight;
          for (const auto& kv : clusterWM) {
            trackPtrToWeight[&(kv.first.track())] = kv.second;
          }
          
          // Rebuild weight map using fitter's TransientTracks
          TransientVertex::TransientTrackToFloatMap newWeightMap;
          int foundCount = 0;
          int notFoundCount = 0;
          edm::LogInfo("SequentialPrimaryVertexFitterAdapter") 
              << "  [DEBUG] Matching fitter tracks to cluster weights: "
              << (printAll ? "*** PRINTING ALL ***" : "first 5 only");
          for (const auto& tt : v.originalTracks()) {
            const reco::Track* tptr = &(tt.track());
            auto it = trackPtrToWeight.find(tptr);
            if (it != trackPtrToWeight.end()) {
              newWeightMap[tt] = it->second;
              if (foundCount < maxPrint) {
                edm::LogInfo("SequentialPrimaryVertexFitterAdapter") 
                    << "    M[" << foundCount << "]: ptr=" << tptr << ", pt=" << tt.track().pt() 
                    << ", w=" << it->second;
              }
              foundCount++;
            } else {
              // Track not in cluster weight map - use default 1.0
              newWeightMap[tt] = 1.0f;
              if (notFoundCount < maxPrint) {
                edm::LogInfo("SequentialPrimaryVertexFitterAdapter") 
                    << "    NOT FOUND[" << notFoundCount << "]: ptr=" << tptr << ", pt=" << tt.track().pt() 
                    << ", using w=1.0";
              }
              notFoundCount++;
            }
          }
          
          v.weightMap(newWeightMap);
          
          edm::LogInfo("SequentialPrimaryVertexFitterAdapter") 
              << "  Weight map rebuilt: " << foundCount << " matched, " 
              << notFoundCount << " not found (used 1.0)";
        }
        pvs.push_back(v);
      }
      clusterIdx++;
    }
    return pvs;
  };

protected:
  // configuration
  const VertexFitter<5>* fitter;  // Kalman or Adaptive
  bool useClusterWeights_;        // true: use cluster weights (GNN), false: use fitter geometric weights
};
#endif
