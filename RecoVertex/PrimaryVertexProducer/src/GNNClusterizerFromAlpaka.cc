/**
 * GNNClusterizerFromAlpaka.cc - Builds TransientVertex from Alpaka SoA outputs
 *
 * This does NOT run inference - it consumes pre-computed outputs from 
 * GNNVertexProducerAlpaka and builds TransientVertex objects.
 */

#include "RecoVertex/PrimaryVertexProducer/interface/GNNClusterizerFromAlpaka.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include <algorithm>
#include <vector>

namespace vertexgnn {

  GNNClusterizerFromAlpaka::GNNClusterizerFromAlpaka(const edm::ParameterSet& conf)
      : existenceThreshold_(conf.getParameter<double>("existenceThreshold")),
        trackAssignmentThreshold_(conf.getParameter<double>("trackAssignmentThreshold")),
        verbose_(conf.getUntrackedParameter<bool>("verbose", false)) {}

  std::vector<TransientVertex> GNNClusterizerFromAlpaka::vertices(
      const std::vector<reco::TransientTrack>& tracks,
      const float* A_flat,
      const float* z_hat,
      const float* t_hat,
      const float* p,
      const float* pi_0,
      const float* pi_1,
      const float* pi_2,
      int N,
      int K) const {
    
    std::vector<TransientVertex> clusters;
    
    if (N == 0 || K == 0) {
      return clusters;
    }

    // =========================================================================
    // PROCESS ASSIGNMENTS (same logic as GNNClusterizer.cc)
    // =========================================================================
    
    // Step 1: For each track, find slot with maximum assignment probability
    std::vector<int> track_to_slot(N, -1);
    std::vector<float> track_max_prob(N, 0.0f);
    
    for (int i = 0; i < N; ++i) {
      float max_prob = -1.0f;
      int best_slot = -1;
      for (int k = 0; k < K; ++k) {
        float prob = A_flat[i * K + k];
        if (prob > max_prob) {
          max_prob = prob;
          best_slot = k;
        }
      }
      track_to_slot[i] = best_slot;
      track_max_prob[i] = max_prob;
    }

    // Step 2: Count tracks per slot
    std::vector<int> slot_track_count(K, 0);
    for (int i = 0; i < N; ++i) {
      int slot = track_to_slot[i];
      if (slot >= 0 && slot < K) {
        slot_track_count[slot]++;
      }
    }

    // Step 3: Identify active slots (existence prob > threshold AND has tracks)
    std::vector<bool> slot_active(K, false);
    int num_active = 0;
    for (int k = 0; k < K; ++k) {
      if (p[k] > existenceThreshold_ && slot_track_count[k] >= 1) {
        slot_active[k] = true;
        num_active++;
      }
    }

    if (verbose_) {
      edm::LogInfo("GNNClusterizerFromAlpaka") << "Active slots: " << num_active 
                                                << " (threshold=" << existenceThreshold_ << ")";
    }

    // Step 4: Build TransientVertex for each active slot
    GlobalError dummyError(0.01, 0, 0.01, 0., 0., 0.01);
    
    for (int k = 0; k < K; ++k) {
      if (!slot_active[k]) continue;

      // Collect tracks assigned to this slot
      std::vector<reco::TransientTrack> cluster_tracks;
      std::vector<float> cluster_weights;
      
      for (int i = 0; i < N; ++i) {
        // Only include track if assigned to this slot AND probability exceeds threshold
        if (track_to_slot[i] == k && track_max_prob[i] >= trackAssignmentThreshold_) {
          cluster_tracks.push_back(tracks[i]);
          cluster_weights.push_back(track_max_prob[i]);
        }
      }

      if (cluster_tracks.empty()) continue;

      // Use GNN's predicted z position
      double z_cluster = static_cast<double>(z_hat[k]);
      GlobalPoint pos(0, 0, z_cluster);
      TransientVertex vertex(pos, dummyError, cluster_tracks, 0);
      
      // Populate weight map with GNN assignment probabilities
      TransientVertex::TransientTrackToFloatMap weightMap;
      for (size_t i = 0; i < cluster_tracks.size(); ++i) {
        weightMap[cluster_tracks[i]] = cluster_weights[i];
      }
      vertex.weightMap(weightMap);
      
      clusters.push_back(vertex);

      if (verbose_) {
        edm::LogInfo("GNNClusterizerFromAlpaka") 
            << "  Vertex k=" << k 
            << ": z=" << z_cluster
            << ", t=" << t_hat[k]
            << ", p=" << p[k]
            << ", ntracks=" << cluster_tracks.size();
      }
    }

    if (verbose_) {
      edm::LogInfo("GNNClusterizerFromAlpaka") << "Final vertex count: " << clusters.size();
    }

    return clusters;
  }

  void GNNClusterizerFromAlpaka::fillPSetDescription(edm::ParameterSetDescription& desc) {
    desc.add<double>("existenceThreshold", 0.5)
        ->setComment("Slot existence probability threshold (0-1)");
    desc.add<double>("trackAssignmentThreshold", 0.5)
        ->setComment("Minimum track assignment probability");
    desc.addUntracked<bool>("verbose", false);
  }

}  // namespace vertexgnn
