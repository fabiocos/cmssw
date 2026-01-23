
#ifndef GNNClusterizer_h
#define GNNClusterizer_h

/**
 * GNNClusterizer - End-to-end Vertex Slot Model based track clustering
 * 
 * This implementation uses a VertexSlotModel ONNX export that directly outputs
 * track-to-vertex assignments and vertex positions, eliminating the need for
 * external clustering algorithms.
 */

#include "RecoVertex/PrimaryVertexProducer/interface/TrackClusterizerInZ.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/stream/EDAnalyzer.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"

#include <memory>
#include <vector>

using namespace cms::Ort;

class GNNClusterizer final : public TrackClusterizerInZ {
public:
  GNNClusterizer(const edm::ParameterSet& conf, const ONNXRuntime* onnxRuntime);
  ~GNNClusterizer() override {}

  std::vector<TransientVertex> vertices(const std::vector<reco::TransientTrack>& tracks) const override;
  std::vector<std::vector<reco::TransientTrack>> clusterize(
      const std::vector<reco::TransientTrack>& tracks) const override;
  static void fillPSetDescription(edm::ParameterSetDescription& desc);
  static std::unique_ptr<ONNXRuntime> initializeGlobalCache(const edm::ParameterSet& conf);
  static void globalEndJob(const ONNXRuntime* cache);

  // Cached outputs from last vertices() call - VertexSlotModel outputs
  const std::vector<float>& lastAssignments() const { return last_assignments_; }  // N*K
  const std::vector<float>& lastZHat() const { return last_z_hat_; }               // K
  const std::vector<float>& lastTHat() const { return last_t_hat_; }               // K
  const std::vector<float>& lastP() const { return last_p_; }                      // K (existence prob)
  const std::vector<float>& lastPi() const { return last_pi_; }                    // N*3 (PID weights)
  int lastNumSlots() const { return last_num_slots_; }
  int lastTrackCount() const { return last_n_tracks_; }
  bool hasLastOutputs() const { return last_valid_; }

private:
  const ONNXRuntime* onnxRuntime_;
  std::string nnVersion_;
  bool verbose_;
  double d0CutOff_;
  double vertexSize_;
  
  // VertexSlotModel-specific parameters
  double existenceThreshold_;
  double trackAssignmentThreshold_;  // Minimum track assignment probability
  int numSlots_;

  // Cached outputs from last inference
  mutable std::vector<float> last_assignments_;  // N*K
  mutable std::vector<float> last_z_hat_;        // K
  mutable std::vector<float> last_t_hat_;        // K
  mutable std::vector<float> last_p_;            // K
  mutable std::vector<float> last_pi_;           // N*3
  mutable int last_num_slots_ = 0;
  mutable int last_n_tracks_ = 0;
  mutable bool last_valid_ = false;
};

#endif
