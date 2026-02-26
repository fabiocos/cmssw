#ifndef RecoVertex_PrimaryVertexProducer_interface_VertexGNNSoA_h
#define RecoVertex_PrimaryVertexProducer_interface_VertexGNNSoA_h

/**
 * SoA layouts for GNN vertex producer using PyTorchAlpaka.
 * 
 * Designed for proper TensorCollection integration with batch=N (tracks):
 * - Input: TrackFeaturesSoA [N tracks × 13 features]
 * - Output: GNNOutputSoA [N tracks × K slots per output]
 * 
 * All outputs have batch=N to enable:
 *   model_.forward(event.queue(), inputs, outputs)
 * 
 * Slot predictions (z, t, p) are replicated per track for TensorCollection
 * compatibility, but values are identical across all tracks.
 * PID weights (pi_0, pi_1, pi_2) are unique per track.
 */

#include <Eigen/Core>
#include <Eigen/Dense>

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace vertexgnn {

  // Fixed number of vertex slots (compile-time constant for Eigen types)
  // NOTE: Must match the TorchScript model's num_slots parameter!
  // v17p1 model uses num_slots=180
  constexpr int kNumSlots = 180;

  // Type aliases for Eigen column types
  using SlotVector = Eigen::Vector<float, kNumSlots>;           // [K] values per track
  using PIDVector = Eigen::Vector<float, 3>;                    // [3] PID weights per track (pion, kaon, proton)

  // =========================================================================
  // INPUT: Track features (13 per track)
  // =========================================================================
  // Layout: [N tracks] with 13 columns → TensorCollection creates [N, 13]
  
  GENERATE_SOA_LAYOUT(TrackFeaturesLayout,
                      SOA_COLUMN(float, vz),
                      SOA_COLUMN(float, dz),
                      SOA_COLUMN(float, pt),
                      SOA_COLUMN(float, eta),
                      SOA_COLUMN(float, mva),
                      SOA_COLUMN(float, pl),
                      SOA_COLUMN(float, t_pi),
                      SOA_COLUMN(float, t_k),
                      SOA_COLUMN(float, t_p),
                      SOA_COLUMN(float, s_pi),
                      SOA_COLUMN(float, s_k),
                      SOA_COLUMN(float, s_p),
                      SOA_COLUMN(float, has_time))

  using TrackFeaturesSoA = TrackFeaturesLayout<>;

  // =========================================================================
  // OUTPUT: Combined GNN outputs (all batch=N for TensorCollection)
  // =========================================================================
  // All outputs indexed by track (N), each storing K slot values.
  // Model outputs 5 tensors: (A[N,K], z[N,K], t[N,K], p[N,K], pi[N,3])
  // 
  // A: Assignment probabilities - unique per (track, slot)
  // z, t, p: Slot predictions - replicated across tracks (same for all tracks)
  // pi: PID weights - per track [N, 3] (pion, kaon, proton)
  
  GENERATE_SOA_LAYOUT(GNNOutputLayout,
                      SOA_EIGEN_COLUMN(SlotVector, A),       // [N, K] assignment probs
                      SOA_EIGEN_COLUMN(SlotVector, z_hat),   // [N, K] z predictions (replicated)
                      SOA_EIGEN_COLUMN(SlotVector, t_hat),   // [N, K] t predictions (replicated)
                      SOA_EIGEN_COLUMN(SlotVector, p),       // [N, K] existence probs (replicated)
                      SOA_EIGEN_COLUMN(PIDVector, pi))       // [N, 3] PID weights per track

  using GNNOutputSoA = GNNOutputLayout<>;

}  // namespace vertexgnn

#endif  // RecoVertex_PrimaryVertexProducer_interface_VertexGNNSoA_h
