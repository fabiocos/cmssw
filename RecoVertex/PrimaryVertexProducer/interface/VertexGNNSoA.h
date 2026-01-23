#ifndef RecoVertex_PrimaryVertexProducer_interface_VertexGNNSoA_h
#define RecoVertex_PrimaryVertexProducer_interface_VertexGNNSoA_h

/**
 * SoA layouts for GNN vertex producer using PyTorchAlpaka.
 * 
 * Input: TrackFeaturesSoA with 13 features per track
 * Output: SlotPredictionsSoA with per-slot predictions (z_hat, t_hat, p, pi[4])
 *         AssignmentSoA for assignment matrix A[N, K]
 */

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace vertexgnn {

  // =========================================================================
  // INPUT: Track features (13 per track)
  // =========================================================================
  // Features must be contiguous in memory for TensorCollection::add to work
  // Order: vz, dz, pt, eta, mva, pl, t_pi, t_k, t_p, s_pi, s_k, s_p, has_time
  
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
  // OUTPUT: Per-slot predictions (K slots)
  // =========================================================================
  // z_hat, t_hat, p: per-slot scalars
  // pi: 4 PID weights per slot (stored as 4 separate columns to avoid Eigen)
  
  GENERATE_SOA_LAYOUT(SlotPredictionsLayout,
                      SOA_COLUMN(float, z_hat),      // predicted z position
                      SOA_COLUMN(float, t_hat),      // predicted t position  
                      SOA_COLUMN(float, p),          // existence probability
                      SOA_COLUMN(float, pi_0),       // PID weight 0
                      SOA_COLUMN(float, pi_1),       // PID weight 1
                      SOA_COLUMN(float, pi_2),       // PID weight 2
                      SOA_COLUMN(float, pi_3))       // PID weight 3

  using SlotPredictionsSoA = SlotPredictionsLayout<>;

  // =========================================================================
  // OUTPUT: Assignment matrix A[N, K]
  // =========================================================================
  // Stored as flat array of size N*K, accessed as A[i*K + k]
  
  GENERATE_SOA_LAYOUT(AssignmentLayout,
                      SOA_COLUMN(float, prob))       // assignment probability

  using AssignmentSoA = AssignmentLayout<>;

}  // namespace vertexgnn

#endif  // RecoVertex_PrimaryVertexProducer_interface_VertexGNNSoA_h
