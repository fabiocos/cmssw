#ifndef RecoVertex_PrimaryVertexProducer_interface_VertexGNNSoA_h
#define RecoVertex_PrimaryVertexProducer_interface_VertexGNNSoA_h

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace vertexgnn {

  // Input SoA: 13 features per track
  // Layout: [N tracks] x [13 features]
  // Features: vz, dz, pt, eta, mva, pl, t_pi, t_k, t_p, s_pi, s_k, s_p, has_time
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

  // Output SoA for per-slot (vertex) predictions
  // Layout: [K slots]
  GENERATE_SOA_LAYOUT(SlotPredictionsLayout,
                      SOA_COLUMN(float, z_hat),      // predicted z position
                      SOA_COLUMN(float, t_hat),      // predicted t position
                      SOA_COLUMN(float, p),          // existence probability
                      SOA_COLUMN(float, pi_0),       // PID weight 0
                      SOA_COLUMN(float, pi_1),       // PID weight 1
                      SOA_COLUMN(float, pi_2),       // PID weight 2
                      SOA_COLUMN(float, pi_3))       // PID weight 3

  using SlotPredictionsSoA = SlotPredictionsLayout<>;

  // Output SoA for assignment matrix
  // Layout: [N tracks] - each track stores its assignment probabilities to K slots
  // Note: For dynamic K, we use a flat representation
  // The assignment matrix A[N,K] is stored as N*K floats
  GENERATE_SOA_LAYOUT(AssignmentLayout,
                      SOA_COLUMN(float, prob))  // Flat assignment probabilities

  using AssignmentSoA = AssignmentLayout<>;

}  // namespace vertexgnn

#endif  // RecoVertex_PrimaryVertexProducer_interface_VertexGNNSoA_h
