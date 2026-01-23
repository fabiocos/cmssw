#ifndef GNNClusterizerFromAlpaka_h
#define GNNClusterizerFromAlpaka_h

/**
 * GNNClusterizerFromAlpaka - Builds TransientVertex objects from Alpaka SoA outputs
 *
 * This clusterizer does NOT run inference - it CONSUMES the already-computed
 * outputs from GNNVertexProducerAlpaka:
 *   - Slot predictions: z_hat, t_hat, p per slot [K]
 *   - PID weights: pi_0, pi_1, pi_2 per track [N]
 *   - Assignment probabilities: A[N,K] track-to-slot
 *
 * Pipeline:
 *   GNNVertexProducerAlpaka (runs PyTorchAlpaka inference)
 *        ↓
 *   GNNClusterizerFromAlpaka (builds vertices from SoA)
 *        ↓
 *   std::vector<TransientVertex>
 */

#include "RecoVertex/PrimaryVertexProducer/interface/TrackClusterizerInZ.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include <vector>

namespace vertexgnn {

  class GNNClusterizerFromAlpaka {
  public:
    GNNClusterizerFromAlpaka(const edm::ParameterSet& conf);

    /**
     * Build vertices from pre-computed GNN outputs
     *
     * @param tracks           Original TransientTrack collection
     * @param A_flat           Assignment matrix A[N*K] (flattened)
     * @param z_hat            Predicted z positions [K]
     * @param t_hat            Predicted t positions [K]
     * @param p                Slot existence probabilities [K]
     * @param pi_0..pi_2       PID weights per track [N] (pion, kaon, proton)
     * @param N                Number of tracks
     * @param K                Number of slots
     */
    std::vector<TransientVertex> vertices(
        const std::vector<reco::TransientTrack>& tracks,
        const float* A_flat,
        const float* z_hat,
        const float* t_hat,
        const float* p,
        const float* pi_0,
        const float* pi_1,
        const float* pi_2,
        int N,
        int K
    ) const;

    static void fillPSetDescription(edm::ParameterSetDescription& desc);

  private:
    double existenceThreshold_;
    double trackAssignmentThreshold_;
    bool verbose_;
  };

}  // namespace vertexgnn

#endif
