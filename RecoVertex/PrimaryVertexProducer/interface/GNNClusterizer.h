

#ifndef GNNClusterizer_h
#define GNNClusterizer_h

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

  // cached outputs from last vertices(...) call
  const std::vector<float>& lastBeta() const { return last_beta_; }              // size N
  const std::vector<float>& lastPhi() const { return last_phi_; }                // size N
  const std::vector<float>& lastPidLogits() const { return last_pid_logits_; }   // flattened N*3
  const std::vector<float>& lastEmbeddings() const { return last_embeddings_; }  // flattened N*D
  const std::vector<float>& lastPCA() const { return last_pca_flat_; }           // flattened N*3
  int lastEmbeddingDim() const { return last_embedding_dim_; }
  int lastTrackCount() const { return last_n_tracks_; }
  bool hasLastOutputs() const { return last_valid_; }

private:
  const ONNXRuntime* onnxRuntime_;
  std::string nnVersion_;
  double nnWorkingPoint_;
  std::string AlgoVersion_;
  bool verbose_;
  double zSep;
  double d0CutOff_;
  double t_beta_;
  double t_d_;
  double eps_;
  int min_cluster_size_;
  std::vector<int> pca_dim_num_;
  std::vector<double> t_d_per_dim_;
  double min_score_thresh_;
  double vertexSize_;
  mutable std::vector<float> last_beta_;        // N
  mutable std::vector<float> last_phi_;         // N
  mutable std::vector<float> last_pid_logits_;  // N*3
  mutable std::vector<float> last_embeddings_;  // N*D
  mutable std::vector<float> last_pca_flat_;    // N*3
  mutable int last_embedding_dim_ = 0;
  mutable int last_n_tracks_ = 0;
  mutable bool last_valid_ = false;
};

class UnionFind {
public:
  UnionFind(int n);
  int find(int x);
  void unite(int x, int y);

private:
  std::vector<int> parent;
};

#endif
