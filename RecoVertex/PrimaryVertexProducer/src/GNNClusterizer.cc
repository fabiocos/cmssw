#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "RecoVertex/PrimaryVertexProducer/interface/GNNClusterizer.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"
#include "RecoVertex/PrimaryVertexProducer/interface/DAClusterizerInZ_vect.h"
#include "RecoVertex/PrimaryVertexProducer/interface/ObjectCondensationClustering.h"
#include "RecoVertex/PrimaryVertexProducer/interface/DBSCANClusterizer.h"
#include "RecoVertex/PrimaryVertexProducer/interface/BetaPhiDBSCANClusterizer.h"
#include "RecoVertex/PrimaryVertexProducer/interface/BetaScoreClusterer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/isFinite.h"

#include "vdt/vdtMath.h"
#include <iostream>
#include <unordered_map>
#include <array>
#include "TPrincipal.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include <cmath>
#include <cassert>
#include <limits>
#include <iomanip>
#include <map>
#include <vector>
#include <math.h>
#include <queue>

#define cputime
#ifdef cputime
#include <chrono>
typedef std::chrono::duration<int, std::micro> microseconds_type;
#endif

using namespace std;
enum class ClusterAlgo {
    Object,
    DBSCAN,
    BetaDBSCAN,
    OCScore
};
std::unordered_map<std::string, ClusterAlgo> algoMap = {
    {"ObjectCond", ClusterAlgo::Object},
    {"dbscan", ClusterAlgo::DBSCAN}, 
    {"beta_phi_dbscan", ClusterAlgo::BetaDBSCAN},
    {"OC_Score", ClusterAlgo::OCScore}
};

float combine_sigma(float sig_h, float sig_tmtd, float eps = 1e-12f) {
    // Step 1: ensure non-negative inputs
    float h = std::max(sig_h, 0.0f);
    float t = std::max(sig_tmtd, 0.0f);

    // Step 2: compute sqrt(h^2 + t^2 + eps)
    float out = std::sqrt(h*h + t*t + eps);

    // Step 3: replace NaN or Inf with large number (1e9)
    if (!std::isfinite(out)) {
        out = 1e9f;
    }

    return out;
}

GNNClusterizer::GNNClusterizer(const edm::ParameterSet& conf, const ONNXRuntime* onnxRuntime)
    : onnxRuntime_(onnxRuntime),
      nnVersion_(conf.getParameter<std::string>("nnVersion")),
      nnWorkingPoint_(conf.getParameter<double>("nnWorkingPoint")),
      AlgoVersion_(conf.getParameter<std::string>("AlgoVersion")),	
      verbose_(conf.getUntrackedParameter<bool>("verbose", true)),
      zSep(conf.getParameter<double>("zSeparation")),
      d0CutOff_(conf.getParameter<double>("d0CutOff")),
      t_beta_(conf.getParameter<double>("t_beta")),
      t_d_(conf.getParameter<double>("t_d")),
      eps_(conf.getParameter<double>("eps")),	
      min_cluster_size_(conf.getParameter<int>("min_cluster_size")),	
      pca_dim_num_(conf.getParameter<std::vector<int>>("pca_dim_num")),	
      t_d_per_dim_(conf.getParameter<std::vector<double>>("t_d_per_dim")),	
      min_score_thresh_(conf.getParameter<double>("min_score_thresh")),	
      vertexSize_(conf.getParameter<double>("vertexSize")) {}

void GNNClusterizer::globalEndJob(const ONNXRuntime* cache) {}

namespace {
  inline double local_exp(double const& inp) { return vdt::fast_exp(inp); }
}  // namespace

std::unique_ptr<ONNXRuntime> GNNClusterizer::initializeGlobalCache(const edm::ParameterSet& conf) {
  return std::make_unique<ONNXRuntime>(conf.getParameter<edm::FileInPath>("onnxModelPath").fullPath());
}

std::vector<TransientVertex> GNNClusterizer::vertices(const std::vector<reco::TransientTrack>& tracks) const {
  std::string clusterAlgo_;

  auto it = algoMap.find(AlgoVersion_);
     if (it != algoMap.end()) {
      clusterAlgo_ = it->first;
    } else {
      throw std::invalid_argument("Invalid algorithm name: " + AlgoVersion_);
    }

  cms::Ort::FloatArrays input_values;     // Stores float inputs for ONNX
  std::vector<TransientVertex> clusters;  // Output cluster vertices
  int N = tracks.size();                  // Number of tracks
  int shapeFeatures = 12;                 // Dim input features 

  // Reset caches for this call
  last_valid_ = false;
  last_n_tracks_ = N;
  last_embedding_dim_ = 0;
  last_beta_.clear();
  last_phi_.clear();
  last_pid_logits_.clear();
  last_embeddings_.clear();
  last_pca_flat_.clear();


  if (N == 0) {
    return clusters;  // Return empty if no tracks
  }

  // Boolean masks
  std::vector<bool> tmtd_ok(N, false);
  std::vector<bool> pi_ok(N, false);
  std::vector<bool> k_ok(N, false);
  std::vector<bool> p_ok(N, false);


  // Prepare input data: features ('x')
  std::vector<float> features;       // 'x'
  std::vector<float> row_splits;

  features.reserve(N * shapeFeatures);

  for (size_t i = 0; i < tracks.size(); ++i) {
    const reco::TransientTrack& ttrack = tracks[i];
    float track_vz = ttrack.track().vz(); // Feature: vz
    float track_dz = ttrack.track().dzError(); // Feature: dz
    float track_pt = ttrack.track().pt(); // Feature: pt
    float track_eta = ttrack.track().eta(); // Feature: eta
    float track_mva = ttrack.MVAquality(); // Feature: mva
    float track_pl = ttrack.pathLength(); // Feature: PL
    float track_t_pi = ttrack.MTDtime() - ttrack.trackTime_pi(); // Feature: t_pi
    float track_t_k = ttrack.MTDtime() - ttrack.trackTime_k(); // Feature: t_k
    float track_t_p = ttrack.MTDtime() - ttrack.trackTime_p(); // Feature: t_p
    float track_s_tmtd = ttrack.MTDtimeErr();
    float track_s_pi = ttrack.sigma_time_pi(); 
    float track_s_k  = ttrack.sigma_time_k(); 
    float track_s_p  = ttrack.sigma_time_p();  

    // Build masks
    tmtd_ok[i] = (track_s_tmtd >= 0.0f);
    pi_ok[i]   = (track_s_pi   >= 0.0f) && tmtd_ok[i];
    k_ok[i]    = (track_s_k    >= 0.0f) && tmtd_ok[i];
    p_ok[i]    = (track_s_p    >= 0.0f) && tmtd_ok[i];

    float feat_t_pi = 0.0f, feat_t_k = 0.0f, feat_t_p = 0.0f;
    float feat_s_pi = 0.2f, feat_s_k = 0.2f, feat_s_p = 0.2f;

    if (pi_ok[i]) {
        feat_t_pi = track_t_pi;
        feat_s_pi = combine_sigma(track_s_pi, track_s_tmtd);
    }
    if (k_ok[i]) {
        feat_t_k = track_t_k;
        feat_s_k = combine_sigma(track_s_k, track_s_tmtd);
    }
    if (p_ok[i]) {
        feat_t_p = track_t_p;
        feat_s_p = combine_sigma(track_s_p, track_s_tmtd);
    }

    features.push_back(track_vz);
    features.push_back(track_dz);
    features.push_back(track_pt);
    features.push_back(track_eta);
    features.push_back(track_mva);
    features.push_back(track_pl);
    features.push_back(feat_t_pi);
    features.push_back(feat_t_k);
    features.push_back(feat_t_p);
    features.push_back(feat_s_pi);
    features.push_back(feat_s_k);
    features.push_back(feat_s_p);
  }

  for (auto &v : features) {
    if (!std::isfinite(v)) v = 1e9f;  // replace NaN/Inf
  }


  row_splits.push_back(0);
  row_splits.push_back(N);

  // Define input names consistent with ONNX model
  std::vector<std::string> input_names = {"x", "row_splits"} ;

  // Add float inputs to input_values
  input_values.emplace_back(features);  // 'x'
  input_values.emplace_back(row_splits);

  // Define input dimensions
  std::vector<std::vector<long int>> input_dims;
  input_dims.push_back({1, N, shapeFeatures});  // 'x': [1, N, shapeFeatures]
  input_dims.push_back({1, 2});                 // 'row_dimensions': [1, N]

#ifdef cputime
  std::chrono::duration<int, std::micro> tcpu_inference(0), tcpu_out(0);
  auto start_inference = std::chrono::high_resolution_clock::now();
#endif


  // Run ONNX inference
  auto output_values = onnxRuntime_->run(input_names, input_values, input_dims);

#ifdef cputime
  auto stop_inference = std::chrono::high_resolution_clock::now();
  tcpu_inference = std::chrono::duration_cast<std::chrono::microseconds>(stop_inference - start_inference);
  edm::LogInfo("PrimaryVertexProducer") << "###TIME inference " << tcpu_inference;
#endif


  std::vector<float>& beta_predictions = output_values[0];
  std::vector<float>& embeddings_flat = output_values[1];
  std::vector<float>& learned_phi = output_values[2];
  std::vector<float>& pid_logits = output_values[3];

  // --- Cache raw outputs (flattened) ---
  last_beta_ = beta_predictions;                 // size N
  last_phi_  = learned_phi;                      // size N
  last_pid_logits_ = pid_logits;                 // size N*3
 
  // Determine embedding dimension
  int embedding_dim = embeddings_flat.size() / N;

  last_embedding_dim_ = embedding_dim;
  last_embeddings_ = embeddings_flat;            // size N*D  

  // Reshape into [N][embedding_dim]
  std::vector<std::vector<double>> embedding_vectors(N, std::vector<double>(embedding_dim));

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < embedding_dim; ++j) {
      embedding_vectors[i][j] = embeddings_flat[i * embedding_dim + j];
    }
  }

    // Reshape PID logits into [N][3] + normalize
  std::vector<std::vector<double>> pid_vectors(N, std::vector<double>(3));

  for (int i = 0; i < N; ++i) {
    double l0 = pid_logits[i * 3 + 0];
    double l1 = pid_logits[i * 3 + 1];
    double l2 = pid_logits[i * 3 + 2];

    double max_l = std::max({l0, l1, l2});
    double e0 = std::exp(l0 - max_l);
    double e1 = std::exp(l1 - max_l);
    double e2 = std::exp(l2 - max_l);
    double sum = e0 + e1 + e2;

    pid_vectors[i][0] = e0 / sum;
    pid_vectors[i][1] = e1 / sum;
    pid_vectors[i][2] = e2 / sum;
  }

  // Make PCA on embeddings
  TPrincipal pca(embedding_dim, "D");

  // Add each embedding as a row
  for (int i = 0; i < N; ++i) {
      Double_t row[embedding_dim];
      for (int j = 0; j < embedding_dim; ++j) {
      row[j] = embedding_vectors[i][j];
      }
      pca.AddRow(row);
  }

  // Compute principal components
  pca.MakePrincipals();

  // Compute PCA first 3 components into last_pca_flat_ (N*3) ----
last_pca_flat_.assign(static_cast<size_t>(N) * 3, std::numeric_limits<float>::quiet_NaN());
for (int i = 0; i < N; ++i) {
  std::vector<Double_t> row(std::max(1, embedding_dim)), proj(std::max(1, embedding_dim));
  const Double_t* orig = pca.GetRow(i);
  for (int j = 0; j < embedding_dim; ++j) row[j] = orig[j];
  pca.X2P(row.data(), proj.data());
  for (int c = 0; c < 3; ++c) {
    float val = (c < embedding_dim) ? static_cast<float>(proj[c])
                                    : std::numeric_limits<float>::quiet_NaN();
    last_pca_flat_[3 * i + c] = val;
  }
}

  // ---- Mark caches valid ----
  last_valid_ = true;


 // Build PCA vectors for clustering from last_pca_flat_
std::vector<std::vector<double>> pca_vec;
pca_vec.reserve(N);
for (int i = 0; i < N; ++i) {
  std::vector<double> v(3, std::numeric_limits<double>::quiet_NaN());
  v[0] = static_cast<double>(last_pca_flat_[3 * i + 0]);
  v[1] = static_cast<double>(last_pca_flat_[3 * i + 1]);
  v[2] = static_cast<double>(last_pca_flat_[3 * i + 2]);
  pca_vec.emplace_back(std::move(v));
}

#ifdef cputime
  auto start_out = std::chrono::high_resolution_clock::now();
#endif

  //double eps = 0.00025;
  double epsilon_ = 1e-6;
  double scaling_factor_ = 2000;
  double c_adapt = 0.5;
  std::vector<int> labels;
  if (clusterAlgo_ == "dbscan") {
     DBSCANClusterizer dbscan(eps_, min_cluster_size_);
     labels = dbscan.run_clustering(embedding_vectors);
  } 
   else if (clusterAlgo_ == "ObjectCond") {    
     ObjectCondensationClustering clustering(t_beta_, t_d_, min_cluster_size_, pca_dim_num_, t_d_per_dim_, true);
     //std::cout<<" t_beta_ "<<t_beta_<<" t_d_ "<<t_d_<<std::endl;
     labels = clustering.cluster(pca_vec,beta_predictions);
  }
    else if (clusterAlgo_ == "beta_phi_dbscan") {
	    //std::cout<<" Algo betadBSCAN t_beta_ "<<t_beta_<<" t_d_ "<<t_d_<<std::endl;
    BetaPhiDBSCANClusterizer beta_phi_clusterer(t_beta_,t_d_, c_adapt, min_cluster_size_);
    labels = beta_phi_clusterer.cluster(embedding_vectors,beta_predictions, learned_phi);
  }
  else {
     BetaScoreClusterer OC_Scoreclusterer(t_beta_, min_score_thresh_, epsilon_, scaling_factor_, min_cluster_size_);
     labels = OC_Scoreclusterer.cluster(embedding_vectors,beta_predictions, learned_phi);	     
  }	  
     
#ifdef cputime
  auto stop_out = std::chrono::high_resolution_clock::now();
  tcpu_out = std::chrono::duration_cast<std::chrono::microseconds>(stop_out - start_out);
  edm::LogInfo("GNNClusterizer") << "###TIME out " << tcpu_out;
#endif

  std::map<int, std::vector<int>> clusters_map;  // Map from cluster id to list of track indices
  //std::cout<<" labels size "<<labels.size()<<std::endl;
  for (size_t i = 0; i < labels.size(); ++i) {
    int label = labels[i];
  //  int idx = selected_nodes[i];  // Original index in tracks
    if (label == -1) {
      continue;
    }
    clusters_map[label].push_back(i);
  }

  GlobalError dummyError(0.01, 0, 0.01, 0., 0., 0.01);

  for (const auto& cluster_pair : clusters_map) {
    const std::vector<int>& cluster_indices = cluster_pair.second;
    std::vector<reco::TransientTrack> cluster_tracks;
    double sumwz = 0.0;
    double sumw = 0.0;
    //double t_tkwt = 1.0;
    for (int idx : cluster_indices) {
      const reco::TransientTrack& track = tracks[idx];
      cluster_tracks.push_back(track);

      //// double weight = 1.0 / (track.track().dzError() * track.track().dzError());
      //if (d0CutOff_ > 0) {
      //  Measurement1D atIP = track.stateAtBeamLine().transverseImpactParameter();  // error contains beamspot
      //  double t_tkwt = 1. / (1. + local_exp(std::pow(atIP.value() / atIP.error(), 2) -
      //                                       std::pow(d0CutOff_, 2)));  // reduce weight for high ip tracks
      //  if (edm::isNotFinite(t_tkwt) || t_tkwt < std::numeric_limits<double>::epsilon()) {
      //    edm::LogWarning("GNN Clusterizer") << "rejected track t_tkwt " << t_tkwt;
      //    continue;  // usually is > 0.99
      //  }
      //}  // d0 cutoff

      //auto const& t_mom = track.stateAtBeamLine().trackStateAtPCA().momentum();
      //reco::BeamSpot beamspot = track.stateAtBeamLine().beamSpot();
      //double t_dz2 = std::pow(track.track().dzError(), 2)  // track errror
      //               +
      //               (std::pow(beamspot.BeamWidthX() * t_mom.x(), 2) + std::pow(beamspot.BeamWidthY() * t_mom.y(), 2)) *
      //                   std::pow(t_mom.z(), 2) / std::pow(t_mom.perp2(), 2)  // beam spot width
      //               + std::pow(vertexSize_, 2);  // intrinsic vertex size, safer for outliers and short lived decays
      //t_dz2 = 1. / t_dz2;
      //if (edm::isNotFinite(t_dz2) || t_dz2 < std::numeric_limits<double>::min()) {
      //  edm::LogWarning("GNN Clusterizer") << "rejected track t_dz2 " << t_dz2;
      //  continue;
      //}

      //double weight = t_tkwt * t_dz2;
      double weight = 1.;
      double z = track.track().vz();
      sumwz += weight * z;
      sumw += weight;
    }
    double z_cluster = sumwz / sumw;
    //if (verbose_) {
    //           std::cout << "z_cluster:" << z_cluster << std::endl;
    //           }
    GlobalPoint pos(0, 0, z_cluster);
    TransientVertex vertex(pos, dummyError, cluster_tracks, 0);
    clusters.push_back(vertex);
  }
  //std::cout << "cluster size:" << clusters.size() << std::endl;
  return clusters;
}

std::vector<std::vector<reco::TransientTrack>> GNNClusterizer::clusterize(
    const std::vector<reco::TransientTrack>& tracks) const {
  return std::vector<std::vector<reco::TransientTrack>>();
}

void GNNClusterizer::fillPSetDescription(edm::ParameterSetDescription& desc) {
  DAClusterizerInZ_vect::fillPSetDescription(desc);
  desc.add<double>("zSeparation", 0.001)->setComment("Epsilon value for DBSCAN clustering");
  desc.add<double>("eps", 0.001)->setComment("Epsilon DBSCAN clustering");
  desc.addUntracked<bool>("verbose", true);
  desc.add<std::string>("AlgoVersion", "dbscan")->setComment("Clustering Algorithm for GravNet");
  desc.add<edm::FileInPath>("onnxModelPath", edm::FileInPath("RecoVertex/PrimaryVertexProducer/data/gravnet_da.onnx"))->setComment("Path to the ONNX model");
  desc.add<double>("t_beta", 0.7)->setComment("Beta threshold for node selection.");
  desc.add<double>("t_d", 1.0)->setComment("Eucl Distance.");
  desc.add<int>("min_cluster_size", 4)->setComment("Min tracks in a cluster");
  desc.add<std::vector<int>>("pca_dim_num", {0,1,2})->setComment("Number PCA components to use for clustering");
  desc.add<std::vector<double>>("t_d_per_dim", {0.05,0.005,0.005})->setComment("Thresholds per PCA component");
  desc.add<double>("min_score_thresh", 1e-2)->setComment("Minimum score threshold");  
  desc.add<std::string>("nnVersion", "gravnet_v1")->setComment("GNN version tag.");
  desc.add<double>("nnWorkingPoint", 0.00)->setComment("Beta threshold for node selection.");
}
