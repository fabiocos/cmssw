/**
 * GNNClusterizer.cc - End-to-end Vertex Slot Model based track clustering
 * 
 * This implementation uses a VertexSlotModel ONNX export that directly outputs
 * track-to-vertex assignments and vertex positions.
 * 
 * ONNX Model Inputs:
 *   - x: (1, N, 13) track features
 * 
 * ONNX Model Outputs:
 *   - A: (N, K) soft assignment probabilities
 *   - z_hat: (1, K) predicted vertex z positions [cm]
 *   - t_hat: (1, K) predicted vertex times [ns]
 *   - p: (1, K) slot existence probabilities
 *   - pi: (N, 3) PID mixture weights
 */

#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "RecoVertex/PrimaryVertexProducer/interface/GNNClusterizer.h"
#include "RecoVertex/VertexPrimitives/interface/VertexException.h"
#include "RecoVertex/PrimaryVertexProducer/interface/DAClusterizerInZ_vect.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/isFinite.h"

#include <iostream>
#include <unordered_map>
#include <array>
#include <cmath>
#include <cassert>
#include <limits>
#include <iomanip>
#include <map>
#include <vector>
#include <algorithm>

#define cputime
#ifdef cputime
#include <chrono>
typedef std::chrono::duration<int, std::micro> microseconds_type;
#endif

using namespace std;

/**
 * Combine timing uncertainties in quadrature
 */
float combine_sigma(float sig_h, float sig_tmtd, float eps = 1e-12f) {
    float h = std::max(sig_h, 0.0f);
    float t = std::max(sig_tmtd, 0.0f);
    float out = std::sqrt(h*h + t*t + eps);
    if (!std::isfinite(out)) {
        out = 1e9f;
    }
    return out;
}

GNNClusterizer::GNNClusterizer(const edm::ParameterSet& conf, const ONNXRuntime* onnxRuntime)
    : onnxRuntime_(onnxRuntime),
      nnVersion_(conf.getParameter<std::string>("nnVersion")),
      verbose_(conf.getUntrackedParameter<bool>("verbose", true)),
      d0CutOff_(conf.getParameter<double>("d0CutOff")),
      vertexSize_(conf.getParameter<double>("vertexSize")),
      existenceThreshold_(conf.getParameter<double>("existenceThreshold")),
      trackAssignmentThreshold_(conf.getParameter<double>("trackAssignmentThreshold")),
      numSlots_(conf.getParameter<int>("numSlots")) {}

void GNNClusterizer::globalEndJob(const ONNXRuntime* cache) {}

std::unique_ptr<ONNXRuntime> GNNClusterizer::initializeGlobalCache(const edm::ParameterSet& conf) {
  std::string backend_str = conf.getParameter<std::string>("onnxBackend");
  cms::Ort::Backend backend = cms::Ort::Backend::cpu;
  if (backend_str == "CUDA") {
    backend = cms::Ort::Backend::cuda;
  }

  auto session_options = ONNXRuntime::defaultSessionOptions(backend);
  edm::LogInfo("GNNClusterizer") << "Initializing ONNXRuntime with backend: " << backend_str;
  return std::make_unique<ONNXRuntime>(conf.getParameter<edm::FileInPath>("onnxModelPath").fullPath(), &session_options);
}

std::vector<TransientVertex> GNNClusterizer::vertices(const std::vector<reco::TransientTrack>& tracks) const {
  cms::Ort::FloatArrays input_values;
  std::vector<TransientVertex> clusters;
  int N = tracks.size();
  int K = numSlots_;
  int shapeFeatures = 13;  // Input feature dimension (matches model.py)

  // Reset caches
  last_valid_ = false;
  last_n_tracks_ = N;
  last_num_slots_ = K;
  last_assignments_.clear();
  last_z_hat_.clear();
  last_t_hat_.clear();
  last_p_.clear();
  last_pi_.clear();

  if (N == 0) {
    return clusters;
  }

  // =========================================================================
  // PREPARE INPUT FEATURES (matching model.py input_dim=13)
  // =========================================================================
  // Features: [z, dz, pt, eta, mva, pl, t_pi, t_k, t_p, s_pi, s_k, s_p, has_time]
  
  std::vector<float> features;
  features.reserve(N * shapeFeatures);

  for (size_t i = 0; i < tracks.size(); ++i) {
    const reco::TransientTrack& ttrack = tracks[i];
    
    // Basic track properties
    float track_vz = ttrack.track().vz();
    float track_dz = ttrack.track().dzError();
    float track_pt = ttrack.track().pt();
    float track_eta = ttrack.track().eta();
    float track_mva = ttrack.MVAquality();
    float track_pl = ttrack.pathLength();

    // Ensure missing MVA/PL are -1.0f
    if (!std::isfinite(track_mva)) track_mva = -1.0f;
    if (!std::isfinite(track_pl))  track_pl  = -1.0f;
    
    // Timing information
    float track_t_pi = ttrack.MTDtime() - ttrack.trackTime_pi();
    float track_t_k = ttrack.MTDtime() - ttrack.trackTime_k();
    float track_t_p = ttrack.MTDtime() - ttrack.trackTime_p();
    float track_s_tmtd = ttrack.MTDtimeErr();
    float track_s_pi = ttrack.sigma_time_pi();
    float track_s_k = ttrack.sigma_time_k();
    float track_s_p = ttrack.sigma_time_p();

    // Build timing masks
    bool tmtd_ok = (track_s_tmtd >= 0.0f);
    bool pi_ok = (track_s_pi >= 0.0f) && tmtd_ok;
    bool k_ok = (track_s_k >= 0.0f) && tmtd_ok;
    bool p_ok = (track_s_p >= 0.0f) && tmtd_ok;

    // Apply masks to timing features
    float feat_t_pi = 0.0f, feat_t_k = 0.0f, feat_t_p = 0.0f;
    float feat_s_pi = 0.2f, feat_s_k = 0.2f, feat_s_p = 0.2f;

    if (pi_ok) {
        feat_t_pi = track_t_pi;
        feat_s_pi = combine_sigma(track_s_pi, track_s_tmtd);
    }
    if (k_ok) {
        feat_t_k = track_t_k;
        feat_s_k = combine_sigma(track_s_k, track_s_tmtd);
    }
    if (p_ok) {
        feat_t_p = track_t_p;
        feat_s_p = combine_sigma(track_s_p, track_s_tmtd);
    }

    // has_time flag
    float has_time = tmtd_ok ? 1.0f : 0.0f;

    // Push all 13 features
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
    features.push_back(has_time);
  }

  // Replace NaN/Inf
  for (auto& v : features) {
    if (!std::isfinite(v)) v = 1e9f;
  }

  // =========================================================================
  // RUN ONNX INFERENCE
  // =========================================================================
  std::vector<std::string> input_names = {"x"};
  input_values.emplace_back(features);

  std::vector<std::vector<long int>> input_dims;
  input_dims.push_back({1, N, shapeFeatures});

#ifdef cputime
  std::chrono::duration<int, std::micro> tcpu_inference(0), tcpu_out(0);
  auto start_inference = std::chrono::high_resolution_clock::now();
#endif

  auto output_values = onnxRuntime_->run(input_names, input_values, input_dims);

#ifdef cputime
  auto stop_inference = std::chrono::high_resolution_clock::now();
  tcpu_inference = std::chrono::duration_cast<std::chrono::microseconds>(stop_inference - start_inference);
  edm::LogInfo("GNNClusterizer") << "###TIME inference " << tcpu_inference.count() << " us";
#endif

  // =========================================================================
  // PARSE ONNX OUTPUTS
  // =========================================================================
  std::vector<float>& A_flat = output_values[0];      // N*K
  std::vector<float>& z_hat_flat = output_values[1];  // K
  std::vector<float>& t_hat_flat = output_values[2];  // K
  std::vector<float>& p_flat = output_values[3];      // K
  std::vector<float>& pi_flat = output_values[4];     // N*3

  // Cache raw outputs
  last_assignments_ = A_flat;
  last_z_hat_ = z_hat_flat;
  last_t_hat_ = t_hat_flat;
  last_p_ = p_flat;
  last_pi_ = pi_flat;
  last_valid_ = true;

  // Infer K from z_hat size
  K = static_cast<int>(z_hat_flat.size());
  last_num_slots_ = K;

  if (verbose_) {
    edm::LogInfo("GNNClusterizer") << "ONNX outputs: N=" << N << ", K=" << K;
  }

#ifdef cputime
  auto start_out = std::chrono::high_resolution_clock::now();
#endif

  // =========================================================================
  // PROCESS ASSIGNMENTS
  // =========================================================================
  
  // Step 1: For each track, find the slot with maximum assignment probability
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

  // Step 3: Identify active slots
  std::vector<bool> slot_active(K, false);
  int num_active = 0;
  for (int k = 0; k < K; ++k) {
    if (p_flat[k] > existenceThreshold_ && slot_track_count[k] >= 1) {
      slot_active[k] = true;
      num_active++;
    }
  }

  if (verbose_) {
    edm::LogInfo("GNNClusterizer") << "Active slots: " << num_active;
  }

  // Step 4: Build vertex clusters with GNN weights
  GlobalError dummyError(0.01, 0, 0.01, 0., 0., 0.01);
  
  for (int k = 0; k < K; ++k) {
    if (!slot_active[k]) continue;

    // Collect tracks assigned to this slot
    std::vector<reco::TransientTrack> cluster_tracks;
    std::vector<float> cluster_weights;
    
    for (int i = 0; i < N; ++i) {
      if (track_to_slot[i] == k) {
        cluster_tracks.push_back(tracks[i]);
        cluster_weights.push_back(track_max_prob[i]);
      }
    }

    if (cluster_tracks.empty()) continue;

    // Use GNN's predicted z position (from z_hat output)
    double z_cluster = static_cast<double>(z_hat_flat[k]);
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
      edm::LogInfo("GNNClusterizer") << "  Vertex k=" << k 
                                     << ": z=" << z_cluster
                                     << ", p=" << p_flat[k]
                                     << ", ntracks=" << cluster_tracks.size();
    }
  }

#ifdef cputime
  auto stop_out = std::chrono::high_resolution_clock::now();
  tcpu_out = std::chrono::duration_cast<std::chrono::microseconds>(stop_out - start_out);
  edm::LogInfo("GNNClusterizer") << "###TIME clustering " << tcpu_out.count() << " us";
#endif

  if (verbose_) {
    edm::LogInfo("GNNClusterizer") << "Final vertex count: " << clusters.size();
  }

  return clusters;
}

std::vector<std::vector<reco::TransientTrack>> GNNClusterizer::clusterize(
    const std::vector<reco::TransientTrack>& tracks) const {
  return std::vector<std::vector<reco::TransientTrack>>();
}

void GNNClusterizer::fillPSetDescription(edm::ParameterSetDescription& desc) {
  DAClusterizerInZ_vect::fillPSetDescription(desc);
  
  // ONNX model configuration
  desc.add<edm::FileInPath>("onnxModelPath", 
      edm::FileInPath("RecoVertex/PrimaryVertexProducer/data/vertex_slot_model.onnx"))
      ->setComment("Path to the VertexSlotModel ONNX export");
  desc.add<std::string>("nnVersion", "vertex_slot_v8p9")
      ->setComment("Model version tag");
  
  // VertexSlotModel-specific parameters
  desc.add<double>("existenceThreshold", 0.5)
      ->setComment("Slot existence probability threshold (0-1)");
  desc.add<double>("trackAssignmentThreshold", 0.0)
      ->setComment("Minimum track assignment probability");
  desc.add<int>("numSlots", 200)
      ->setComment("Number of slots K in the model");
  desc.add<std::string>("onnxBackend", "CPU")
      ->setComment("ONNX execution backend: CPU, CUDA");

  desc.addUntracked<bool>("verbose", false);
}
