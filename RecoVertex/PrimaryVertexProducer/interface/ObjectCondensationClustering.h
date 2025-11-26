#ifndef OBJECT_CONDENSATION_CLUSTERING_HPP
#define OBJECT_CONDENSATION_CLUSTERING_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

class ObjectCondensationClustering {
public:
    // require_all_dims=true  -> a hit is assigned if it is within threshold on all selected PCA axes
    // require_all_dims=false -> a hit is assigned if it is within threshold on any selected PCA axis
    ObjectCondensationClustering(double t_beta = 0.1,
		                 double t_d = 1.0, 
				 int min_cluster_size = 4,
				 std::vector<int> use_dims = {},                  // ex. {0} for PC1, {0,1} for (PC1, PC2) , {0,1,2} for (PC1,PC2,PC3)
				 std::vector<double> t_d_per_dim = {},            // per axis thresholds ex. {0.25,0.15} for (PC1,PC2)
				 bool require_all_dims = true)
        : t_beta_(t_beta), 
	  t_d_global_(t_d), 
	  min_cluster_size_(min_cluster_size),
          use_dims_(std::move(use_dims)),
	  t_d_per_dim_(std::move(t_d_per_dim)),
	  require_all_dims_(require_all_dims) {}

    // N x D PCA coordinates (or embeddings); beta of size N
    std::vector<int> cluster(const std::vector<std::vector<double>>& embeddings,
        		     const std::vector<float>& beta);

private:
    // legacy L2 distance
    static double euclidean_distance(const std::vector<double>& a, const std::vector<double>& b) {
        double dist2 = 0.0;
        const size_t D = std::min(a.size(), b.size());
        for (size_t i = 0; i < D; ++i) {
            const double d = a[i] - b[i];
            dist2 += d * d;
        }
        return std::sqrt(dist2);
    }

    // fall back to global t_d if missing
    void ensure_thresholds_(size_t ndims_selected) {
        if (t_d_per_dim_.size() < ndims_selected) {
            t_d_per_dim_.resize(ndims_selected, t_d_global_);
        }
    }
    
    // compare delta on each selected axis to its threshold
    bool within_per_axis_box_(const std::vector<double>& a,
                              const std::vector<double>& b,
                              const size_t D_input) const {
        if (use_dims_.empty()) {
            // fallback to global euclidean distance
            return euclidean_distance(a, b) < t_d_global_;
        }

        int n_within = 0;
        for (size_t j = 0; j < use_dims_.size(); ++j) {
            const int d = use_dims_[j];
            if (d < 0 || static_cast<size_t>(d) >= D_input) continue;
            const double thr = (j < t_d_per_dim_.size() ? t_d_per_dim_[j] : t_d_global_);
            const double diff = std::abs(a[d] - b[d]);
            if (diff < thr) ++n_within;
        }
        return require_all_dims_ ? (n_within == static_cast<int>(use_dims_.size()))
                                 : (n_within >= 1);
    }

    double t_beta_; 
    double t_d_global_; 
    int min_cluster_size_;

    std::vector<int>    use_dims_;
    std::vector<double> t_d_per_dim_;
    bool require_all_dims_;
};

// main function
std::vector<int> ObjectCondensationClustering::cluster(
    const std::vector<std::vector<double>>& embeddings,  // N x D (PCA coord) 
    const std::vector<float>& beta                       // N
) {
    const size_t N = embeddings.size();
    if (N == 0) return {};
    const size_t D = embeddings[0].size();

    ensure_thresholds_(use_dims_.size());

    // 1) Select seed candidates with beta > t_beta_
    std::vector<int> seed_indices;
    seed_indices.reserve(N);
    for (size_t i = 0; i < N; ++i) {
        if (beta[i] > t_beta_) seed_indices.push_back(static_cast<int>(i));
    }
    if (seed_indices.empty()) return std::vector<int>(N, -1);

    // 2) Sort by beta descending
    std::sort(seed_indices.begin(), seed_indices.end(),
              [&beta](int a, int b){ return beta[a] > beta[b]; });

    // 3) suppression of seeds
    std::vector<int> chosen_seeds;
    chosen_seeds.reserve(seed_indices.size());
    for (int idx : seed_indices) {
        const auto& e_cand = embeddings[idx];
        bool keep = true;
        for (int s_idx : chosen_seeds) {
            const auto& e_seed = embeddings[s_idx];
            // If already-chosen seed, drop it
            if (within_per_axis_box_(e_cand, e_seed, D)) {
                keep = false;
                break;
            }
        }
        if (keep) chosen_seeds.push_back(idx);
    }

    if (chosen_seeds.empty()) return std::vector<int>(N, -1);

    // 4) Assign hits to clusters
    std::vector<int> cluster_labels(N, -1);
    for (size_t cid = 0; cid < chosen_seeds.size(); ++cid) {
        const int seed_idx = chosen_seeds[cid];
        const auto& e_seed = embeddings[seed_idx];

        for (size_t i = 0; i < N; ++i) {
            if (cluster_labels[i] != -1) continue;
            const auto& e_i = embeddings[i];
            if (within_per_axis_box_(e_i, e_seed, D)) {
                cluster_labels[i] = static_cast<int>(cid);
            }
        }
    }

    // 5) Filter small clusters
    std::unordered_map<int, int> cluster_counts;
    for (int lab : cluster_labels) if (lab != -1) ++cluster_counts[lab];

    std::unordered_set<int> valid;
    for (const auto& kv : cluster_counts) {
        if (kv.second >= min_cluster_size_) valid.insert(kv.first);
    }

    for (size_t i = 0; i < N; ++i) {
        if (cluster_labels[i] != -1 && valid.find(cluster_labels[i]) == valid.end()) {
            cluster_labels[i] = -1;
        }
    }
    return cluster_labels;
}

#endif // OBJECT_CONDENSATION_CLUSTERING_HPP

