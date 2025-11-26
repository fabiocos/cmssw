#ifndef BETA_SCORE_CLUSTERER_HPP
#define BETA_SCORE_CLUSTERER_HPP

#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

class BetaScoreClusterer {
public:
    BetaScoreClusterer(double t_beta = 0.5,
                       double min_score_thresh = 1e-2,
                       double epsilon = 1e-6,
                       double scaling_factor = 2000,
                       int min_cluster_size = 4)
        : t_beta(t_beta),
          min_score_thresh_(min_score_thresh),
          epsilon_(epsilon),
          scaling_factor_(scaling_factor),
          min_cluster_size(min_cluster_size) {}

    std::vector<int> cluster(const std::vector<std::vector<double>>& embeddings,
                             const std::vector<float>& beta,
                             const std::vector<float>& phi) {
        int N = embeddings.size();
        std::vector<int> labels(N, -1);

        // Step 1: Find valid seeds
        std::vector<int> seed_indices;
        for (int i = 0; i < N; ++i) {
            if (beta[i] >= t_beta) {
                seed_indices.push_back(i);
            }
        }
        if (seed_indices.empty()) return labels;
        // Mapping from seed idx to cluster id
        std::unordered_map<int, int> seed_to_cluster;
	for (std::size_t i = 0; i < seed_indices.size(); ++i) {
         seed_to_cluster[seed_indices[i]] = static_cast<int>(i); // compact ID
          }
        // Step 2: Score-based assignment
        for (int i = 0; i < N; ++i) {
            float best_score = -1;
            int best_seed_idx = -1;

            for (size_t j = 0; j < seed_indices.size(); ++j) {
                int seed_idx = seed_indices[j];
                float dist_sq = euclidean_dist(embeddings[i], embeddings[seed_idx]);
                float phi_sq = phi[seed_idx] * phi[seed_idx];
                float score = beta[seed_idx] / ((dist_sq / phi_sq) * scaling_factor_ + epsilon_);
                //std::cout<<" score "<<score<<std::endl;
                if (score > best_score) {
                    best_score = score;
                    best_seed_idx = seed_idx;
                }
            }

            if (best_score >= min_score_thresh_) {
std::cout<<" best_score "<<best_score<<" inside label "<<seed_to_cluster[best_seed_idx]<<" best_seed_id "<<best_seed_idx<<std::endl;
               // labels[i] = seed_to_cluster[seed_indices[best_seed_idx]];
                labels[i] = seed_to_cluster[best_seed_idx];
            }
        }

        // Step 3: Remove small clusters
        std::unordered_map<int, int> cluster_counts;
        for (int label : labels) {
            if (label != -1) cluster_counts[label]++;
        }

        std::unordered_set<int> invalid_clusters;
        for (const auto& kv : cluster_counts) {
           if (kv.second < min_cluster_size) {
               invalid_clusters.insert(kv.first);
            }
       }
        for (int i = 0; i < N; ++i) {
         if (invalid_clusters.count(labels[i])) {
           labels[i] = -1;
         }
       }
       // Step 4: Remap valid cluster IDs to compact range [0, K-1]
std::unordered_map<int, int> remap_ids;
int new_cluster_id = 0;
std::cout<<" min_cluster_size "<<min_cluster_size<<std::endl;
for (const auto& kv : cluster_counts) {
    if (kv.second >= min_cluster_size) {
        remap_ids[kv.first] = new_cluster_id++;
    }
}

// Step 5: Apply remapping
for (int& label : labels) {
    if (label != -1 && remap_ids.count(label)) {
        label = remap_ids[label];
    } else {
        label = -1;  // discard anything not in remap (just to be safe)
    }
}

        std::cout<<" labels "<<labels.size()<<std::endl;
        return labels;
    }

private:
    double t_beta;
    double min_score_thresh_;
    double epsilon_;
    double scaling_factor_;
    int min_cluster_size;

    float euclidean_dist(const std::vector<double>& a, const std::vector<double>& b) {
        float sum = 0;
        for (size_t i = 0; i < a.size(); ++i) {
            float diff = a[i] - b[i];
            sum += diff * diff;
        }
        return sum;
    }
};

#endif // BETA_SCORE_CLUSTERER_HPP

