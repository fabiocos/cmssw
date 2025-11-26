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
    ObjectCondensationClustering(double t_beta = 0.1, double t_d = 1.0, int min_cluster_size = 4)
        : t_beta(t_beta), t_d(t_d), min_cluster_size(min_cluster_size) {}

    std::vector<int> cluster(
        const std::vector<std::vector<double>>& embeddings,  // N x D matrix of embeddings
        const std::vector<float>& beta,	// N vector of beta values
	const std::vector<float>& phi
    );

private:
    double euclidean_distance(const std::vector<double>& a, const std::vector<double>& b);

    double t_beta;  // Threshold for selecting condensation points
    double t_d;     // Distance threshold for clustering
    int min_cluster_size;  // Minimum cluster size for validity
};

double ObjectCondensationClustering::euclidean_distance(const std::vector<double>& a, const std::vector<double>& b) {
    double dist = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        dist += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return std::sqrt(dist);
}

// Main function to perform object condensation clustering
std::vector<int> ObjectCondensationClustering::cluster(
    const std::vector<std::vector<double>>& embeddings,  // N x D matrix of embeddings
    const std::vector<float>& beta,                      // N vector of beta values
    const std::vector<float>& phi
) {
    size_t N = embeddings.size();
    size_t D = embeddings[0].size();
    float min_phi_value = 1e-3; 
    // 1. Select condensation points (seeds) with beta > t_beta
    std::vector<int> seed_indices;
    for (size_t i = 0; i < N; ++i) {
        if (beta[i] > t_beta) {
            seed_indices.push_back(i);
        }
    }

    if (seed_indices.empty()) {
        // No seeds found; all hits are noise
        return std::vector<int>(N, -1);
    }

    // 2. Sort seeds by descending beta
    std::sort(seed_indices.begin(), seed_indices.end(), [&beta](int a, int b) {
        return beta[a] > beta[b];
    });

    // 3. Select seeds ensuring each is at least t_d away from previously selected seeds
    std::vector<int> chosen_seeds;
    for (int idx : seed_indices) {
        const std::vector<double>& emb_candidate = embeddings[idx];
        bool keep_it = true;
        for (int s_idx : chosen_seeds) {
            double phi_i = std::max(min_phi_value, phi[idx]);       // candidate seed
            double phi_j = std::max(min_phi_value, phi[s_idx]); 		
            double dist = euclidean_distance(emb_candidate, embeddings[s_idx]);
	    double scale = std::max(min_phi_value, phi[s_idx]);
           // if (dist < (t_d*scale)) {
              if (dist < t_d * std::min(phi_i, phi_j)) { 
	       keep_it = false;
                break;
            }
        }
        if (keep_it) {
            chosen_seeds.push_back(idx);
        }
    }

    // 4. Assign hits to clusters based on distance to seeds
    std::vector<int> cluster_labels(N, -1);  // Initialize all as noise
    for (size_t i = 0; i < N; ++i) {
    int best_seed = -1;
    double best_dist = std::numeric_limits<double>::infinity();
    float phi_i = std::max(phi[i], min_phi_value);
    double radius_i = t_d * phi_i;

    for (size_t cid = 0; cid < chosen_seeds.size(); ++cid) {
        int seed_idx = chosen_seeds[cid];
        double dist = euclidean_distance(embeddings[i], embeddings[seed_idx]);

        if (dist < radius_i && dist < best_dist) {
            best_dist = dist;
            best_seed = cid;
        }
    }
    if (best_seed != -1) {
        cluster_labels[i] = best_seed;
        }
    }
    /*
    for (size_t cid = 0; cid < chosen_seeds.size(); ++cid) {
        int seed_idx = chosen_seeds[cid];
	float seed_phi = std::max(phi[seed_idx], min_phi_value);
        for (size_t i = 0; i < N; ++i) {
          if (cluster_labels[i] == -1) {		
            //float scaled_phi = std::max(phi[i], min_phi_value); 		
            double dist = euclidean_distance(embeddings[i], embeddings[seed_idx]);
            if (dist < (t_d*seed_phi)) {
                cluster_labels[i] = cid;
            }
        }
      }
    }*/
    // 5. Discard clusters with fewer than min_cluster_size tracks
    std::unordered_map<int, int> cluster_counts;
    for (int label : cluster_labels) {
        if (label != -1) {
            ++cluster_counts[label];
        }
    }

    std::vector<int> valid_clusters;
    for (const auto& entry : cluster_counts) {
        if (entry.second >= min_cluster_size) {
            valid_clusters.push_back(entry.first);
        }
    }

    // Create a mask for hits belonging to valid clusters
    std::unordered_set<int> valid_clusters_set(valid_clusters.begin(), valid_clusters.end());
    std::vector<int> cluster_labels_filtered = cluster_labels;
    for (size_t i = 0; i < N; ++i) {
        if (valid_clusters_set.find(cluster_labels[i]) == valid_clusters_set.end()) {
            cluster_labels_filtered[i] = -1;  // Set to noise
        }
    }

    return cluster_labels_filtered;
}

#endif // OBJECT_CONDENSATION_CLUSTERING_HPP

