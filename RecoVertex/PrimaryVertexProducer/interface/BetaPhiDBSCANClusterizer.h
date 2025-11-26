#ifndef BETAPHI_DBSCAN_CLUSTERIZER_HPP
#define BETAPHI_DBSCAN_CLUSTERIZER_HPP

#include <vector>
#include <cmath>
#include <queue>
#include <algorithm>

class BetaPhiDBSCANClusterizer {
public:
    // Constructor with your naming convention
    BetaPhiDBSCANClusterizer(double t_beta = 0.1, double c_reach = 0.7, double c_adapt = 1.0, int min_cluster_size = 4)
        : t_beta(t_beta), c_reach(c_reach), c_adapt(c_adapt), min_cluster_size(min_cluster_size) {}

    std::vector<int> cluster(
        const std::vector<std::vector<double>>& embeddings,
        const std::vector<float>& beta,
        const std::vector<float>& phi
    ) {
        int N = embeddings.size();
        std::vector<int> labels(N, -1);
        if (N < min_cluster_size) return labels;

        // 1. Select seed points based on beta threshold
        std::vector<int> seeds;
        for (int i = 0; i < N; ++i) {
            if (beta[i] > t_beta)
                seeds.push_back(i);
        }

        // 2. Sort seeds by beta descending
        std::sort(seeds.begin(), seeds.end(), [&](int a, int b) {
            return beta[a] > beta[b];
        });
        std::unordered_set<int> seed_set(seeds.begin(), seeds.end());
        int cluster_id = 0;
        for (int seed_idx : seeds) {
            if (labels[seed_idx] != -1) continue;

            //std::vector<int> init_neighbors = get_adaptive_neighbors(seed_idx, embeddings, phi, beta);
            //if ((int)init_neighbors.size() < min_cluster_size) continue;
            double max_reach = c_reach * phi[seed_idx];

            std::queue<int> q;
            q.push(seed_idx);
            labels[seed_idx] = cluster_id;
            std::unordered_set<int> cluster_pts = { seed_idx };

            while (!q.empty()) {
                int idx = q.front(); q.pop();
                std::vector<int> neighbors = get_adaptive_neighbors(idx, embeddings, phi, beta, seed_set);
                neighbors.erase(std::remove_if(neighbors.begin(), neighbors.end(), [&](int nb) { return seed_set.count(nb); } ),neighbors.end());
                std::cout<<" Neighbours size "<<neighbors.size()<<std::endl;

                if ((int)neighbors.size() < min_cluster_size) continue;
                int valid_neighbors = 0;
                for (int nb : neighbors) {
                    if (labels[nb] != -1) continue;
                    double dist_to_seed = euclidean_dist(embeddings[nb], embeddings[seed_idx]);
                    if (dist_to_seed > max_reach) continue;
		    ++valid_neighbors;
                    labels[nb] = cluster_id;
                    cluster_pts.insert(nb);
                    q.push(nb);
 
                }
            }

            if ((int)cluster_pts.size() < min_cluster_size) {
                for (int p : cluster_pts) labels[p] = -1;
            } else {
                cluster_id++;
            }
        }

        return labels;
    }

private:
    double t_beta;          // Seed beta threshold
    double c_reach;         // Global reachability
    double c_adapt;         // Adaptive neighborhood scaling
    int min_cluster_size;

    double euclidean_dist(const std::vector<double>& a, const std::vector<double>& b) const {
        double sum = 0.0;
        for (size_t i = 0; i < a.size(); ++i)
            sum += (a[i] - b[i]) * (a[i] - b[i]);
        return std::sqrt(sum);
    }
    std::vector<int> get_adaptive_neighbors(
        int idx,
        const std::vector<std::vector<double>>& x,
        const std::vector<float>& phi,
	const std::vector<float>& beta, const std::unordered_set<int>& seed_set
    ) const {
        std::vector<int> neighbors;
        int N = x.size();
        for (int j = 0; j < N; ++j) {
            if (j == idx) continue;
            double dist = euclidean_dist(x[idx], x[j]);
            double threshold = c_adapt * (phi[j]);
	    if (dist >= threshold)
              continue;
	     if (seed_set.count(j))
               continue;

               neighbors.push_back(j);

	      double ratio = dist / threshold;
        }
        return neighbors;
    }
};


#endif // BETAPHI_DBSCAN_CLUSTERIZER_HPP

