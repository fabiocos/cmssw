#ifndef DBSCANCLUSTERIZER_H
#define DBSCANCLUSTERIZER_H

#include <vector>
#include <queue>
#include <cmath>
#include <iostream>

class DBSCANClusterizer {
public:
    DBSCANClusterizer(double eps, size_t min_pts)
        : eps_squared(eps * eps), min_pts(min_pts) {}

    std::vector<int> run_clustering(const std::vector<std::vector<double>>& points) {
        const size_t n = points.size();
        std::vector<int> labels(n, -1);
        std::vector<bool> visited(n, false);
        std::vector<bool> is_core(n, false);
        int cluster_id = 0;

        // Identify core points
        for (size_t i = 0; i < n; ++i) {
            if (find_neighbors(points, i).size() >= min_pts) {
                is_core[i] = true;
            }
        }

        // Expand clusters
        for (size_t point_idx = 0; point_idx < n; ++point_idx) {
            if (visited[point_idx] || !is_core[point_idx]) continue;

            cluster_id++;
            labels[point_idx] = cluster_id;
            visited[point_idx] = true;

            std::queue<size_t> seed_queue;
            seed_queue.push(point_idx);

            while (!seed_queue.empty()) {
                size_t current = seed_queue.front();
                seed_queue.pop();

                std::vector<size_t> neighbors = find_neighbors(points, current);
                for (size_t neighbor : neighbors) {
                    if (!visited[neighbor]) {
                        visited[neighbor] = true;
                        labels[neighbor] = cluster_id;
                        if (is_core[neighbor]) {
                            seed_queue.push(neighbor);
                        }
                    }
                }
            }
        }

        return labels;
    }

private:
    double eps_squared;
    size_t min_pts;

    double euclidean_distance(const std::vector<double>& a, const std::vector<double>& b) {
        double dist_squared = 0.0;
        for (size_t i = 0; i < a.size(); ++i) {
            double diff = a[i] - b[i];
            dist_squared += diff * diff;
        }
        return dist_squared;  // Squared distance to avoid sqrt
    }

    std::vector<size_t> find_neighbors(const std::vector<std::vector<double>>& points, size_t point_idx) {
        std::vector<size_t> neighbors;
        for (size_t j = 0; j < points.size(); ++j) {
            if (euclidean_distance(points[point_idx], points[j]) <= eps_squared) {
                neighbors.push_back(j);
            }
        }
        return neighbors;
    }
};

#endif // DBSCANCLUSTERIZER_H

