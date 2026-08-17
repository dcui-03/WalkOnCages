#include "polynet.hpp"

#include <Eigen/Core>
#include <vector>
#include <utility>
#include <limits>
#include <queue>
#include <algorithm>

// Closest-point BVH over edges
namespace Polynet {
    // Compute bounding volume hierarchy
    int polynet::computeBVH() {
        BVH.clear();

        std::vector<int> activeEdges;
        activeEdges.reserve(E.size());
        for (int e = 0; e < E.size(); e++) {
            if (E[e].active) {
                activeEdges.push_back(e);
            }
        }
        if (activeEdges.empty()) {
            return -1;
        }

        int root = buildBVHNode(activeEdges, 0);
        if (root < 0) {
            return -1;
        }
        max_depth = 0;
        for (int i = 0; i < BVH.size(); i++) {
            max_depth = std::max(BVH[i].depth, max_depth);
        }
        return 1;
    }

    int polynet::buildBVHNode(const std::vector<int>& edges, int depth) {
        if (edges.empty()) {
            return -1;
        }
        int max_depth = 16;       // Hard stop so we don't recurse too much
        int leafEdgeCount = 50;   // Maximum number of edges a leaf node can have
        int nodeIdx = static_cast<int>(BVH.size());
        BVH.emplace_back();

        BVH[nodeIdx].depth = depth;
        BVH[nodeIdx].leaf = false;
        BVH[nodeIdx].edges.clear();
        BVH[nodeIdx].children.clear();

        // 1. Compute tight AABB around all edges in this node
        bool found = false;
        Eigen::Vector3d minV;
        Eigen::Vector3d maxV;
        for (int e_idx = 0; e_idx < edges.size(); e_idx++) {
            int e = edges[e_idx];
            if (e < 0 || e >= E.size() || !E[e].active) {
                continue;
            }
            int he = E[e].he;
            Eigen::Vector3d p0 = V[HE[HE[he].twin].dest].new_pos;
            Eigen::Vector3d p1 = V[HE[he].dest].new_pos;
            Eigen::Vector3d pts[2] = {p0, p1};
            for (int p_idx = 0; p_idx < 2; p_idx++) {
                const Eigen::Vector3d& p = pts[p_idx];
                if (!found) {
                    minV = p;
                    maxV = p;
                    found = true;
                } else {
                    minV = minV.cwiseMin(p);
                    maxV = maxV.cwiseMax(p);
                }
            }
        }
        if (!found) {
            return -1;
        }
        BVH[nodeIdx].bdyVerts = {minV, maxV};

        // 2. Leaf stopping criteria
        if (depth >= max_depth || edges.size() <= leafEdgeCount) {
            BVH[nodeIdx].leaf = true;
            BVH[nodeIdx].edges = edges;
            return nodeIdx;
        }

        // 3. Choose split axis using longest AABB extent
        Eigen::Vector3d extent = maxV - minV;
        int axis = 0;
        extent.maxCoeff(&axis);
        // Degenerate box: cannot split meaningfully
        if (extent[axis] <= 1e-12) {
            BVH[nodeIdx].leaf = true;
            BVH[nodeIdx].edges = edges;
            return nodeIdx;
        }

        // 4. Compute edge centroids along selected axis
        std::vector<std::pair<double, int>> centroidEdgePairs;
        centroidEdgePairs.reserve(edges.size());
        for (int e_idx = 0; e_idx < edges.size(); e_idx++) {
            int e = edges[e_idx];
            if (e < 0 || e >= E.size() || !E[e].active) {
                continue;
            }
            int he = E[e].he;
            Eigen::Vector3d p0 = V[HE[HE[he].twin].dest].new_pos;
            Eigen::Vector3d p1 = V[HE[he].dest].new_pos;
            Eigen::Vector3d centroid = (p0 + p1) / 2.0;
            centroidEdgePairs.push_back(std::make_pair(centroid[axis], e));
        }
        // We've hit the leaf edge count. Terminate at this node
        if (centroidEdgePairs.size() <= leafEdgeCount) {
            BVH[nodeIdx].leaf = true;
            for (int i = 0; i < centroidEdgePairs.size(); i++) {
                BVH[nodeIdx].edges.push_back(centroidEdgePairs[i].second);
            }
            return nodeIdx;
        }

        // 5. Sort edges by centroid coordinate
        std::sort(centroidEdgePairs.begin(), centroidEdgePairs.end());
        // 6. Split edge list in half
        int mid = centroidEdgePairs.size() / 2;
        // Safety: Degenerate case
        if (mid <= 0 || mid >= centroidEdgePairs.size()) {
            BVH[nodeIdx].leaf = true;
            for (int i = 0; i < centroidEdgePairs.size(); i++) {
                BVH[nodeIdx].edges.push_back(centroidEdgePairs[i].second);
            }
            return nodeIdx;
        }

        // Distinguish the left and right side edges
        std::vector<int> leftEdges;
        std::vector<int> rightEdges;
        leftEdges.reserve(mid);
        rightEdges.reserve(centroidEdgePairs.size() - mid);
        for (int i = 0; i < centroidEdgePairs.size(); i++) {
            if (i < mid) {
                leftEdges.push_back(centroidEdgePairs[i].second);
            } else {
                rightEdges.push_back(centroidEdgePairs[i].second);
            }
        }

        // 7. Recurse to get children
        int leftChild = buildBVHNode(leftEdges, depth + 1);
        int rightChild = buildBVHNode(rightEdges, depth + 1);
        if (leftChild < 0 || rightChild < 0) {  // If either child encounters an error, then return as leaf
            BVH[nodeIdx].leaf = true;
            BVH[nodeIdx].edges = edges;
            BVH[nodeIdx].children.clear();
            return nodeIdx;
        }
        BVH[nodeIdx].children.push_back(leftChild);
        BVH[nodeIdx].children.push_back(rightChild);

        // Clear any non-leaf nodes' edges for storage (we will never need them)
        BVH[nodeIdx].edges.clear();
        BVH[nodeIdx].leaf = false;
        return nodeIdx;
    }

    // Squared distance from a point to an AABB
    double polynet::pointAABBDist2(const Eigen::Vector3d& p, int box) const {
        if (box < 0 || box >= BVH.size()) {
            return std::numeric_limits<double>::infinity();
        }
        const Eigen::Vector3d& bmin = BVH[box].bdyVerts.first;
        const Eigen::Vector3d& bmax = BVH[box].bdyVerts.second;

        double d2 = 0.0;
        for (int k = 0; k < 3; k++) {
            if (p[k] < bmin[k]) {
                double d = bmin[k] - p[k];
                d2 += d * d;
            } else if (p[k] > bmax[k]) {
                double d = p[k] - bmax[k];
                d2 += d * d;
            }
        }
        return d2;
    }

    // Find the closest point on the polyline to p, snapping to a vertex within snapTol
    int polynet::closestPoint(const Eigen::Vector3d& p, polyBindData& bind, bool snap, double snapTol) const {
        if (BVH.empty()) {
            return -1;
        }

        double bestDist2 = std::numeric_limits<double>::infinity();
        bind = polyBindData();

        std::priority_queue<BVHQueueEntry> q;
        q.push({0, pointAABBDist2(p, 0)});
        while (!q.empty()) {
            BVHQueueEntry curr = q.top();
            q.pop();
            if (curr.dist2 >= bestDist2) {
                break;
            }

            const AABB& box = BVH[curr.box];
            if (box.leaf) {
                for (int i = 0; i < box.edges.size(); i++) {
                    int e = box.edges[i];
                    int he = E[e].he;
                    int v0 = HE[HE[he].twin].dest;
                    int v1 = HE[he].dest;
                    const Eigen::Vector3d& p0 = V[v0].new_pos;
                    const Eigen::Vector3d& p1 = V[v1].new_pos;

                    Eigen::Vector3d vec = p1 - p0;
                    double denom = vec.squaredNorm();
                    double t = denom <= 1e-16 ? 0.0 : (p - p0).dot(vec) / denom;
                    t = std::max(0.0, std::min(1.0, t));
                    Eigen::Vector3d candidatePos = p0 + t * vec;

                    double d2 = (p - candidatePos).squaredNorm();
                    if (d2 < bestDist2) {
                        bestDist2 = d2;
                        bind.pos = candidatePos;
                        if (snap && (candidatePos - p0).squaredNorm() <= snapTol * snapTol) {
                            bind.elType = 0;
                            bind.elIdx = v0;
                            bind.t = 0.0;
                        } else if (snap && (candidatePos - p1).squaredNorm() <= snapTol * snapTol) {
                            bind.elType = 0;
                            bind.elIdx = v1;
                            bind.t = 1.0;
                        } else {
                            bind.elType = 1;
                            bind.elIdx = e;
                            bind.t = t;
                        }
                    }
                }
            } else {
                for (int i = 0; i < box.children.size(); i++) {
                    int child = box.children[i];
                    double childDist2 = pointAABBDist2(p, child);
                    if (childDist2 < bestDist2) {
                        q.push({child, childDist2});
                    }
                }
            }
        }

        if (bind.elIdx < 0) {
            return -1;
        }
        return 1;
    }
}   // namespace Polynet
