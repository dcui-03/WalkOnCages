// mesh_types.hpp
#pragma once

#include "utils/utils.hpp"
#include <Eigen/Core>
#include <vector>

// File with basic structs used by mesh class

namespace Mesh {
    // Bounding Volume Hierarchy for mesh
    struct AABB {
        int depth = 0;
        bool leaf = false;
        std::vector<int> children;
        std::vector<int> faces;                                 // Faces in the bvh
        std::pair<Eigen::Vector3d, Eigen::Vector3d> bdyVerts;   // Bounding box verts
    };

    // Used in projection AABB priority queue
    struct BVHQueueEntry {
        int box;
        double dist2;

        bool operator<(const BVHQueueEntry& other) const {
            return dist2 > other.dist2; // reversed for std::priority_queue min-heap behavior
        }
    };

    struct Vert {
        Eigen::Vector3d pos;
        Eigen::Vector3d n = Eigen::Vector3d::Zero();
        double vArea = 0.0;   // barycentric dual area. To be computed only when necessary
        int he = -1;   // one outgoing halfedge, or -1 if isolated
        bool active = true;     // For safety, say if the component is active (ignore for now)
    };

    struct HalfEdge {
        int dest = -1;

        int twin = -1;
        int next = -1;
        int prev = -1;

        int edge = -1;           // Index to the edge
        int face = -1;           // Index to the face

        bool active = true;     // For safety, say if the component is active (ignore for now)
        bool boundary = false;  // This is necessary for cutmesh face reinitialization

        // Other indices for cut mesh
        int dCN_idx = -1;   // -1 if not connected, dCN HE index otherwise
    };

    // Only store one halfedge for each edge
    // No real need to store the normal
    struct Edge {
        int he = -1;
        Eigen::Vector3d n = Eigen::Vector3d::Zero();
        bool active = true;     // For safety, say if the component is active (ignore for now)
    };

    struct Face {
        int he = -1;
        Eigen::Vector3d n = Eigen::Vector3d::Zero();  // Face normal
        double fArea = 0.0;   // Face area. To be computed only when necessary
        bool active = true;     // For safety, say if the component is active (ignore for now)
    };
}   // namespace Mesh