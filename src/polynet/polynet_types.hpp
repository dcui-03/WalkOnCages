// polynet_types.hpp
#pragma once

#include <Eigen/Core>
#include <vector>
#include <utility>

// File with basic structs used by the polyline network class

namespace Polynet {

    // Projection data onto the corresponding mesh for each vertex
    struct projData {
        int elType = -1;
        int elIdx = -1;
        Eigen::VectorXd coords;
        Eigen::Vector3d projVec = Eigen::Vector3d::Zero();
    };

    // Bounding Volume Hierarchy over edges
    struct AABB {
        int depth = 0;
        bool leaf = false;
        std::vector<int> children;
        std::vector<int> edges;                                 // Edges in the bvh
        std::pair<Eigen::Vector3d, Eigen::Vector3d> bdyVerts;   // Bounding box verts
    };

    // Used in closest-point AABB priority queue
    struct BVHQueueEntry {
        int box;
        double dist2;

        bool operator<(const BVHQueueEntry& other) const {
            return dist2 > other.dist2; // reversed for std::priority_queue min-heap behavior
        }
    };

    // Result of a closest-point query against the polyline
    struct polyBindData {
        int elType = -1;    // 0 = vertex, 1 = edge
        int elIdx = -1;
        double t = -1.0;    // Local edge parameter (0/1 if snapped from a vertex)
        Eigen::Vector3d pos;
    };

    struct Vert {
        Eigen::Vector3d pos;
        Eigen::Vector3d n = Eigen::Vector3d::Zero();    // Init to zero, since most vertices will not receive an initial normal
        std::vector<int> adjHE;    // For high-valence vertices, stores CCW outgoing HE's. For all others, stores a single outgoing halfedge
        bool active = true;     // For safety, say if the component is active (ignore for now)

        // Vertex type by valence (1 = endpoint, 2 = interior, 3 = intersection), set by assignVertType
        int vType = -1;

        // Robustness for future work
        projData proj;

        // Weights
        bool fixed_w = false;
        double w = 1.0;

        // Runtime variables
        Eigen::Vector3d new_pos;
    };

    // NOTE: A halfedge's next/prev can be -1 if this is the end of a curve
    // BUT if entering a high valence vertex, its next should be the CCW outgoing halfedge (and vice versa for prev)
    struct HalfEdge {
        int dest = -1;
        int twin = -1;
        int next = -1;
        int prev = -1;
        int edge = -1;
        bool active = true;     // For safety, say if the component is active (ignore for now)
    };

    // Edge in a curve
    struct Edge {
        int he = -1;
        int curve = -1; // Index to curve
        bool active = true;
    };

    // Polyline curve (a chain of edges between two non-interior vertices, or a closed loop)
    struct Curve {
        int he_start = -1;        // Starting outgoing halfedge
        int he_end = -1;        // Ending outgoing halfedge
        int start = -1;    // One of the endpoint vertices
        int end = -1;      // The other endpoint vertex
        bool active = true;     // For safety, say if the component is active (ignore for now)
    };
}   // namespace Polynet
