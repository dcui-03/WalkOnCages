// dcurvenet_types.hpp
#pragma once

#include "dcurvenet_def.hpp"
#include <Eigen/Core>
#include <vector>
#include <utility>

// File with basic structs used by discrete curvenet class

namespace DCurvenet {

    // Projection data onto the corresponding mesh for each vertex
    struct projData {
        int elType = -1;
        int elIdx = -1;
        Eigen::VectorXd coords;
        Eigen::Vector3d projVec = Eigen::Vector3d::Zero();
    };

    struct Vert {
        Eigen::Vector3d pos;
        Eigen::Vector3d n = Eigen::Vector3d::Zero();    // Init to zero, since most vertices will not receive an initial normal
        std::vector<int> adjHE;    // For control vertices, stores CCW outgoing HE's. For all others, stores a single outgoing halfedge
        bool active = true;     // For safety, say if the component is active (ignore for now)
        
        // cn_idx is redundant due to initialization, but better to be safe
        int cn_idx = -1;    // curvenet index if coincident with a control vertex
        int cn_type = -1;   // curvenet vertex type if coincident with a control vertex

        // Robustness for future work
        projData proj;

        // Weights
        double w = 1.0;

        // Runtime variables
        Eigen::Vector3d new_pos;
    };

    // NOTE: A halfedge's next/prev can be -1 if this is the end of a spline
    // BUT if entering a high valence vertex, its next should be the CCW outgoing halfedge (and vice versa for prev)
    struct HalfEdge {
        int dest = -1;
        int twin = -1;
        int next = -1;
        int prev = -1;
        int edge = -1;
        bool active = true;     // For safety, say if the component is active (ignore for now)

        // Deformation data, including the scaled frames (old and new)
        heDeformData defData;
    };

    // Edge in a curve
    struct Edge {
        int he = -1;
        int curve = -1; // Index to curve
        bool active = true;
    };

    // Discrete curve
    struct Curve {
        int he_start = -1;        // Starting outgoing halfedge
        int he_end = -1;        // Ending outgoing halfedge
        int start = -1;    // One of the endpoint vertices
        int end = -1;      // The other endpoint vertex
        bool active = true;     // For safety, say if the component is active (ignore for now)

        // Other info
        int cn_idx = -1;    // Curvenet index of parent curve
    };
}   // namespace DCurvenet