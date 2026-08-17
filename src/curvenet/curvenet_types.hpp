// curvenet_types.hpp
#pragma once

#include <Eigen/Core>
#include <vector>

// File with basic structs used by mesh class

namespace Curvenet {

    // Result of a closest-point query against the curve network
    struct cnBindData {
        int s = -1;                     // Spline index
        double t = -1.0;                // Parameter value on that spline
        Eigen::Vector3d pos;            // Point on the spline at t
    };

    // Projection data onto the corresponding mesh for tangents and handles
    struct projData {
        int elType = -1;
        int elIdx = -1;
        Eigen::VectorXd coords;
        Eigen::Vector3d projVec = Eigen::Vector3d::Zero();
        // Local coordinate frame at projection point
        Eigen::Matrix3d projFrame = Eigen::Matrix3d::Identity();
    };

    // Control vertices
    // NOTE: Each control stores an (ordered) list of outgoing halfedges
    //       This is in case a spline starts and ends at the same control
    struct Control {
        Eigen::Vector3d pos;
        Eigen::Vector3d n = Eigen::Vector3d::Zero();
        std::vector<int> adjHE;        // Outgoing Halfedge list
        bool active = true;     // For safety, say if the component is active (ignore for now)
        bool sorted = false;    // Safety flag. True when outgoing halfedges are sorted
        int cType = 0;      // Control point type (1 = anchor, 2 = loop, 3 = intersection)

        // Weight, if specified
        bool fixed_w = true;
        double w = 1.0;

        // Runtime info
        Eigen::Vector3d new_pos;
    };

    // NOTE: Halfedge iteration indices (next, prev) prioritize easy iteration over their associated controls
    //       As a result, they may not be ideal for standard halfedge traversal tasks.
    struct HalfEdge {
        int origin = -1;  // To control
        int s = -1;     // associated spline
        // Iterators: Note that if we are at an endpoint, leave entry as -1
        int twin = -1;
        int next = -1;
        int prev = -1;
        Eigen::Vector3d rest_tan;   // Tangent vector (defined in global coordinates, NOT relative to control)

        bool active = true;     // For safety, say if the component is active (ignore for now)
        
        // Runtime info
        Eigen::Vector3d tan;
    };

    // Splines for easy iteration
    struct CubicSpline {
        int he = -1;        // "First" halfedge describing the canonical direction
        int curve = -1;     // Curve index the spline belongs to
        int num_samples = -1;   // Number of samples to take during discretization
        bool active = true;     // For safety, say if the component is active (ignore for now)
    };

    // Curves are chains of Splines
    struct Curve {
        // Just store the chain of spline indices
        std::vector<int> splines;
        bool active = true;
    };
}   // namespace Curvenet