// dcurvenet_def.hpp
#pragma once

#include <Eigen/Core>
#include <vector>
#include <utility>

namespace DCurvenet {

    // Deformation
    struct vertDeformData {
        // Runtime variables
        Eigen::Vector3d new_pos;
    };

    struct scaledFrame {
        Eigen::Vector3d tangent;
        Eigen::Vector3d binormal;
        Eigen::Vector3d normal;
        double l, w, h; // length, width, and height
    };

    struct heDeformData {
        // Note, we can save these each separately for easy access.
        // Paper provides an easy method for computing def grad using components rather than matrices
        scaledFrame restFrame;
        // Runtime info: Altered frame needed for def grad computation
        scaledFrame newFrame;
        Eigen::Matrix3d defGrad;
    };

    struct curveDeformData {
        // Corner normals
        std::pair<Eigen::Vector3d, Eigen::Vector3d> N_pos;  // first is start, second is end
        std::pair<Eigen::Vector3d, Eigen::Vector3d> N_neg;  // first is start, second is end
        // Corner widths
        std::pair<double, double> W_pos;    // first is start, second is end
        std::pair<double, double> W_neg;    // first is start, second is end
    };
}   // namespace DCurvenet