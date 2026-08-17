// dcurvenet_color.hpp
#pragma once

#include <Eigen/Core>

namespace Polynet {

    // Deformation
    struct vertColorData {
        // Runtime variables
        Eigen::Vector4d rgba;
    };

    struct heColorData {
        Eigen::Vector4d rgba;
    };

    struct faceDeformData {
    };
}   // namespace Polynet