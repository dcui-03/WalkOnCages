// cagedeformer_types.hpp
#pragma once

#include <Eigen/Core>
#include <vector>
#include <utility>

// File with basic structs used by mesh class

namespace CageDeformer {

    struct Sample {
        // Basis of relevant vertices to this sample evaluated at its location
        Eigen::Vector3d pos;
        std::vector<std::pair<int, double>> bases;
        double weight;
    };

    struct QueryVert {
        Eigen::Vector3d pos;
        std::vector<Sample> samples;
    };

}   // namespace CageDeformer
