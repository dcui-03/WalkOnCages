// cage.hpp
#pragma once

#include "utils/utils.hpp"
#include <Eigen/Core>
#include <vector>
#include <array>
#include <map>
#include <random>

namespace Cage {

// Generic cage object
class cage {
    public:
        // Init using pointer to some object (generic)
        virtual ~cage() = default;

        // Function for retrieving verices as a matrix
        virtual int matrixVerts(Eigen::MatrixXd& Verts) const = 0;

        // Function for querying closest point
        virtual int closestPoint(const Eigen::Vector3d& p, Utils::projData& hit) const = 0;

        // Function for computing basis
        virtual std::vector<std::pair<int, double>> computeBasis(const Utils::projData& proj) const = 0;

        // Debug colors, one row per cage vertex (same order as matrixVerts), values in [0,1]
        virtual int computeColors(Eigen::MatrixXd& Colors) const = 0;

        // Sample a random unit direction from q_pos, chosen so a raycast along it is likely (mesh) or guaranteed (curve/polyline) to hit the cage
        virtual int sampleDirection(const Eigen::Vector3d& q_pos, std::mt19937& gen, Eigen::Vector3d& direc) const = 0;

        // Cast a ray against the cage, collecting all valid intersections (nearest first)
        virtual int raycast(const Eigen::Vector3d& origin, const Eigen::Vector3d& direc, std::vector<Utils::projData>& hits, double tol) const = 0;

        // Get bounding box diagonal length
        virtual double bboxDiag() const = 0;
    protected:
        // All inherited funcs are public
    private:
        // Point to object of whatever type
};

}   // namespace Cage
