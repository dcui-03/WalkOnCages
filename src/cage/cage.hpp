// cage.hpp
#pragma once

#include <Eigen/Core>
#include <vector>
#include <array>
#include <map>

namespace Cage {

// Generic cage object
class cage {
    public:
        // Init using pointer to some object (generic)
        virtual ~cage() = default;

        // Function for retrieving verices as a matrix
        virtual int matrixVerts(Eigen::MatrixXd& Verts) const = 0;

        // Function for querying closest point
        virtual int closestPoint(const Eigen::Vector3d& p, int& elType, int& elIdx, Eigen::Vector3d& proj, Eigen::VectorXd& coords) const = 0;

        // Function for computing basis
        virtual std::vector<std::pair<int, double>> computeBasis(const int& elType, const int& elIdx, const Eigen::VectorXd& coords) const = 0;

        // Debug colors, one row per cage vertex (same order as matrixVerts), values in [0,1]
        virtual int computeColors(Eigen::MatrixXd& Colors) const = 0;

        // Get bounding box diagonal length
        virtual double bboxDiag() const = 0;
    protected:
        // All inherited funcs are public
    private:
        // Point to object of whatever type
};

}   // namespace Cage
