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

        // TODO: Function for retrieving verices as a matrix
        virtual int matrixVerts(Eigen::MatrixXd& Verts);

        // TODO: Function for querying closest point
        virtual int closestPoint(int& elType, int& elIdx, Eigen::Vector3d& proj, Eigen::VectorXd& coords);

        // TODO: Function for computing basis
        virtual std::vector<std::pair<int, double>> computeBasis(const int& elType, const int& elIdx, const Eigen::VectorXd& coords);
    protected:
        // All inherited funcs are public
    private:
        // Point to object of whatever type
};

}   // namespace Cage