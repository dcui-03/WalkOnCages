#include "meshcage.hpp"
#include <Eigen/Core>
#include <vector>
#include <array>
#include <map>

namespace Cage {

// Init using mesh object
meshcage::meshcage(Mesh::mesh* M) (M: M);

// Function for retrieving verices as a matrix
int meshcage::matrixVerts(Eigen::MatrixXd& Verts) {
    M->vertsAsMatrix(Verts);
    return 1;
}

// Function for querying closest point
int meshcage::closestPoint(const Eigen::Vector3d& p, int& elType, int& elIdx, Eigen::Vector3d& proj, Eigen::VectorXd& coords) {
    Mesh::meshBindData bindData;
    int success = M->computeVBinding(p, bindData);
    if (success != 1) {
        return -1;
    }
    elType = bindData.elType;
    elIdx = bindData.elIdx;
    proj = bindData.proj;
    coords = bindData.coords;
    return 1;
}

// Function for computing basis of mesh element
std::vector<std::pair<int, double>> meshcage::computeBasis(const int& elType, const int& elIdx, const Eigen::VectorXd& coords) {
    std::vector<std::pair<int, double>> bases;
    M->evaluateBasis(elType, elIdx, coords, bases);
    return bases;
}

}   // namespace Cage