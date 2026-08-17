#include "curvenetcage.hpp"
#include <Eigen/Core>
#include <vector>
#include <array>
#include <map>

namespace Cage {

// Init using curvenet object
curvenetcage::curvenetcage(Curvenet::curvenet* CN) (CN: CN) {
    // Create dCN for closest point BVH
    dCN = Polynet::dcurvenet(CN, nullptr);
}

// Function for retrieving verices as a matrix
int curvnetcage::matrixVerts(Eigen::MatrixXd& Verts) {
    CN->vertsAsMatrix(Verts);
    return 1;
}

// Function for querying closest point
int curvenetcage::closestPoint(const Eigen::Vector3d& p, int& elType, int& elIdx, Eigen::Vector3d& proj, Eigen::VectorXd& coords) {
    Curvenet::cnBindData bindData;
    int success = CN->closestPoint(p, &dCN, bindData);
    if (success != 1) {
        return -1;
    }
    elType = 0;
    elIdx = bindData.s;
    proj = bindData.proj;
    coords.resize(1);
    coords[0] = bindData.t;
    return 1;
}

// Function for computing basis of mesh element
std::vector<std::pair<int, double>> curvenetcage::computeBasis(const int& elType, const int& elIdx, const Eigen::VectorXd& coords) {
    std::vector<std::pair<int, double>> bases;
    CN->evaluateBasis(elIdx, coords[0], bases);
    return bases;
}

}   // namespace Cage