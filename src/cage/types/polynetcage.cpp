#include "polynetcage.hpp"
#include <Eigen/Core>
#include <vector>
#include <array>
#include <map>

namespace Cage {

// Init using polynet object
polynetcage::polynetcage(Polynet::polynet* PN) (PN: PN);

// Function for retrieving verices as a matrix
int polynetcage::matrixVerts(Eigen::MatrixXd& Verts) {
    PN->vertsAsMatrix(Verts);
    return 1;
}

// Function for querying closest point
int polynetcage::closestPoint(const Eigen::Vector3d& p, int& elType, int& elIdx, Eigen::Vector3d& proj, Eigen::VectorXd& coords) {
    Polynet::polyBindData bindData;
    int success = CN->closestPoint(p, bindData);
    if (success != 1) {
        return -1;
    }
    elType = bindData.elType;
    elIdx = bindData.elIdx;
    proj = bindData.pos;
    coords.resize(1);
    coords[0] = bindData.t;
    return 1;
}

// Function for computing basis of mesh element
std::vector<std::pair<int, double>> polynetcage::computeBasis(const int& elType, const int& elIdx, const Eigen::VectorXd& coords) {
    std::vector<std::pair<int, double>> bases;
    PN->evaluateBasis(elType, elIdx, coords[0], bases);
    return bases;
}

}   // namespace Cage