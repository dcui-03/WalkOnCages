#include "meshcage.hpp"
#include <Eigen/Core>
#include <vector>
#include <array>
#include <map>

namespace Cage {

// Init using mesh object
meshcage::meshcage(Mesh::mesh* M): M(M) {

}

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
    Mesh::vertProjData proj{elType, elIdx};
    M->evaluateBasis(proj, coords, bases);
    return bases;
}

int meshcage::computeColors(Eigen::MatrixXd& Colors) {
    Eigen::MatrixXd Verts;
    M->vertsAsMatrix(Verts);
    Colors.resize(Verts.rows(), 3);
    for (int v = 0; v < Verts.rows(); v++) {
        Colors.row(v) = ((M->getVNormal(v) + Eigen::Vector3d::Ones()) * 0.5).transpose();
    }
    return 1;
}

double meshcage::bboxDiag() const {
    return M->getBBoxDiag();
}

}   // namespace Cage