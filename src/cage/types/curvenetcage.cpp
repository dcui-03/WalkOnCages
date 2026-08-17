#include "curvenetcage.hpp"
#include <Eigen/Core>
#include <vector>
#include <array>
#include <map>

namespace Cage {

// Init using curvenet object
curvenetcage::curvenetcage(Curvenet::curvenet* CN): CN(CN) {
    // Create dCN for closest point BVH
    dCN = Polynet::dcurvenet(CN, nullptr);
}

// Function for retrieving verices as a matrix
int curvenetcage::matrixVerts(Eigen::MatrixXd& Verts) const {
    return CN->CTasMatrix(Verts);
}

// Function for querying closest point
int curvenetcage::closestPoint(const Eigen::Vector3d& p, int& elType, int& elIdx, Eigen::Vector3d& proj, Eigen::VectorXd& coords) const {
    Curvenet::cnBindData bindData;
    int success = CN->closestPoint(p, &dCN, bindData);
    if (success != 1) {
        return -1;
    }
    elType = 0;
    elIdx = bindData.s;
    proj = bindData.pos;
    coords.resize(1);
    coords[0] = bindData.t;
    return 1;
}

// Function for computing basis of mesh element
std::vector<std::pair<int, double>> curvenetcage::computeBasis(const int& elType, const int& elIdx, const Eigen::VectorXd& coords) const {
    std::vector<std::pair<int, double>> bases;
    CN->evaluateBasis(elIdx, coords[0], bases);
    return bases;
}

int curvenetcage::computeColors(Eigen::MatrixXd& Colors) const {
    Eigen::MatrixXd Verts;
    CN->CTasMatrix(Verts);
    Eigen::Vector3d centroid = Verts.colwise().mean();
    Colors.resize(Verts.rows(), 3);
    for (int v = 0; v < Verts.rows(); v++) {
        Eigen::Vector3d dir = Verts.row(v).transpose() - centroid;
        if (dir.squaredNorm() > 1e-20) {
            dir.normalize();
        }
        Colors.row(v) = ((dir + Eigen::Vector3d::Ones()) * 0.5).transpose();
    }
    return 1;
}

double curvenetcage::bboxDiag() const {
    return CN->getBBoxDiag();
}

}   // namespace Cage