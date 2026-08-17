#include "cage.hpp"

namespace Cage {

int cage::matrixVerts(Eigen::MatrixXd& Verts) const {
    return -1;
}

int cage::closestPoint(const Eigen::Vector3d& p, int& elType, int& elIdx, Eigen::Vector3d& proj, Eigen::VectorXd& coords) const {
    return -1;
}

std::vector<std::pair<int, double>> cage::computeBasis(const int& elType, const int& elIdx, const Eigen::VectorXd& coords) const {
    return {};
}

int cage::computeColors(Eigen::MatrixXd& Colors) const {
    return -1;
}

double cage::bboxDiag() const {
    return -1;
}

}   // namespace Cage
