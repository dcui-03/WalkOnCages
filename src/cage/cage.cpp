#include "cage.hpp"

namespace Cage {

int cage::matrixVerts(Eigen::MatrixXd& Verts) const {
    return -1;
}

int cage::closestPoint(const Eigen::Vector3d& p, Utils::projData& hit) const {
    return -1;
}

std::vector<std::pair<int, double>> cage::computeBasis(const Utils::projData& proj) const {
    return {};
}

int cage::computeColors(Eigen::MatrixXd& Colors) const {
    return -1;
}

int cage::sampleDirection(const Eigen::Vector3d& q_pos, std::mt19937& gen, Eigen::Vector3d& direc) const {
    return -1;
}

int cage::raycast(const Eigen::Vector3d& origin, const Eigen::Vector3d& direc, std::vector<Utils::projData>& hits, double tol) const {
    return -1;
}

double cage::bboxDiag() const {
    return -1;
}

}   // namespace Cage
