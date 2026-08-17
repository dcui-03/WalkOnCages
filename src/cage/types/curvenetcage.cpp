#include "curvenetcage.hpp"
#include <Eigen/Core>
#include <vector>
#include <array>
#include <map>
#include <algorithm>
#include <random>

namespace Cage {

// Init using curvenet object
curvenetcage::curvenetcage(Curvenet::curvenet* CN): CN(CN) {
    // Create dCN for closest point BVH
    dCN = Polynet::dcurvenet(CN, nullptr);
    computeArclenTable();
}

// Concatenate dCN's edge lengths into a global cumulative table so we can sample off it
void curvenetcage::computeArclenTable() {
    int num_edges = dCN.numEdges();
    totalArclen.assign(num_edges + 1, 0.0);
    for (int e = 0; e < num_edges; e++) {
        auto [v0, v1] = dCN.edgeEndpoints(e);
        totalArclen[e + 1] = totalArclen[e] + (v1 - v0).norm();
    }
}

// Pick a dCN edge weighted by length, then place a random point uniformly along it to get our random sample
int curvenetcage::sampleArclengthPoint(std::mt19937& gen, Eigen::Vector3d& pos) const {
    if (totalArclen.size() < 2 || totalArclen.back() <= 1e-16) {
        return -1;
    }
    std::uniform_real_distribution<double> dist(0.0, totalArclen.back());
    double target = dist(gen);
    int e = static_cast<int>(std::upper_bound(totalArclen.begin(), totalArclen.end(), target) - totalArclen.begin()) - 1;
    e = std::clamp(e, 0, static_cast<int>(totalArclen.size()) - 2);

    auto [v0, v1] = dCN.edgeEndpoints(e);
    double segLen = totalArclen[e + 1] - totalArclen[e];
    double alpha = segLen > 1e-16 ? std::clamp((target - totalArclen[e]) / segLen, 0.0, 1.0) : 0.0;
    pos = v0 + alpha * (v1 - v0);
    return 1;
}

// Function for retrieving verices as a matrix
int curvenetcage::matrixVerts(Eigen::MatrixXd& Verts) const {
    return CN->CTasMatrix(Verts);
}

// Function for querying closest point
int curvenetcage::closestPoint(const Eigen::Vector3d& p, Utils::projData& hit) const {
    int success = CN->closestPoint(p, &dCN, hit);
    if (success != 1) {
        return -1;
    }
    hit.elType = 0;
    return 1;
}

// Function for computing basis of mesh element
std::vector<std::pair<int, double>> curvenetcage::computeBasis(const Utils::projData& proj) const {
    std::vector<std::pair<int, double>> bases;
    CN->evaluateBasis(proj.elIdx, proj.coords[0], bases);
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

// Direction toward a random point on the curve network, arclength-weighted
int curvenetcage::sampleDirection(const Eigen::Vector3d& q_pos, std::mt19937& gen, Eigen::Vector3d& direc) const {
    Eigen::Vector3d pos;
    if (sampleArclengthPoint(gen, pos) != 1) {
        return -1;
    }
    direc = pos - q_pos;
    double len = direc.norm();
    if (len <= 1e-12) {
        return -1;
    }
    direc /= len;
    return 1;
}

// Cast a ray against the curve network, refined onto the true splines
int curvenetcage::raycast(const Eigen::Vector3d& origin, const Eigen::Vector3d& direc, std::vector<Utils::projData>& hits, double tol) const {
    return CN->raycast(origin, direc, &dCN, hits, tol);
}

double curvenetcage::bboxDiag() const {
    return CN->getBBoxDiag();
}

}   // namespace Cage