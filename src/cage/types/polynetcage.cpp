#include "polynetcage.hpp"
#include <Eigen/Core>
#include <vector>
#include <array>
#include <map>
#include <algorithm>
#include <random>

namespace Cage {

// Init using polynet object
polynetcage::polynetcage(Polynet::polynet* PN): PN(PN) {
    computeArclenTable();
}

// Concatenate each edge's length into a global cumulative table
void polynetcage::computeArclenTable() {
    int num_edges = PN->numEdges();
    totalArclen.assign(num_edges + 1, 0.0);
    for (int e = 0; e < num_edges; e++) {
        auto [v0, v1] = PN->edgeEndpoints(e);
        totalArclen[e + 1] = totalArclen[e] + (v1 - v0).norm();
    }
}

// Pick an edge weighted by length, then place a random point uniformly along it
int polynetcage::sampleArclengthPoint(std::mt19937& gen, Eigen::Vector3d& pos) const {
    if (totalArclen.size() < 2 || totalArclen.back() <= 1e-16) {
        return -1;
    }
    std::uniform_real_distribution<double> dist(0.0, totalArclen.back());
    double target = dist(gen);
    int e = static_cast<int>(std::upper_bound(totalArclen.begin(), totalArclen.end(), target) - totalArclen.begin()) - 1;
    e = std::clamp(e, 0, static_cast<int>(totalArclen.size()) - 2);

    auto [v0, v1] = PN->edgeEndpoints(e);
    double segLen = totalArclen[e + 1] - totalArclen[e];
    // Safety clamp to endpoints
    double alpha = segLen > 1e-16 ? std::clamp((target - totalArclen[e]) / segLen, 0.0, 1.0) : 0.0;
    pos = v0 + alpha * (v1 - v0);
    return 1;
}

// Function for retrieving verices as a matrix
int polynetcage::matrixVerts(Eigen::MatrixXd& Verts) const {
    PN->vertsAsMatrix(Verts);
    return 1;
}

// Function for querying closest point
int polynetcage::closestPoint(const Eigen::Vector3d& p, Utils::projData& hit) const {
    return PN->closestPoint(p, hit);
}

// Function for computing basis of mesh element
std::vector<std::pair<int, double>> polynetcage::computeBasis(const Utils::projData& proj) const {
    std::vector<std::pair<int, double>> bases;
    PN->evaluateBasis(proj.elType, proj.elIdx, proj.coords[0], bases);
    return bases;
}

int polynetcage::computeColors(Eigen::MatrixXd& Colors) const {
    Eigen::MatrixXd Verts;
    PN->vertsAsMatrix(Verts);
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

// Direction toward a random point on the polyline, arclength-weighted
int polynetcage::sampleDirection(const Eigen::Vector3d& q_pos, std::mt19937& gen, Eigen::Vector3d& direc) const {
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

// Cast a ray against the polyline
int polynetcage::raycast(const Eigen::Vector3d& origin, const Eigen::Vector3d& direc, std::vector<Utils::projData>& hits, double tol) const {
    return PN->raycast(origin, direc, hits, tol);
}

double polynetcage::bboxDiag() const {
    return PN->getBBoxDiag();
}

}   // namespace Cage