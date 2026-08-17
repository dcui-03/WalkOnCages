#include "meshcage.hpp"
#include "utils/wos.hpp"
#include <Eigen/Core>
#include <vector>
#include <array>
#include <map>

namespace Cage {

// Init using mesh object
meshcage::meshcage(Mesh::mesh* M): M(M) {

}

// Function for retrieving verices as a matrix
int meshcage::matrixVerts(Eigen::MatrixXd& Verts) const {
    M->vertsAsMatrix(Verts);
    return 1;
}

// Function for querying closest point
int meshcage::closestPoint(const Eigen::Vector3d& p, Utils::projData& hit) const {
    Utils::frameData bindData;
    int success = M->computeVBinding(p, bindData);
    if (success != 1) {
        return -1;
    }
    hit = bindData.proj;
    return 1;
}

// Function for computing basis of mesh element
std::vector<std::pair<int, double>> meshcage::computeBasis(const Utils::projData& proj) const {
    std::vector<std::pair<int, double>> bases;
    M->evaluateBasis(proj, proj.coords, bases);
    return bases;
}

int meshcage::computeColors(Eigen::MatrixXd& Colors) const {
    Eigen::MatrixXd Verts;
    M->vertsAsMatrix(Verts);
    Colors.resize(Verts.rows(), 3);
    for (int v = 0; v < Verts.rows(); v++) {
        Colors.row(v) = ((M->getVNormal(v) + Eigen::Vector3d::Ones()) * 0.5).transpose();
    }
    return 1;
}

// Uniform random direction on the sphere, independent of q_pos
int meshcage::sampleDirection(const Eigen::Vector3d& q_pos, std::mt19937& gen, Eigen::Vector3d& direc) const {
    return WoS::generateNewDirection(direc, gen);
}

// Cast a ray against the mesh
int meshcage::raycast(const Eigen::Vector3d& origin, const Eigen::Vector3d& direc, std::vector<Utils::projData>& hits, double tol) const {
    return M->raycast(origin, direc, hits);
}

double meshcage::bboxDiag() const {
    return M->getBBoxDiag();
}

}   // namespace Cage