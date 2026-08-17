#include "cagedeformer.hpp"

#include "cage/cage.hpp"
#include "utils/utils.hpp"
#include <Eigen/Core>
#include <vector>
#include <glm/vec3.hpp>

// Polyscope viewer-formatting reader for debug colors.

namespace CageDeformer {

int cagedeformer::colorsPolyscopeFormat(std::vector<glm::vec3>& colors) {
    if (def_cage == nullptr) {
        return -1;
    }
    Eigen::MatrixXd sourceColors;
    if (def_cage->computeColors(sourceColors) != 1) {
        return -1;
    }
    Eigen::MatrixXd resultColors;
    if (applyColor(sourceColors, resultColors) != 1) {
        return -1;
    }
    colors.resize(resultColors.rows());
    for (int q = 0; q < resultColors.rows(); q++) {
        colors[q] = Utils::eigenToGLM(resultColors.row(q).transpose());
    }
    return 1;
}

int cagedeformer::cageColorsPolyscopeFormat(std::vector<glm::vec3>& colors) {
    if (def_cage == nullptr) {
        return -1;
    }
    Eigen::MatrixXd sourceColors;
    if (def_cage->computeColors(sourceColors) != 1) {
        return -1;
    }
    colors.resize(sourceColors.rows());
    for (int v = 0; v < sourceColors.rows(); v++) {
        colors[v] = Utils::eigenToGLM(sourceColors.row(v).transpose());
    }
    return 1;
}

int cagedeformer::colorGradientsPolyscopeFormat(std::vector<glm::vec3>& gradR, std::vector<glm::vec3>& gradG, std::vector<glm::vec3>& gradB) {
    if (def_cage == nullptr) {
        return -1;
    }
    Eigen::MatrixXd sourceColors;
    if (def_cage->computeColors(sourceColors) != 1) {
        return -1;
    }
    Eigen::MatrixXd GradX, GradY, GradZ;
    if (applyColorGradient(sourceColors, GradX, GradY, GradZ) != 1) {
        return -1;
    }
    int n = GradX.rows();
    gradR.resize(n);
    gradG.resize(n);
    gradB.resize(n);
    for (int q = 0; q < n; q++) {
        gradR[q] = Utils::eigenToGLM(Eigen::Vector3d(GradX(q, 0), GradY(q, 0), GradZ(q, 0)));
        gradG[q] = Utils::eigenToGLM(Eigen::Vector3d(GradX(q, 1), GradY(q, 1), GradZ(q, 1)));
        gradB[q] = Utils::eigenToGLM(Eigen::Vector3d(GradX(q, 2), GradY(q, 2), GradZ(q, 2)));
    }
    return 1;
}

}   // namespace CageDeformer
