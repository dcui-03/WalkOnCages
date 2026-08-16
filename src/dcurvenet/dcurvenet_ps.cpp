#include "dcurvenet.hpp"

#include "utils/utils.hpp"
#include <Eigen/Core>
#include <vector>
#include <utility>
#include <glm/glm.hpp>
#include <glm/vec3.hpp>
#include <cmath>
#include <iostream>

namespace DCurvenet {

int dcurvenet::polyscopeFormat(Eigen::MatrixXd& Verts, 
                            std::vector<std::array<int, 2>>& Edges, 
                            std::vector<glm::vec3>& posEdgeTangents,
                            std::vector<glm::vec3>& posEdgeBinormals,
                            std::vector<glm::vec3>& posEdgeNormals,
                            std::vector<glm::vec3>& negEdgeTangents,
                            std::vector<glm::vec3>& negEdgeBinormals,
                            std::vector<glm::vec3>& negEdgeNormals, 
                            std::vector<double>& weights) const {
    Verts.resize(V.size(), 3);
    weights.resize(V.size());
    Edges.resize(E.size());
    posEdgeTangents.resize(E.size());
    posEdgeBinormals.resize(E.size());
    posEdgeNormals.resize(E.size());
    negEdgeTangents.resize(E.size());
    negEdgeBinormals.resize(E.size());
    negEdgeNormals.resize(E.size());

    for (int v = 0; v < V.size(); v++) {
        Verts.row(v) = V[v].new_pos.transpose();
        weights[v] = V[v].w;
    }

    for (int e = 0; e < E.size(); e++) {
        int he_pos = E[e].he;
        int he_neg = HE[he_pos].twin;
        Edges[e] = std::array<int, 2>({HE[he_neg].dest, HE[he_pos].dest});
        // Positive
        if (std::isfinite(HE[he_pos].defData.newFrame.tangent.norm()) && std::isfinite(HE[he_pos].defData.newFrame.l)) {
            posEdgeTangents[e] = Utils::eigenToGLM(HE[he_pos].defData.newFrame.l * HE[he_pos].defData.newFrame.tangent);
        } else {
            if (!std::isfinite(HE[he_pos].defData.newFrame.tangent.norm())) {
                std::cout << "Bad Tangent Vector on HE " << he_pos << std::endl;
            }
            if (!std::isfinite(HE[he_pos].defData.newFrame.l)) {
                std::cout << "Bad Length on HE " << he_pos << std::endl;
            }
            posEdgeTangents[e] = Utils::eigenToGLM(Eigen::Vector3d::Zero());
        }
        if (std::isfinite(HE[he_pos].defData.newFrame.binormal.norm()) && std::isfinite(HE[he_pos].defData.newFrame.w)) {
            posEdgeBinormals[e] = Utils::eigenToGLM(HE[he_pos].defData.newFrame.w * HE[he_pos].defData.newFrame.binormal);
        } else {
            if (!std::isfinite(HE[he_pos].defData.newFrame.binormal.norm())) {
                std::cout << "Bad Binormal Vector on HE " << he_pos << std::endl;
            }
            if (!std::isfinite(HE[he_pos].defData.newFrame.w)) {
                std::cout << "Bad Width on HE " << he_pos << std::endl;
            }
            posEdgeBinormals[e] = Utils::eigenToGLM(Eigen::Vector3d::Zero());
        }
        
        if (std::isfinite(HE[he_pos].defData.newFrame.normal.norm()) && std::isfinite(HE[he_pos].defData.newFrame.h)) {
            posEdgeNormals[e] = Utils::eigenToGLM(HE[he_pos].defData.newFrame.h * HE[he_pos].defData.newFrame.normal);
        } else {
            if (!std::isfinite(HE[he_pos].defData.newFrame.normal.norm())) {
                std::cout << "Bad Normal Vector on HE " << he_pos << std::endl;
            }
            if (!std::isfinite(HE[he_pos].defData.newFrame.h)) {
                std::cout << "Bad Height on HE " << he_pos << std::endl;
            }
            posEdgeNormals[e] = Utils::eigenToGLM(Eigen::Vector3d::Zero());
        }


        // Negative
        if (std::isfinite(HE[he_neg].defData.newFrame.tangent.norm()) && std::isfinite(HE[he_neg].defData.newFrame.l)) {
            negEdgeTangents[e] = Utils::eigenToGLM(HE[he_neg].defData.newFrame.l * HE[he_neg].defData.newFrame.tangent);
        } else {
            if (!std::isfinite(HE[he_neg].defData.newFrame.binormal.norm())) {
                std::cout << "Bad Tangent Vector on HE " << he_neg << std::endl;
            }
            if (!std::isfinite(HE[he_neg].defData.newFrame.w)) {
                std::cout << "Bad Length on HE " << he_neg << std::endl;
            }
            negEdgeTangents[e] = Utils::eigenToGLM(Eigen::Vector3d::Zero());
        }
        
        if (std::isfinite(HE[he_neg].defData.newFrame.binormal.norm()) && std::isfinite(HE[he_neg].defData.newFrame.w)) {
            negEdgeBinormals[e] = Utils::eigenToGLM(HE[he_neg].defData.newFrame.w * HE[he_neg].defData.newFrame.binormal);
        } else {
            if (!std::isfinite(HE[he_neg].defData.newFrame.binormal.norm())) {
                std::cout << "Bad Binormal Vector on HE " << he_neg << std::endl;
            }
            if (!std::isfinite(HE[he_neg].defData.newFrame.w)) {
                std::cout << "Bad Width on HE " << he_neg << std::endl;
            }
            negEdgeBinormals[e] = Utils::eigenToGLM(Eigen::Vector3d::Zero());
        }
        
        if (std::isfinite(HE[he_neg].defData.newFrame.normal.norm()) && std::isfinite(HE[he_neg].defData.newFrame.h)) {
            negEdgeNormals[e] = Utils::eigenToGLM(HE[he_neg].defData.newFrame.h * HE[he_neg].defData.newFrame.normal);
        } else {
            if (!std::isfinite(HE[he_neg].defData.newFrame.normal.norm())) {
                std::cout << "Bad Normal Vector on HE " << he_neg << std::endl;
            }
            if (std::isfinite(HE[he_neg].defData.newFrame.h)) {
                std::cout << "Bad Height on HE " << he_neg << std::endl;
            }
            negEdgeNormals[e] = Utils::eigenToGLM(Eigen::Vector3d::Zero());
        }
    }
    return 1;
}


}   // namespace DCurvenet