#include "dcurvenet.hpp"

#include "curvenet/curvenet.hpp"
#include "utils/utils.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <cmath>
#include <algorithm>

namespace DCurvenet {

    // Move a vertex
    int dcurvenet::moveVert(int v, Eigen::Vector3d new_pos) {
        if (v >= V.size()) {
            return -1;
        }
        V[v].new_pos = new_pos;
        return 1;
    }

    bool dcurvenet::isPositiveHalfedge(int he) const {
        return E[HE[he].edge].he == he;
    }

    // Compute the deformation gradient of a halfedge
    Eigen::Matrix3d dcurvenet::computeHEDefGrad(int he) {
        const scaledFrame& heRestFrame = HE[he].defData.restFrame;
        const scaledFrame& heNewFrame = HE[he].defData.newFrame;
        Eigen::Matrix3d F = (heNewFrame.l / heRestFrame.l) * (heNewFrame.tangent * heRestFrame.tangent.transpose()) + 
                            (heNewFrame.w / heRestFrame.w) * (heNewFrame.binormal * heRestFrame.binormal.transpose()) + 
                            (heNewFrame.h / heRestFrame.h) * (heNewFrame.normal * heRestFrame.normal.transpose());
        return F;
    }

    int dcurvenet::computeDefGradAll() {
        for (int he = 0; he < HE.size(); he++) {
            HE[he].defData.defGrad = computeHEDefGrad(he);
        }
        return 1;
    }

    // Accumulate rotation matrices, starting from an initial halfedge and tracing forward until we hit the goal vertex
    double dcurvenet::accumulateRotations(int start_he, int end_v, std::vector<Eigen::Matrix3d>& rots, std::vector<double>& lens) {
        rots.clear();
        lens.clear();
        int he_curr = start_he;
        int he_prev = start_he;
        if (start_he < 0) {
            return 0.0;
        }
        Eigen::Vector3d tan_prev = HE[he_curr].defData.newFrame.tangent;
        Eigen::Matrix3d curr_rot = Eigen::Matrix3d::Identity();
        double curr_len = 0.0;
        // Traverse halfedges until we hit the end vertex
        do {
            Eigen::Vector3d tan_curr = HE[he_curr].defData.newFrame.tangent;
            curr_rot = Utils::computeRotation(tan_prev, tan_curr) * curr_rot;   // Left multiply to accumulate rotations
            rots.push_back(curr_rot);
            curr_len += HE[he_curr].defData.newFrame.l;
            lens.push_back(curr_len);

            he_prev = he_curr;
            he_curr = HE[he_curr].next;
            tan_prev = tan_curr;
        } while ((HE[he_prev].dest != end_v) && (he_curr != -1));
        return curr_len;
    }

    // Total torsion
    double dcurvenet::computeTorsion(Eigen::Vector3d n_1, Eigen::Vector3d n_k, Eigen::Matrix3d Om_k, Eigen::Vector3d t_k) {
        Eigen::Vector3d twist = Om_k * n_1;
        double y = twist.dot(n_k.cross(t_k));
        double x = twist.dot(n_k);
        return std::atan2(y, x);
    }

}   // namespace DCurvenet