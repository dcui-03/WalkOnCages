#include "utils.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <igl/point_mesh_squared_distance.h>
#include <glm/vec3.hpp>
#include <limits>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <iostream>

namespace Utils {

// HELPERS FOR CONVERSION/COPYING

// GLM::vec3 to Eigen::Vector3d converter
Eigen::Vector3d glmToEigen(const glm::vec3 input) {
    Eigen::Vector3d output;
    output(0) = static_cast<double>(input.x);
    output(1) = static_cast<double>(input.y);
    output(2) = static_cast<double>(input.z);
    return output;
}

// Eigen::Vector3d to GLM::vec3 converter
glm::vec3 eigenToGLM(const Eigen::Vector3d input) {
    glm::vec3 output;
    output.x = static_cast<float>(input(0));
    output.y = static_cast<float>(input(1));
    output.z = static_cast<float>(input(2));
    return output;
}

// Convert an eigen matrix with 3 columns to a std::vector
void EigM3toStdV(const Eigen::MatrixXd& mat, std::vector<Eigen::Vector3d>& vec) {
    // Assert that matrix columns = 3
    if (mat.cols() != 3) {
        return;
    }
    vec.resize(mat.rows());
    for (int v = 0; v < mat.rows(); v++) {
        vec[v] = mat.row(v).transpose();
    }
    return;
}

// Entire mesh conversion routine Eigen to GLM
void meshConversionEigentoGLM(const std::vector<Eigen::Vector3d>& Eig, std::vector<glm::vec3>& GLM) {
    GLM.clear();
    GLM.resize(Eig.size());
    for (int v = 0; v < Eig.size(); v++) {
        GLM[v] = eigenToGLM(Eig[v]);
    }
    return;
}

// Entire mesh conversion routine GLM to Eigen
void meshConversionGLMtoEigen(std::vector<Eigen::Vector3d>& Eig, const std::vector<glm::vec3>& GLM) {
    Eig.clear();
    Eig.resize(GLM.size());
    for (int v = 0; v < GLM.size(); v++) {
        Eig[v] = glmToEigen(GLM[v]);
    }
    return;
}

// Copy positions and connectivity into a copied container
void copyPositions(const std::vector<Eigen::Vector3d>& V_old, std::vector<Eigen::Vector3d>& V_new) {
    V_new.clear();
    V_new.resize(V_old.size());
    for (int v = 0; v < V_old.size(); v++) {
        Eigen::Vector3d new_v = {V_old[v](0), V_old[v](1), V_old[v](2)};
        V_new[v] = new_v;
    }
    return;
}

void copyConnectivity(const std::vector<std::vector<int>>& T_old, std::vector<std::vector<int>>& T_new) {
    T_new.clear();
    T_new.resize(T_old.size());
    for (int f = 0; f < T_old.size(); f++) {
        std::vector<int> f_idxs;
        for (int v = 0; v < T_old[f].size(); v++) {
            f_idxs.push_back(T_old[f][v]);
        }
        T_new[f] = f_idxs;
    }
    return;
}

// SORTING
void doubleListIdxSort(std::vector<double>& ref_List, std::vector<int>& idx_List) {
    if (ref_List.size() != idx_List.size()) {
        return;
    }

    const int n = ref_List.size();

    for (int i = 0; i < n - 1; ++i) {
        bool swapped = false;

        for (int j = 0; j < n - i - 1; ++j) {
            if (ref_List[j] > ref_List[j + 1]) {
                std::swap(ref_List[j], ref_List[j + 1]);
                std::swap(idx_List[j], idx_List[j + 1]);
                swapped = true;
            }
        }

        if (!swapped) {
            break;
        }
    }
}

// Insert an integer entry in a list between two specified values
bool insertIdxBetweenPair(std::vector<int>& idxList, int a, int b, int new_idx) {
    for (int i = 0; i < idxList.size(); ++i) {
        int j = (i + 1) % idxList.size();
        if (idxList[i] == a && idxList[j] == b) {
            idxList.insert(idxList.begin() + j, new_idx);
            return true;
        }
    }
    return false;
}

// Flattens an Eigen::Matrix3d into a 9x1 row vector
// NOTE: Does so column-wise!
Eigen::VectorXd flattenMatrix3d(const Eigen::Matrix3d& F) {
    return F.reshaped();
}

// Compresses a 9x1 Eigen::VectorXd into an Eigen::Matrix3d
// Assumes column-wise storage
Eigen::Matrix3d compressVector9d(const Eigen::VectorXd& f) {
    assert(f.size() == 9);
    return Eigen::Map<const Eigen::Matrix3d>(f.data());
}

// MESH HELPERS
std::pair<int, int> undirectedKey(int a, int b) {
    return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
}

int edgeDirRelativeToKey(int a, int b) {
    return (a < b) ? +1 : -1;
}

bool orientFacesConsistently(std::vector<std::vector<int>>& F_List) {
    struct FaceEdgeUse {
        int face = -1;
        int localEdge = -1;
        int dir = 0; // +1 if stored as min->max, -1 if max->min
    };
    using EdgeKey = std::pair<int, int>;

    std::map<EdgeKey, std::vector<FaceEdgeUse>> edgeUses;

    // 1. Build undirected edge -> incident face uses.
    for (int f = 0; f < F_List.size(); f++) {
        const auto& face = F_List[f];
        int n = face.size();

        if (n < 3) {
            return false;
        }

        for (int i = 0; i < n; i++) {
            int a = face[i];
            int b = face[(i + 1) % n];

            if (a == b) {
                return false;
            }

            EdgeKey key = undirectedKey(a, b);
            edgeUses[key].push_back(FaceEdgeUse{
                f,
                i,
                edgeDirRelativeToKey(a, b)
            });
        }
    }

    // 2. Reject nonmanifold edges for this mesh structure.
    for (const auto& kv : edgeUses) {
        if (kv.second.size() > 2) {
            return false;
        }
    }

    // 3. Build face adjacency with "same direction?" relation.
    std::vector<std::vector<std::pair<int, bool>>> faceAdj(F_List.size());

    for (const auto& kv : edgeUses) {
        const auto& uses = kv.second;

        if (uses.size() != 2) {
            continue; // boundary edge
        }

        const FaceEdgeUse& a = uses[0];
        const FaceEdgeUse& b = uses[1];

        // If two faces use the shared undirected edge in the same direction,
        // one of them must be flipped.
        bool sameDir = (a.dir == b.dir);

        faceAdj[a.face].push_back({b.face, sameDir});
        faceAdj[b.face].push_back({a.face, sameDir});
    }

    // 4. BFS assign flip parity per connected component.
    // flip[f] == 0 means keep original orientation.
    // flip[f] == 1 means reverse this face.
    std::vector<int> flip(F_List.size(), -1);

    for (int root = 0; root < F_List.size(); root++) {
        if (flip[root] != -1) {
            continue;
        }

        flip[root] = 0;
        std::queue<int> q;
        q.push(root);

        while (!q.empty()) {
            int f = q.front();
            q.pop();

            for (auto [g, sameDir] : faceAdj[f]) {
                // If sameDir is true, neighbor must have opposite flip parity.
                // If sameDir is false, neighbor must have same flip parity.
                int requiredFlip = flip[f] ^ static_cast<int>(sameDir);

                if (flip[g] == -1) {
                    flip[g] = requiredFlip;
                    q.push(g);
                } else if (flip[g] != requiredFlip) {
                    // Contradiction: non-orientable or inconsistent connectivity.
                    return false;
                }
            }
        }
    }

    // 5. Apply flips.
    for (int f = 0; f < F_List.size(); f++) {
        if (flip[f]) {
            std::reverse(F_List[f].begin(), F_List[f].end());
        }
    }

    return true;
}

// GEOMETRY HELPERS

// Rotation-variant SVD
// Compute the Rotation Variant SVD
// Input an empty U, Sigma, V
void rotationVariantSVD(Eigen::Matrix3d& mat, Eigen::Matrix3d& U, Eigen::Vector3d& Sigma, Eigen::Matrix3d& V) {
    // Compute SVD to get Sigma, U and V
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(mat, Eigen::ComputeFullU | Eigen::ComputeFullV);
    U = svd.matrixU();
    Sigma = svd.singularValues();
    V = svd.matrixV();

    // Compute L and remove reflections from U and V
    Eigen::Matrix3d L;
    L.setIdentity();
    L(2,2) = (U*V.transpose()).determinant();
    Sigma(2) = Sigma(2) * L(2,2); // To keep it a vector (taken from HOBAK)

    double u_det = U.determinant();
    double v_det = V.determinant();
    if (u_det < 0 && v_det > 0) {
        U = U * L;
    } else if (u_det > 0 && v_det < 0) {
        V = V * L;
    }

    return;
}

// Compute Polar Decomposition given rotation variant SVD
// Returns a vector containing R and then S
void polarDecomposition(Eigen::Matrix3d& mat, Eigen::Matrix3d& R, Eigen::Matrix3d& S) {
    Eigen::Matrix3d U, V; 
    Eigen::Vector3d Sigma;
    rotationVariantSVD(mat, U, Sigma, V);

    // Put together R and S from the inputs
    R = U * V.transpose();
    S = V * Sigma.asDiagonal() * V.transpose();
    return;
}

// Compute 3D signed angle between two vectors
double signedAngle(const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, const Eigen::Vector3d& axis, bool positive) {
    const double eps = 1e-12;

    if (v0.squaredNorm() <= eps || v1.squaredNorm() <= eps || axis.squaredNorm() <= eps) {
        return 0.0;
    }

    Eigen::Vector3d a = v0.normalized();
    Eigen::Vector3d b = v1.normalized();
    Eigen::Vector3d n = axis.normalized();

    double sinTheta = n.dot(a.cross(b));
    double cosTheta = std::clamp(a.dot(b), -1.0, 1.0);

    double sAngle = std::atan2(sinTheta, cosTheta);

    if (positive && sAngle < 0.0) {
        sAngle += 2.0 * M_PI;
    }

    return sAngle;
}

// Find the closest point to a triangle
Eigen::Vector3d triangleClosestPoint(const std::vector<Eigen::Vector3d> triVerts, const Eigen::Vector3d p) {
    const double eps = 1e-8;
    const Eigen::Vector3d& a = triVerts[0];
    const Eigen::Vector3d& b = triVerts[1];
    const Eigen::Vector3d& c = triVerts[2];

    const Eigen::Vector3d ab = b - a;
    const Eigen::Vector3d ac = c - a;
    const Eigen::Vector3d ap = p - a;

    const double d1 = ab.dot(ap);
    const double d2 = ac.dot(ap);

    // Vertex region outside A
    if (d1 <= 0.0 && d2 <= 0.0) {
        return a;
    }

    const Eigen::Vector3d bp = p - b;
    const double d3 = ab.dot(bp);
    const double d4 = ac.dot(bp);

    // Vertex region outside B
    if (d3 >= 0.0 && d4 <= d3) {
        return b;
    }

    // Edge region AB
    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        const double denom = d1 - d3;
        if (std::abs(denom) <= eps) {
            return Utils::closestPointOnSegment3D(p, a, b, true);
        }

        const double t = d1 / denom;
        return a + t * ab;
    }

    const Eigen::Vector3d cp = p - c;
    const double d5 = ab.dot(cp);
    const double d6 = ac.dot(cp);

    // Vertex region outside C
    if (d6 >= 0.0 && d5 <= d6) {
        return c;
    }

    // Edge region AC
    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        const double denom = d2 - d6;
        if (std::abs(denom) <= eps) {
            return Utils::closestPointOnSegment3D(p, a, c, true);
        }
        const double t = d2 / denom;
        return a + t * ac;
    }

    // Edge region BC
    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        const double denom = (d4 - d3) + (d5 - d6);
        if (std::abs(denom) <= eps) {
            return Utils::closestPointOnSegment3D(p, b, c, true);
        }
        const double t = (d4 - d3) / denom;
        return b + t * (c - b);
    }

    // Inside face region
    const double denom = va + vb + vc;
    if (std::abs(denom) <= eps) {
        // Degenerate triangle fallback: closest point among three edges.
        Eigen::Vector3d pab = Utils::closestPointOnSegment3D(p, a, b, true);
        Eigen::Vector3d pac = Utils::closestPointOnSegment3D(p, a, c, true);
        Eigen::Vector3d pbc = Utils::closestPointOnSegment3D(p, b, c, true);

        double dab = (p - pab).squaredNorm();
        double dac = (p - pac).squaredNorm();
        double dbc = (p - pbc).squaredNorm();

        if (dab <= dac && dab <= dbc) return pab;
        if (dac <= dab && dac <= dbc) return pac;
        return pbc;
    }

    const double invDenom = 1.0 / denom;
    const double vBary = vb * invDenom;
    const double wBary = vc * invDenom;
    return a + vBary * ab + wBary * ac;
}

// Overload with returned element type
Eigen::Vector3d triangleClosestPoint(const std::vector<Eigen::Vector3d>& triVerts, const Eigen::Vector3d& p, 
                                         int& projType, int& projIdx, double snapTol) {
    const double eps = 1e-8;
    projType = -1;
    projIdx = -1;

    if (triVerts.size() != 3) {
        return Eigen::Vector3d::Zero();
    }

    const Eigen::Vector3d& a = triVerts[0];
    const Eigen::Vector3d& b = triVerts[1];
    const Eigen::Vector3d& c = triVerts[2];

    const Eigen::Vector3d ab = b - a;
    const Eigen::Vector3d ac = c - a;
    const Eigen::Vector3d ap = p - a;

    const double d1 = ab.dot(ap);
    const double d2 = ac.dot(ap);

    Eigen::Vector3d cp;

    // Vertex region outside A
    if (d1 <= 0.0 && d2 <= 0.0) {
        cp = a;
        projType = 0;
        projIdx = 0;
    } else {
        const Eigen::Vector3d bp = p - b;
        const double d3 = ab.dot(bp);
        const double d4 = ac.dot(bp);

        // Vertex region outside B
        if (d3 >= 0.0 && d4 <= d3) {
            cp = b;
            projType = 0;
            projIdx = 1;
        } else {
            const double vc = d1 * d4 - d3 * d2;

            // Edge region AB
            if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
                const double denom = d1 - d3;
                double t = 0.0;
                if (std::abs(denom) > eps) {
                    t = d1 / denom;
                }
                cp = a + t * ab;
                projType = 1;
                projIdx = 0; // edge 0 -> 1
            } else {
                const Eigen::Vector3d cpv = p - c;
                const double d5 = ab.dot(cpv);
                const double d6 = ac.dot(cpv);
                // Vertex region outside C
                if (d6 >= 0.0 && d5 <= d6) {
                    cp = c;
                    projType = 0;
                    projIdx = 2;
                } else {
                    const double vb = d5 * d2 - d1 * d6;
                    // Edge region AC
                    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
                        const double denom = d2 - d6;
                        double t = 0.0;
                        if (std::abs(denom) > eps) {
                            t = d2 / denom;
                        }
                        cp = a + t * ac;
                        projType = 1;
                        projIdx = 2; // edge 2 -> 0
                    } else {
                        const double va = d3 * d6 - d5 * d4;
                        // Edge region BC
                        if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
                            const double denom = (d4 - d3) + (d5 - d6);
                            double t = 0.0;
                            if (std::abs(denom) > eps) {
                                t = (d4 - d3) / denom;
                            }
                            cp = b + t * (c - b);
                            projType = 1;
                            projIdx = 1; // edge 1 -> 2
                        } else {
                            const double denom = va + vb + vc;
                            if (std::abs(denom) <= eps) {
                                // Degenerate triangle fallback.
                                Eigen::Vector3d pab = closestPointOnSegment3D(p, a, b, true);
                                Eigen::Vector3d pbc = closestPointOnSegment3D(p, b, c, true);
                                Eigen::Vector3d pca = closestPointOnSegment3D(p, c, a, true);

                                double dab = (p - pab).squaredNorm();
                                double dbc = (p - pbc).squaredNorm();
                                double dca = (p - pca).squaredNorm();

                                if (dab <= dbc && dab <= dca) {
                                    cp = pab;
                                    projType = 1;
                                    projIdx = 0;
                                } else if (dbc <= dab && dbc <= dca) {
                                    cp = pbc;
                                    projType = 1;
                                    projIdx = 1;
                                } else {
                                    cp = pca;
                                    projType = 1;
                                    projIdx = 2;
                                }
                            } else {
                                const double invDenom = 1.0 / denom;
                                const double vBary = vb * invDenom;
                                const double wBary = vc * invDenom;
                                cp = a + vBary * ab + wBary * ac;
                                projType = 2;
                                projIdx = -1;
                            }
                        }
                    }
                }
            }
        }
    }

    // Optional snapping override.
    if (snapTol > 0.0) {
        for (int i = 0; i < 3; i++) {
            if ((cp - triVerts[i]).norm() <= snapTol) {
                projType = 0;
                projIdx = i;
                return triVerts[i];
            }
        }

        for (int i = 0; i < 3; i++) {
            int j = (i + 1) % 3;
            Eigen::Vector3d ecp = closestPointOnSegment3D(cp, triVerts[i], triVerts[j], true);

            if ((cp - ecp).norm() <= snapTol) {
                projType = 1;
                projIdx = i;
                return ecp;
            }
        }
    }

    return cp;
}

// Computes the closest point to a bilinear patch
// NOTE: Currently uses a "cheap" convergence check by thresholding u, v
int bilinearPatchClosestPoint(const std::vector<Eigen::Vector3d>& patchVerts, const Eigen::Vector3d& p, double& u, double& v,double eps, int max_iter) {
    if (patchVerts.size() != 4) {
        return -1;
    }
    // 1. Find closest point on segments
    Eigen::Vector3d new_guess;
    double closest = std::numeric_limits<double>::infinity();
    int closest_idx = 0;
    for (int e = 0; e < 4; e++) {
        Eigen::Vector3d closestOnEdge = closestPointOnSegment3D(p, patchVerts[e], patchVerts[(e+1)%4]);
        double dist = (p - closestOnEdge).squaredNorm();
        if (dist < closest) {
            closest = dist;
            new_guess = closestOnEdge;
            closest_idx = e;
        }
    }
    double edge_len = (patchVerts[(closest_idx+1)%4] - patchVerts[closest_idx]).norm();
    double seg_len = (new_guess - patchVerts[closest_idx]).norm();
    if (edge_len <= eps) {  // Minor guard against bad behavior
        edge_len = 1.0;
        seg_len = 0.0;
    }
    double t = seg_len / edge_len;
    if (closest_idx == 0) {
        u = t;
        v = 0.0;
    } else if (closest_idx == 1) {
        u = 1.0;
        v = t;
    } else if (closest_idx == 2) {
        u = 1.0 - t;
        v = 1.0;
    } else {
        u = 0.0;
        v = 1.0 - t;
    }
    // 3. Given this start guess, optimize
    // double curr_loss = (new_guess - p).squaredNorm();
    Eigen::Vector3d a = patchVerts[1] - patchVerts[0];
    Eigen::Vector3d b = patchVerts[2] - patchVerts[1];
    Eigen::Vector3d c = patchVerts[3] - patchVerts[2];
    Eigen::Vector3d d = patchVerts[0] - patchVerts[3];
    int counter = 0;
    for (int i = 0; i < max_iter; i++) {
        double old_u = u;
        double old_v = v;
        // Update u with v fixed.
        Eigen::Vector3d A = (1.0 - v) * patchVerts[0] + v * patchVerts[3];
        Eigen::Vector3d B = (1.0 - v) * (a) - v * (c);
        double denom = B.squaredNorm();
        if (denom > 1e-16) {
            u = std::clamp((p - A).dot(B) / denom, 0.0, 1.0);
        }

        // Update v with u fixed.
        Eigen::Vector3d C = (1.0 - u) * patchVerts[0] + u * patchVerts[1];
        Eigen::Vector3d D = u * (b) - (1.0 - u) * (d);
        denom = D.squaredNorm();
        if (denom > 1e-16) {
            v = std::clamp((p - C).dot(D) / denom, 0.0, 1.0);
        }
        counter++;
        // If we have stagnated, then return
        if (std::abs(u - old_u) <= eps && std::abs(v - old_v) <= eps) {
            break;
        }
    }
    // std::cout << "Converged in " << counter << " iterations: (" << u << ", " << v << ")" << std::endl;
    return 1;
}

Eigen::Vector3d bilinearPatchClosestPoint(const std::vector<Eigen::Vector3d>& patchVerts, const Eigen::Vector3d& p, double eps, int max_iter) {
    double u, v;
    if (bilinearPatchClosestPoint(patchVerts, p, u, v, eps, max_iter) != -1) {
        return Eigen::Vector3d::Zero();
    }
    return bilinearPatch(patchVerts, u, v);
}

// Overload with returned element type
Eigen::Vector3d bilinearPatchClosestPoint(const std::vector<Eigen::Vector3d>& patchVerts, const Eigen::Vector3d& p,
                                              int& projType, int& projIdx, double snapTol, double eps, int max_iter) {
    projType = -1;
    projIdx = -1;

    if (patchVerts.size() != 4) {
        return Eigen::Vector3d::Zero();
    }

    Eigen::Vector3d new_guess;
    double u = 0.0;
    double v = 0.0;

    double closest = std::numeric_limits<double>::infinity();
    int closest_idx = 0;

    for (int e = 0; e < 4; e++) {
        Eigen::Vector3d closestOnEdge = closestPointOnSegment3D(
            p,
            patchVerts[e],
            patchVerts[(e + 1) % 4],
            true
        );
        double dist = (p - closestOnEdge).squaredNorm();
        if (dist < closest) {
            closest = dist;
            new_guess = closestOnEdge;
            closest_idx = e;
        }
    }

    double edge_len = (patchVerts[(closest_idx + 1) % 4] - patchVerts[closest_idx]).norm();
    double seg_len = (new_guess - patchVerts[closest_idx]).norm();

    if (edge_len <= eps) {
        edge_len = 1.0;
        seg_len = 0.0;
    }

    double t = seg_len / edge_len;

    if (closest_idx == 0) {
        u = t;
        v = 0.0;
    } else if (closest_idx == 1) {
        u = 1.0;
        v = t;
    } else if (closest_idx == 2) {
        u = 1.0 - t;
        v = 1.0;
    } else {
        u = 0.0;
        v = 1.0 - t;
    }

    Eigen::Vector3d a = patchVerts[1] - patchVerts[0];
    Eigen::Vector3d b = patchVerts[2] - patchVerts[1];
    Eigen::Vector3d c = patchVerts[3] - patchVerts[2];
    Eigen::Vector3d d = patchVerts[0] - patchVerts[3];

    for (int iter = 0; iter < max_iter; iter++) {
        double old_u = u;
        double old_v = v;

        Eigen::Vector3d A = (1.0 - v) * patchVerts[0] + v * patchVerts[3];
        Eigen::Vector3d B = (1.0 - v) * a - v * c;

        double denom = B.squaredNorm();
        if (denom > 1e-16) {
            u = std::clamp((p - A).dot(B) / denom, 0.0, 1.0);
        }

        Eigen::Vector3d C = (1.0 - u) * patchVerts[0] + u * patchVerts[1];
        Eigen::Vector3d D = u * b - (1.0 - u) * d;

        denom = D.squaredNorm();
        if (denom > 1e-16) {
            v = std::clamp((p - C).dot(D) / denom, 0.0, 1.0);
        }

        if (std::abs(u - old_u) <= eps && std::abs(v - old_v) <= eps) {
            break;
        }
    }

    Eigen::Vector3d cp = bilinearPatch(patchVerts, u, v);

    double featureTol = snapTol > 0.0 ? snapTol : eps;

    // Vertex classification.
    if (snapTol > 0.0) {
        for (int i = 0; i < 4; i++) {
            if ((cp - patchVerts[i]).norm() <= snapTol) {
                projType = 0;
                projIdx = i;
                return patchVerts[i];
            }
        }
    }

    // Parametric vertex classification, useful even when snapTol == 0.
    if (u <= featureTol && v <= featureTol) {
        projType = 0;
        projIdx = 0;
    } else if (u >= 1.0 - featureTol && v <= featureTol) {
        projType = 0;
        projIdx = 1;
    } else if (u >= 1.0 - featureTol && v >= 1.0 - featureTol) {
        projType = 0;
        projIdx = 2;
    } else if (u <= featureTol && v >= 1.0 - featureTol) {
        projType = 0;
        projIdx = 3;
    } else if (v <= featureTol) {
        projType = 1;
        projIdx = 0;
    } else if (u >= 1.0 - featureTol) {
        projType = 1;
        projIdx = 1;
    } else if (v >= 1.0 - featureTol) {
        projType = 1;
        projIdx = 2;
    } else if (u <= featureTol) {
        projType = 1;
        projIdx = 3;
    } else {
        projType = 2;
        projIdx = -1;
    }

    // If snapping to edge is requested, replace cp with exact segment projection.
    if (snapTol > 0.0 && projType == 1) {
        int i = projIdx;
        int j = (i + 1) % 4;
        cp = closestPointOnSegment3D(cp, patchVerts[i], patchVerts[j], true);
    }

    return cp;
}

// Evaluate bilinear patch point at specified u, v
// NOTE: I am not guarding values outside the range 0, 1
Eigen::Vector3d bilinearPatch(const std::vector<Eigen::Vector3d>& patchVerts, double u, double v) {
    return (1 - v) * ((1 - u) * patchVerts[0] + u * patchVerts[1]) + v * ((1 - u) * patchVerts[3] + u * patchVerts[2]);
}

Eigen::Vector3d polygonClosestPointNewell(const std::vector<Eigen::Vector3d>& polyVerts, const Eigen::Vector3d& p,
                                          const Eigen::Vector3d& polyNormal, int& projType, int& projIdx, double snapTol) {
    projType = -1;
    projIdx = -1;
    const int n = polyVerts.size();

    if (n < 3) {
        return Eigen::Vector3d::Zero();
    }

    Eigen::Vector3d barycenter = Eigen::Vector3d::Zero();
    for (int i = 0; i < n; i++) {
        barycenter += polyVerts[i];
    }
    barycenter /= static_cast<double>(n);

    Eigen::Vector3d t1;
    Eigen::Vector3d t2;
    buildPlaneBasis(polyNormal, t1, t2);

    Eigen::Vector3d pProj3D = projectPointOntoPlane(polyNormal, barycenter, p);
    Eigen::Vector2d pProj2D = convertTo2D(pProj3D, barycenter, t1, t2);
    std::vector<Eigen::Vector2d> poly2D(n);

    for (int i = 0; i < n; i++) {
        Eigen::Vector3d viProj3D = projectPointOntoPlane(polyNormal, barycenter, polyVerts[i]);
        poly2D[i] = convertTo2D(viProj3D, barycenter, t1, t2);
    }
    Eigen::Vector3d cp3D;
    Eigen::Vector2d cp2D;

    bool inside = pointInPolygon2D(pProj2D, poly2D);
    if (inside) {
        cp2D = pProj2D;
        cp3D = revertTo3D(cp2D, barycenter, t1, t2);
        projType = 2;
        projIdx = -1;
    } else {
        double bestDist2 = std::numeric_limits<double>::infinity();
        for (int i = 0; i < n; i++) {
            int j = (i + 1) % n;
            Eigen::Vector3d edgeCP = closestPointOnSegment3D(p, polyVerts[i], polyVerts[j], true);
            double d2 = (p - edgeCP).squaredNorm();

            if (d2 < bestDist2) {
                bestDist2 = d2;
                cp3D = edgeCP;
                projType = 1;
                projIdx = i;
            }
        }
        if (snapTol > 0.0) {
            for (int i = 0; i < n; i++) {
                if ((cp3D - polyVerts[i]).norm() <= snapTol) {
                    projType = 0;
                    projIdx = i;
                    return polyVerts[i];
                }
            }
        }
        return cp3D;
    }

    // Snap inside-face projection to local vertices/edges.
    if (snapTol > 0.0) {
        for (int i = 0; i < n; i++) {
            if ((cp2D - poly2D[i]).norm() <= snapTol) {
                projType = 0;
                projIdx = i;
                return polyVerts[i];
            }
        }

        for (int i = 0; i < n; i++) {
            int j = (i + 1) % n;
            Eigen::Vector2d edgeCP2D = closestPointOnSegment2D(cp2D, poly2D[i], poly2D[j], true);
            double d = (cp2D - edgeCP2D).norm();
            if (d <= snapTol) {
                Eigen::Vector2d e = poly2D[j] - poly2D[i];
                double denom = e.squaredNorm();
                double t = 0.0;
                if (denom > 1e-20) {
                    t = (edgeCP2D - poly2D[i]).dot(e) / denom;
                    t = std::clamp(t, 0.0, 1.0);
                }
                projType = 1;
                projIdx = i;
                return polyVerts[i] + t * (polyVerts[j] - polyVerts[i]);
            }
        }
    }
    return cp3D;
}


// Converts a 3D direction into an angle in the tangent plane spanned by t1, t2.
// Returns false if the projected direction is degenerate.
bool directionAngleInPlane(
    const Eigen::Vector3d& origin,
    const Eigen::Vector3d& target,
    const Eigen::Vector3d& normal,
    const Eigen::Vector3d& t1,
    const Eigen::Vector3d& t2,
    double& theta
) {
    const double eps = 1e-12;

    Eigen::Vector3d n = normal;
    if (n.squaredNorm() <= eps) {
        return false;
    }
    n.normalize();

    Eigen::Vector3d d = target - origin;
    d = d - d.dot(n) * n;

    if (d.squaredNorm() <= eps) {
        return false;
    }

    const double x = d.dot(t1);
    const double y = d.dot(t2);

    theta = std::atan2(y, x);
    if (theta < 0.0) {
        theta += 2.0 * M_PI;
    }

    return true;
}

// Returns true if two angular values are effectively the same direction.
bool anglesCoincident(double a, double b, double eps) {
    double diff = std::abs(a - b);
    diff = std::min(diff, 2.0 * M_PI - diff);
    return diff <= eps;
}

// Given two unit vectors, compute the rotation from one to the other
// Rotation formulation taken from https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
Eigen::Matrix3d computeRotation(const Eigen::Vector3d& u, const Eigen::Vector3d& v) {
    const double eps = 1e-8;

    if (u.norm() <= eps || v.norm() <= eps) {
        return Eigen::Matrix3d::Identity();
    }

    Eigen::Vector3d a = u.normalized();
    Eigen::Vector3d b = v.normalized();

    double cosUV = std::clamp(a.dot(b), -1.0, 1.0);

    if (cosUV >= 1.0 - eps) {
        return Eigen::Matrix3d::Identity();
    }

    if (cosUV <= -1.0 + eps) {
        Eigen::Vector3d axis = a.unitOrthogonal();
        return Eigen::AngleAxisd(M_PI, axis).toRotationMatrix();
    }

    Eigen::Vector3d axis = a.cross(b).normalized();
    double theta = std::acos(cosUV);

    return Eigen::AngleAxisd(theta, axis).toRotationMatrix();
}

// Overload
Eigen::Matrix3d computeRotation(const Eigen::Vector3d& axis, const double& theta) {
    const double eps = 1e-8;

    if (axis.norm() <= eps) {
        return Eigen::Matrix3d::Identity();
    }

    return Eigen::AngleAxisd(theta, axis.normalized()).toRotationMatrix();
}

// Find basis vectors for a planar region (ex. tangent plane)
// Build plane basis given only n, and unitialized t1, t2
void buildPlaneBasis(const Eigen::Vector3d& n, Eigen::Vector3d& t1, Eigen::Vector3d& t2) {
    const double eps = 1e-12;
    if (n.squaredNorm() <= eps) {
        t1 = Eigen::Vector3d::UnitX();
        t2 = Eigen::Vector3d::UnitY();
        return;
    }

    Eigen::Vector3d n_norm = n.normalized();
    if (std::abs(n(0)) < 0.9)
        t1 = n_norm.cross(Eigen::Vector3d::UnitX()).normalized();
    else
        t1 = n_norm.cross(Eigen::Vector3d::UnitY()).normalized();
    t2 = n_norm.cross(t1); // already unit
    return;
}

// Given a point on a plane basis and the plane basis, convert to 2D planar point
Eigen::Vector2d convertTo2D(const Eigen::Vector3d& p, const Eigen::Vector3d& origin, const Eigen::Vector3d& t1, const Eigen::Vector3d& t2) {
    Eigen::Vector3d vec = p - origin;
    return Eigen::Vector2d(vec.dot(t1), vec.dot(t2));
}

// Given a 2D planar point and the plane basis, revert to its 3D counterapart
Eigen::Vector3d revertTo3D(const Eigen::Vector2d& p, const Eigen::Vector3d& origin, const Eigen::Vector3d& t1, const Eigen::Vector3d& t2) {
    return origin + p(0) * t1 + p(1) * t2;
}

// Projects a vector onto the tangent plane of a normal vector.
// If the resulting projection is near-degenerate, then return a random unit tangent vector.
double projectVectorOntoTangentPlane(const Eigen::Vector3d& normal, const Eigen::Vector3d& vec, Eigen::Vector3d& proj, double scale) {
    const double eps = 1e-12;
    // Ensure normal is unit
    Eigen::Vector3d n = normal.normalized();

    // Project onto tangent plane (i.e., subtract out projection onto normal vector)
    proj = vec - vec.dot(n) * n;
    double len = proj.norm();

    // Under degeneracy, sample a random unit vector in the tangent plane using some
    // constructed tangent plane basis
    if (len < eps) {
        // Build orthonormal tangent basis
        Eigen::Vector3d t1;
        Eigen::Vector3d t2;
        buildPlaneBasis(normal, t1, t2);

        // Pick vector by sampling a random angle in [0, 2pi)
        static std::mt19937 gen(std::random_device{}());
        static std::uniform_real_distribution<double> angle_dist(0.0, 2.0 * M_PI);
        // Construct proj using basis vectors
        double theta = angle_dist(gen);
        proj = std::cos(theta) * t1 + std::sin(theta) * t2;
        return -1;
    }
    // Normalize, scale, and return length
    proj /= len;
    proj *= scale;
    return len;
}

// Projects a point onto the tangent plane of a normal given a center 
Eigen::Vector3d projectPointOntoPlane(const Eigen::Vector3d& normal, const Eigen::Vector3d& center, const Eigen::Vector3d& p) {
    const double eps = 1e-12;
    if (normal.squaredNorm() <= eps) {
        return p;
    }
    Eigen::Vector3d n = normal.normalized();
    return p - (p - center).dot(n) * n;
}

// Check if a 2D point is in a 2D polygon
// To do this, we do raycasting to the segment
bool pointInPolygon2D(const Eigen::Vector2d& p, const std::vector<Eigen::Vector2d>& poly) {
    bool inside = false;
    int n = poly.size();

    for (int v0 = 0; v0 < n; v0++) {
        int v1 = (v0 + n - 1) % n;
        const Eigen::Vector2d& p0 = poly[v0];
        const Eigen::Vector2d& p1 = poly[v1];
        // If 
        bool intersect = ((p0(1) > p(1)) != (p1(1) > p(1))) &&
                         (p(0) < (p1(0) - p0(0)) * (p(1) - p0(1)) / (p1(1) - p0(1)) + p0(0));

        if (intersect) {
            inside = !inside;
        }
    }
    return inside;
}

// Get the closest point on a segment in 2D, where the endpoints are defined
// To do this, project onto parameterized segment and snap t to [0, 1]
Eigen::Vector2d closestPointOnSegment2D(const Eigen::Vector2d& p, const Eigen::Vector2d& v0, const Eigen::Vector2d& v1, bool clip) {
    Eigen::Vector2d vec = v1 - v0;
    double denom = vec.squaredNorm();
    if (denom < 1e-16) {
        return v0;
    }
    double t = (p - v0).dot(vec) / denom;
    if (clip) {
        t = std::max(0.0, std::min(1.0, t));
    }
    return v0 + t * vec;
}

double cross2D(const Eigen::Vector2d& a, const Eigen::Vector2d& b) {
    return a.x() * b.y() - a.y() * b.x();
}

// TODO: Fix this function so that it properly sets u
bool raycastToSegment2D(const Eigen::Vector2d& p, const Eigen::Vector2d& direc, const Eigen::Vector2d& v0, const Eigen::Vector2d& v1,
                        double& t, double& u) {
    double eps = 1e-12;
    const Eigen::Vector2d seg = v1 - v0;
    Eigen::Vector2d direc_norm = direc.normalized();
    const double denom = cross2D(direc_norm, seg);

    // Parallel or nearly parallel
    if (std::abs(denom) < eps) {
        return false;
    }

    const Eigen::Vector2d rhs = v0 - p;

    t = cross2D(rhs, seg) / denom;
    u = cross2D(rhs, direc_norm) / denom;

    // Ray constraint and segment constraint
    if (t < -eps) {
        return false;
    }
    if (u < -eps || u > 1.0 + eps) {
        return false;
    }

    // Clamp tiny numerical drift
    if (t < 0.0) {
        t = 0.0;
    }
    if (u < 0.0) {
        u = 0.0;
    } else if (u > 1.0) {
        u = 1.0;
    }
    return true;
}

// Get the closest point on a segment in 3D, where the endpoints are defined
// To do this, project onto parameterized segment and snap t to [0, 1]
Eigen::Vector3d closestPointOnSegment3D(const Eigen::Vector3d& p, const Eigen::Vector3d& v0, const Eigen::Vector3d& v1, bool clip) {
    Eigen::Vector3d vec = v1 - v0;
    double denom = vec.squaredNorm();
    if (denom < 1e-16)
        return v0;
    double t = (p - v0).dot(vec) / denom;
    if (clip) {
        t = std::max(0.0, std::min(1.0, t));
    }
    return v0 + t * vec;
}

// MEAN VALUE COORDINATES
// Compute angle between any two 2D vectors given the four endpoints
// Vectors are computed as p1 - p0, p3 - p2
double vectorAngle(const Eigen::Vector2d& p0, const Eigen::Vector2d& p1, const Eigen::Vector2d& p2, const Eigen::Vector2d& p3) {
    Eigen::Vector2d v0 = p1 - p0;
    Eigen::Vector2d v1 = p3 - p2;
    return std::atan2(v0(0)*v1(1) - v0(1)*v1(0), v0(0)*v1(0) + v0(1)*v1(1));
}

// Compute the sign of a double value
double computeSign(const double& value) {
    if (value > 0.0) {
        return 1.0;
    } else if (value < 0.0) {
        return -1.0;
    }
    return 0.0;
}

void meanValueCoordinates(const Eigen::Vector2d& target, const std::vector<Eigen::Vector2d>& cage, Eigen::VectorXd& weights) {
    weights.setZero();

    double W = 0.0;
    std::vector<double> beta(cage.size());
    std::vector<double> gamma(cage.size());
    std::vector<double> s(cage.size());
    std::vector<double> r(cage.size());
    for (int v = 0; v < cage.size(); v++) {
        beta[v] = vectorAngle(cage[v], cage[(v+1)%cage.size()], target, cage[v]);
        gamma[v] = vectorAngle(cage[(v+1)%cage.size()], cage[v], cage[(v+1)%cage.size()], target);
        s[v] = beta[v] + gamma[v];
        r[v] = (cage[v] - target).norm();
    }
    std::vector<double> w(cage.size());
    for (int v = 0; v < cage.size(); v++) {
        int v_m1 = (v + cage.size() - 1)%cage.size();
        double alpha_m1p1 = vectorAngle(target, cage[v_m1], target, cage[(v+1)%cage.size()]);
        double s_m1p1 = M_PI * (computeSign(s[v_m1]) + computeSign(s[v])) - s[v_m1] - s[v];
        if (computeSign(alpha_m1p1) != computeSign(s_m1p1)) {    // NOTE: Potential numerical issue here
            alpha_m1p1 *= -1.0;
        }
        w[v] = r[v_m1] * std::sin(alpha_m1p1/2.0);
        for (int u = 0; u < cage.size(); u++) {
            if ((u != v_m1) && (u != v)) {
                w[v] *= r[u] * std::sin(std::abs(s[v])/2.0);
            }
        }
        W += w[v];
    }
    // Normalize to get final weights
    for (int v = 0; v < cage.size(); v++) {
        weights[v] = w[v]/W;
    }

    return;
}

} // namespace Utils