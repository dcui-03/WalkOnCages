#include "mesh.hpp"

#include "../utils/decUtils.hpp"
#include "../utils/utils.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <limits>
#include <queue>
#include <utility>
#include <algorithm>
#include <iostream>

// Utility functions for mesh (projection, etc.)

namespace Mesh {

// Project a vertex onto the mesh. If multiple, just picks the one with smaller index.
// Also returns the element type that was landed on.
// For non-planar faces, I am just going to fit a Newell plane using the barycenter and vector area + a barycentric height interpolation
int mesh::computeVProjection(const Eigen::Vector3d& v, Eigen::Vector3d& proj, int& elIdx, bool snap, bool fast) const {
    vertProjData projData = computeVProjection(v, proj, snap, fast);
    elIdx = projData.elIdx;
    return projData.elType;
}
vertProjData mesh::computeVProjection(const Eigen::Vector3d& v,
                                      Eigen::Vector3d& proj,
                                      bool snap,
                                      bool fast) const {
    vertProjData projData({-1, -1});

    if (active_f == 0) {
        return projData;
    }

    // 1. Closest point search
    int success = -1;
    // If BVH is not empty, then find the closest point using the BVH
    if (!BVH.empty()) {
        success = closestFaceBVH(v, proj, projData, snap);
    }

    // Fallback to brute force if BVH is unavailable or failed
    if (success != 1) {
        double minDist2 = std::numeric_limits<double>::infinity();

        for (int f = 0; f < F.size(); f++) {
            if (!F[f].active) {
                continue;
            }

            Eigen::Vector3d candidateProj;
            vertProjData candidateData({2, f});

            if (closestPointOnFace(f, v, candidateProj, candidateData, snap) != 1) {
                continue;
            }

            double dist2 = (v - candidateProj).squaredNorm();

            if (dist2 < minDist2) {
                minDist2 = dist2;
                proj = candidateProj;
                projData = candidateData;
            }
        }
    }
    if (projData.elIdx < 0) {   // Invalid
        return projData;
    }

    // 2. If fast or snapped to edge/vert, return immediately
    if (fast || projData.elType != 2) {
        return projData;
    }

    // 3. Slow MVC-based lifting for non-planar 5+ sided faces
    std::vector<int> fVerts = faceAdjVertIdxs(projData.elIdx);
    int fSize = fVerts.size();
    std::vector<Eigen::Vector3d> fVertsPos = adjVerts(fVerts);
    Eigen::Vector3d fN = F[projData.elIdx].n;
    // If we hit a triangle or bilinear patch, then we have an exact hit. No need to lift
    if (fSize <= 4) {
        return projData;
    }

    Eigen::VectorXd fHeight = computeFaceHeight(projData.elIdx);
    // Check if we are on a planar face
    bool planar = true;
    for (int i = 0; i < fSize; i++) {
        if (std::abs(fHeight(i)) >= 1e-6) {
            planar = false;
            break;
        }
    }
    // If we are on a planar face, then we have an "exact" answer
    if (planar) {
        return projData;
    }

    // Otherwise, do lifting
    Eigen::Vector3d barycenter = DECUtils::computeBarycenter(fVertsPos);
    Eigen::Vector3d t1;
    Eigen::Vector3d t2;
    Utils::buildPlaneBasis(fN, t1, t2);
    Eigen::Vector2d v_proj2d = Utils::convertTo2D(proj, barycenter, t1, t2);
    std::vector<Eigen::Vector2d> fVerts2D(fSize);
    for (int fv = 0; fv < fSize; fv++) {
        Eigen::Vector3d fv_proj3D = Utils::projectPointOntoPlane(fN, barycenter, fVertsPos[fv]);
        fVerts2D[fv] = Utils::convertTo2D(fv_proj3D, barycenter, t1, t2);
    }
    Eigen::VectorXd MVCWeights(fSize);
    Utils::meanValueCoordinates(v_proj2d, fVerts2D, MVCWeights);
    double h = MVCWeights.dot(fHeight);
    proj += h * fN;

    return projData;
}

// Squared distance from a point to an AABB
double mesh::pointAABBDist2(const Eigen::Vector3d& p, int box) const {
    if (box < 0 || box >= BVH.size()) {
        return std::numeric_limits<double>::infinity();
    }
    const Eigen::Vector3d& bmin = BVH[box].bdyVerts.first;
    const Eigen::Vector3d& bmax = BVH[box].bdyVerts.second;

    double d2 = 0.0;
    // Find the closest distance to the bounding box by testing against each axis
    for (int k = 0; k < 3; k++) {
        if (p[k] < bmin[k]) {
            double d = bmin[k] - p[k];
            d2 += d * d;
        } else if (p[k] > bmax[k]) {
            double d = p[k] - bmax[k];
            d2 += d * d;
        }
    }
    return d2;
}

int mesh::closestPointOnFace(int f, const Eigen::Vector3d& v, Eigen::Vector3d& v_proj,
                             vertProjData& candidateProj, bool snap) const {
    if (f < 0 || f >= F.size() || !F[f].active) {
        return -1;
    }
    std::vector<int> fVerts = faceAdjVertIdxs(f);
    const int fSize = static_cast<int>(fVerts.size());
    if (fSize < 3) {    // Invalid face
        return -1;
    }

    std::vector<Eigen::Vector3d> fVertsPos = adjVerts(fVerts);
    double snapTol = snap ? 1e-5 * bboxDiag : 1e-6 * meanE;

    int projType, projIdx;

    if (fSize == 3) {
        v_proj = Utils::triangleClosestPoint(fVertsPos, v, projType, projIdx, snapTol);
    } else if (fSize == 4) {
        v_proj = Utils::bilinearPatchClosestPoint(fVertsPos, v, projType, projIdx, snapTol);
    } else {
        v_proj = Utils::polygonClosestPointNewell(fVertsPos, v, F[f].n, projType, projIdx, snapTol);
    }

    if (projType == 0) {
        int local_v = projIdx;
        if (local_v < 0 || local_v >= fSize) {
            return -1;
        }
        candidateProj.elType = 0;
        candidateProj.elIdx = fVerts[local_v];
        return 1;
    }

    if (projType == 1) {
        int local_e = projIdx;
        if (local_e < 0 || local_e >= fSize) {
            return -1;
        }
        int local_next = (local_e + 1) % fSize;
        auto it = vertPairToHE.find({fVerts[local_e], fVerts[local_next]});
        if (it == vertPairToHE.end()) {
            return -1;
        }
        candidateProj.elType = 1;
        candidateProj.elIdx = HE[it->second].edge;
        return 1;
    }

    if (projType == 2) {
        candidateProj.elType = 2;
        candidateProj.elIdx = f;
        return 1;
    }

    return -1;
}

// Closest face test for BVH 
int mesh::closestFaceBVH(const Eigen::Vector3d& v, Eigen::Vector3d& proj, vertProjData& projData, bool snap) const {
    if (BVH.empty()) {
        return -1;
    }
    // Initialize closest distance and projData
    double bestDist2 = std::numeric_limits<double>::infinity();
    projData = {-1, -1};
    // Init priority queue with best nodes
    std::priority_queue<BVHQueueEntry> q;
    q.push({0, pointAABBDist2(v, 0)});
    // Pop bbox from queue until empty (i.e., found closest face)
    while (!q.empty()) {
        BVHQueueEntry curr = q.top();
        q.pop();
        // Since queue is ordered, if the closest remaining box cannot improve,
        // no later box can improve either
        if (curr.dist2 >= bestDist2) {
            break;
        }

        const AABB& box = BVH[curr.box];
        if (box.leaf) {     // If we are at a leaf, then we need to test distances to the participating faces
            for (int i = 0; i < box.faces.size(); i++) {
                int f = box.faces[i];
                Eigen::Vector3d candidateProj;
                vertProjData candidateData({2, f});
                if (closestPointOnFace(f, v, candidateProj, candidateData, snap) != 1) {
                    continue;
                }
                double d2 = (v - candidateProj).squaredNorm();
                if (d2 < bestDist2) {
                    bestDist2 = d2;
                    proj = candidateProj;
                    projData = candidateData;
                }
            }
        } else {    // Otherwise, trace down to lower depth bbox
            for (int i = 0; i < box.children.size(); i++) {
                int child = box.children[i];
                double childDist2 = pointAABBDist2(v, child);
                if (childDist2 < bestDist2) {   // If we have a closer candidate, then push onto queue
                    q.push({child, childDist2});
                }
            }
        }
    }

    if (projData.elIdx < 0) {   // Invalid projData
        return -1;
    }

    return 1;
}

// Alternative to vertex projection that also computes additional bind data
int mesh::computeVBinding(const Eigen::Vector3d& p, meshBindData& bind, bool snap, bool fast) const {
    bind = meshBindData();
    Eigen::Vector3d proj;
    vertProjData projData = computeVProjection(p, proj, snap, fast);

    if (projData.elIdx < 0) {
        return -1;
    }

    bind.elType = projData.elType;
    bind.elIdx = projData.elIdx;
    bind.proj = proj;
    bind.offset = p - proj;

    if (computeBindCoords(bind.elType, bind.elIdx, proj, bind.coords) != 1) {
        return -1;
    }
    if (computeBindFrame(bind.elType, bind.elIdx, bind.restFrame) != 1) {
        return -1;
    }

    return 1;
}

int mesh::computeBindCoords(int elType, int elIdx, const Eigen::Vector3d& proj, Eigen::VectorXd& coords) const {
    coords.resize(0);
    const double eps = 1e-12;

    if (elType == 2) {  // Landed on face; need MVC
        int f = elIdx;
        if (f < 0 || f >= F.size() || !F[f].active) {
            return -1;
        }
        std::vector<int> fVerts = faceAdjVertIdxs(f);
        int fSize = fVerts.size();
        if (fSize < 3) {    // Invalid face
            return -1;
        }
        std::vector<Eigen::Vector3d> fVertsPos = adjVerts(fVerts);
        coords.resize(fSize);
        coords.setZero();

        if (fSize == 3) {   // Triangle. Get Barycentric coords
            const Eigen::Vector3d& a = fVertsPos[0];
            const Eigen::Vector3d& b = fVertsPos[1];
            const Eigen::Vector3d& c = fVertsPos[2];

            Eigen::Vector3d v0 = b - a;
            Eigen::Vector3d v1 = c - a;
            Eigen::Vector3d v2 = proj - a;

            double d00 = v0.dot(v0);
            double d01 = v0.dot(v1);
            double d11 = v1.dot(v1);
            double d20 = v2.dot(v0);
            double d21 = v2.dot(v1);

            double denom = d00 * d11 - d01 * d01;
            if (std::abs(denom) <= eps) {
                return -1;
            }
            double v = (d11 * d20 - d01 * d21) / denom;
            double w = (d00 * d21 - d01 * d20) / denom;
            double u = 1.0 - v - w;
            coords(0) = u;
            coords(1) = v;
            coords(2) = w;
            return 1;
        } else if (fSize == 4) {
            double u, v;
            Utils::bilinearPatchClosestPoint(fVertsPos, proj, u, v, 0.0);
            coords.resize(2);
            coords(0) = u;
            coords(1) = v;
            return 1;
        }

        // 5+ polygons: MVC on the face/Newell plane.
        Eigen::Vector3d n = F[f].n;
        if (n.squaredNorm() <= eps) {
            return -1;
        }
        n.normalize();

        Eigen::Vector3d barycenter = DECUtils::computeBarycenter(fVertsPos);

        Eigen::Vector3d t0;
        Eigen::Vector3d t1;
        Utils::buildPlaneBasis(n, t0, t1);

        std::vector<Eigen::Vector2d> fVerts2D(fSize);

        for (int i = 0; i < fSize; i++) {
            Eigen::Vector3d vi_proj =
                Utils::projectPointOntoPlane(n, barycenter, fVertsPos[i]);
            fVerts2D[i] = Utils::convertTo2D(vi_proj, barycenter, t0, t1);
        }

        Eigen::Vector3d projPlane = Utils::projectPointOntoPlane(n, barycenter, proj);
        Eigen::Vector2d proj2D = Utils::convertTo2D(projPlane, barycenter, t0, t1);

        Utils::meanValueCoordinates(proj2D, fVerts2D, coords);

        return 1;
    }

    if (elType == 1) {  // Edge case. Simply compute t-value
        int e = elIdx;
        if (e < 0 || e >= E.size() || !E[e].active) {
            return -1;
        }
        int he = E[e].he;
        int v1 = HE[he].dest;
        int v0 = HE[HE[he].twin].dest;

        Eigen::Vector3d a = V[v0].pos;
        Eigen::Vector3d b = V[v1].pos;
        Eigen::Vector3d ab = b - a;

        double denom = ab.squaredNorm();
        if (denom <= eps) {
            return -1;
        }

        double t = (proj - a).dot(ab) / denom;
        t = std::clamp(t, 0.0, 1.0);    // Safety
        coords.resize(1);
        coords(0) = t;
        return 1;
    }

    if (elType == 0) {  // Vertex. Do nothing
        int v = elIdx;
        if (v < 0 || v >= V.size() || !V[v].active) {
            return -1;
        }
        coords.resize(0);
        return 1;
    }
    return -1;
}

// Compute the frame at the bind point
int mesh::computeBindFrame(int elType, int elIdx, Eigen::Matrix3d& frame) const {
    const double eps = 1e-12;
    frame.setIdentity();
    Eigen::Vector3d t;
    Eigen::Vector3d b;
    Eigen::Vector3d n;

    if (elType == 2) {  // Face
        int f = elIdx;
        int he = F[f].he;
        n = F[f].n;
        t = V[HE[he].dest].pos - V[HE[HE[he].twin].dest].pos;
    } else if (elType == 1) {   // Edge
        int e = elIdx;
        int he = E[e].he;
        n = E[e].n;
        t = V[HE[he].dest].pos - V[HE[HE[he].twin].dest].pos;
    } else if (elType == 0) {
        int v = elIdx;
        n = V[v].n;
        int he = V[v].he;

        if (he >= 0 && he < HE.size()) {
            t = V[HE[he].dest].pos - V[v].pos;
        }

        // Project onto tangent plane
        t = t - t.dot(n) * n;
        // Fallback if we get a degenerate tangent
        if (t.squaredNorm() <= eps) {
            std::vector<int> adjHE = vertAdjHEs(v);
            bool found = false;
            for (int i = 0; i < adjHE.size(); i++) {
                int currHE = adjHE[i];
                Eigen::Vector3d candidate = V[HE[currHE].dest].pos - V[v].pos;
                candidate = candidate - candidate.dot(n) * n;
                if (candidate.squaredNorm() > eps) {
                    t = candidate;
                    found = true;
                    break;
                }
            }
            if (!found) {
                return -1;
            }
        }
    } else {
        return -1;
    }
    if (n.squaredNorm() <= eps) {
        return -1;
    }
    n.normalize();
    t = t - t.dot(n) * n;
    if (t.squaredNorm() <= eps) {
        Utils::buildPlaneBasis(n, t, b);
    } else {
        t.normalize();
        b = n.cross(t);
    }
    frame.col(0) = t;
    frame.col(1) = b;
    frame.col(2) = n;
    return 1;
}

// Recover a point location from stored coordinates
int mesh::recoverCoords(int elType, int elIdx, const Eigen::VectorXd& coords, Eigen::Vector3d& p, bool fast) {
    if (elType == 0) {  // verts
        if (elIdx < 0 || elIdx >= V.size()) {
            return -1;
        } else {
            p = V[elIdx].pos;
        }
    } else if (elType == 1) {  // edges
        if ((elIdx < 0 || elIdx >= E.size()) || 
            coords.size() != 1 || (coords(0) < 0.0 || coords(1) > 1.0)) {
            return -1;
        } else {
            Eigen::Vector3d v1 = V[HE[E[elIdx].he].dest].pos;
            Eigen::Vector3d v0 = V[HE[HE[E[elIdx].he].twin].dest].pos;
            p = v0 + coords(0) * (v1 - v0);
        }
    } else if (elType == 2) {
        if (elIdx < 0 || elIdx >= F.size()) {
            return -1;
        }
        std::vector<Eigen::Vector3d> adjVerts = faceAdjVerts(elIdx);
        if (coords.size() != adjVerts.size() || adjVerts.size() <= 2) {
            return -1;
        }
        // Special exception for quads
        if (adjVerts.size() == 4) {
            p = Utils::bilinearPatch(adjVerts, coords(0), coords(1));
        } else {    // Classic coordinates
            // Stack
            Eigen::MatrixXd adjVertsStack = (DECUtils::posOp(adjVerts)).transpose();
            p = adjVertsStack * coords;
        }
    } else {
        return -1;
    }
    return 1;
}

}   // namespace Mesh