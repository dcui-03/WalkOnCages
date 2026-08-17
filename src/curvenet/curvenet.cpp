#define _USE_MATH_DEFINES
#include "curvenet.hpp"

#include "utils/utils.hpp"
#include "mesh/mesh.hpp"
#include "dcurvenet/dcurvenet.hpp"
#include <Eigen/Core>
#include <cmath>
#include <array>
#include <algorithm>
#include <utility>
#include <map>
#include <stdexcept>


namespace Curvenet {
    // Initialize from an existing list of controls, splines
    curvenet::curvenet(const std::vector<Eigen::Vector3d>& Controls, const std::vector<Eigen::Vector3d>& Tangents, const std::vector<std::array<int, 4>>& Splines, const Mesh::mesh* M, int alpha): alpha(alpha) {
        for (int c = 0; c < Controls.size(); c++) {
            int new_c = addControl(Controls[c]);
            inputCtoC[c] = new_c;
        }
        for (int s = 0; s < Splines.size(); s++) {
            std::array<int, 4> S = Splines[s];
            if (S[1] < 0 || S[1] >= Tangents.size() || S[2] < 0 || S[2] >= Tangents.size() ||
                S[0] < 0 || S[0] >= Controls.size() || S[3] < 0 || S[3] >= Controls.size()) {
                throw std::runtime_error("Failed to initialize curve network.");
            }
            int c0 = inputCtoC.at(S[0]);
            int c1 = inputCtoC.at(S[3]);
            std::pair<int, int> he = addSpline(c0, c1, Tangents[S[1]], Tangents[S[2]], M);
            inputTtoHE[S[1]] = he.first;
            inputTtoHE[S[2]] = he.second;
        }
        // Get projection data
        if (M) {
            ctrlProjDataFromMesh(*M);
            tanProjDataFromMesh(*M);
            sortAdjHEAll();
        }
        assignCtrlTypeAll();
        if (traceCurves() == -1) {
            throw std::runtime_error("Failed to trace curve network.");
        }
        computeBBoxDiag();
        return;
    }

    // Update curve network with new positions
    // No topological changes, so just update positions
    void curvenet::updateCurveNet(std::vector<Eigen::Vector3d> Controls, std::vector<Eigen::Vector3d> Tangents) {
        for (auto& idxPair : inputCtoC) {
            C[idxPair.second].new_pos = Controls[idxPair.first];
        }
        for (auto& idxPair : inputTtoHE) {
            HE[idxPair.second].tan = Tangents[idxPair.first];
        }
        return;
    }

    // Assign weight to a control
    int curvenet::assignWeight(int c, bool fixed_w, double w) {
        if (c < 0 || c >= C.size() || w < 0.0 || w > 1.0) {
            return -1;
        }
        C[c].w = w;
        C[c].fixed_w = fixed_w;
        if (fixed_w == false) {
            C[c].w = 1.0;
        }
        return 1;
    }
    void curvenet::resetWeights() {
        for (int c = 0; c < C.size(); c++) {
            C[c].w = 1.0;
            C[c].fixed_w = true;
        }
        return;
    }
    // empty initializer
    curvenet::curvenet() {
        C.clear();
        HE.clear();
        S.clear();
        vertData.clear();
        tanData.clear();
    }

    // Create new control
    int curvenet::addControl(Eigen::Vector3d pos) {
        int c = C.size();
        C.emplace_back();
        C[c].pos = pos;
        C[c].new_pos = pos;
        vertData.emplace_back();
        return c;
    }
    // Add a spline to the spline list given indices of the points
    std::pair<int, int> curvenet::addSpline(int start, int end, Eigen::Vector3d t0, Eigen::Vector3d t1, const Mesh::mesh* M) {
        // Create 2 new halfedges and a new spline
        int he0 = HE.size();
        int he1 = he0+1;
        HE.emplace_back();
        HE.emplace_back();
        tanData.emplace_back();
        tanData.emplace_back();
        int s = S.size();
        S.emplace_back();

        // Re-wire
        HE[he0].twin = he1;
        HE[he1].twin = he0;
        HE[he0].origin = start;
        HE[he1].origin = end;
        HE[he0].tan = t0;
        HE[he1].tan = t1;
        HE[he0].rest_tan = t0;
        HE[he1].rest_tan = t1;
        HE[he0].s = s;
        HE[he1].s = s;
        S[s].he = he0;
        // Determine sampling
        // NOTE: Setting sampling for the estimate to 75 for now
        S[s].num_samples = computeNumSamples(arclenEst(s, 75), M);
        // Insert spline into vertex list
        C[start].adjHE.push_back(he0);
        C[end].adjHE.push_back(he1);
        return std::make_pair(he0, he1);
    }

    // Change control normal
    int curvenet::editControlN(int c, Eigen::Vector3d normal) {
        // Catch degenerate cases
        if ((!C[c].active) || normal.squaredNorm() == 0.0) {
            return -1;
        }
        C[c].n = normal.normalized();   // Always normalize for safety
        return 1;
    }

    // Compute normals for each vertex by projecting onto a mesh
    int curvenet::ctrlProjDataFromMesh(const Mesh::mesh& m) {
        #pragma omp parallel for
        for (int c = 0; c < C.size(); c++) {
            if (!C[c].active) {
                continue;
            }
            Utils::frameData bindData;
            if (m.computeVBinding(C[c].pos, bindData) != 1) {
                throw std::runtime_error("curvenet::ctrlProjDataFromMesh(): invalid bind data");
            }
            Eigen::Vector3d n = m.getNormal(bindData.proj.elType, bindData.proj.elIdx);
            // Copy over data
            editControlN(c, n);
            vertData[c] = bindData;
        }
        setMesh = true;
        return 1;
    }

    // Compute normals for each vertex by projecting onto a mesh
    int curvenet::tanProjDataFromMesh(const Mesh::mesh& m) {
        #pragma omp parallel for
        for (int he = 0; he < HE.size(); he++) {
            if (!HE[he].active) {
                continue;
            }
            Utils::frameData bindData;
            if (m.computeVBinding(HE[he].rest_tan, bindData) != 1) {
                throw std::runtime_error("curvenet::tanProjDataFromMesh(): invalid bind data");
            }
            // Copy over data
            tanData[he] = bindData;
        }
        setMesh = true;
        return 1;
    }

    // Sort adjacent tangent vectors to a control point
    int curvenet::sortAdjHE(int c) {
        if ((!C[c].active) || (C[c].n == Eigen::Vector3d::Zero())) {
            return -1;
        }

        std::vector<int> adjHE = C[c].adjHE;
        std::vector<Eigen::Vector3d> adjT = ctrlAdjTans(c);

        // Compute a local angle for each adjacent halfedge
        std::vector<double> adjTheta;
        Eigen::Vector3d t1, t2;
        Utils::buildPlaneBasis(C[c].n, t1, t2);
        for (int t = 0; t < adjT.size(); t++) {
            double theta;
            bool success = Utils::directionAngleInPlane(C[c].pos, adjT[t], C[c].n, t1, t2, theta);
            // Insert theta into the local adjacency list
            if (!success) { // If degenerate, just give it a big angle so that it gets sorted to the end
                theta = 3.0*M_PI;
            }
            adjTheta.push_back(theta);
        }

        // Sort by angle
        Utils::doubleListIdxSort(adjTheta, adjHE);
        // If a curve connects to itself and is the only one, then do not rewire.
        if (adjTheta.size() == 2 && (HE[adjHE[0]].twin == adjHE[1])) {
            C[c].adjHE = adjHE;
            C[c].sorted = true;
            return 1;   // success
        }
        // Do the necessary re-wiring
        for (int he_idx = 0; he_idx < adjHE.size(); he_idx++) {
            int he = adjHE[he_idx];
            int he_next = adjHE[(he_idx + 1)%adjHE.size()];
            HE[he].prev = HE[he_next].twin;
            HE[HE[he_next].twin].next = he;
        }
        // Assign back to control
        C[c].adjHE = adjHE;
        C[c].sorted = true;
        return 1;   // success
    }
    int curvenet::sortAdjHEAll() {
        if (!setMesh) {
            return -1;
        }
        for (int c = 0; c < C.size(); c++) {
            if (C[c].active && !C[c].sorted) {
                if (sortAdjHE(c) == -1) {
                    return -1;
                }
            }
        }
        return 1;
    }
    
    // Compute what kind of vertex each control is using the valence of splines
    int curvenet::assignCtrlType(int c) {
        std::vector<int> adjHE = C[c].adjHE;
        // NOTE: Treat self-loops as not anchors
        C[c].cType = std::min(static_cast<int>(adjHE.size()), 3);
        return C[c].cType;
    }
    int curvenet::assignCtrlTypeAll() {
        for (int c = 0; c < C.size(); c++) {
            if (C[c].active) {
                assignCtrlType(c);
            }
        }
        return 1;
    }
    // Trace spline curves
    int curvenet::traceCurves() {
        Crv.clear();
        std::vector<bool> splineFound(S.size(), false);
        // Search from intersection and anchor controls
        for (int c = 0; c < C.size(); c++) {
            int cType = C[c].cType;
            if (cType == 2) {  // Skip all loops for now
                    continue;
            }
            const std::vector<int> adjHE = C[c].adjHE;
            // Start from each outgoing halfedge and trace until we hit a stop point
            for (int he = 0; he < adjHE.size(); he++) {
                int s0 = HE[adjHE[he]].s;
                if (splineFound[s0]) {   // Skip any splines we've already seen
                    continue;
                }

                // Otherwise, need to do a spline trace
                // Create new curve object
                int crv = Crv.size();
                Crv.emplace_back();
                bool curveEnd = false;
                int counter = 0;
                std::vector<int> splineList;
                int curr_he = adjHE[he];
                // Trace splines until we hit an intersection or an anchor
                do {
                    int s = HE[curr_he].s;
                    // If current spline is not oriented with the curve, then flip it.
                    if (S[s].he != curr_he) {
                        S[s].he = curr_he;
                    }
                    int curr_end = HE[HE[curr_he].twin].origin;
                    // Add spline to list
                    splineList.push_back(s);
                    S[s].curve = crv;
                    splineFound[s] = true;
                    counter++;
                    // Go to next spline
                    if (C[curr_end].cType == 1 || C[curr_end].cType == 3) {
                        curveEnd = true;
                    } else {    // Must be 2 outgoing HE's from this one. Pick the one we haven't gone to yet
                        int next_he = nextHEFromControl(curr_he, curr_end);
                        if (next_he < 0) {
                            return -1;
                        }
                        curr_he = next_he;
                    }
                } while (!curveEnd && counter < S.size());
                if (!curveEnd) {    // Error check
                    return -1;
                }

                Crv[crv].splines = splineList;
            }
        }

        // Start from remaining splines, which must form closed loops
        for (int c = 0; c < C.size(); c++) {
            int cType = C[c].cType;
            if (cType == 1 || cType == 3) {  // Skip anchors and intersections
                    continue;
            }
            const std::vector<int> adjHE = C[c].adjHE;
            // Start from each outgoing halfedge and trace until we hit the start point again
            for (int he = 0; he < adjHE.size(); he++) {
                int s0 = HE[adjHE[he]].s;
                if (splineFound[s0]) {   // Skip any splines we've already seen
                    continue;
                }

                // Otherwise, need to do a spline trace
                // Create new curve object
                int crv = Crv.size();
                Crv.emplace_back();
                bool curveEnd = false;
                int counter = 0;
                std::vector<int> splineList;
                int curr_he = adjHE[he];
                // Trace splines until we hit the start again
                do {
                    int s = HE[curr_he].s;
                    // If current spline is not oriented with the curve, then flip it.
                    if (S[s].he != curr_he) {
                        S[s].he = curr_he;
                    }
                    int curr_end = HE[HE[curr_he].twin].origin;
                    // Add spline to list
                    splineList.push_back(s);
                    S[s].curve = crv;
                    splineFound[s] = true;
                    counter++;
                    // Go to next spline
                    if (curr_end == c) {    // Hit the start point again
                        curveEnd = true;
                    } else {    // Must be 2 outgoing HE's from this one. Pick the one we haven't gone to yet
                        int next_he = nextHEFromControl(curr_he, curr_end);
                        if (next_he < 0) {
                            return -1;
                        }
                        curr_he = next_he;
                    }
                } while (!curveEnd && counter < S.size());
                if (!curveEnd) {    // Error check
                    return -1;
                }

                Crv[crv].splines = splineList;
            }
        }
        return 1;
    }

    // Helper for computeCurves
    // Compute the next halfedge from a degree 2 control
    int curvenet::nextHEFromControl(int curr_he, int curr_end) {
        if (curr_end < 0 || curr_end >= C.size()) {
            return -1;
        }
        if (C[curr_end].cType != 2) {
            return -1;
        }
        const std::vector<int>& adjHE = C[curr_end].adjHE;
        if (adjHE.size() != 2) {
            return -1;
        }
        int back_he = HE[curr_he].twin;
        if (adjHE[0] == back_he) {
            return adjHE[1];
        }
        if (adjHE[1] == back_he) {
            return adjHE[0];
        }
        // The halfedge we arrived on does not actually end at this control according to the control's adjacency list
        return -1;
    }

    // Controls and tangents as a matrix
    int curvenet::CTasMatrix(Eigen::MatrixXd& Verts) {
        Verts.resize(C.size() + HE.size(), 3);

        for (int c = 0; c < C.size(); c++) {
            Verts.row(c) = C[c].new_pos.transpose();
        }
        for (int t = 0; t < HE.size(); t++) {
            Verts.row(C.size() + t) = HE[t].tan.transpose();
        }
        return 1;
    }

    // Find the closest point on the curve network to p, via dCN's polyline BVH
    int curvenet::closestPoint(const Eigen::Vector3d& p, const Polynet::dcurvenet* dCN, Utils::projData& bind, bool snap, double snapTol) const {
        if (dCN == nullptr) {
            return -1;
        }
        if (dCN->closestPoint(p, bind, snap, snapTol) != 1) {
            return -1;
        }
        // Optimize t based on the locally guess
        bind.coords(0) = optimizeT(bind.coords(0), bind.elIdx, p);
        // Compute the explicit position
        bind.pos = tSampleBezier(bind.elIdx, bind.coords(0));
        return 1;
    }

    // Newton iterations to refine an initial guess t-value to get true closest point
    double curvenet::optimizeT(double t, int s, const Eigen::Vector3d& p, int max_iter) const {
        t = std::clamp(t, 0.0, 1.0);
        double curr_t = t;
        for (int i = 0; i < max_iter; i++) {
            Eigen::Vector3d B = tSampleBezier(s, t);
            Eigen::Vector3d B_ = tBezier_first(s, t);
            Eigen::Vector3d B__ = tBezier_second(s, t);
            Eigen::Vector3d diff = B - p;
            double grad = diff.dot(B_);
            double hessian = B_.squaredNorm() + diff.dot(B__);

            // refine t
            double dt = grad / hessian;
            curr_t -= dt;
            curr_t = std::clamp(curr_t, 0.0, 1.0);
            if (dt <= 1e-6) {  // convergence check
                return curr_t;
            }
        }
        return curr_t;
    }

    // Cast a ray against dCN's polyline BVH, then Newton-refine each hit onto its true spline
    int curvenet::raycast(const Eigen::Vector3d& origin, const Eigen::Vector3d& direc, const Polynet::dcurvenet* dCN, std::vector<Utils::projData>& hits, double tol) const {
        if (dCN == nullptr) {
            return -1;
        }
        if (dCN->raycast(origin, direc, hits, tol) != 1) {
            return -1;
        }
        for (Utils::projData& hit : hits) {
            int s = hit.elIdx;
            hit.coords(0) = optimizeT(hit.coords(0), s, hit.pos);
            hit.pos = tSampleBezier(s, hit.coords(0));
        }
        return 1;
    }

    // Evaluate basis functions on a spline
    int curvenet::evaluateBasis(int s, double t, std::vector<std::pair<int, double>>& basis) {
        if (s < 0 || s >= S.size()) {
            return -1;
        }
        basis.clear();
        t = std::clamp(t, 0.0, 1.0);
        const double eps = 1e-9;
        int t0 = S[s].he;
        int c0 = HE[t0].origin;
        int t1 = HE[t0].twin;
        int c1 = HE[t1].origin;
        if (t <= eps) { // Start of spline
            basis.push_back({c0, 1.0});
        } else if (t >= 1.0 - eps) {  // End of spline
            basis.push_back({c1, 1.0});
        } else {    // Somewhere in the middle
            basis.resize(4);
            // Fill in with correct bases
            // NOTE: The interior "tangent" vertices have a global C + HE indexing here...
            basis[0] = {c0, std::pow(1.0 - t, 3)};
            basis[1] = {t0 + C.size(), 3.0 * std::pow(1.0 - t, 2) * t};
            basis[2] = {t1 + C.size(), 3.0 * (1.0 - t) * t * t};
            basis[3] = {c1, std::pow(t, 3)};
        }
        return 1;
    }

    // Compute the length of the diagonal of the bounding box
    void curvenet::computeBBoxDiag() {
        bool found = false;
        Eigen::Vector3d minV;
        Eigen::Vector3d maxV;
        for (int c = 0; c < C.size(); c++) {
            if (!C[c].active) {
                continue;
            }
            const Eigen::Vector3d& p = C[c].new_pos;
            if (!found) {
                minV = p;
                maxV = p;
                found = true;
            } else {
                minV = minV.cwiseMin(p);
                maxV = maxV.cwiseMax(p);
            }
        }
        bboxDiag = found ? (maxV - minV).norm() : 0.0;
        return;
    }

    double curvenet::getBBoxDiag() const {
        return bboxDiag;
    }
}   // namespace Curvenet