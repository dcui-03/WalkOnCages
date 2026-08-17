#include "pscurvenet.hpp"

#include "pscurvenet_types.hpp"
#include "utils/utils.hpp"
#include <Eigen/Core>
#include <vector>
#include <array>
#include <map>

namespace psCurvenet {

// Constructor takes four points [start, tangent 1, tangent 2, end], and associated normals
// NOTE: Constructor assumes you already have no duplicates in your inputs
// Empty constructor
pscurvenet::pscurvenet() {
    resetCurvenet();
}

// --------- UPDATE CURVENET -----------
// Hard reset everything
void pscurvenet::resetCurvenet() {
    C.clear();
    S.clear();
    psTangentToS.clear();
    recomputeMap = false;
}

// Update control position
void pscurvenet::updateControlPos(int c, Eigen::Vector3d new_pos, bool project) {
    if (c >= C.size()) {
        return;
    }
    // Recompute tangent direction
    Eigen::Vector3d oldPos = C[c].pos;
    C[c].pos = new_pos;
    for (int s = 0; s < S.size(); s++) {
        if (S[s].start == c) {
            Eigen::Vector3d new_t0 = new_pos + (S[s].t0 - oldPos);
            updateTangentPos(s, true, new_t0, project);
        } else if (S[s].end == c) {
            Eigen::Vector3d new_t1 = new_pos + (S[s].t1 - oldPos);
            updateTangentPos(s, false, new_t1, project);
        }
    }
    return;
}

void pscurvenet::updateControlNormal(int c, const Eigen::Vector3d& new_normal, bool rotation, bool project) {
    if (c >= C.size()) {
        return;
    }
    const Eigen::Vector3d& old_normal = C[c].n;
    // Chain rotations
    Eigen::Matrix3d rot = Utils::computeRotation(C[c].n, new_normal);
    C[c].n = new_normal.normalized();
    // Recompute tangent directions
    for (int s = 0; s < S.size(); s++) {
        if (S[s].start == c) {
            if (rotation) {
                rotateTangentPos(s, true, rot);
            } else {
                updateTangentPos(s, true, S[s].t0, project);
            }
        } else if (S[s].end == c) {
            if (rotation) {
                rotateTangentPos(s, false, rot);
            } else {
                updateTangentPos(s, false, S[s].t1, project);
            }
        }
    }
    return;
}

// Update tangent pos given tangent network index
bool pscurvenet::updateTangentPos(int psT_idx, const Eigen::Vector3d& new_pos, bool project) {
    return updateTangentPos(psTangentToS[psT_idx].first, psTangentToS[psT_idx].second, new_pos, project);
}

// Update tangent position
bool pscurvenet::updateTangentPos(int s, bool t0, const Eigen::Vector3d& new_pos, bool project) {
    if (s >= S.size()) {
        return false;
    }
    double eps = 1e-8;
    // Prevent tangent movement if we are too close to an endpoint
    if (((C[S[s].start].pos - new_pos).norm() <= eps) || ((C[S[s].end].pos - new_pos).norm() <= eps)) {
        return false;
    }

    if (project) {
        Eigen::Vector3d new_direc;
        if (t0) {
            new_direc = new_pos - C[S[s].start].pos;
            Eigen::Vector3d proj_direc = new_direc;
            Utils::projectVectorOntoTangentPlane(C[S[s].start].n, new_direc, proj_direc, new_direc.norm());
            S[s].t0 = C[S[s].start].pos + proj_direc;
        } else {
            new_direc = new_pos - C[S[s].end].pos;
            Eigen::Vector3d proj_direc = new_direc;
            Utils::projectVectorOntoTangentPlane(C[S[s].end].n, new_direc, proj_direc, new_direc.norm());
            S[s].t1 = C[S[s].end].pos + proj_direc;
        }
    } else {
        if (t0) {
            S[s].t0 = new_pos;
        } else {
            S[s].t1 = new_pos;
        }
    }

    return true;
}

bool pscurvenet::rotateTangentPos(int s, bool t0, Eigen::Matrix3d rotation) {
    if (s >= S.size()) {
        return false;
    }
    double eps = 1e-6;
    Eigen::Vector3d pos;
    if (t0) {
        pos = S[s].t0;
    } else {
        pos = S[s].t1;
    }
    Eigen::Vector3d new_direc;
    if (t0) {
        new_direc = rotation * (pos - C[S[s].start].pos);
        S[s].t0 = C[S[s].start].pos + new_direc;
    } else {
        new_direc = rotation * (pos - C[S[s].end].pos);
        S[s].t1 = C[S[s].end].pos + new_direc;
    }

    return true;
}

// Add a control and return its index. offset shifts pos along normal before storing
int pscurvenet::addControl(Eigen::Vector3d pos, Eigen::Vector3d normal, double offset) {
    recomputeMap = true;
    int c = C.size();
    C.emplace_back();
    C[c].pos = pos + offset * normal.normalized();
    C[c].n = normal;
    return c;
}

// Add a spline given only the start and end. Estimate t0 and t1 from these
int pscurvenet::addSpline(int c0, int c1, double init_factor) {
    init_factor = std::max(1.0, init_factor);
    double eps = 1e-6;
    Eigen::Vector3d t0 = C[c1].pos - C[c0].pos;
    Eigen::Vector3d t1 = -1 * t0;
    double dist = t0.norm();
    if (c0 == c1 || dist <= eps) {
        Utils::buildPlaneBasis(C[c0].n, t0, t1);
        return addSpline(c0, C[c0].pos + t0, C[c1].pos + t1, c1);
    }
    // Project each onto local tangent plane
    Eigen::Vector3d t0_proj;
    Eigen::Vector3d t1_proj;
    double t0_norm = Utils::projectVectorOntoTangentPlane(C[c0].n, t0, t0_proj, dist / init_factor);
    double t1_norm = Utils::projectVectorOntoTangentPlane(C[c1].n, t1, t1_proj, dist / init_factor);
    // If projection is ill-posed, then just pick a random orth direc
    if (t0_norm <= eps) {
        Eigen::Vector3d t0_temp;
        Utils::buildPlaneBasis(C[c0].n, t0_proj, t0_temp);
        t0_proj *= dist / init_factor;
    }
    if (t1_norm <= eps) {
        Eigen::Vector3d t1_temp;
        Utils::buildPlaneBasis(C[c0].n, t1_proj, t1_temp);
        // Prefer that t1 is roughly along the same direc as t0
        if (t0_proj.dot(t1_proj) < 0) {
            t1_proj *= -1.0;
        }
        t1_proj *= dist / init_factor;
    }
    // Add in the new spline
    return addSpline(c0, C[c0].pos + t0_proj, C[c1].pos + t1_proj, c1);
}
// Add a spline and return its index
int pscurvenet::addSpline(int c0, Eigen::Vector3d t0_pos, Eigen::Vector3d t1_pos, int c1) {
    if (c0 >= C.size() || c1 >= C.size()) {
        return -1;
    } else if (!C[c0].active || !C[c1].active) {
        return -1;
    }
    int s = S.size();
    S.emplace_back();
    S[s].start = c0;
    S[s].end = c1;
    S[s].t0 = t0_pos;
    S[s].t1 = t1_pos;
    recomputeMap = true;
    return s;
}
// Remove a control.
int pscurvenet::removeControl(int c) {
    if (c >= C.size()) {
        return -1;
    }
    // Remove control from list
    C.erase(C.begin() + c);
    // All controls c and above are now relabeled
    // Hard reset the splines
    std::vector<Spline> newS;
    newS.reserve(S.size());
    for (int s = 0; s < static_cast<int>(S.size()); ++s) {
        int start = S[s].start;
        int end   = S[s].end;
        // Drop any splines that contained the control
        if (start == c || end == c) {
            continue;
        }
        // Make a copy and decrement any invalid indices
        Spline spline = S[s];
        if (spline.start > c) {
            spline.start -= 1;
        }
        if (spline.end > c) {
            spline.end -= 1;
        }
        newS.push_back(spline);
    }
    // Swap the old for the new
    S.swap(newS);
    recomputeMap = true;
    constructTangentMap();
    return 1;
}

int pscurvenet::removeSplineByTangent(int ps_TIdx) {
    constructTangentMap();
    auto it = psTangentToS.find(ps_TIdx);
    if (it == psTangentToS.end()) {
        return -1;
    }
    return removeSpline(it->second.first);
}

// Remove a spline
int pscurvenet::removeSpline(int s) {
    if (s >= S.size()) {
        return -1;
    }
    // Remove spline from list
    S.erase(S.begin() + s);
    recomputeMap = true;
    constructTangentMap();
    return 1;
}

// Clean up loose controls
int pscurvenet::cleanupControls() {
    // First, identify any loose controls
    std::vector<bool> activeC(C.size(), false);
    for (int s = 0; s < S.size(); s++) {
        activeC[S[s].start] = true;
        activeC[S[s].end] = true;
    }
    std::vector<int> looseC;
    for (int c = 0; c < C.size(); c++) {
        if (!activeC[c]) {
            looseC.push_back(c);
        }
    }
    // Remove loose controls
    for (int c = 0; c < looseC.size(); c++) {
        removeControl(looseC[c]);
        for (int c_new = c; c_new < looseC.size(); c_new++) {
            looseC[c_new]--;
        }
    }
    recomputeMap = true;
    return looseC.size();
}


// --------- GETTERS + POLYSCOPE CONVERSION -----------
// Get the normal at a control
Eigen::Vector3d pscurvenet::getNormal(int c) {
    return C[c].n;
}
// Control positions but returned as an Eigen::MatrixXd
void pscurvenet::cPosAsMatrix(Eigen::MatrixXd& cPos) {
    cPos.resize(C.size(), 3);
    for (int c = 0; c < C.size(); c++) {
        cPos.row(c) = C[c].pos.transpose();
    }
    return;
}

// Tangent positions but returned as an Eigen::MatrixXd
void pscurvenet::tPosAsMatrix(Eigen::MatrixXd& tPos) {
    tPos.resize(2*S.size(), 3);
    for (int s = 0; s < S.size(); s++) {
        tPos.row(2*s) = S[s].t0.transpose();
        tPos.row(2*s + 1) = S[s].t1.transpose();
    }
    constructTangentMap();
    return;
}
// Curve network as a chain of discrete splines
void pscurvenet::cnAsCurveNetwork(Eigen::MatrixXd& verts, std::vector<std::array<int, 2>>& connectivity) {
    std::vector<std::vector<Eigen::Vector3d>> dSplines(S.size());
    int num_verts = C.size();
    int num_edges = 0;
    // First, compute all splines as polylines
    for (int s = 0; s < S.size(); s++) {
        int num_samples = 30;   // Low sampling rate for polyscope
        dSplines[s] = sampleBezierNaive(s, num_samples);
        num_verts += num_samples - 2;   // Exclude start and end
        num_edges += num_samples - 1;
    }
    verts.resize(num_verts, 3);
    connectivity.resize(num_edges);
    // Push back all controls
    for (int c = 0; c < C.size(); c++) {
        verts.row(c) = C[c].pos.transpose();
    }
    // Push back all splines
    int currVert = C.size();
    int currEdge = 0;
    for (int s = 0; s < dSplines.size(); s++) {
        int start = S[s].start;
        int end = S[s].end;
        int currPrev = start;
        for (int v = 1; v < dSplines[s].size() - 1; v++) {
            int newVert = currVert;
            verts.row(newVert) = dSplines[s][v].transpose();

            connectivity[currEdge] = {currPrev, newVert};
            currPrev = newVert;

            currVert++;
            currEdge++;
        }
        // Add in the final edge
        connectivity[currEdge] = {currPrev, end};
        currEdge++;
    }
    return;
}
// Curve network that is just from controls to tangents as a curvenetwork
void pscurvenet::tansAsCurveNetwork(Eigen::MatrixXd& verts, std::vector<std::array<int, 2>>& connectivity) {
    verts.resize(C.size() + 2 * S.size(), 3);
    connectivity.resize(2 * S.size());
    // First, copy in all the controls
    for (int c = 0; c < C.size(); c++) {
        verts.row(c) = C[c].pos.transpose();
    }
    // Iterate over splines and add in verts/edges where needed
    for (int s = 0; s < S.size(); s++) {
        int currVert = C.size() + 2 * s;
        int currEdge = 2 * s;
        verts.row(currVert) = S[s].t0.transpose();
        verts.row(currVert + 1) = S[s].t1.transpose();
        connectivity[currEdge] = {S[s].start, currVert};
        connectivity[currEdge + 1] = {S[s].end, currVert + 1};
    }
    return;
}

void pscurvenet::cnAsStdVector(std::vector<Eigen::Vector3d>& controls, std::vector<Eigen::Vector3d>& tangents, std::vector<std::array<int, 4>>& splines) {
    controls.resize(C.size());
    tangents.resize(2 * S.size());
    splines.resize(S.size());

    for (int c = 0; c < C.size(); c++) {
        controls[c] = C[c].pos;
    }

    for (int s = 0; s < S.size(); s++) {
        tangents[2*s] = S[s].t0;
        tangents[2*s+1] = S[s].t1;
        splines[s] = std::array<int, 4>({S[s].start, 2*s, 2*s+1, S[s].end});
    }
    return;
}

// --------- SAMPLING -----------
// Sample a bezier curve at time t
Eigen::Vector3d pscurvenet::tSampleBezier(const Eigen::Vector3d& c0, const Eigen::Vector3d& c1, const Eigen::Vector3d& c2, const Eigen::Vector3d& c3, double t) {
    Eigen::Vector3d sample = std::pow(1 - t, 3) * c0 +
                        3 * std::pow(1 - t, 2) * t * c1 +
                        3 * (1 - t) * std::pow(t, 2) * c2 +
                        std::pow(t, 3) * c3;
    return sample;
}

Eigen::Vector3d pscurvenet::tSampleBezier(int s, double t) {
    return tSampleBezier(C[S[s].start].pos, S[s].t0, S[s].t1, C[S[s].end].pos, t);
}
// NOTE: This is a naive, fast sampler that uniformly samples t's. Re-implement if desired
// Returns n_samples points on the curve, including the endpoints
std::vector<Eigen::Vector3d> pscurvenet::sampleBezierNaive(int s, int n_samples) {
    if (n_samples < 2) {
        return {};
    }
    std::vector<Eigen::Vector3d> samples(n_samples);
    double h = 1.0/(n_samples - 1);
    double t = 0.0;

    samples[0] = C[S[s].start].pos;  // Start
    for (int i = 1; i < n_samples - 1; i++) {
        t += h;
        t = std::min(1.0, t);
        samples[i] = tSampleBezier(s, t);
    }
    samples[n_samples-1] = C[S[s].end].pos;   // End

    return samples;
}

// --------- ITERATORS -----------
// Get the adjacent splines to a control vertex
std::vector<int> pscurvenet::ctrlAdjSplines(int c) const {
    std::vector<int> adjS;
    if (c >= C.size()) {
        return adjS;
    }
    for (int s = 0; s < S.size(); s++) {
        if (S[s].start == c || S[s].end == c) {
            adjS.push_back(s);
        }
    }
    return adjS;
}

// Returns vertices adjacent to a spline
std::pair<int, int> pscurvenet::splineAdjCtrls(int s) const {
    if (s >= S.size()) {
        return std::make_pair(-1, -1);
    }
    return std::make_pair(S[s].start, S[s].end);
}

void pscurvenet::constructTangentMap() {
    if (!recomputeMap) {
        return;
    }
    psTangentToS.clear();
    int tanIdx = 0;
    for (int s = 0; s < S.size(); s++) {
        psTangentToS[tanIdx] = std::make_pair(s, true);
        psTangentToS[++tanIdx] = std::make_pair(s, false);
        tanIdx++;
    }
    recomputeMap = false;
    return;
}

}   // namespace psCurvenet
