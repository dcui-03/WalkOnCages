#include "dcurvenet.hpp"

#include "curvenet/curvenet.hpp"
#include "mesh/mesh.hpp"
#include <Eigen/Core>
#include <vector>
#include <algorithm>
#include <stdexcept>

namespace Polynet {
    // Takes the original curvenet and discretizes it
    dcurvenet::dcurvenet(Curvenet::curvenet* CN, const Mesh::mesh* M, bool sampleNaive): CN(CN), sampleNaive(sampleNaive) {
        const std::vector<Curvenet::Control>& cnCtrl = CN->controls();
        int num_curves = CN->numCurves();
        // Defensive reset
        V.clear();
        HE.clear();
        E.clear();
        C.clear();
        vertData.clear();
        curveCNIdx.clear();
        edgeSpline.clear();
        edgeSplineT.clear();
        // 1. First copy in the control points
        int num_controls = cnCtrl.size();
        std::vector<int> ctrlVerts(num_controls);
        for (int c = 0; c < num_controls; c++) {
            int new_v = addVert(cnCtrl[c], c);
            inputCtoV[c] = new_v;
            ctrlVerts[c] = new_v;
        }
        // 2. Initialize new dCN copies of the CN curves
        for (int crv = 0; crv < num_curves; crv++) {
            int c = addCurve(crv, sampleNaive);
            inputCrvToC[crv] = c;
        }

        // 3. Iterate over every control vertex that is an intersection/anchor and rewire its outgoing/incoming halfedges
        // NOTE: This is not really necessary. We already have the outgoing order of halfedges per control vertex, and we will never need to traverse betwen curves
        for (int v_idx = 0; v_idx < ctrlVerts.size(); v_idx++) {
            int v = ctrlVerts[v_idx];
            rewireVertAdjHE(v, CN->setMesh);
        }
        // Compute projection data
        if (M) {
            computeProjData(M);
            setMesh = true;
        }
        if (computeBVH() != 1) {
            throw std::runtime_error("Failed to build BVH.");
        }
    }

    // Update the new frames on all halfedges
    void dcurvenet::updateDiscCurveNet() {
        const std::vector<Curvenet::Control>& cnCtrl = CN->C;
        const std::vector<Curvenet::CubicSpline>& cnSpline = CN->S;
        const std::vector<Curvenet::Curve>& cnCurve = CN->Crv;
        // 1. Copy new control positions to their corresponding dvert
        for (const auto& idxPair : inputCtoV) {
            V[idxPair.second].new_pos = cnCtrl[idxPair.first].new_pos;
        }
        // 2. Iterate over curves and recompute
        // TODO: Can we parallelize this?
        // #pragma omp parallel for
        for (const auto& idxPair : inputCrvToC) {
            std::vector<int> splines = cnCurve[idxPair.first].splines;
            int c = idxPair.second;
            int he_curr = C[c].he_start;
            int v_prev = C[c].start;
            for (int s_idx = 0; s_idx < splines.size(); s_idx++) {
                int s = splines[s_idx];
                int n_samples = cnSpline[s].num_samples;
                // Assume n_samples will always be > 3
                const std::vector<Eigen::Vector3d> samples = sampleNaive ? CN->sampleBezierNaive(s, n_samples) : CN->unifSample(s, n_samples);
                // Exploit the fact that this matches the halfedge direction that the curve was constructed from
                for (int i = 1; i < samples.size(); i++) {
                    // for numerical reasons, only copy in non-controls
                    if (i != samples.size() - 1) {
                        V[HE[he_curr].dest].new_pos = samples[i];
                    }
                    int v_dest = HE[he_curr].dest;
                    v_prev = v_dest;
                    he_curr = HE[he_curr].next;
                }
            }
        }
        return;
    }

    dcurvenet::dcurvenet() {

    }

    // Add a new vertex that matches an existing control
    int dcurvenet::addVert(Curvenet::Control ctrl, int ctrl_idx) {
        int v = polynet::addVert(ctrl.pos, ctrl.n);
        V[v].adjHE.resize(ctrl.adjHE.size(), -1);
        vertData.emplace_back();
        vertData[v].cn_idx = ctrl_idx;
        vertData[v].cn_type = ctrl.cType;
        return v;
    }
    // Add a vertex given its parameters
    int dcurvenet::addVert(Eigen::Vector3d new_pos, Eigen::Vector3d new_n) {
        int v = polynet::addVert(new_pos, new_n);
        vertData.emplace_back();
        return v;
    }
    // Add a curve that matches an input curvenet curve
    int dcurvenet::addCurve(int crv, bool sampleNaive) {
        const std::vector<Curvenet::HalfEdge>& cnHE = CN->HE;
        const std::vector<Curvenet::CubicSpline>& cnSpline = CN->S;
        const std::vector<Curvenet::Curve>& cnCrv = CN->Crv;
        if (crv < 0 || crv >= cnCrv.size()) {
            return -1;
        }
        const std::vector<int>& crvSplines = cnCrv[crv].splines;
        if (crvSplines.empty()) {
            return -1;
        }
        int c = C.size();
        C.emplace_back();
        curveCNIdx.push_back(crv);

        // Track the previous positive halfedge and previous negative halfedge
        // so that each newly added segment is wired into the curve chain.
        int prev_he0 = -1;
        int next_he1 = -1;
        // Current origin vertex in the discrete curvenet
        int temp_origin = -1;

        for (int s_idx = 0; s_idx < crvSplines.size(); s_idx++) {
            int s = crvSplines[s_idx];
            if (s < 0 || s >= cnSpline.size()) {
                return -1;
            }
            // Get the CN halfedges for this spline
            int cn_he = cnSpline[s].he;
            int cn_he_twin = cnHE[cn_he].twin;
            // Get the start and end CN control indices
            int cn_start = cnHE[cn_he].origin;
            int cn_end   = cnHE[cn_he_twin].origin;
            // Get the start and end dCN vertices
            int s_start = inputCtoV.at(cn_start);
            int s_end   = inputCtoV.at(cn_end);
            double temp_t = 0.0;   // t of temp_origin on this spline (always 0 at the spline's own start)
            if (s_idx == 0) {
                C[c].start = s_start;
                temp_origin = s_start;
            }
            if (s_idx == static_cast<int>(crvSplines.size()) - 1) {
                C[c].end = s_end;
            }

            // Sample the spline
            int n_samples = cnSpline[s].num_samples;
            std::vector<Eigen::Vector3d> samples;
            std::vector<double> sampleT;
            if (sampleNaive) {
                samples = CN->sampleBezierNaive(s, n_samples);
                sampleT.resize(n_samples);
                double h = 1.0 / (n_samples - 1);
                for (int i = 0; i < n_samples; i++) {
                    sampleT[i] = std::min(1.0, i * h);
                }
            } else {
                sampleT = CN->unifSampleT(s, n_samples);
                samples.resize(sampleT.size());
                for (int i = 0; i < static_cast<int>(sampleT.size()); i++) {
                    samples[i] = CN->tSampleBezier(s, sampleT[i]);
                }
            }
            // Get the local spline index for the start and end vertices
            std::vector<int> startLocalSplineIdx = CN->controlLocalSplineIdx(cn_start, s);
            std::vector<int> endLocalSplineIdx   = CN->controlLocalSplineIdx(cn_end, s);

            // Self-loop: parent curvenet has two local halfedges at the same control
            if (cn_start == cn_end && startLocalSplineIdx.size() == 2) {
                endLocalSplineIdx[0] = startLocalSplineIdx[1];
            }
            // Insert the discretized spline (i = 0 is the start vert)
            for (int i = 1; i < n_samples; i++) {
                int v = -1;

                if (i == n_samples - 1) {   // Last sample is the end vert
                    v = s_end;
                } else {    // Interior sample
                    v = addVert(samples[i]);
                }
                // Add a new edge
                int new_edge = addEdge(temp_origin, v, prev_he0, next_he1, c);
                edgeSpline.push_back(s);
                edgeSplineT.push_back({temp_t, sampleT[i]});

                int he0 = E[new_edge].he;       // positive/canonical direction
                int he1 = HE[he0].twin;         // negative/opposite direction

                // Assign adjacent halfedge(s) for each prev vert
                if (vertData[temp_origin].cn_idx < 0) {    // Interior vert
                    V[temp_origin].adjHE.clear();
                    V[temp_origin].adjHE.push_back(he0);
                } else if (i == 1) {    // Previous vert must be the spline start
                    V[temp_origin].adjHE[startLocalSplineIdx[0]] = he0;
                    if (s_idx == 0) {   // If we are starting the curve, make it the curve starting halfedge
                        C[c].he_start = he0;
                    }
                }
                // Assign outgoing halfedge for the end vert when the spline ends
                if (i == n_samples - 1) {
                    V[v].adjHE[endLocalSplineIdx[0]] = he1;
                    if (s_idx == crvSplines.size() - 1) {  // We are ending the curve, so add the outgoing halfedge to the end vert
                        C[c].he_end = he1;
                    }
                }

                prev_he0 = he0;
                next_he1 = he1;
                temp_origin = v;
                temp_t = sampleT[i];
            }
        }

        return c;
    }

    // Propagate weights to rest of curvenet using a Laplacian
    int dcurvenet::propagateWeights() {
        // Sync fixed/value info from the source curvenet's controls into the base's generic fields
        for (int v = 0; v < V.size(); v++) {
            int v_type = vertData[v].cn_type;
            int v_idx = vertData[v].cn_idx;
            if (v_type != -1 && CN->C[v_idx].fixed_w) {
                V[v].fixed_w = true;
                V[v].w = CN->C[v_idx].w;
            } else {
                V[v].fixed_w = false;
            }
        }
        return polynet::propagateWeights();
    }

    bool dcurvenet::hasOrderedConnectivity() const {
        return CN->setMesh;
    }

    // Find the closest point on the discretized curve network, then recover which curvenet
    // spline/t-value it corresponds to
    int dcurvenet::closestPoint(const Eigen::Vector3d& p, Curvenet::cnBindData& bind, bool snap, double snapTol) const {
        polyBindData pBind;
        if (polynet::closestPoint(p, pBind, snap, snapTol) != 1) {
            return -1;
        }

        int e = -1;
        double local_t = 0.0;
        if (pBind.elType == 1) {   // Edge
            e = pBind.elIdx;
            local_t = pBind.t;
        } else if (pBind.elType == 0) {   // Vertex: borrow any incident edge
            int v = pBind.elIdx;
            if (V[v].adjHE.empty()) {
                return -1;
            }
            int he = V[v].adjHE[0];
            e = HE[he].edge;
            int origin_v = HE[HE[E[e].he].twin].dest;
            local_t = (v == origin_v) ? 0.0 : 1.0;
        } else {
            return -1;
        }
        if (e < 0 || e >= edgeSpline.size()) {
            return -1;
        }

        bind = Curvenet::cnBindData();
        bind.s = edgeSpline[e];
        bind.t = (1.0 - local_t) * edgeSplineT[e].first + local_t * edgeSplineT[e].second;
        bind.pos = pBind.pos;
        return 1;
    }
}   // namespace Polynet
