#include "dcurvenet.hpp"

#include "curvenet/curvenet.hpp"
#include "mesh/mesh.hpp"
#include "utils/utils.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <vector>
#include <cmath>
#include <utility>

namespace DCurvenet {
    // Takes the original curvenet and discretizes it
    dcurvenet::dcurvenet(Curvenet::curvenet* CN, Mesh::mesh* M): CN(CN), M(M) {
        const std::vector<Curvenet::Control>& cnCtrl = CN->controls();
        int num_curves = CN->numCurves();
        // Defensive reset
        V.clear();
        HE.clear();
        E.clear();
        C.clear();
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
            int c = addCurve(crv);
            inputCrvToC[crv] = c;
        }

        // 3. Iterate over every control vertex that is an intersection/anchor and rewire its outgoing/incoming halfedges
        // NOTE: This is not really necessary. We already have the outgoing order of halfedges per control vertex, and we will never need to traverse betweem curves
        for (int v_idx = 0; v_idx < ctrlVerts.size(); v_idx++) {
            int v = ctrlVerts[v_idx];
            rewireVertAdjHE(v);
        }
        // Compute projection data
        computeProjData();

        // 4. Compute all corner normals and widths
        std::vector<curveDeformData> curveData;
        allCornerNormalsAndWidths(curveData);
        // 5. Transport normals and widths along all splines
        transportNormalsAndWidths(curveData);
        // 6. Compute scaled frames on all splines
        computeScaledFrames();

        // 7. Cleanup by copying realtime frames to rest frames
        copyFramesToRest();
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
                const std::vector<Eigen::Vector3d> samples = CN->unifSample(s, n_samples);
                // Exploit the fact that this matches the halfedge direction that the curve was constructed from
                for (int i = 1; i < samples.size(); i++) {
                    // for numerical reasons, only copy in non-controls
                    if (i != samples.size() - 1) {
                        V[HE[he_curr].dest].new_pos = samples[i];
                    }
                    int v_dest = HE[he_curr].dest;
                    Eigen::Vector3d tangent = V[v_dest].new_pos - V[v_prev].new_pos;
                    HE[he_curr].defData.newFrame.tangent = tangent.normalized();
                    HE[he_curr].defData.newFrame.l = tangent.norm();

                    int he_twin = HE[he_curr].twin;
                    HE[he_twin].defData.newFrame.tangent = -tangent.normalized();
                    HE[he_twin].defData.newFrame.l = tangent.norm();

                    v_prev = v_dest;
                    he_curr = HE[he_curr].next;
                }
            }
        }
        std::vector<curveDeformData> curveData;
        // Compute corner normals and widths
        allCornerNormalsAndWidths(curveData);
        // Transport all normals and widths
        transportNormalsAndWidths(curveData);
        // Compute scaled frames on all splines
        computeScaledFrames();
        // Check to make sure new frames are valid
        validateFrames();
        // Compute the deformation gradients
        computeDefGradAll();
        return;
    }

    dcurvenet::dcurvenet() {

    }

    // Add a new vertex that matches an existing control
    int dcurvenet::addVert(Curvenet::Control ctrl, int ctrl_idx) {
        int v = V.size();
        V.emplace_back();
        V[v].n = ctrl.n;
        V[v].pos = ctrl.pos;
        V[v].new_pos = ctrl.pos;
        V[v].cn_idx = ctrl_idx;
        V[v].cn_type = ctrl.cType;
        V[v].adjHE.resize(ctrl.adjHE.size(), -1);
        return v;
    }
    // Add a vertex given its parameters
    int dcurvenet::addVert(Eigen::Vector3d new_pos, Eigen::Vector3d new_n, int ctrl_idx, int ctrl_type, int adjSize) {
        int v = V.size();
        V.emplace_back();
        V[v].n = new_n;
        V[v].pos = new_pos;
        V[v].new_pos = new_pos;
        V[v].cn_idx = ctrl_idx;
        V[v].cn_type = ctrl_type;
        V[v].adjHE.resize(adjSize, -1);
        return v;
    }
    // Add a new edge in and return its halfedges
    int dcurvenet::addEdge(int origin, int dest, int prev_he0, int next_he1, int c) {
        if (origin < 0 || origin >= V.size() || dest < 0 || dest >= V.size()) {
            return -1;
        }
        int e = E.size();
        E.emplace_back();
        int he0 = HE.size();
        int he1 = he0+1;
        HE.emplace_back();
        HE.emplace_back();
        
        // Rewire
        E[e].he = he0;
        E[e].curve = c;
        HE[he0].dest = dest;
        HE[he1].dest = origin;
        Eigen::Vector3d edgeVec = V[dest].pos - V[origin].pos;
        HE[he0].twin = he1;
        HE[he1].twin = he0;
        HE[he0].prev = prev_he0;
        HE[he1].next = next_he1;
        HE[he0].edge = e;
        HE[he1].edge = e;
        if (prev_he0 != -1) {
            HE[prev_he0].next = he0;
        }
        if (next_he1 != -1) {
            HE[next_he1].prev = he1;
        }
        HE[he0].defData.newFrame.tangent = edgeVec.normalized();
        HE[he1].defData.newFrame.tangent = -1 * edgeVec.normalized();
        HE[he0].defData.newFrame.l = edgeVec.norm();
        HE[he1].defData.newFrame.l = edgeVec.norm();

        return e;
    }
    // Add a curve that matches an input curvenet curve
    int dcurvenet::addCurve(int crv) {
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
        C[c].cn_idx = crv;

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
            if (s_idx == 0) {
                C[c].start = s_start;
                temp_origin = s_start;
            }
            if (s_idx == static_cast<int>(crvSplines.size()) - 1) {
                C[c].end = s_end;
            }

            // Sample the spline
            int n_samples = cnSpline[s].num_samples;
            const std::vector<Eigen::Vector3d> samples = CN->unifSample(s, n_samples);
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

                int he0 = E[new_edge].he;       // positive/canonical direction
                int he1 = HE[he0].twin;         // negative/opposite direction

                // Assign adjacent halfedge(s) for each prev vert
                if (V[temp_origin].cn_idx < 0) {    // Interior vert
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
            }
        }

        return c;
    }
    // Rewire incoming and outgoing halfedges of intersection and anchor vertices
    int dcurvenet::rewireVertAdjHE(int v) {
        const std::vector<int> adjHE = V[v].adjHE;
        for (int he = 0; he < adjHE.size(); he++) {
            int he0 = adjHE[he];
            int he1 = adjHE[(he+1)%adjHE.size()];
            HE[he0].prev = HE[he1].twin;
            HE[HE[he1].twin].next = he0;
        }
        return 1;
    }


    // For intersections, computes their corner normals. For non-controls, this method does nothing (return -1)
    int dcurvenet::vertCornerNormalsWidths(int v, std::vector<curveDeformData>& curveData) {
        if (!V[v].active) {
            return -1;
        }
        // Only original curvenet controls have corner data
        if (V[v].cn_idx < 0 || V[v].cn_type < 1) {
            return -1;
        }
        double eps = 1e-12;
        const std::vector<int> adjHE = V[v].adjHE;
        std::vector<Eigen::Vector3d> adjNormals(adjHE.size());
        // Explictly handle anchors and closed curves
        if (V[v].cn_type < 3) {
            // Check if we are on the start of a curve
            int he0 = adjHE[0];
            int c = E[HE[he0].edge].curve;
            if ((C[c].start != v) && V[v].cn_type == 2) {   // valence 2 control that is not the start of the curve
                return -1;
            }
            for (int he_idx = 0; he_idx < adjHE.size(); he_idx++) {
                int he = adjHE[he_idx];
                const Eigen::Vector3d& tan = HE[he].defData.newFrame.tangent;
                // Gram-schmidt to ensure the normal is orthogonal to each tangent
                Eigen::Vector3d n = V[v].n - V[v].n.dot(tan) * tan;
                if (n.norm() <= eps) {  // Extremely rare case where normal == tangent or normal == -tangent
                    // Just pick a random direction orthogonal to the tangent 
                    if (std::abs(tan(0)) < 0.9) {
                        n = tan.cross(Eigen::Vector3d::UnitX());
                    } else {
                        n = tan.cross(Eigen::Vector3d::UnitY());
                    }
                }
                double l = n.norm();
                n.normalize();
                int next_he = adjHE[(he_idx+1)%adjHE.size()];
                double w = HE[he].defData.newFrame.l + l * (HE[next_he].defData.newFrame.l - HE[he].defData.newFrame.l);
                if (isPositiveHalfedge(he)) {
                    curveData[c].N_pos.first = n;
                    curveData[c].W_pos.first = w;
                    curveData[c].N_neg.first = n;
                    curveData[c].W_neg.first = w;
                } else {
                    curveData[c].N_pos.second = n;
                    curveData[c].W_pos.second = w;
                    curveData[c].N_neg.second = n;
                    curveData[c].W_neg.second = w;
                }
            }
            return 1;
        }

        // Otherwise, we are at an intersection and need to compute normals explicitly
        std::vector<bool> skipList(adjNormals.size(), false);
        for (int he = 0; he < adjHE.size(); he++) {
            int he0 = adjHE[he];
            int he1 = adjHE[(he+1)%adjHE.size()];
            int c0 = E[HE[he0].edge].curve;
            int c1 = E[HE[he1].edge].curve;
            Eigen::Vector3d cornerNormal;
            // First check if we are parallel. If so skip for now
            double dotProdTest = (HE[he0].defData.newFrame.tangent.normalized()).dot(HE[he1].defData.newFrame.tangent.normalized());
            if (dotProdTest <= -1.0+eps) {
                skipList[he] = true;
                continue;
            } else if (dotProdTest >= 1 - eps) {    // Two outgoing HEs are coincident
                const Eigen::Vector3d& tan = HE[he0].defData.newFrame.tangent;
                cornerNormal = V[v].n - V[v].n.dot(tan) * tan;
                if (cornerNormal.norm() <= eps) {  // Extremely rare case where normal == tangent or normal == -tangent
                    // Just pick a random direction orthogonal to the tangent 
                    if (std::abs(tan(0)) < 0.9) {
                        cornerNormal = tan.cross(Eigen::Vector3d::UnitX());
                    } else {
                        cornerNormal = tan.cross(Eigen::Vector3d::UnitY());
                    }
                }
            } else {
                // Compute corner normal as usual
                cornerNormal = HE[he0].defData.newFrame.tangent.cross(HE[he1].defData.newFrame.tangent);
            }
            // Flip normals if the signed angle was > 180 (ex., simplified check against the vertex normal)
            if (cornerNormal.normalized().dot(V[v].n) < 0.0) {
                cornerNormal *= -1.0;
            }

            // Assign corner normals to curves
            // First figure out if this is the start halfedge of the curve
            if (isPositiveHalfedge(he0)) {
                curveData[c0].N_pos.first = cornerNormal.normalized();
            } else {
                curveData[c0].N_neg.second = cornerNormal.normalized();
            }
            if (isPositiveHalfedge(he1)) {
                curveData[c1].N_neg.first = cornerNormal.normalized();
            } else {
                curveData[c1].N_pos.second = cornerNormal.normalized();
            }
            adjNormals[he] = cornerNormal;
        }
        // Loop over any degenerate cases (straight angles)
        for (int he = 0; he < skipList.size(); he++) {
            if (!skipList[he]) {
                continue;
            }
            int he0 = adjHE[he];
            int he_m1_local = (he + adjHE.size() - 1)%adjHE.size();
            int he1_local = (he + 1)%adjHE.size();
            int he1 = adjHE[he1_local];
            // Get adjacent curves
            int c0 = E[HE[he0].edge].curve;
            int c1 = E[HE[he1].edge].curve;

            // Average the adjacent vectors
            Eigen::Vector3d cornerNormal = adjNormals[he_m1_local] + adjNormals[he1_local];
            // Extremely unlikely, but just in case, put in a safeguard...
            if (cornerNormal.norm() <= eps) {
                Eigen::Vector3d tan = HE[he0].defData.newFrame.tangent;
                cornerNormal = V[v].n - V[v].n.dot(tan) * tan;
                if (cornerNormal.norm() <= eps) {
                    // Just pick a random direction orthogonal to the tangent 
                    if (std::abs(tan(0)) < 0.9) {
                        cornerNormal = tan.cross(Eigen::Vector3d::UnitX());
                    } else {
                        cornerNormal = tan.cross(Eigen::Vector3d::UnitY());
                    }
                }
            }

            // Assign corner normals to curves
            if (isPositiveHalfedge(he0)) {
                curveData[c0].N_pos.first = cornerNormal.normalized();
            } else {
                curveData[c0].N_neg.second = cornerNormal.normalized();
            }
            if (isPositiveHalfedge(he1)) {
                curveData[c1].N_neg.first = cornerNormal.normalized();
            } else {
                curveData[c1].N_pos.second = cornerNormal.normalized();
            }
            adjNormals[he] = cornerNormal;
        }
        // Compute widths
        // Use the computed normal norm list to calculate widths
        for (int he = 0; he < adjHE.size(); he++) {
            // Loop's halfedge indices
            int he0_local = he;
            int he1_local = (he + 1) % adjHE.size();
            int he_m1_local = (he - 1 + adjHE.size()) % adjHE.size();
            // Get the global halfedge index
            int he0 = adjHE[he0_local];
            int he1 = adjHE[he1_local];
            int he_m1 = adjHE[he_m1_local];
            // Get the curve and the corresponding lengths
            int c0 = E[HE[he0].edge].curve;
            double he0_len = HE[he0].defData.newFrame.l;
            double he_m1_len = HE[he_m1].defData.newFrame.l;
            double he1_len = HE[he1].defData.newFrame.l;

            // Compute corner widths
            double cornerWidth0 = he0_len + adjNormals[he0_local].norm() * (he1_len - he0_len);
            double cornerWidth1 = he0_len + adjNormals[he_m1_local].norm() * (he_m1_len - he0_len);

            // Assign corner normals to curves
            // First figure out if this is the start halfedge of the curve
            if (isPositiveHalfedge(he0)) {
                curveData[c0].W_pos.first = cornerWidth0;
                curveData[c0].W_neg.first = cornerWidth1;
            } else {
                curveData[c0].W_neg.second = cornerWidth0;
                curveData[c0].W_pos.second = cornerWidth1;
            }
        }
        return 1;
    }

    // Corner normals on only intersections
    int dcurvenet::allCornerNormalsAndWidths(std::vector<curveDeformData>& curveData) {
        curveData.clear();
        curveData.resize(C.size());
        // TODO: Parallelize? NOTE: This may not be safe, since we operate on curveData simultaneously
        #pragma omp parallel for
        for (int v = 0; v < V.size(); v++) {
            vertCornerNormalsWidths(v, curveData);
        }
        return 1;
    }

    // Transport corner normals and widths from the two ends of a curve
    int dcurvenet::transportNWOnCurve(int c, const curveDeformData& cData) {
        if (!C[c].active) {
            return -1;
        }
        // First, query what kinds of endpoints we have
        int start = C[c].start;
        int end = C[c].end;
        int he_start = C[c].he_start;
        int he_end = C[c].he_end;
        int start_type = V[C[c].start].cn_type;
        int end_type = V[C[c].end].cn_type;

        // Case 1: Both endpoints are intersections or both are anchors
        if (start_type == end_type) {
            // POSITIVE SIDE
            std::vector<Eigen::Matrix3d> rots;
            std::vector<double> lens;
            Eigen::Vector3d start_n = cData.N_pos.first;
            Eigen::Vector3d end_n = cData.N_pos.second;
            double start_w = cData.W_pos.first;
            double end_w = cData.W_pos.second;
            // 1. For the pos side, first trace until the end vertex, accumulating rotations alpha values
            double total_len = accumulateRotations(C[c].he_start, end, rots, lens);
            // 2. Compute the torsion angle theta
            double torsion = computeTorsion(start_n, end_n, rots[rots.size()-1], -1*HE[C[c].he_end].defData.newFrame.tangent);
            // 3. Propagate rotations and widths to halfedges
            int he_curr = he_start;
            for (int he = 0; he < rots.size(); he++) {
                double beta = lens[he]/total_len;
                Eigen::Matrix3d he_torsion = Utils::computeRotation(HE[he_curr].defData.newFrame.tangent, beta*torsion);
                HE[he_curr].defData.newFrame.normal = he_torsion * rots[he] * start_n;
                HE[he_curr].defData.newFrame.w = (1-beta)*start_w + (beta)*end_w;
                he_curr = HE[he_curr].next;
            }

            // NEGATIVE SIDE
            // The same thing but backwards (so we can trace using the same function of next halfedges)
            start_n = cData.N_neg.second;
            end_n = cData.N_neg.first;
            start_w = cData.W_neg.second;
            end_w = cData.W_neg.first;
            total_len = accumulateRotations(C[c].he_end, start, rots, lens);
            torsion = computeTorsion(start_n, end_n, rots[rots.size()-1], -1*HE[C[c].he_start].defData.newFrame.tangent);
            he_curr = he_end;
            for (int he = 0; he < rots.size(); he++) {
                double beta = lens[he]/total_len;
                Eigen::Matrix3d he_torsion = Utils::computeRotation(HE[he_curr].defData.newFrame.tangent, beta*torsion);
                HE[he_curr].defData.newFrame.normal = he_torsion * rots[he] * start_n;
                HE[he_curr].defData.newFrame.w = (1-beta)*start_w + (beta)*end_w;
                he_curr = HE[he_curr].next;
            }
            return 1;
        }
        // Case 2: Start is an intersection, end is an anchor
        else if ((start_type == 3) && (end_type == 1)) {
            // POSITIVE SIDE
            std::vector<Eigen::Matrix3d> rots;
            std::vector<double> lens;
            Eigen::Vector3d start_n = cData.N_pos.first;
            double start_w = cData.W_pos.first;
            // 1. For the pos side, first trace until the end vertex, accumulating rotations and alpha values
            double total_len = accumulateRotations(C[c].he_start, end, rots, lens);
            // 3. Propagate rotations and widths to halfedges
            int he_curr = he_start;
            for (int he = 0; he < rots.size(); he++) {
                HE[he_curr].defData.newFrame.normal = rots[he] * start_n;
                HE[he_curr].defData.newFrame.w = start_w;
                he_curr = HE[he_curr].next;
            }

            // NEGATIVE SIDE
            start_n = cData.N_neg.first;
            start_w = cData.W_neg.first;
            he_curr = HE[he_start].twin;
            // Do not re-initialize rotations, since we only trace from interesection to anchor
            for (int he = 0; he < rots.size(); he++) {
                HE[he_curr].defData.newFrame.normal = rots[he] * start_n;
                HE[he_curr].defData.newFrame.w = start_w;
                he_curr = HE[he_curr].prev;
            }
            return 1;
        }
        // Case 3: Start is an anchor, end is an intersection
        else if ((start_type == 1) && (end_type == 3)) {
            // Same as case 2 but trace backwards instead of forwards
            // NEGATIVE SIDE
            std::vector<Eigen::Matrix3d> rots;
            std::vector<double> lens;
            Eigen::Vector3d start_n = cData.N_neg.second;
            double start_w = cData.W_neg.second;
            double total_len = accumulateRotations(C[c].he_end, start, rots, lens);
            int he_curr = he_end;
            for (int he = 0; he < rots.size(); he++) {
                HE[he_curr].defData.newFrame.normal = rots[he] * start_n;
                HE[he_curr].defData.newFrame.w = start_w;
                he_curr = HE[he_curr].next;
            }

            start_n = cData.N_pos.second;
            start_w = cData.W_pos.second;
            he_curr = HE[he_end].twin;
            // Do not re-initialize rotations, since we only trace from interesection to anchor
            for (int he = 0; he < rots.size(); he++) {
                HE[he_curr].defData.newFrame.normal = rots[he] * start_n;
                HE[he_curr].defData.newFrame.w = start_w;
                he_curr = HE[he_curr].prev;
            }
            return 1;
        }
        return -1;
    }
    // Transport normals for all curves
    int dcurvenet::transportNormalsAndWidths(const std::vector<curveDeformData>& curveData) {
        // TODO: Parallelize?
        #pragma omp parallel for
        for (int c = 0; c < C.size(); c++) {
            transportNWOnCurve(c, curveData[c]);
        }
        return 1;
    }

    int dcurvenet::computeScaledFrames() {
        // TODO: Parallelize?
        #pragma omp parallel for
        for (int he = 0; he < HE.size(); he++) {
            if (HE[he].active) {
                // Re-orthogonalize normals for safety
                HE[he].defData.newFrame.normal = (HE[he].defData.newFrame.normal - 
                                                    HE[he].defData.newFrame.normal.dot(HE[he].defData.newFrame.tangent) * 
                                                    HE[he].defData.newFrame.tangent).normalized();
                // Positive side
                HE[he].defData.newFrame.binormal = (HE[he].defData.newFrame.tangent.cross(HE[he].defData.newFrame.normal)).normalized();
                HE[he].defData.newFrame.h = std::sqrt(std::abs(HE[he].defData.newFrame.l * HE[he].defData.newFrame.w));
            }
        }
        return 1;
    }

    // Compute projection data for this dcurvenet point
    int dcurvenet::computeProjData() {
        for (int v = 0; v < V.size(); v++) {
            if (!V[v].active) {
                continue;
            }
            Mesh::meshBindData bindData;
            if (M->computeVBinding(V[v].pos, bindData) != 1) {
                return -1;
            }
            V[v].proj.coords = bindData.coords;
            V[v].proj.elType = bindData.elType;
            V[v].proj.elIdx = bindData.elIdx;
            V[v].proj.projVec = bindData.offset;
        }
        return 1;
    }

    // Make sure frames don't 
    int dcurvenet::validateFrames() {
        double eps = 1e-12;
        for (int he = 0; he < HE.size(); he++) {
            const scaledFrame& newFrame = HE[he].defData.newFrame;
            if (newFrame.l <= eps || newFrame.w <= eps || newFrame.h <= eps) {
                throw std::runtime_error("validateFrames(): degenerate frame");
            }
            if (!std::isfinite(newFrame.l) || !std::isfinite(newFrame.w) || !std::isfinite(newFrame.h)) {
                throw std::runtime_error("validateFrames(): non-finite scaled frame scale");
            }
        }
        return 1;
    }
    // Copy scaled frame data to new local variables
    int dcurvenet::copyFrameToRest(int he) {
        if (he >= HE.size()) {
            return -1;
        }
        HE[he].defData.restFrame.tangent = HE[he].defData.newFrame.tangent;
        HE[he].defData.restFrame.binormal = HE[he].defData.newFrame.binormal;
        HE[he].defData.restFrame.normal = HE[he].defData.newFrame.normal;

        HE[he].defData.restFrame.l = HE[he].defData.newFrame.l;
        HE[he].defData.restFrame.w = HE[he].defData.newFrame.w;
        HE[he].defData.restFrame.h = HE[he].defData.newFrame.h;
        return 1;
    }
    int dcurvenet::copyFramesToRest() {
        for (int he = 0; he < HE.size(); he++) {
            copyFrameToRest(he);
        }
        return 1;
    }

    // Propagate weights to rest of curvenet using Laplacian
    int dcurvenet::propagateWeights() {
        using T = Eigen::Triplet<double>;
        std::vector<T> tripletList;
        tripletList.reserve(E.size() * 4);  // Conservative overestimate
        // First, build Laplacian system
        Eigen::SparseMatrix<double> cnL(V.size(), V.size());
        // Eigen::VectorXd M(V.size());
        Eigen::VectorXd f(V.size());
        f.setZero();
        // M.setZero();
        std::vector<bool> fixed(V.size(), false);
        int num_fixed = 0;
        // Process fixed verts first
        for (int v = 0; v < V.size(); v++) {
            int v_type = V[v].cn_type;
            int v_idx = V[v].cn_idx;

            if (v_type != -1 && CN->C[v_idx].fixed_w) {
                fixed[v] = true;
                V[v].w = CN->C[v_idx].w;
                f[v] = V[v].w;
                tripletList.push_back(T(v, v, 1.0));
                num_fixed++;
            }
        }
        if (num_fixed == 0) {
            for (int v = 0; v < V.size(); v++) {
                V[v].w = 1.0;
            }
            return 1;
        }
        for (int e = 0; e < E.size(); e++) {
            int v0 = HE[E[e].he].dest;
            int v1 = HE[HE[E[e].he].twin].dest;

            int v0_idx = V[v0].cn_idx;
            int v1_idx = V[v1].cn_idx;
            // If both are fixed, then skip
            if (fixed[v0] && fixed[v1]) {   // Unlikely case, but check anyways
                continue;
            }
            double e_len = std::max((V[v1].pos - V[v0].pos).norm(), 1e-8);
            double weight = 1/e_len;
            bool f0 = fixed[v0];
            bool f1 = fixed[v1];
            // Add to diagonal entries
            if (!f0 && !f1) {
                tripletList.push_back(T(v0, v0,  weight));
                tripletList.push_back(T(v1, v1,  weight));
                tripletList.push_back(T(v0, v1, -weight));
                tripletList.push_back(T(v1, v0, -weight));
            } else if (!f0 && f1) {
                tripletList.push_back(T(v0, v0, weight));
                f[v0] += weight * V[v1].w;
            } else if (f0 && !f1) {
                tripletList.push_back(T(v1, v1, weight));
                f[v1] += weight * V[v0].w;
            }
            // M[v0] += e_len;
            // M[v1] += e_len;
        }
        cnL.setFromTriplets(tripletList.begin(), tripletList.end());
        // M *= 0.5;
        // Eigen::VectorXd RHS = M.asDiagonal() * f;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> factorL;
        factorL.analyzePattern(cnL);
        factorL.factorize(cnL);
        Eigen::VectorXd new_weights = factorL.solve(f);
        
        // Redistribute weights
        for (int v = 0; v < V.size(); v++) {
            V[v].w = std::clamp(new_weights[v], 0.0, 1.0);
        }

        return 1;
    }
}   // namespace DCurvenet