#include "polynet.hpp"

#include "mesh/mesh.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <array>
#include <algorithm>
#include <stdexcept>

namespace Polynet {
    // Build directly from a raw vertex/edge list and trace its own curves
    polynet::polynet(const std::vector<Eigen::Vector3d>& V_list, const std::vector<std::array<int, 2>>& E_list, const Mesh::mesh* M) {
        V.clear();
        HE.clear();
        E.clear();
        C.clear();
        // 1. Add all verts
        for (int v = 0; v < V_list.size(); v++) {
            addVert(V_list[v]);
        }
        // 2. Add all edges
        for (int e = 0; e < E_list.size(); e++) {
            int v0 = E_list[e][0];
            int v1 = E_list[e][1];
            int new_e = addEdge(v0, v1);
            if (new_e < 0) {
                continue;
            }
            int he0 = E[new_e].he;
            int he1 = HE[he0].twin;
            V[v0].adjHE.push_back(he0);
            V[v1].adjHE.push_back(he1);
        }
        // 3. Rewire every vertex's fan. No normal to sort a high-valence fan by, so
        //    prev/next are left ambiguous (-1) there; valence 1/2 is always unambiguous.
        for (int v = 0; v < V.size(); v++) {
            rewireVertAdjHE(v, false);
        }
        // 4. Trace curves by valence
        assignVertTypeAll();
        traceCurves();
        // Compute projection data
        if (M) {
            computeProjData(M);
            setMesh = true;
        }
        if (computeBVH() != 1) {
            throw std::runtime_error("Failed to build BVH.");
        }
    }

    // empty initializer
    polynet::polynet() {

    }

    // Add a vertex given its parameters
    int polynet::addVert(Eigen::Vector3d new_pos, Eigen::Vector3d new_n) {
        int v = V.size();
        V.emplace_back();
        V[v].n = new_n;
        V[v].pos = new_pos;
        V[v].new_pos = new_pos;
        return v;
    }
    // Add a new edge in and return its halfedges
    int polynet::addEdge(int origin, int dest, int prev_he0, int next_he1, int c) {
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

        return e;
    }
    // Rewire incoming and outgoing halfedges of a vertex
    int polynet::rewireVertAdjHE(int v, bool ordered) {
        const std::vector<int> adjHE = V[v].adjHE;
        if (adjHE.size() > 2 && !ordered) {
            // High-valence vertex with no defined ordering: leave prev/next as -1 (end of curve)
            return 1;
        }
        for (int he = 0; he < adjHE.size(); he++) {
            int he0 = adjHE[he];
            int he1 = adjHE[(he+1)%adjHE.size()];
            HE[he0].prev = HE[he1].twin;
            HE[HE[he1].twin].next = he0;
        }
        return 1;
    }

    // Compute projection data for this polynet's points
    int polynet::computeProjData(const Mesh::mesh* M) {
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

    // Propagate weights to rest of network using a Laplacian
    // Fixed vertices are those with V[v].fixed_w == true, using V[v].w as the boundary value
    int polynet::propagateWeights() {
        using T = Eigen::Triplet<double>;
        std::vector<T> tripletList;
        tripletList.reserve(E.size() * 4);  // Conservative overestimate
        // First, build Laplacian system
        Eigen::SparseMatrix<double> cnL(V.size(), V.size());
        Eigen::VectorXd f(V.size());
        f.setZero();
        std::vector<bool> fixed(V.size(), false);
        int num_fixed = 0;
        // Process fixed verts first
        for (int v = 0; v < V.size(); v++) {
            if (V[v].fixed_w) {
                fixed[v] = true;
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
        }
        cnL.setFromTriplets(tripletList.begin(), tripletList.end());
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

    int polynet::vertsAsMatrix(Eigen::MatrixXd& Verts) {
        Verts.resize(V.size(), 3);
        for (int v = 0; v < V.size(); v++) {
            Verts.row(v) = V[v].pos.transpose();
        }
        return 1;
    }

    // Evaluate basis function
    int polynet::evaluateBasis(int elType, int elIdx, double t, std::vector<std::pair<int, double>>& basis) {
        basis.clear();
        if (elType == 0) {
            if (elIdx > V.size() || elIdx < 0) {
                return -1;
            }
            basis.push_back({elIdx, 1.0});
        } else if (elType == 1) {
            if (elIdx < 0 || elIdx > E.size()) {
                return -1;
            }
            int he0 = E[elType].he;
            int v1 = HE[he0].dest;
            int v0 = HE[HE[he0].twin].dest;
            basis.push_back({v0, 1.0 - t});
            basis.push_back({v1, t});
        } else {
            return -1;
        }
        return 1;
    }
}   // namespace Polynet
