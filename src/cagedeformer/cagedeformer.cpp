#include "cagedeformer.hpp"

#include "mesh/mesh.hpp"
#include "curvenet/curvenet.hpp"
#include "dcurvenet/dcurvenet.hpp"
#include "utils/decUtils.hpp"
#include "utils/utils.hpp"
#include "utils/wos.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <chrono>
#include <iostream>


namespace CageDeformer {
    int cagedeformer::computeCoordinates(const Eigen::MatrixXd& query_pos, int coordType, int num_samples, int max_samples) {
        // First make sure we have sufficient conditions
        if (cage_CN == nullptr && cage_M == nullptr) {
            return -1;
        }
        // Check the coordinate type
        if (coordType < 0 || coordType > 2) {
            return -1;
        }
        // Check the query positions
        if ((query_pos.rows() < 1) || (query_pos.cols() != 3)) {
            return -1;
        }


        // Same initial setup for all coords. Only differ by weights and sampling criteria
        // First, get all query points and create a query for each
        std::vector<Query> queries(query_pos.rows());
        // Get all cage controls' positions
        // TODO
        Eigen::MatrixXd cage_pos;
        
        // Set up the coordinate matrix
        coords.resize(queries.size(), cage_pos.rows());
        coords.setZero();

        // For each query point...
        #pragma omp parallel for
        for (int q = 0; q < queries.size(); q++) {
            int num_success;
            for (int s = 0; s < max_samples; s++) {
                // Launch k WoS walks, until we get at least t samples
                // TODO: make adaptable/generalize to more cage types
                Mesh::meshBindData hitData;
                // TODO: tune the correct walk parameters
                int success = WoS::WalkOnSpheres_Mesh(queries[q].pos, cage_M, hitData, 0, 60, 1e-3 * M->bboxDiag);
                if (success != 1) { // Failed, continue
                    continue;
                }
                // TODO: not every coordType wants to do a traditional WoS, some want to do raycasting/visibility... may need to modify/rethink?
                // For each sample, compute its relevant verts and basis coords
                Sample y;
                y.pos = hitData.proj;
                if (hitData.elType == 0) {  // A single vertex
                    y.bases.push_back({hitData.elIdx, 1.0});
                } else if (hitData.elType == 1) {
                    int he0 = cage_M->E[hitData.elIdx].he;
                    int v1 = cage_M->HE[he0].dest;
                    int v0 = cage_M->V[cage_M->HE[he].twin].dest;
                    y.bases.push_back({v0, hitData.coords[0]});
                    y.bases.push_back({v1, 1.0 - hitData.coords[0]});
                } else if (hitData.elType == 2) {
                    std::vector<int> fVerts = cage_M->faceAdjVertIdxs(hitData.elIdx);
                    for (int v = 0; v < fVerts.size(); v++) {
                        y.bases.push_back({v, hitData.coords[v]});
                    }
                }
                // Get the correct corresponding weight for the sample
                if (coordType == 0) {
                    y.weight = 1.0;
                } else if (coordType == 1) {    // TODO
                    y.weight = 1.0
                } else if (coordType == 2) {    // TODO
                    y.weight = 1.0;
                }
                // Append sample to the query point's estimates
                queries[q].samples.push_back(y);
            }
            
            // After recording all estimates, compute M -> M_inv once
            Eigen::Matrix4d M;
            M.setZero();
            for (int s = 0; s < queries[q].samples.size(); s++) {
                const Sample& y = queries[q].samples[s];
                Eigen::Vector4d homog_y = {y.pos[0], y.pos[1], y.pos[2], 1};
                M += queries[q].weight * (homog_y * homog_y.transpose());
            }
            // Invert M
            Eigen::Matrix4d M_inv = M.inverse();

            std::map<int, Eigen::Vector4d> m_v;   // m_v evaluated 
            // For each v which contributed, compute m_v
            for (int s = 0; s < queries[q].samples.size(); s++) {
                const Sample& y = queries[q].samples[s];
                Eigen::Vector4d homog_y = {y.pos[0], y.pos[1], y.pos[2], 1};
                // For each contributing vert, gather all its influence
                for (int v = 0; v < y.bases.size(); v++) {
                    if (auto it = m_v.find(y.bases[v].first); it != m_v.end()) {    // New entry to list
                        m_v[y.bases.first] = y.weight * y.bases[v].second * homog_y;
                    } else {
                        m_v[y.bases.first] += y.weight * y.bases[v].second * homog_y;
                    }
                }
            }

            // Compute the alpha value for each v
            Eigen::Vector4d homog_q = {queries[q].pos[0], queries[q].pos[1], queries[q].pos[2], 1};
            for (const auto& [v, m_vx] : m_v) {
                double alpha = homog_q.transpose() * M_inv * m_vx;
                // Fill in the i-th row of the coords matrix with these values
                coords(q, v) = alpha;
            }
        }

        
        return 1;
    }

    // The user supplies the correct smoothign operator for the query object
    int cagedeformer::applySmoothing(const Eigen::SparseMatrix<double>& smoothingOp) {
        // 1. Check that the dimensions of the smoothing operator are fine
        if (smoothingOp.rows() != coords.rows() || smoothingOp.cols() != coords.rows()) {
            return -1;
        }
        // 2. For each column in the coord matrix, apply the smoothing operator (we can do this all at once actually!)
        coords = smoothingOp * coords;
        return 1;
    }

    // Apply deformation
    // First, query the new vertex locations, then apply the coords operator
    int cagedeformer::applyDeformation(Eigen::MatrixXd& query_pos) {
        // Get the relevant coords
        if (cage_CN == nullptr && cage_M == nullptr) {
            return -1;
        }

        // TODO: Get CN or Mesh coords
        Eigen::MatrixXd cage_pos;
        if (cage_pos.rows() != coords.cols()) {
            return -1;
        }

        // Apply the coordinate operator
        query_pos = coords * cage_pos;
        return 1;
    }

    // Compute harmonic coordinates
    int cagedeformer::computeHarmonicCoordinates() {
        double w_k = 1.0;   // Harmonic weights implicitly weight via poisson kernel
    }
}   // namespace ProfileMover