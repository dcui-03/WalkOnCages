#include "cagedeformer.hpp"

#include "cage/cage.hpp"
#include "query/query.hpp"
#include "utils/decUtils.hpp"
#include "utils/utils.hpp"
#include "utils/wos.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <chrono>
#include <iostream>


namespace CageDeformer {
    cagedeformer::cagedeformer(Cage::cage* C, Query::query* Q): def_cage(C), def_query(Q);

    int cagedeformer::computeCoordinates(int coordType, int num_samples, int max_samples) {
        // First make sure we have sufficient conditions
        if (def_cage == nullptr && def_query == nullptr) {
            return -1;
        }
        // Check the coordinate type
        if (coordType < 0 || coordType > 2) {
            return -1;
        }
        // Check the query positions
        Eigen::MatrixXd query_pos;
        def_query->matrixVerts(query_pos);
        if ((query_pos.rows() < 1) || (query_pos.cols() != 3)) {
            return -1;
        }
        // Same initial setup for all coords. Only differ by weights and sampling criteria
        // First, get all query points and create a query for each
        std::vector<QueryVert> queries(query_pos.rows());

        // Get all cage controls' positions
        Eigen::MatrixXd cage_pos;
        def_cage->matrixVerts(cage_pos);
        if ((cage_pos.rows() < 1) || (cage_pos.cols() != 3)) {
            return -1;
        }
        // Set up the coordinate matrix
        coords.resize(queries.size(), cage_pos.rows());
        coords.setZero();

        // For each query point...
        #pragma omp parallel for
        for (int q = 0; q < queries.size(); q++) {
            int num_success;
            for (int s = 0; s < max_samples; s++) {
                // Launch k WoS walks, until we get at least t samples
                // TODO: tune the correct walk parameters
                // TODO: Figure out how to replace hit data...
                int success = WoS::WalkOnSpheres(queries[q].pos, def_cage, TODO, 0, 60, 1e-3 * M->bboxDiag);
                if (success != 1) { // Failed, continue
                    continue;
                }
                // TODO: not every coordType wants to do a traditional WoS, some want to do raycasting/visibility... may need to modify/rethink?
                // For each sample, compute its relevant verts and basis coords
                Sample y;
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

            std::map<int, Eigen::Vector4d> m_v;   // m_v is evaluated at vertices
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

    // Get the smoothing operator from the user
    int cagedeformer::applySmoothing() {
        // For each column in the coord matrix, apply the smoothing operator (we can do this all at once actually!)
        Eigen::MatrixXd result;
        if (def_query->applySmoothing(coords, result) != 1) {
            return -1;
        }
        coords = result;    // Apply changes
        return 1;
    }

    // Apply deformation
    // First, query the new vertex locations, then apply the coords operator
    int cagedeformer::applyDeformation(Eigen::MatrixXd& query_pos) {
        Eigen::MatrixXd cage_pos = def_cage->matrixVerts();
        // Apply the coordinate operator
        query_pos = coords * cage_pos;
        return 1;
    }

    // Compute harmonic coordinates
    int cagedeformer::computeHarmonicCoordinates() {
        double w_k = 1.0;   // Harmonic weights implicitly weight via poisson kernel
    }
}   // namespace ProfileMover