#include "cagedeformer.hpp"

#include "cage/cage.hpp"
#include "query/query.hpp"
#include "cagedeformer_types.hpp"
#include "utils/decUtils.hpp"
#include "utils/utils.hpp"
#include "utils/wos.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/QR>
#include <chrono>
#include <iostream>
#include <random>


namespace CageDeformer {
    cagedeformer::cagedeformer(): def_cage(nullptr), def_query(nullptr) {

    }

    void cagedeformer::applyCage(Cage::cage* C) {
        def_cage = C;
        return;
    }
    void cagedeformer::applyQuery(Query::query* Q) {
        def_query = Q;
        return;
    }

    int cagedeformer::computeCoordinates(int coordType, int num_samples, int max_samples) {
        // First make sure we have sufficient conditions
        if (def_cage == nullptr || def_query == nullptr) {
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
        for (int q = 0; q < queries.size(); q++) {
            queries[q].pos = query_pos.row(q).transpose();
        }

        // Get all cage controls' positions
        Eigen::MatrixXd cage_pos;
        def_cage->matrixVerts(cage_pos);
        if ((cage_pos.rows() < 1) || (cage_pos.cols() != 3)) {
            return -1;
        }
        // Set up the coordinate matrix and its per-component gradient
        coords.resize(queries.size(), cage_pos.rows());
        coords.setZero();
        gradX.resize(queries.size(), cage_pos.rows());
        gradX.setZero();
        gradY.resize(queries.size(), cage_pos.rows());
        gradY.setZero();
        gradZ.resize(queries.size(), cage_pos.rows());
        gradZ.setZero();

        // Tracks whether each query point got at least one successful sample
        std::vector<char> query_ok(queries.size(), 0);

        // For each query point...
        #pragma omp parallel for
        for (int q = 0; q < queries.size(); q++) {
            // One persistent generator per thread; seeded once, reused across queries/walks
            thread_local std::mt19937 gen{std::random_device{}()};

            int num_success = -1;
            if (coordType == 0) {
                num_success = computeHarmonicCoordinates(queries[q].pos, num_samples, max_samples, gen, queries[q].samples);
            } else if (coordType == 1) {
                num_success = computeMVCoordinates(queries[q].pos, num_samples, max_samples, gen, queries[q].samples);
            } else if (coordType == 2) {
                num_success = computePositiveMVCoordinates(queries[q].pos, num_samples, max_samples, gen, queries[q].samples);
            }
            query_ok[q] = (num_success > 0) ? 1 : 0;

            solveAlpha(q, queries[q].pos, queries[q].samples);
        }

        // Fail loudly if any query point never got a single successful sample
        for (int q = 0; q < query_ok.size(); q++) {
            if (!query_ok[q]) {
                return -1;
            }
        }
        applySmoothing();    // Apply smoothing
        return 1;
    }

    // Harmonic coordinates: WoS random walks, boundary hit weighted uniformly
    int cagedeformer::computeHarmonicCoordinates(const Eigen::Vector3d& q_pos, int num_samples, int max_samples, std::mt19937& gen, std::vector<Sample>& samples) {
        int num_success = 0;
        double eps = 1e-3 * def_cage->bboxDiag();
        for (int s = 0; s < max_samples; s++) {
            // Launch a WoS walk, until we get at least num_samples samples.
            // The first hop is stratified across attempts to reduce clustering; the rest of the walk is genuinely random.
            int elType, elIdx;
            Eigen::Vector3d proj;
            Eigen::VectorXd wos_coords;
            if (def_cage->closestPoint(q_pos, elType, elIdx, proj, wos_coords) != 1) {
                continue;
            }
            double d0 = (proj - q_pos).norm();
            int success;
            if (d0 <= eps) {
                success = 1;
            } else {
                Eigen::Vector3d firstDirec;
                WoS::stratifySamples(s, firstDirec);
                Eigen::Vector3d p1 = q_pos + d0 * firstDirec;
                success = WoS::WalkOnSpheres(p1, def_cage, elType, elIdx, proj, wos_coords, gen, 1, 60, eps);
            }
            if (success != 1) { // Failed, continue
                continue;
            }
            Sample y;
            y.pos = proj;
            y.bases = def_cage->computeBasis(elType, elIdx, wos_coords);
            y.weight = 1.0;    // Harmonic weights implicitly weight via poisson kernel
            samples.push_back(y);
            num_success++;
            if (num_success >= num_samples) {
                break;
            }
        }
        return num_success;
    }

    // MVC: not yet implemented, needs cage raycasting
    int cagedeformer::computeMVCoordinates(const Eigen::Vector3d& q_pos, int num_samples, int max_samples, std::mt19937& gen, std::vector<Sample>& samples) {
        return -1;
    }

    // Positive MVC: not yet implemented, needs cage raycasting
    int cagedeformer::computePositiveMVCoordinates(const Eigen::Vector3d& q_pos, int num_samples, int max_samples, std::mt19937& gen, std::vector<Sample>& samples) {
        return -1;
    }

    // Solve for alpha_v using M and m_v
    void cagedeformer::solveAlpha(int q, const Eigen::Vector3d& q_pos, const std::vector<Sample>& samples) {
        Eigen::Matrix4d M;
        M.setZero();
        for (const Sample& y : samples) {
            Eigen::Vector4d homog_y = {y.pos[0], y.pos[1], y.pos[2], 1};
            M += y.weight * (homog_y * homog_y.transpose());
        }
        // Rank-deficiency-safe solve (samples can be near-coplanar, e.g. near a flat cage face)
        Eigen::CompleteOrthogonalDecomposition<Eigen::Matrix4d> M_cod(M);

        std::map<int, Eigen::Vector4d> m_v;   // m_v is evaluated at vertices
        // For each v which contributed, compute m_v
        for (const Sample& y : samples) {
            Eigen::Vector4d homog_y = {y.pos[0], y.pos[1], y.pos[2], 1};
            // For each contributing vert, gather all its influence
            for (int v = 0; v < y.bases.size(); v++) {
                if (auto it = m_v.find(y.bases[v].first); it != m_v.end()) {    // Existing entry, accumulate
                    m_v[y.bases[v].first] += y.weight * y.bases[v].second * homog_y;
                } else {    // New entry to list
                    m_v[y.bases[v].first] = y.weight * y.bases[v].second * homog_y;
                }
            }
        }

        // Compute the alpha value for each v, along with its gradient (the linear part of u_x)
        Eigen::Vector4d homog_q = {q_pos[0], q_pos[1], q_pos[2], 1};
        for (const auto& [v, m_vx] : m_v) {
            Eigen::Vector4d u_x = M_cod.solve(m_vx);
            coords(q, v) = homog_q.transpose() * u_x;
            gradX(q, v) = u_x(0);
            gradY(q, v) = u_x(1);
            gradZ(q, v) = u_x(2);
        }
        return;
    }

    // Get the smoothing operator from the user, and apply it to coords and each gradient component
    int cagedeformer::applySmoothing() {
        Eigen::MatrixXd result;
        if (def_query->applySmoothing(coords, result) != 1) {
            return -1;
        }
        coords = result;

        if (def_query->applySmoothing(gradX, result) != 1) {
            return -1;
        }
        gradX = result;

        if (def_query->applySmoothing(gradY, result) != 1) {
            return -1;
        }
        gradY = result;

        if (def_query->applySmoothing(gradZ, result) != 1) {
            return -1;
        }
        gradZ = result;

        return 1;
    }

    // Apply deformation
    // First, query the new vertex locations, then apply the coords operator
    int cagedeformer::applyDeformation(Eigen::MatrixXd& query_pos) {
        Eigen::MatrixXd cage_pos;
        def_cage->matrixVerts(cage_pos);
        // Apply the coordinate operator
        query_pos = coords * cage_pos;
        return 1;
    }

    // Interpolate arbitrary per-cage-vertex values
    int cagedeformer::applyColor(const Eigen::MatrixXd& Colors, Eigen::MatrixXd& Result) {
        if (Colors.rows() != coords.cols()) {
            return -1;
        }
        Result = coords * Colors;
        return 1;
    }

    // Spatial gradient of arbitrary per-cage-vertex values, using our pre-computed gradient operators
    int cagedeformer::applyColorGradient(const Eigen::MatrixXd& Colors, Eigen::MatrixXd& GradX, Eigen::MatrixXd& GradY, Eigen::MatrixXd& GradZ) {
        if (Colors.rows() != gradX.cols()) {
            return -1;
        }
        GradX = gradX * Colors;
        GradY = gradY * Colors;
        GradZ = gradZ * Colors;
        return 1;
    }
}   // namespace CageDeformer
