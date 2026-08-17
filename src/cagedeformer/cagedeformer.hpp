// cagedeformer.hpp
#pragma once

#include "cage/cage.hpp"
#include "query/query.hpp"
#include "cagedeformer_types.hpp"

#include <vector>
#include <array>
#include <map>
#include <random>
#include <stdexcept>
#include <Eigen/Core>
#include <glm/vec3.hpp>

namespace CageDeformer {

class cagedeformer {
    public:
        // Constructor, which first builds the mesh
        cagedeformer();

        // Apply a boundary cage
        void applyCage(Cage::cage* C);
        void applyQuery(Query::query* Q);

        // Compute the actual coordinates
        // coordType: 0 is harmonic, 1 is MVC, 2 is positive-MVC
        int computeCoordinates(int coordType = 0, int num_samples = 20, int max_samples = 10000);
        // The user supplies the correct smoothing operator for the query object
        // Assuming for now that the operator is sparse (unlikely to be dense)
        int applySmoothing();

        // Apply deformation
        // First, query the new vertex locations, then apply the coords operator
        int applyDeformation(Eigen::MatrixXd& query_pos);

        // Interpolate arbitrary per-cage-vertex values (e.g. debug colors) using the same coords operator
        int applyColor(const Eigen::MatrixXd& Colors, Eigen::MatrixXd& Result);
        // Debug colors, per query vertex, ready for a polyscope vertex color quantity (cagedeformer_ps.cpp)
        int colorsPolyscopeFormat(std::vector<glm::vec3>& colors);
        // The cage's own (uninterpolated) colors, ready for a polyscope point cloud color quantity
        int cageColorsPolyscopeFormat(std::vector<glm::vec3>& colors);

        // Get spatial gradients using stored values (for colors specifically)
        int applyColorGradient(const Eigen::MatrixXd& Colors, Eigen::MatrixXd& GradX, Eigen::MatrixXd& GradY, Eigen::MatrixXd& GradZ);
        // Visualize color gradients in polyscope format
        int colorGradientsPolyscopeFormat(std::vector<glm::vec3>& gradR, std::vector<glm::vec3>& gradG, std::vector<glm::vec3>& gradB);
    protected:
        // No class inheritance
    private:
        // Diff coordinate types. Each fills samples with this query's estimates and
        // returns the number of successful samples (-1 on hard failure)
        int computeHarmonicCoordinates(const Eigen::Vector3d& q_pos, int num_samples, int max_samples, std::mt19937& gen, std::vector<Sample>& samples);
        int computeMVCoordinates(const Eigen::Vector3d& q_pos, int num_samples, int max_samples, std::mt19937& gen, std::vector<Sample>& samples);
        int computePositiveMVCoordinates(const Eigen::Vector3d& q_pos, int num_samples, int max_samples, std::mt19937& gen, std::vector<Sample>& samples);

        // Shared moment-fit solve, coordType-agnostic. Writes row q of coords
        void solveAlpha(int q, const Eigen::Vector3d& q_pos, const std::vector<Sample>& samples);

        // Coordinates for each vertex as an n x m matrix (i.e., each row is a query; each col is a cage vert)
        Eigen::MatrixXd coords;
        // Gradient of each coordinate w.r.t. the query position, same shape as coords, split by component
        Eigen::MatrixXd gradX, gradY, gradZ;

        // Store current coordinate type
        int coordType = 0;
        // Pointers to various needs
        Cage::cage* def_cage;
        Query::query* def_query;
};

}   // namespace CageDeformer