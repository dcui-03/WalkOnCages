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

        // Compute the actual coordinates. Smooths automatically
        // coordType: 0 is harmonic, 1 is MVC, 2 is positive-MVC
        int computeCoordinates(int coordType = 0, int num_samples = 20, int max_samples = 10000);

        // Apply deformation
        // First, query the new vertex locations, then apply the coords operator
        int applyDeformation(Eigen::MatrixXd& query_pos);

        // Interpolate colors from the cage to the query
        int applyColor(const Eigen::MatrixXd& Colors, Eigen::MatrixXd& Result);
        // Debug colors on query in polyscope format
        int colorsPolyscopeFormat(std::vector<glm::vec3>& colors);
        // Visualize cage's colors in polyscope format
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

        // Ask the cage for a set of random samples (via raycasting)
        int computeRandomSamples(const Eigen::Vector3d& q_pos, std::mt19937& gen, std::vector<Utils::projData>& hits);

        // Solve for gradX/gradY/gradZ/beta (the components of u_x) at query q
        void solveAlpha(int q, const std::vector<Sample>& samples);

        // Smooths gradX/gradY/gradZ/beta
        int applySmoothing(int num_samples);

        // Coordinates for each vertex as an n x m matrix (i.e., each row is a query; each col is a cage vert)
        Eigen::MatrixXd coords;
        // Gradient of each coordinate w.r.t. the query position (first three coords of u_x)
        Eigen::MatrixXd gradX, gradY, gradZ;
        // Homogeneous coordinate of u_x
        Eigen::MatrixXd beta;

        // Store current coordinate type
        int coordType = 0;
        // Pointers to various needs
        Cage::cage* def_cage;
        Query::query* def_query;
};

}   // namespace CageDeformer