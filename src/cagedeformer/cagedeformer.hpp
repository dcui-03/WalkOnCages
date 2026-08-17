// cagedeformer.hpp
#pragma once

#include "cage/cage.hpp"
#include "query/query.hpp"

#include <vector>
#include <array>
#include <map>
#include <stdexcept>
#include <Eigen/Core>

namespace CageDeformer {

class cagedeformer {
    public:
        // Constructor, which first builds the mesh
        cagedeformer();

        // Apply a boundary cage
        void applyCage(Cage::cage C);
        void applyQuery(Query::query Q);

        // Compute the actual coordinates
        // coordType: 0 is harmonic, 1 is MVC, 2 is positive-MVC
        int computeCoordinates(int coordType = 0, int num_samples = 20, int max_samples = 1000);
        // The user supplies the correct smoothing operator for the query object
        // Assuming for now that the operator is sparse (unlikely to be dense)
        int applySmoothing();

        // Apply deformation
        // First, query the new vertex locations, then apply the coords operator
        int applyDeformation(Eigen::MatrixXd& query_pos);
    protected:
        // No class inheritance
    private:
        // Diff coordinate types
        int computeHarmonicCoordinates();
        int computeMVCoordinates();
        int computePositiveMVCoordinates();
        // Apply Smoothing based on the query's smoothing operator
        // Returns -1 if no available smoothing operator
        int smoothCoordinates();

        // Coordinates for each vertex as an n x m matrix (i.e., each row is a query; each col is a cage vert)
        Eigen::MatrixXd coords;

        // Store current coordinate type
        int coordType = 0;
        // Pointers to various needs
        Cage::cage* def_cage;
        Query::query* def_query; 
};

}   // namespace ProfileMover