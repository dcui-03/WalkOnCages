// cagedeformer.hpp
#pragma once

#include "curvenet/curvenet.hpp"
#include "dcurvenet/dcurvenet.hpp"
#include "mesh/mesh.hpp"
#include "cutmesh/cutmesh.hpp"
#include "cagedeformer_types.hpp"
#include "utils/decUtils.hpp"
#include <vector>
#include <array>
#include <map>
#include <stdexcept>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

namespace CageDeformer {

class cagedeformer {
    public:
        // Constructor, which first builds the mesh
        cagedeformer();

        // Apply a boundary cage
        void applyCage(Curvenet::curvenet* cage_CN);
        void applyCage(Mesh::mesh* mesh_CN);

        // Compute the actual coordinates
        // coordType: 0 is harmonic, 1 is MVC, 2 is positive-MVC
        int computeCoordinates(Eigen::MatrixXd query_pos, int coordType = 0, int num_samples = 20, int max_samples = 1000);
        // The user supplies the correct smoothing operator for the query object
        // Assuming for now that the operator is sparse (unlikely to be dense)
        int applySmoothing(const Eigen::SparseMatrix<double>& smoothingOp);

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

        // Coordinates for each vertex as an n x m matrix (i.e., each row is the)
        Eigen::MatrixXd coords;

        // Store current coordinate type
        int coordType = 0;
        // Pointers to various needs
        Curvenet::curvenet* cage_CN;
        Mesh::mesh* cage_M;    
};

}   // namespace ProfileMover