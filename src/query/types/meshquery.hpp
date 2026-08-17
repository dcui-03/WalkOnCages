// meshquery.hpp
#pragma once

#include "query/query.hpp"
#include "mesh/mesh.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>

namespace Query {

// Object to query barycentric coordinates of
class meshquery : public query {
    public:
        // Init with pointer to mesh
        meshquery(Mesh::mesh* M);

        // Get relevant vertices as a matrix
        virtual int matrixVerts(Eigen::MatrixXd& Verts) override;

        // Apply smoothing
        virtual int applySmoothing(const Eigen::MatrixXd& Input, Eigen::MatrixXd& Result) override;

        // Get a smoothing operator
        virtual int computeSmoothingOp(int num_samples) override;
    protected:
    private:
        // Pointer to mesh object
        Mesh::mesh* M;
        // Smoothing operators
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> AtL_inv;
        Eigen::VectorXd A;
        double t;
};

}   // namespace Query
