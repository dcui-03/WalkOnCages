#include "meshquery.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <cmath>
#include <iostream>

namespace Query {

// Init with pointer to mesh
meshquery::meshquery(Mesh::mesh* M);

// Get relevant vertices as a matrix
int meshquery::matrixVerts(Eigen::MatrixXd& Verts) {
    int success = computeSmoothingOp();
    if (success != 1) {
        std::cout << "No available smoothing operator, or construction failed. Using without smoothing." << std::endl;
    }
};

// Get relevant vertices as a matrix
int meshquery::matrixVerts(Eigen::MatrixXd& Verts) {
    M->vertsAsMatrix(Verts);
    return 1;
}

int meshquery::applySmoothing(const Eigen::MatrixXd& Input, Eigen::MatrixXd& Result) {
    if (Input.rows() != A.size()) { // Size mismatch
        return -1;
    }
    // Compute smoothing
    Result = AtL_inv.solve(A.asDiagonal() * Input);
    return 1;
}

// Get a smoothing operator
int meshquery::computeSmoothingOp() {
    Eigen::SparseMatrix<double> L;
    int L_success = M->computeLaplacian(L);
    int A_success = M->computeMass(A);

    if (L_success == -1 || A_success == -1) {
        return -1;
    }
    // Compute the t
    double t = std::max(1.0, 5.0 - std::log10(double(M->getNumActiveV()))) * M->getSquaredMeanE();
    Eigen::SparseMatrix<double> LHS = A + t * L;
    AtL_inv.analyzePattern(LHS);
    AtL_inv.factorize(LHS);
    smoothAvailable = true;
    return 1;
};

}   // namespace Query