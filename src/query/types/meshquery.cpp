#include "meshquery.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <cmath>
#include <iostream>

namespace Query {

// Init with pointer to mesh
meshquery::meshquery(Mesh::mesh* M): M(M) {
    computeSmoothingOp();
}

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
    double scale_factor = 5.0 - std::log10(double(M->getNumActiveV()));
    t = std::max(1.0, scale_factor) * M->getSquaredMeanE();
    // t = 1e-2 * M->getSquaredMeanE();
    std::cout << "Smoothing timestep: " << t << std::endl;
    std::cout << "Mean Squared Edge Length: " << M->getSquaredMeanE() << std::endl;
    std::cout << "Scaling Factor: " << scale_factor << std::endl;
    Eigen::SparseMatrix<double> LHS = Eigen::SparseMatrix<double>(A.asDiagonal()) + t * L;
    AtL_inv.analyzePattern(LHS);
    AtL_inv.factorize(LHS);
    smoothAvailable = true;
    return 1;
}

}   // namespace Query
