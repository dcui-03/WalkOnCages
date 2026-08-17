// decUtils.hpp
#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
// TODO: HEADER file for whatever HE mesh data struct we use

namespace DECUtils {
    typedef Eigen::Triplet<double> T;
    // For input vertex lists, assume that they are correctly order CCW and closed
    // Polygonal Mesh Laplacian operator, which assemble's face-based Laplacians
    Eigen::SparseMatrix<double> LaplacianOp(const std::vector<Eigen::Vector3d>& V,
                                            const std::vector<std::vector<int>>& F,
                                            double lambda = 1.0);

    // All functions below are FACE-BASED
    // HIGHER OPERATORS
    // Face Laplacian operator
    Eigen::MatrixXd faceLaplacianOp(const std::vector<Eigen::Vector3d>& f, double lambda = 1.0);
    Eigen::MatrixXd faceLaplacianOp(const Eigen::MatrixXd& D, const Eigen::MatrixXd& M);

    // Face Divergence operator (d0*)
    // TODO: Check this
    Eigen::MatrixXd divOp(const std::vector<Eigen::Vector3d>& f, double lambda = 1.0);
    Eigen::MatrixXd divOp(const Eigen::MatrixXd& D, const Eigen::MatrixXd& M);

    // Face Curl operator (aka d1)
    Eigen::VectorXd curlOp(int n);

    // Face-based Metric on 1-forms (*)
    Eigen::MatrixXd metricOp(const std::vector<Eigen::Vector3d>& f, double lambda = 1.0);
    Eigen::MatrixXd metricOp(double area, const Eigen::MatrixXd& V, const Eigen::MatrixXd& U, double lambda = 1.0);  // Compute proj. on-the-fly or by calling projOp()


    // MUSICAL ISOMORPHISMS
    // Face Sharp operator
    // I.e., generate a face gradient vector for function given 1-forms
    Eigen::MatrixXd sharpOp(const std::vector<Eigen::Vector3d>& f);
    Eigen::MatrixXd sharpOp(double area, Eigen::Vector3d& fNormal, Eigen::MatrixXd& B, Eigen::Vector3d& barycenter);

    // Face Flat operator
    // I.e., generate a 1-form from face-based gradient
    Eigen::MatrixXd flatOp(const std::vector<Eigen::Vector3d>& f);
    Eigen::MatrixXd flatOp(Eigen::Vector3d& fNormal, Eigen::MatrixXd& E);

    // Face Projection operator
    Eigen::MatrixXd projOp(const std::vector<Eigen::Vector3d>& f);
    Eigen::MatrixXd projOp(const Eigen::MatrixXd& V, const Eigen::MatrixXd& U);


    // FUNDAMENTAL OPERATORS
    // Get face Position operator
    Eigen::MatrixXd posOp(const std::vector<Eigen::Vector3d>& f);

    // Get face difference operator (aka d0)
    Eigen::SparseMatrix<double> diffOp(int n);

    // Get face average operator
    Eigen::SparseMatrix<double> avgOp(int n);

    // Get face midpoint operator
    Eigen::MatrixXd midpntOp(const std::vector<Eigen::Vector3d>& f);
    Eigen::MatrixXd midpntOp(const Eigen::SparseMatrix<double> A, const Eigen::MatrixXd& X);

    // Get face edge operator
    Eigen::MatrixXd edgeOp(const std::vector<Eigen::Vector3d>& f);
    Eigen::MatrixXd edgeOp(const Eigen::MatrixXd& D, const Eigen::MatrixXd& X);

    // Face Gradient operator
    Eigen::MatrixXd gradOp(const std::vector<Eigen::Vector3d>& f);
    Eigen::MatrixXd gradOp(double area, Eigen::Vector3d& fNormal, Eigen::MatrixXd& E, Eigen::MatrixXd& A);


    // FACE-BASED VALUES
    // Helper function to get vector area and stores normalized in fNormal
    // Returns magnitude of vector area
    double vectorArea(const std::vector<Eigen::Vector3d>& f, Eigen::Vector3d& fNormal);

    // compute barycenter (c)
    Eigen::Vector3d computeBarycenter(const std::vector<Eigen::Vector3d>& f);
    Eigen::Vector3d computeBarycenter(const Eigen::MatrixXd& X);

    // compute skew-symmetric cross product matrix from a vector
    Eigen::Matrix3d vectorCrossProd(const Eigen::Vector3d& vec);

    // VERTEX-BASED VALUES
    // Area of barycentric dual for a vertex
    // TODO: Add mesh data struct as input
    double vertexArea(int v_idx);
} // namespace DECUtils