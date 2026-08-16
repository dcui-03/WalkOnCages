#include "decUtils.hpp"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <mesh/mesh.hpp>
#include <vector>

namespace DECUtils {
    // For input vertex lists, assume that they are correctly order CCW and closed
    // NOTE: Since we are doing this face-based, we don't need half edges
    // Polygonal Mesh Laplacian operator, which assemble's face-based Laplacians
    Eigen::SparseMatrix<double> LaplacianOp(const std::vector<Eigen::Vector3d>& V,
                                            const std::vector<std::vector<int>>& F,
                                            double lambda) {
        int num_faces = F.size();
        int n = V.size();
        Eigen::SparseMatrix<double> L(n, n);
        // Get number of entries
        int num_entries = 0;
        for (int f_idx = 0; f_idx < num_faces; f_idx++) {
            num_entries += F[f_idx].size()*F[f_idx].size();
        }
        std::vector<T> tripletList;
        tripletList.reserve(num_entries);
        // Construct by face and assemble
        for (int f_idx = 0; f_idx < num_faces; f_idx++) {
            // compute per-face Laplacian
            std::vector<Eigen::Vector3d> f(F[f_idx].size());
            for (int v_idx = 0; v_idx < F[f_idx].size(); v_idx++) {
                f[v_idx] = V[F[f_idx][v_idx]];
            }
            Eigen::MatrixXd faceL = faceLaplacianOp(f, lambda);
            for (int i = 0; i < faceL.rows(); i++) {
                for (int j = 0; j < faceL.cols(); j++) {
                    tripletList.push_back(T(F[f_idx][i], F[f_idx][j], faceL(i, j)));
                }
            }
        }
        L.setFromTriplets(tripletList.begin(), tripletList.end());
        return L;
    }


    // All functions below are FACE-BASED
    // HIGHER OPERATORS
    // Face Laplacian operator
    Eigen::MatrixXd faceLaplacianOp(const std::vector<Eigen::Vector3d>& f, double lambda) {
        Eigen::SparseMatrix<double> D = diffOp(f.size());
        Eigen::MatrixXd M = metricOp(f, lambda);
        return D.transpose() * M * D;
    }
    Eigen::MatrixXd faceLaplacianOp(const Eigen::MatrixXd& D, const Eigen::MatrixXd& M) {
        return D.transpose() * M * D;
    }

    // Face Divergence operator (d0*)
    Eigen::MatrixXd divOp(const std::vector<Eigen::Vector3d>& f, double lambda) {
        Eigen::SparseMatrix<double> D = diffOp(f.size());
        Eigen::MatrixXd M = metricOp(f, lambda);
        return D.transpose() * M;
    }
    Eigen::MatrixXd divOp(const Eigen::MatrixXd& D, const Eigen::MatrixXd& M) {
        return D.transpose() * M;
    }

    // Face Curl operator (aka d1)
    Eigen::VectorXd curlOp(int n) {
        return Eigen::VectorXd::Ones(n);
    }

    // Face-based Metric on 1-forms (*)
    Eigen::MatrixXd metricOp(const std::vector<Eigen::Vector3d>& f, double lambda) {
        Eigen::Vector3d fNormal;
        double area = vectorArea(f, fNormal);
        Eigen::MatrixXd U = sharpOp(f);
        Eigen::MatrixXd V = flatOp(f);
        Eigen::MatrixXd P = projOp(V, U);
        return area * (U.transpose() * U) + lambda * (P.transpose() * P);
    }
    Eigen::MatrixXd metricOp(double area, const Eigen::MatrixXd& V, const Eigen::MatrixXd& U, double lambda) {
        Eigen::MatrixXd P = projOp(V, U);
        return area * (U.transpose() * U) + lambda * (P.transpose() * P);
    }


    // MUSICAL ISOMORPHISMS
    // Face Sharp operator
    // I.e., generate a face gradient vector for function given 1-forms
    Eigen::MatrixXd sharpOp(const std::vector<Eigen::Vector3d>& f) {
        int n = f.size();
        Eigen::Vector3d fNormal;
        double area = vectorArea(f, fNormal);
        Eigen::MatrixXd B = midpntOp(f);
        Eigen::Vector3d barycenter = computeBarycenter(f);
        double hodge_star;
        if (area <= 1e-8) { // Safety
            hodge_star = 0.0;
        } else {
            hodge_star = 1/area;
        }
        return hodge_star * vectorCrossProd(fNormal) * (B.transpose() - barycenter * curlOp(n).transpose());
    }
    Eigen::MatrixXd sharpOp(double area, Eigen::Vector3d& fNormal, Eigen::MatrixXd& B, Eigen::Vector3d& barycenter) {
        int n = B.rows();
        double hodge_star;
        if (area <= 1e-8) { // Safety
            hodge_star = 0.0;
        } else {
            hodge_star = 1/area;
        }
        return hodge_star * vectorCrossProd(fNormal) * (B.transpose() - barycenter * curlOp(n).transpose());
    }

    // Face Flat operator
    // I.e., generate a 1-form from face-based gradient
    Eigen::MatrixXd flatOp(const std::vector<Eigen::Vector3d>& f) {
        Eigen::Vector3d fNormal;
        double area = vectorArea(f, fNormal);
        Eigen::MatrixXd E = edgeOp(f);
        return E * (Eigen::MatrixXd::Identity(3, 3) - fNormal * fNormal.transpose());
    }
    Eigen::MatrixXd flatOp(Eigen::Vector3d& fNormal, Eigen::MatrixXd& E) {
        return E * (Eigen::MatrixXd::Identity(3, 3) - fNormal * fNormal.transpose());
    }

    // Face Projection operator
    Eigen::MatrixXd projOp(const std::vector<Eigen::Vector3d>& f) {
        int n = f.size();
        Eigen::MatrixXd U = sharpOp(f);
        Eigen::MatrixXd V = flatOp(f);
        return Eigen::MatrixXd::Identity(n, n) - (V * U);
    }
    Eigen::MatrixXd projOp(const Eigen::MatrixXd& V, const Eigen::MatrixXd& U) {
        int n = V.rows();
        return Eigen::MatrixXd::Identity(n, n) - (V * U);
    }


    // FUNDAMENTAL OPERATORS
    // Get face Position operator
    Eigen::MatrixXd posOp(const std::vector<Eigen::Vector3d>& f) {
        int n = f.size();
        Eigen::MatrixXd pos(n, 3);
        for (int v = 0; v < n; v++) {
            pos.row(v) = f[v].transpose();
        }
        return pos;
    }

    // Get face difference operator (aka d0)
    Eigen::SparseMatrix<double> diffOp(int n) {
        Eigen::SparseMatrix<double> diff(n, n);
        std::vector<T> tripletList;
        tripletList.reserve(n * n);
        for (int v = 0; v < n; v++) {
            tripletList.push_back(T(v, v, -1.0));
            tripletList.push_back(T(v, (v+1)%n, 1.0));
        }
        diff.setFromTriplets(tripletList.begin(), tripletList.end());
        return diff;
    }

    // Get face average operator
    Eigen::SparseMatrix<double> avgOp(int n) {
        Eigen::SparseMatrix<double> avg(n, n);
        std::vector<T> tripletList;
        tripletList.reserve(n * n);
        for (int v = 0; v < n; v++) {
            tripletList.push_back(T(v, v, 0.5));
            tripletList.push_back(T(v, (v+1)%n, 0.5));
        }
        avg.setFromTriplets(tripletList.begin(), tripletList.end());
        return avg;
    }

    // Get face midpoint operator
    Eigen::MatrixXd midpntOp(const std::vector<Eigen::Vector3d>& f) {
        Eigen::MatrixXd X = posOp(f);
        Eigen::SparseMatrix<double> A = avgOp(f.size());
        Eigen::MatrixXd B = A * X;
        return B;
    }
    Eigen::MatrixXd midpntOp(const Eigen::SparseMatrix<double>& A, const Eigen::MatrixXd& X) {
        return A * X;
    }

    // Get face edge operator
    Eigen::MatrixXd edgeOp(const std::vector<Eigen::Vector3d>& f) {
        Eigen::MatrixXd X = posOp(f);
        Eigen::SparseMatrix<double> D = diffOp(f.size());
        Eigen::MatrixXd E = D * X;
        return E;
    }
    Eigen::MatrixXd edgeOp(const Eigen::SparseMatrix<double>& D, const Eigen::MatrixXd& X) {
        return D * X;
    }

    // Face Gradient operator
    Eigen::MatrixXd gradOp(const std::vector<Eigen::Vector3d>& f) {
        Eigen::SparseMatrix<double> A = avgOp(f.size());
        Eigen::MatrixXd E = edgeOp(f);
        Eigen::Vector3d fNormal;
        double area = vectorArea(f, fNormal);

        return (-1/area) * vectorCrossProd(fNormal) * E.transpose() * A;
    }
    Eigen::MatrixXd gradOp(double area, Eigen::Vector3d& fNormal, Eigen::MatrixXd& E, Eigen::MatrixXd& A) {
        return (-1/area) * vectorCrossProd(fNormal) * E.transpose() * A;
    }


    // FACE-BASED VALUES
    // Helper function to get vector area and stores normalized in fNormal
    // Returns magnitude of vector area
    double vectorArea(const std::vector<Eigen::Vector3d>& f, Eigen::Vector3d& fNormal) {
        fNormal = Eigen::Vector3d::Zero();
        int n = f.size();
        for (int v = 0; v < n; v++) {
            fNormal += f[v].cross(f[(v+1)%n]);
        }
        fNormal *= 0.5;
        double area = fNormal.norm();
        if (area > 1e-8) {
            fNormal /= area;
        } else {
            fNormal = Eigen::Vector3d::Zero();
        }
        return area;
    }

    // compute barycenter (c)
    Eigen::Vector3d computeBarycenter(const std::vector<Eigen::Vector3d>& f) {
        int n = f.size();
        Eigen::MatrixXd X = posOp(f);
        return (1.0/n) * (X.transpose() * Eigen::VectorXd::Ones(n));
    }
    Eigen::Vector3d computeBarycenter(const Eigen::MatrixXd& X) {
        int n = X.rows();
        return (1.0/n) * (X.transpose() * Eigen::VectorXd::Ones(n));
    }

    // compute skew-symmetric cross product matrix from a vector
    Eigen::Matrix3d vectorCrossProd(const Eigen::Vector3d& vec) {
        Eigen::Matrix3d vecCP = Eigen::Matrix3d::Zero();
        vecCP(0, 1) = -1 * vec[2];
        vecCP(0, 2) = vec[1];

        vecCP(1, 0) = vec[2];
        vecCP(1, 2) = -1 * vec[0];

        vecCP(2, 0) = -1 * vec[1];
        vecCP(2, 1) = vec[0];
        return vecCP;
    }

    // VERTEX-BASED VALUES
    // Area of barycentric dual for a vertex
    // TODO: Add mesh data struct as input
    double vertexArea(int v_idx) {
        double area = 0.0;
        return area;
    }
} // namespace DECUtils