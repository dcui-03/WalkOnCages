// query.hpp
#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>

namespace Query {

// Object to query barycentric coordinates of
class query {
    public:
        // Init with pointer to some object
        virtual ~query() = default;

        // Get relevant vertices as a matrix
        virtual int matrixVerts(Eigen::MatrixXd& Verts) = 0;

        // Apply smoothing
        virtual int applySmoothing(const Eigen::MatrixXd& Input, Eigen::MatrixXd& Result) = 0;

        // (Re)compute the smoothing operator; num_samples scales its strength
        virtual int computeSmoothingOp(int num_samples) = 0;
    protected:
        bool smoothAvailable = false;
    private:
        // Pointer to query object

        // Store smoothing operators for easy reuse
};

}   // namespace Query
