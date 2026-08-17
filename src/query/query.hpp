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
        virtual int matrixVerts(Eigen::MatrixXd& Verts);

        // Apply smoothing
        int applySmoothing(const Eigen::VectorXd& Input, Eigen::VectorXd& Result);
    protected:
        // Get a smoothing operator
        virtual int computeSmoothingOp();

        bool smoothAvailable = false;
    private:
        // Pointer to query object

        // Store smoothing operators for easy reuse
};

}   // namespace Query