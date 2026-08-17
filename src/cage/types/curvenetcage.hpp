// curvenetcage.hpp
#pragma once

#include "cage/cage.hpp"
#include "curvenet/curvenet.hpp"
#include "dcurvenet/dcurvenet.hpp"
#include <Eigen/Core>
#include <vector>
#include <array>
#include <map>

namespace Cage {

// curvenet cage object
class curvenetcage : public cage {
    public:
        // Init using curvneet object
        curvenetcage(Curvenet::curvenet* CN);

        // Function for retrieving verices as a matrix
        virtual int matrixVerts(Eigen::MatrixXd& Verts) override;

        // Function for querying closest point
        virtual int closestPoint(const Eigen::Vector3d& p, int& elType, int& elIdx, Eigen::Vector3d& proj, Eigen::VectorXd& coords) override;

        // Function for computing basis of mesh element
        virtual std::vector<std::pair<int, double>> computeBasis(const int& elType, const int& elIdx, const Eigen::VectorXd& coords) override;
    protected:

    private:
        // Point to curvenet
        Curvenet::curvenet* CN;
        Polynet::dcurvenet dCN:
};

}   // namespace Cage