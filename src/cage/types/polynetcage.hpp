// polynetcage.hpp
#pragma once

#include "cage/cage.hpp"
#include "polynet/polynet.hpp"
#include <Eigen/Core>
#include <vector>
#include <array>
#include <map>

namespace Cage {

// polynet cage object
class polynetcage : public cage {
    public:
        // Init using curvneet object
        polynetcage(Polynet::polynet* PN);

        // Function for retrieving verices as a matrix
        virtual int matrixVerts(Eigen::MatrixXd& Verts) override;

        // Function for querying closest point
        virtual int closestPoint(const Eigen::Vector3d& p, int& elType, int& elIdx, Eigen::Vector3d& proj, Eigen::VectorXd& coords) override;

        // Function for computing basis of mesh element
        virtual std::vector<std::pair<int, double>> computeBasis(const int& elType, const int& elIdx, const Eigen::VectorXd& coords) override;

        // Debug colors: direction from the centroid of all verts, remapped to [0,1]
        virtual int computeColors(Eigen::MatrixXd& Colors) override;

        // Get bounding box diagonal length
        virtual double bboxDiag() const override;
    protected:

    private:
        // Point to curvenet
        Polynet::polynet* PN;
};

}   // namespace Cage