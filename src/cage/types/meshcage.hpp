// meshcage.hpp
#pragma once

#include "cage/cage.hpp"
#include "mesh/mesh.hpp"
#include <Eigen/Core>
#include <vector>
#include <array>
#include <map>

namespace Cage {

// mesh cage object
class meshcage : public cage {
    public:
        // Init using mesh object
        meshcage(Mesh::mesh* M);

        // Function for retrieving verices as a matrix
        virtual int matrixVerts(Eigen::MatrixXd& Verts) const override;

        // Function for querying closest point
        virtual int closestPoint(const Eigen::Vector3d& p, int& elType, int& elIdx, Eigen::Vector3d& proj, Eigen::VectorXd& coords) const override;

        // Function for computing basis of mesh element
        virtual std::vector<std::pair<int, double>> computeBasis(const int& elType, const int& elIdx, const Eigen::VectorXd& coords) const override;

        // Debug colors: vertex normal remapped to [0,1]
        virtual int computeColors(Eigen::MatrixXd& Colors) const override;

        // Get bounding box diagonal length
        virtual double bboxDiag() const override;
    protected:

    private:
        // Point to mesh
        Mesh::mesh* M;
};

}   // namespace Cage