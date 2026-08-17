// meshcage.hpp
#pragma once

#include "cage/cage.hpp"
#include "mesh/mesh.hpp"
#include <Eigen/Core>
#include <vector>
#include <array>
#include <map>
#include <random>

namespace Cage {

// mesh cage object
class meshcage : public cage {
    public:
        // Init using mesh object
        meshcage(Mesh::mesh* M);

        // Function for retrieving verices as a matrix
        virtual int matrixVerts(Eigen::MatrixXd& Verts) const override;

        // Function for querying closest point
        virtual int closestPoint(const Eigen::Vector3d& p, Utils::projData& hit) const override;

        // Function for computing basis of mesh element
        virtual std::vector<std::pair<int, double>> computeBasis(const Utils::projData& proj) const override;

        // Debug colors: vertex normal remapped to [0,1]
        virtual int computeColors(Eigen::MatrixXd& Colors) const override;

        // Uniform random direction on the sphere
        virtual int sampleDirection(const Eigen::Vector3d& q_pos, std::mt19937& gen, Eigen::Vector3d& direc) const override;

        // Cast a ray against the mesh
        virtual int raycast(const Eigen::Vector3d& origin, const Eigen::Vector3d& direc, std::vector<Utils::projData>& hits, double tol) const override;

        // Get bounding box diagonal length
        virtual double bboxDiag() const override;
    protected:

    private:
        // Point to mesh
        Mesh::mesh* M;
};

}   // namespace Cage