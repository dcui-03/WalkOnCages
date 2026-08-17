// polynetcage.hpp
#pragma once

#include "cage/cage.hpp"
#include "polynet/polynet.hpp"
#include <Eigen/Core>
#include <vector>
#include <array>
#include <map>
#include <random>

namespace Cage {

// polynet cage object
class polynetcage : public cage {
    public:
        // Init using curvneet object
        polynetcage(Polynet::polynet* PN);

        // Function for retrieving verices as a matrix
        virtual int matrixVerts(Eigen::MatrixXd& Verts) const override;

        // Function for querying closest point
        virtual int closestPoint(const Eigen::Vector3d& p, Utils::projData& hit) const override;

        // Function for computing basis of mesh element
        virtual std::vector<std::pair<int, double>> computeBasis(const Utils::projData& proj) const override;

        // Debug colors: direction from the centroid of all verts, remapped to [0,1]
        virtual int computeColors(Eigen::MatrixXd& Colors) const override;

        // Direction toward a random point on the polyline. Guarantees a raycast along it hits at least the targeted edge
        virtual int sampleDirection(const Eigen::Vector3d& q_pos, std::mt19937& gen, Eigen::Vector3d& direc) const override;

        // Cast a ray against the polyline
        virtual int raycast(const Eigen::Vector3d& origin, const Eigen::Vector3d& direc, std::vector<Utils::projData>& hits, double tol) const override;

        // Get bounding box diagonal length
        virtual double bboxDiag() const override;
    protected:

    private:
        // Concatenated per-edge lengths as global separators
        std::vector<double> totalArclen;
        // Build totalArclen from PN's edges
        void computeArclenTable();
        // Draw a uniformly-random-by-arclength point on the polyline
        int sampleArclengthPoint(std::mt19937& gen, Eigen::Vector3d& pos) const;

        // Point to curvenet
        Polynet::polynet* PN;
};

}   // namespace Cage