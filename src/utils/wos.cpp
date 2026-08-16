#include "wos.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <../curvenet/curvenet.hpp>
#include <../mesh/mesh.hpp>
#include <vector>
#include <random>
#include <cmath>


// Functions for Walk on Spheres sampling

namespace WoS {
    // Regular Walk on Spheres given a point and a mesh boundary
    int WalkOnSpheres_Mesh(const Eigen::Vector3d& p, const Mesh::mesh* M, Mesh::meshBindData& hitData, const int iter = 0, const int max_iter = 20, double eps = 1e-6) {
        if (iter >= max_iter) {
            return -1;
        }
        // First, find the closest point to the mesh and its coordinates
        M->computeVBinding(p, hitData);
        // Compute the distance
        double d = (bindData.proj - p).norm();
        // If we are too close, terminate
        if (d <= eps) {
            return 1;
        }
        // Generate random sample
        Eigen::Vector3d newDirec;
        generateNewDirection(newDirec);
        // Compute the next sample
        Eigen::Vector3d p_next = p + d * newDirec;
        // Recurse
        return WalkOnSpheres_Mesh(p_next, M, hitData, iter+1, max_iter, eps);
    }

    // Walk on Spheres on a curve network
    int WalkOnSpheres_Curvenet(const Eigen::Vector3d& p, const Curvenet::curvenet* CN, const int iter = 0, const int max_iter = 50) {
        if (iter >= max_iter) {
            return -1;
        }
        // First, find the closest point to the curvenet and its coordinates
        // TODO

        // Compute the distance
        double d = (bindData.proj - p).norm();
        // If we are too close, terminate
        if (d <= eps) {
            return 1;
        }
        // Generate random sample
        Eigen::Vector3d newDirec;
        generateNewDirection(newDirec);
        // Compute the next sample
        Eigen::Vector3d p_next = p + d * newDirec;
        // Recurse
        return WalkOnSpheres_Mesh(p_next, M, hitData, iter+1, max_iter, eps);
    }

    // Walk on Sphere on a polyline network
    // WalkOnSphere_Polynet(const Eigen::Vector3d& p,  Curvenet::curvnet* CN, const int iter = 0, const int max_iter = 20);

    // Generate a new random walk direction
    int generateNewDirection(Eigen::Vector3d& newDirec) {
        newDirec.setZero();

        std::random_device rand;
        std::mt19937 gen(rand());
        std::uniform_real_distribution<double> theta(0.0, 2 * M_PI);
        std::uniform_real_distribution<double> z(-1.0, 1.0);
        double r = sqrt(pow(1 - z, 2));

        newDirec = {r*std::cos(theta), r*std::sin(theta), z};
        return 1;
    }
} // namespace WoS