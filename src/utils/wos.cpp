#include "wos.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "cage/cage.hpp"
#include <../mesh/mesh.hpp>
#include <vector>
#include <random>
#include <cmath>


// Functions for Walk on Spheres sampling
namespace WoS {
    // Regular Walk on Spheres given a point and a mesh boundary
    int WalkOnSpheres(const Eigen::Vector3d& p, const Cage::cage* C, 
                           int& elType, int& elIdx, Eigen::Vector3d& proj, Eigen::VectorXd& coords, 
                           const int iter = 0, const int max_iter = 20, double eps = 1e-6) {
        if (iter >= max_iter) {
            return -1;
        }
        // First, find the closest point and its coordinates
        C->closestPoint(elType, elIdx, proj, coords);
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
        return WalkOnSpheres(p_next, M, hitData, iter+1, max_iter, eps);
    }

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