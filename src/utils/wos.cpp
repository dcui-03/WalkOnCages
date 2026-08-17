#define _USE_MATH_DEFINES
#include "wos.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "cage/cage.hpp"
#include "mesh/mesh.hpp"
#include <vector>
#include <random>
#include <cmath>


namespace {
    // Radical inverse (digit-reversal) of index in the given base: the building block of a Halton low-discrepancy sequence
    double haltonRadicalInverse(int index, int base) {
        double result = 0.0;
        double f = 1.0 / base;
        while (index > 0) {
            result += f * (index % base);
            index /= base;
            f /= base;
        }
        return result;
    }
}

// Functions for Walk on Spheres sampling
namespace WoS {
    // Regular Walk on Spheres given a point and a mesh boundary
    int WalkOnSpheres(const Eigen::Vector3d& p, const Cage::cage* C,
                           Utils::projData& hit,
                           std::mt19937& gen,
                           const int iter, const int max_iter, double eps) {
        if (iter >= max_iter) {
            return -1;
        }
        // First, find the closest point and its coordinates
        if (C->closestPoint(p, hit) != 1) {
            return -1;
        }
        // Compute the distance
        double d = (hit.pos - p).norm();
        // If we are too close, terminate
        if (d <= eps) {
            return 1;
        }
        // Generate random sample
        Eigen::Vector3d newDirec;
        generateNewDirection(newDirec, gen);
        // Compute the next sample
        Eigen::Vector3d p_next = p + d * newDirec;
        // Recurse
        return WalkOnSpheres(p_next, C, hit, gen, iter+1, max_iter, eps);
    }

    // Generate a new random walk direction
    int generateNewDirection(Eigen::Vector3d& newDirec, std::mt19937& gen) {
        newDirec.setZero();

        std::uniform_real_distribution<double> theta_dist(0.0, 2 * M_PI);
        std::uniform_real_distribution<double> z_dist(-1.0, 1.0);
        double theta = theta_dist(gen);
        double z = z_dist(gen);
        double r = std::sqrt(std::max(0.0, 1.0 - z * z));

        newDirec = {r*std::cos(theta), r*std::sin(theta), z};
        return 1;
    }

    // Deterministic, well-spread direction for the index-th sample
    int stratifySamples(int index, Eigen::Vector3d& newDirec) {
        double u = haltonRadicalInverse(index + 1, 2);
        double v = haltonRadicalInverse(index + 1, 3);
        double theta = 2 * M_PI * u;
        double z = 2.0 * v - 1.0;
        double r = std::sqrt(std::max(0.0, 1.0 - z * z));

        newDirec = {r*std::cos(theta), r*std::sin(theta), z};
        return 1;
    }
} // namespace WoS
