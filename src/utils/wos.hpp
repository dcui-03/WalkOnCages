// wos.hpp
#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include "cage/cage.hpp"
#include "mesh/mesh.hpp"
#include <vector>
#include <random>


// Functions for Walk on Spheres sampling

namespace WoS {
    // Regular Walk on Spheres given a point and a boundary
    // gen is caller-owned (e.g. one persistent generator per thread), never seeded internally
    int WalkOnSpheres(const Eigen::Vector3d& p, const Cage::cage* C,
                           Utils::projData& hit,
                           std::mt19937& gen,
                           const int iter = 0, const int max_iter = 20, double eps = 1e-6);

    int generateNewDirection(Eigen::Vector3d& newDirec, std::mt19937& gen);

    // Deterministic, well-spread direction for the index-th sample (a 2D Halton sequence,
    // bases 2 and 3, mapped onto the sphere). Used to reduce clustering across a query point's
    // first-hop directions vs. drawing them i.i.d. at random
    int stratifySamples(int index, Eigen::Vector3d& newDirec);
} // namespace WoS