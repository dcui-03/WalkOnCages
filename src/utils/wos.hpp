// wos.hpp
#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include "cage/cage.hpp"
#include <mesh/mesh.hpp>
#include <vector>


// Functions for Walk on Spheres sampling

namespace WoS {
    // Regular Walk on Spheres given a point and a boundary
    int WalkOnSpheres(const Eigen::Vector3d& p, const Cage::cage* C, 
                           int& elType, int& elIdx, Eigen::Vector3d& proj, Eigen::VectorXd& coords,
                           const int iter = 0, const int max_iter = 20, double eps = 1e-6);

    int generateNewDirection(Eigen::Vector3d& newDirec);
} // namespace WoS