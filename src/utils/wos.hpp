// wos.hpp
#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <curvenet/curvenet.hpp>
#include <mesh/mesh.hpp>
#include <vector>


// Functions for Walk on Spheres sampling

namespace WoS {
    // Regular Walk on Spheres given a point and a mesh boundary
    int WalkOnSpheres_Mesh(const Eigen::Vector3d& p, const Mesh::mesh* M, Mesh::meshBindData& hitData, const int iter = 0, const int max_iter = 20, double eps = 1e-6);

    // Walk on Spheres on a curve network
    int WalkOnSpheres_Curvenet(const Eigen::Vector3d& p, const Curvenet::curvenet* CN, const int iter = 0, const int max_iter = 50, double eps = 1e-7);

    // Walk on Sphere on a polyline network
    // WalkOnSphere_Polynet(const Eigen::Vector3d& p,  Curvenet::curvnet* CN, const int iter = 0, const int max_iter = 20);

    int generateNewDirection(Eigen::Vector3d& newDirec);
} // namespace WoS