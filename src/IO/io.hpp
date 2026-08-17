// IO.hpp
#pragma once

#include <vector>
#include <array>
#include <Eigen/Dense>
#include <vector>
#include <string>

namespace IO {

    static int parseOBJVertexIndex(const std::string& token, int numVerticesSoFar);

    bool readOBJ(
        const std::string& filename,
        std::vector<Eigen::Vector3d>& V,
        std::vector<std::vector<int>>& F
    );

    bool readOBJ(
        const std::string& filename,
        Eigen::MatrixXd& V,
        std::vector<std::vector<int>>& F
    );

    // Extracts the unique undirected edges of a face list, for wireframe display
    void facesToWireframe(
        const std::vector<std::vector<int>>& F,
        std::vector<std::array<int, 2>>& E
    );

} // namespace IO