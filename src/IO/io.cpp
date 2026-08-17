#include "io.hpp"

#include <vector>
#include <array>
#include <set>
#include <utility>
#include <Eigen/Dense>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

namespace IO {

    static int parseOBJVertexIndex(const std::string& token, int numVerticesSoFar) {
        if (token.empty()) {
            throw std::runtime_error("Empty face token in OBJ.");
        }

        // Extract substring before first '/'
        auto slashPos = token.find('/');
        std::string vStr = (slashPos == std::string::npos) ? token : token.substr(0, slashPos);

        if (vStr.empty()) {
            throw std::runtime_error("Malformed OBJ face token: " + token);
        }

        int idx = std::stoi(vStr);

        // OBJ indices are 1-based.
        // Negative indices are relative to the end: -1 means last defined vertex.
        if (idx > 0) {
            return idx - 1; // convert to 0-based
        } else if (idx < 0) {
            return numVerticesSoFar + idx; // since idx is negative
        } else {
            throw std::runtime_error("OBJ index 0 is invalid.");
        }
    }

    // Reads OBJ into std::vector<Eigen::Vector3d> and polygon face list
    bool readOBJ(
        const std::string& filename,
        std::vector<Eigen::Vector3d>& V,
        std::vector<std::vector<int>>& F
    ) {
        V.clear();
        F.clear();

        std::ifstream in(filename);
        if (!in) {
            std::cerr << "Failed to open OBJ file: " << filename << std::endl;
            return false;
        }

        std::string line;
        while (std::getline(in, line)) {
            if (line.empty()) continue;

            std::istringstream iss(line);
            std::string tag;
            iss >> tag;

            if (tag.empty() || tag[0] == '#') {
                continue;
            }

            if (tag == "v") {
                double x, y, z;
                if (!(iss >> x >> y >> z)) {
                    std::cerr << "Malformed vertex line: " << line << std::endl;
                    return false;
                }
                V.emplace_back(x, y, z);
            }
            else if (tag == "f") {
                std::vector<int> face;
                std::string token;

                while (iss >> token) {
                    int idx;
                    try {
                        idx = parseOBJVertexIndex(token, static_cast<int>(V.size()));
                    } catch (const std::exception& e) {
                        std::cerr << "Error parsing face token \"" << token
                                << "\": " << e.what() << std::endl;
                        return false;
                    }

                    if (idx < 0 || idx >= static_cast<int>(V.size())) {
                        std::cerr << "Face index out of range in line: " << line << std::endl;
                        return false;
                    }

                    face.push_back(idx);
                }

                if (face.size() < 3) {
                    std::cerr << "Face has fewer than 3 vertices: " << line << std::endl;
                    return false;
                }

                F.push_back(std::move(face));
            }

            // Ignore everything else: vt, vn, usemtl, mtllib, g, o, ...
        }

        return true;
    }

    // Overload: reads OBJ into Eigen::MatrixXd and polygon face list
    bool readOBJ(
        const std::string& filename,
        Eigen::MatrixXd& V,
        std::vector<std::vector<int>>& F
    ) {
        std::vector<Eigen::Vector3d> Vvec;
        if (!readOBJ(filename, Vvec, F)) {
            V.resize(0, 3);
            return false;
        }

        V.resize(static_cast<Eigen::Index>(Vvec.size()), 3);
        for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(Vvec.size()); ++i) {
            V.row(i) = Vvec[i].transpose();
        }

        return true;
    }

    // Extracts the unique undirected edges of a face list, for wireframe display
    void facesToWireframe(
        const std::vector<std::vector<int>>& F,
        std::vector<std::array<int, 2>>& E
    ) {
        E.clear();
        std::set<std::pair<int, int>> seen;
        for (const std::vector<int>& face : F) {
            int n = static_cast<int>(face.size());
            for (int i = 0; i < n; i++) {
                int a = face[i];
                int b = face[(i + 1) % n];
                std::pair<int, int> key = (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
                if (seen.insert(key).second) {
                    E.push_back({a, b});
                }
            }
        }
        return;
    }
} // namespace IO