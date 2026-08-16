// pscurvenet_io.cpp

#include "pscurvenet.hpp"

#include <Eigen/Core>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <array>

namespace psCurvenet {

namespace {

bool readNextDataLine(std::istream& in, std::string& line) {
    while (std::getline(in, line)) {
        // Strip comments.
        const std::size_t hash = line.find('#');
        if (hash != std::string::npos) {
            line = line.substr(0, hash);
        }

        // Skip whitespace-only lines.
        std::istringstream ss(line);
        std::string first;
        if (ss >> first) {
            return true;
        }
    }

    return false;
}

bool parseVec3Line(
    const std::string& line,
    const std::string& expectedTag,
    Eigen::Vector3d& v
) {
    std::istringstream ss(line);

    std::string tag;
    double x, y, z;

    if (!(ss >> tag >> x >> y >> z)) {
        return false;
    }

    if (tag != expectedTag) {
        return false;
    }

    v = Eigen::Vector3d(x, y, z);
    return true;
}

bool parseCountLine(
    const std::string& line,
    const std::string& expectedTag,
    int& count
) {
    std::istringstream ss(line);

    std::string tag;
    if (!(ss >> tag >> count)) {
        return false;
    }

    if (tag != expectedTag) {
        return false;
    }

    return count >= 0;
}

bool parseSplineLine(
    const std::string& line,
    std::array<int, 4>& s
) {
    std::istringstream ss(line);

    if (!(ss >> s[0] >> s[1] >> s[2] >> s[3])) {
        return false;
    }

    return true;
}

} // anonymous namespace


int pscurvenet::saveCurvenet(const std::string& filepath) const {
    std::ofstream out(filepath);

    if (!out.is_open()) {
        std::cout << "pscurvenet::saveCurvenet(): failed to open file: "
                  << filepath << std::endl;
        return -1;
    }

    out << std::setprecision(17);

    // Metadata header. Leave mostly blank for now, but keep OBJ-like comments.
    out << "# curvenet file\n";
    out << "#\n";

    // Controls.
    out << "c " << C.size() << "\n";

    for (int i = 0; i < static_cast<int>(C.size()); i++) {
        out << "v "
            << C[i].pos(0) << " "
            << C[i].pos(1) << " "
            << C[i].pos(2) << "\n";
    }

    // One normal per control, in matching order.
    for (int i = 0; i < static_cast<int>(C.size()); i++) {
        out << "n "
            << C[i].n(0) << " "
            << C[i].n(1) << " "
            << C[i].n(2) << "\n";
    }

    // Tangents. We intentionally write two tangent points per spline.
    // Tangent index 2*s is spline s's t0.
    // Tangent index 2*s + 1 is spline s's t1.
    const int numTangents = 2 * static_cast<int>(S.size());

    out << "t " << numTangents << "\n";

    for (int s = 0; s < static_cast<int>(S.size()); s++) {
        out << "v "
            << S[s].t0(0) << " "
            << S[s].t0(1) << " "
            << S[s].t0(2) << "\n";

        out << "v "
            << S[s].t1(0) << " "
            << S[s].t1(1) << " "
            << S[s].t1(2) << "\n";
    }

    // Splines.
    out << "s " << S.size() << "\n";

    for (int s = 0; s < static_cast<int>(S.size()); s++) {
        const int t0Idx = 2 * s;
        const int t1Idx = 2 * s + 1;

        out << S[s].start << " "
            << t0Idx << " "
            << t1Idx << " "
            << S[s].end << "\n";
    }

    if (!out.good()) {
        std::cout << "pscurvenet::saveCurvenet(): write failed for file: "
                  << filepath << std::endl;
        return -1;
    }

    return 1;
}


int pscurvenet::loadCurvenet(const std::string& filepath) {
    std::ifstream in(filepath);

    if (!in.is_open()) {
        std::cout << "pscurvenet::loadCurvenet(): failed to open file: "
                  << filepath << std::endl;
        return -1;
    }

    std::string line;

    // Controls count.
    int numControls = 0;
    if (!readNextDataLine(in, line) || !parseCountLine(line, "c", numControls)) {
        std::cout << "pscurvenet::loadCurvenet(): expected line `c <num controls>`."
                  << std::endl;
        return -1;
    }

    std::vector<Control> newC(numControls);

    // Control positions.
    for (int c = 0; c < numControls; c++) {
        Eigen::Vector3d p;

        if (!readNextDataLine(in, line) || !parseVec3Line(line, "v", p)) {
            std::cout << "pscurvenet::loadCurvenet(): expected control position line `v x y z` at control "
                      << c << "." << std::endl;
            return -1;
        }

        newC[c].pos = p;
        newC[c].active = true;
    }

    // Control normals.
    for (int c = 0; c < numControls; c++) {
        Eigen::Vector3d n;

        if (!readNextDataLine(in, line) || !parseVec3Line(line, "n", n)) {
            std::cout << "pscurvenet::loadCurvenet(): expected normal line `n x y z` at control "
                      << c << "." << std::endl;
            return -1;
        }

        if (n.squaredNorm() > 1e-20) {
            newC[c].n = n.normalized();
        } else {
            newC[c].n = Eigen::Vector3d::UnitZ();
        }
    }

    // Tangents count.
    int numTangents = 0;
    if (!readNextDataLine(in, line) || !parseCountLine(line, "t", numTangents)) {
        std::cout << "pscurvenet::loadCurvenet(): expected line `t <num tangents>`."
                  << std::endl;
        return -1;
    }

    std::vector<Eigen::Vector3d> tangents(numTangents);

    for (int t = 0; t < numTangents; t++) {
        Eigen::Vector3d p;

        if (!readNextDataLine(in, line) || !parseVec3Line(line, "v", p)) {
            std::cout << "pscurvenet::loadCurvenet(): expected tangent position line `v x y z` at tangent "
                      << t << "." << std::endl;
            return -1;
        }

        tangents[t] = p;
    }

    // Splines count.
    int numSplines = 0;
    if (!readNextDataLine(in, line) || !parseCountLine(line, "s", numSplines)) {
        std::cout << "pscurvenet::loadCurvenet(): expected line `s <num splines>`."
                  << std::endl;
        return -1;
    }

    std::vector<Spline> newS(numSplines);

    for (int s = 0; s < numSplines; s++) {
        std::array<int, 4> idx;

        if (!readNextDataLine(in, line) || !parseSplineLine(line, idx)) {
            std::cout << "pscurvenet::loadCurvenet(): expected spline line "
                      << "`startControl t0 t1 endControl` at spline "
                      << s << "." << std::endl;
            return -1;
        }

        const int c0 = idx[0];
        const int t0 = idx[1];
        const int t1 = idx[2];
        const int c1 = idx[3];

        if (c0 < 0 || c0 >= numControls ||
            c1 < 0 || c1 >= numControls ||
            t0 < 0 || t0 >= numTangents ||
            t1 < 0 || t1 >= numTangents) {
            std::cout << "pscurvenet::loadCurvenet(): invalid spline indices at spline "
                      << s << ": "
                      << c0 << " " << t0 << " " << t1 << " " << c1
                      << std::endl;
            return -1;
        }

        newS[s].start = c0;
        newS[s].end = c1;
        newS[s].t0 = tangents[t0];
        newS[s].t1 = tangents[t1];
        newS[s].active = true;
    }

    // Only commit after the entire file has parsed successfully.
    C.swap(newC);
    S.swap(newS);

    recomputeMap = true;
    constructTangentMap();

    return 1;
}

} // namespace psCurvenet