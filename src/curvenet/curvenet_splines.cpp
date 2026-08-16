#include "curvenet_splines.hpp"

#include "utils/utils.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <algorithm>

// File with basic structs used by mesh class

namespace Curvenet {
/*
Spline::Spline() {
}
void Spline::setHE(int new_he) {
    if (new_he < 0) {
        he = -1;
    } else {
        he = new_he;
    }
    return;
}

void Spline::setCurve(int new_crv) {
    if (new_crv < 0) {
        curve = -1;
    } else {
        curve = new_crv;
    }
    return;
}

bool Spline::setNumSamples(int n) {
    if (n < 3) {
        return false;
    }
    num_samples = n;
    return true;
}
// Compute num samples using a ratio b/w arclen. estimate and a parameter alpha
int Spline::compute_num_samples(double arclen, double meanE, int alpha) {
    if (meanE <= 1e-16) return 3;
    return std::max(3, static_cast<int>(std::ceil(alpha * arclen / meanE)));
}

CubicBezier::CubicBezier() {
}
Eigen::Vector3d tSample(double t);
std::vector<Eigen::Vector3d> sampleNaive(int num_samples = 50);
double arclenEst(int num_samples = 50);
std::vector<Eigen::Vector3d> unifSamples(int num_samples) {

}


CatmullRom::CatmullRom() {   // Constructor
}
Eigen::Vector3d tSample(double t);
std::vector<Eigen::Vector3d> sampleNaive(int num_samples = 50);
double arclenEst(int num_samples = 50);
std::vector<Eigen::Vector3d> unifSamples(int num_samples);  
*/
}   // namespace Curvenet