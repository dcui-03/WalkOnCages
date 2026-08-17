#include "curvenet.hpp"

#include <Eigen/Core>
#include <cmath>
#include <algorithm>

namespace Curvenet {
    // Sample a bezier curve at time t
    Eigen::Vector3d curvenet::tSampleBezier(const Eigen::Vector3d& c0, const Eigen::Vector3d& c1, const Eigen::Vector3d& c2, const Eigen::Vector3d& c3, double t) const {
        Eigen::Vector3d sample = std::pow(1 - t, 3) * c0 +
                        3 * std::pow(1 - t, 2) * t * c1 +
                        3 * (1 - t) * std::pow(t, 2) * c2 +
                        std::pow(t, 3) * c3;
        return sample;
    }
    // First derivative of Bezier curve
    Eigen::Vector3d curvenet::tBezier_first(const Eigen::Vector3d& c0, const Eigen::Vector3d& c1, const Eigen::Vector3d& c2, const Eigen::Vector3d& c3, double t) const {
        Eigen::Vector3d deriv = std::pow(1.0 - t, 2) * (c1 - c0) + 
                        2 * (1.0 - t) * t * (c2 - c1) +
                        t * t * (c3 - c2);
        return 3 * deriv;
    }
    // Second derivative of Bezier curve
    Eigen::Vector3d curvenet::tBezier_second(const Eigen::Vector3d& c0, const Eigen::Vector3d& c1, const Eigen::Vector3d& c2, const Eigen::Vector3d& c3, double t) const {
        Eigen::Vector3d deriv = (1.0 - t) * (c2 - 2 * c1 + c0) + 
                        t * (c3 - 2 * c2 + c1);
        return 6 * deriv;
    }
    Eigen::Vector3d curvenet::tSampleBezier(int s, double t) const {
        Eigen::Vector3d c0, c1, c2, c3;
        splineCtrlPts(s, c0, c1, c2, c3);
        return tSampleBezier(c0, c1, c2, c3, t);
    }

    Eigen::Vector3d curvenet::tBezier_first(int s, double t) const {
        Eigen::Vector3d c0, c1, c2, c3;
        splineCtrlPts(s, c0, c1, c2, c3);
        return tBezier_first(c0, c1, c2, c3, t);
    }
    Eigen::Vector3d curvenet::tBezier_second(int s, double t) const {
        Eigen::Vector3d c0, c1, c2, c3;
        splineCtrlPts(s, c0, c1, c2, c3);
        return tBezier_second(c0, c1, c2, c3, t);
    }
    // Get the 4 control points of a spline
    void curvenet::splineCtrlPts(int s, Eigen::Vector3d& c0, Eigen::Vector3d& c1, Eigen::Vector3d& c2, Eigen::Vector3d& c3) const {
        int he = S[s].he;
        c0 = C[HE[he].origin].new_pos;
        c1 = HE[he].tan;
        c2 = HE[HE[he].twin].tan;
        c3 = C[HE[HE[he].twin].origin].new_pos;
    }

    // NOTE: This is a naive, fast sampler that uniformly samples t's. Re-implement if desired
    // Returns n_samples points on the curve, including the endpoints
    std::vector<Eigen::Vector3d> curvenet::sampleBezierNaive(int s, int n_samples) const {
        if (n_samples < 2) {
            return {};
        }
        std::vector<Eigen::Vector3d> samples(n_samples);
        int he = S[s].he;
        double h = 1.0/(n_samples - 1);
        double t = 0.0;

        samples[0] = C[HE[he].origin].new_pos;  // Start
        for (int i = 1; i < n_samples - 1; i++) {
            t += h;
            t = std::min(1.0, t);
            samples[i] = tSampleBezier(s, t);
        }
        samples[n_samples-1] = C[HE[HE[he].twin].origin].new_pos;   // End

        return samples;
    }

    // Estimate the arclength
    double curvenet::arclenEst(const std::vector<Eigen::Vector3d>& samples) const {
        if (samples.size() < 2) {
            return 0.0;
        }
        double length = 0.0;
        for (int i = 0; i < samples.size() - 1; i++) {
            length += (samples[i+1] - samples[i]).norm();
        }
        return length;
    }
    double curvenet::arclenEst(int s, int n_samples) const {
        std::vector<Eigen::Vector3d> samples = sampleBezierNaive(s, n_samples);
        return arclenEst(samples);
    }

    // Compute number of samples to take on a spline given a user parameter alpha
    int curvenet::computeNumSamples(double arclen, const Mesh::mesh* M) {
        if (!M) {
            // No mesh available: fall back to a fixed arclength-based heuristic
            return std::max(3, static_cast<int>(std::ceil(arclen / defaultSample)));
        }
        double meanE = M->getMeanE();
        if (meanE <= 1e-16) return 3;
        return std::max(3, static_cast<int>(std::ceil(alpha * arclen / meanE)));
    }

    // Uniformly sample based on arclength estimator and returns the length of the returned curve
    std::vector<Eigen::Vector3d> curvenet::unifSample(int s, int n_samples) const {
        std::vector<double> ts = unifSampleT(s, n_samples);
        std::vector<Eigen::Vector3d> samples(ts.size());
        for (int i = 0; i < static_cast<int>(ts.size()); ++i) {
            samples[i] = tSampleBezier(s, ts[i]);
        }
        return samples;
    }

    // Same as unifSample, but returns the t-values instead of the sampled points
    std::vector<double> curvenet::unifSampleT(int s, int n_samples) const {
        if (n_samples < 2) {
            return {};
        }

        std::vector<Eigen::Vector3d> regSamples = sampleBezierNaive(s, 50);
        std::vector<double> regT(regSamples.size());
        for (int i = 0; i < static_cast<int>(regT.size()); ++i) {
            regT[i] = static_cast<double>(i) / static_cast<double>(regT.size() - 1);
        }
        std::vector<double> unifT(n_samples);

        // Build cumulative arclength list
        std::vector<double> cumLen(regSamples.size(), 0.0);
        for (int i = 1; i < static_cast<int>(regSamples.size()); ++i) {
            cumLen[i] = cumLen[i - 1] + (regSamples[i] - regSamples[i - 1]).norm();
        }
        const double totalLen = cumLen.back();

        // Degenerate curve fallback
        if (totalLen <= 1e-16) {
            // If curve is ~0 length, just return the start
            for (int i = 0; i < n_samples; ++i) {
                unifT[i] = regT.front();
            }
            return unifT;
        }
        unifT[0] = regT.front();
        unifT[n_samples - 1] = regT.back();

        int seg = 1;
        for (int i = 1; i < n_samples - 1; ++i) {
            double targetLen = totalLen * static_cast<double>(i) / static_cast<double>(n_samples - 1);
            // Check if we overflow this segment
            while ((seg < cumLen.size() - 1) && (cumLen[seg] < targetLen)) {
                ++seg;
            }

            double segLen = cumLen[seg] - cumLen[seg - 1];
            if (segLen <= 1e-16) {
                unifT[i] = regT[seg];
                continue;
            }
            // Increment by computing the percentage difference
            double alpha = (targetLen - cumLen[seg - 1]) / segLen;
            alpha = std::max(0.0, std::min(1.0, alpha));
            unifT[i] = (1.0 - alpha) * regT[seg - 1] + alpha * regT[seg];
        }

        return unifT;
    }

}   // namespace Curvenet