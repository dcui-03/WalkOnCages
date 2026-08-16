#include "curvenet.hpp"

#include "utils/utils.hpp"
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
    Eigen::Vector3d curvenet::tSampleBezier(int s, double t) const {
        Eigen::Vector3d c0, c1, c2, c3;
        int he = S[s].he;
        c0 = C[HE[he].origin].new_pos;
        c1 = HE[he].tan;
        c2 = HE[HE[he].twin].tan;
        c3 = C[HE[HE[he].twin].origin].new_pos;
        return tSampleBezier(c0, c1, c2, c3, t);
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
    int curvenet::computeNumSamples(double arclen) {
        if (meanE <= 1e-16) return 3;
        return std::max(3, static_cast<int>(std::ceil(alpha * arclen / meanE)));
    }

    // Uniformly sample based on arclength estimator and returns the length of the returned curve
    std::vector<Eigen::Vector3d> curvenet::unifSample(int s, int n_samples) const {
        if (n_samples < 2) {
            return {};
        }

        std::vector<Eigen::Vector3d> regSamples = sampleBezierNaive(s, 50);
        std::vector<Eigen::Vector3d> unifSamples(n_samples);

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
                unifSamples[i] = regSamples.front();
            }
            return unifSamples;
        }
        unifSamples[0] = regSamples.front();
        unifSamples[n_samples - 1] = regSamples.back();

        int seg = 1;
        for (int i = 1; i < n_samples - 1; ++i) {
            double targetLen = totalLen * static_cast<double>(i) / static_cast<double>(n_samples - 1);
            // Check if we overflow this segment
            while ((seg < cumLen.size() - 1) && (cumLen[seg] < targetLen)) {
                ++seg;
            }

            double segLen = cumLen[seg] - cumLen[seg - 1];
            if (segLen <= 1e-16) {
                unifSamples[i] = regSamples[seg];
                continue;
            }
            // Increment by computing the percentage difference
            double alpha = (targetLen - cumLen[seg - 1]) / segLen;
            alpha = std::max(0.0, std::min(1.0, alpha));
            unifSamples[i] = (1.0 - alpha) * regSamples[seg - 1] + alpha * regSamples[seg];
        }

        return unifSamples;
    }

}   // namespace Curvenet