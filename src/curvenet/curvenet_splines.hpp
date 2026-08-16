// curvenet_types.hpp
#pragma once

#include <Eigen/Core>
#include <vector>

// File with basic structs used by mesh class

namespace Curvenet {

    class Spline {  // Generic spline class
        public:
            Spline();

            void setHE(int new_he);
            void setCurve(int new_crv);
            bool setNumSamples(int n);

            // Sparse attributes
            int he = -1;        // One "halfedge" of the spline
            int curve = -1;     // The participating curve
            int num_samples = -1;       // Number of samples to take on the spline
            bool active = true;
        protected:
            // Compute num samples using a ratio b/w arclen. estimate and a parameter alpha
            virtual int compute_num_samples(double arclen, double meanE, int alpha = 5);
            // Sample at a parameter t
            virtual Eigen::Vector3d tSample(double t);
            // Naive sampling via uniform t-selection
            virtual std::vector<Eigen::Vector3d> sampleNaive(int num_samples = 50);
            // Estimate arclength via sampling
            virtual double arclenEst(int num_samples = 50);
            // Pseudo-uniformly sample the arclength discretization
            virtual std::vector<Eigen::Vector3d> unifSamples(int num_samples);
    };

    class CubicBezier : public Spline { // Cubic Bezier class
        public:
            CubicBezier();  // Constructor
        
        protected:
            Eigen::Vector3d tSample(double t) override;
            std::vector<Eigen::Vector3d> sampleNaive(int num_samples = 50) override;
            double arclenEst(int num_samples = 50) override;
            std::vector<Eigen::Vector3d> unifSamples(int num_samples) override;  
    };

    class CatmullRom : public Spline {  // Catmull Rom class
        public:
            CatmullRom();   // Constructor

        protected:
            Eigen::Vector3d tSample(double t) override;
            std::vector<Eigen::Vector3d> sampleNaive(int num_samples = 50) override;
            double arclenEst(int num_samples = 50) override;
            std::vector<Eigen::Vector3d> unifSamples(int num_samples) override;  
    };
}   // namespace Curvenet