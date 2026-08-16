// pscurvenet.hpp
#pragma once

#include "pscurvenet_types.hpp"
#include <Eigen/Core>
#include <vector>
#include <array>
#include <map>
#include <string>

namespace psCurvenet {

// Editable curvenet class with various operations for test editing
// NOTE: This curvenet representation is NOT associated with profilemover class and is for Polyscope testing purposes only
class pscurvenet {
    public:
        // Empty constructor
        pscurvenet();

        // --------- UPDATE CURVENET -----------
        // Hard reset everything
        void resetCurvenet();
        // Update control position
        void updateControlPos(int c, Eigen::Vector3d new_pos, bool project = true);
        void updateControlNormal(int c, const Eigen::Vector3d& new_normal, bool rotation = true, bool project = true);
        // Rotate tangent using some rotation matrix
        bool rotateTangentPos(int psT_idx, Eigen::Matrix3d rotation);
        bool rotateTangentPos(int s, bool t0, Eigen::Matrix3d rotation);
        // Update tangent position, or returns false
        bool updateTangentPos(int psT_idx, const Eigen::Vector3d& new_pos, bool project = true);
        bool updateTangentPos(int s, bool t0, const Eigen::Vector3d& new_pos, bool project = true);
        // Add a control and return its index
        int addControl(Eigen::Vector3d pos, Eigen::Vector3d normal = Eigen::Vector3d::UnitZ());
        // Add a spline given only the start and end. Estimate t0 and t1 from these
        int addSpline(int c0, int c1, double init_factor = 2.5);
        // Add a spline and return its index
        int addSpline(int c0, Eigen::Vector3d t0_pos, Eigen::Vector3d t1_pos, int c1);
        // Remove a control
        int removeControl(int c);
        // Remove a spline given two end controls
        // NOTE: If two controls are connected by multiple splines, only removes the first one
        int removeSpline(int c0, int c1);
        // Remove a spline
        int removeSpline(int s);
        // Remove spline by polyscope tangent index
        int removeSplineByTangent(int ps_TIdx);
        // Clean up all loose controls
        int cleanupControls();

        // --------- IO -----------
        int saveCurvenet(const std::string& filepath) const;
        int loadCurvenet(const std::string& filepath);

        // --------- COMPLEX OPERATIONS -----------
        // Merge two splines about a specified control point
        int mergeSplines(int s0, int s1, int c);

        // --------- GETTERS + POLYSCOPE CONVERSION -----------
        Eigen::Vector3d getNormal(int c);
        // Control positions but returned as an Eigen::MatrixXd
        void cPosAsMatrix(Eigen::MatrixXd& cPos);
        // Tangent positions but returned as an Eigen::MatrixXd
        void tPosAsMatrix(Eigen::MatrixXd& tPos);
        // Curve network as a chain of discrete splines
        void cnAsCurveNetwork(Eigen::MatrixXd& verts, std::vector<std::array<int, 2>>& connectivity);
        void tansAsCurveNetwork(Eigen::MatrixXd& verts, std::vector<std::array<int, 2>>& connectivity);

        // Convert to vector of control + tangent positions, and spline indices
        void cnAsStdVector(std::vector<Eigen::Vector3d>& controls, std::vector<Eigen::Vector3d>& tangents, std::vector<std::array<int, 4>>& splines);

        // --------- SAMPLING -----------
        // Sample a bezier curve at time t
        Eigen::Vector3d tSampleBezier(const Eigen::Vector3d& c0, const Eigen::Vector3d& c1, const Eigen::Vector3d& c2, const Eigen::Vector3d& c3, double t);
        Eigen::Vector3d tSampleBezier(int s, double t);
        // NOTE: This is a naive, fast sampler that uniformly samples t's. Re-implement if desired
        // Returns n_samples points on the curve, including the endpoints
        std::vector<Eigen::Vector3d> sampleBezierNaive(int s, int n_samples = 30);
    protected:
        // No class inheritance
    private:
        // --------- ITERATORS -----------
        // Get the adjacent splines to a control vertex
        std::vector<int> ctrlAdjSplines(int c) const;

        // Returns vertices adjacent to a spline
        std::pair<int, int> splineAdjCtrls(int s) const;

        // Reconstruct the map from tangents
        void constructTangentMap();

        // Store attributes as maps, where the key is the index
        std::vector<Control> C;
        std::vector<Spline> S;
        
        bool recomputeMap = false;
        std::map<int, std::pair<int, bool>> psTangentToS;    // Maps polyscope tangent to spline and t0/t1
};

}   // namespace psCurvenet