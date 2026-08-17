// curvenet.hpp
#pragma once

#include "curvenet_types.hpp"
#include "mesh/mesh.hpp"
#include <Eigen/Core>
#include <vector>
#include <array>
#include <map>

namespace Polynet {
    class dcurvenet;
}

namespace Curvenet {

// Curvenet class consisting of cubic bezier splines
// NOTE: This class is NOT meant to be an editable curve network. It is purely static and has the purpose of supporting a profilemover object.
class curvenet {
    public:
        // Constructor takes four points [start, tangent 1, tangent 2, end], and associated normals
        // NOTE: Constructor assumes you already have no duplicates in your inputs
        curvenet(const std::vector<Eigen::Vector3d>& Controls, const std::vector<Eigen::Vector3d>& Tangents, const std::vector<std::array<int, 4>>& Splines, const Mesh::mesh* M = nullptr, int alpha = 5);
        // Empty constructor
        curvenet();

        // --------- UPDATE CURVENET -----------
        void updateCurveNet(std::vector<Eigen::Vector3d> Controls, std::vector<Eigen::Vector3d> Tangents);
        // Assign weight to a control
        int assignWeight(int c, bool fixed_w = true, double w = 1.0);
        void resetWeights();

        // --------- GETTERS -----------
        const int numCurves() const { return Crv.size(); }
        // TODO: REMOVE the below getters and use friend classes instead
        const std::vector<Control>& controls() const { return C; }
        std::vector<int> controlLocalSplineIdx(int c, int s) const;
        // Get a combined list of controls and tangents as a matrix
        int CTasMatrix(Eigen::MatrixXd& Verts);

        // --------- SAMPLING -----------
        // Sample a bezier curve at time t
        Eigen::Vector3d tSampleBezier(const Eigen::Vector3d& c0, const Eigen::Vector3d& c1, const Eigen::Vector3d& c2, const Eigen::Vector3d& c3, double t) const;
        Eigen::Vector3d tSampleBezier(int s, double t) const;
        // Get first derivative of bezier curve at a sample point
        Eigen::Vector3d tBezier_first(const Eigen::Vector3d& c0, const Eigen::Vector3d& c1, const Eigen::Vector3d& c2, const Eigen::Vector3d& c3, double t) const;
        Eigen::Vector3d tBezier_first(int s, double t) const;
        // Get second derivative of bezier curve at a sample point
        Eigen::Vector3d tBezier_second(const Eigen::Vector3d& c0, const Eigen::Vector3d& c1, const Eigen::Vector3d& c2, const Eigen::Vector3d& c3, double t) const;
        Eigen::Vector3d tBezier_second(int s, double t) const;
        // Evaluate basis functions on a spline given a t-value
        int evaluateBasis(int s, double t, std::vector<std::pair<int, double>>& basis);
        // NOTE: This is a naive, fast sampler that uniformly samples t's. Re-implement if desired
        // Returns n_samples points on the curve, including the endpoints
        std::vector<Eigen::Vector3d> sampleBezierNaive(int s, int n_samples = 50) const;
        // Estimate the arclength
        double arclenEst(int s, int n_samples = 50) const;
        double arclenEst(const std::vector<Eigen::Vector3d>& samples) const;
        // Uniformly sample based on arclength estimator
        // Takes a user parameter alpha which helps control sampling
        // Returns the length of the computed curve
        std::vector<Eigen::Vector3d> unifSample(int s, int n_samples = 50) const;
        // Same as unifSample, but returns the t-values instead of the sampled points
        std::vector<double> unifSampleT(int s, int n_samples = 50) const;
        // Compute number of samples to take on a spline given a user parameter alpha
        // Falls back to a fixed arclength-based heuristic if no mesh is given
        int computeNumSamples(double arclen, const Mesh::mesh* M);

        // --------- EDITING -----------
        // Exposed normal augmentation
        // NOTE: Normals should only applied to vertices once, by projection onto the mesh
        int editControlN(int c, Eigen::Vector3d normal);

        // --------- OTHER -----------
        int ctrlProjDataFromMesh(const Mesh::mesh& m);
        int tanProjDataFromMesh(const Mesh::mesh& m);
        int sortAdjHEAll();
        int assignCtrlTypeAll();

        // --------- CLOSEST POINT -----------
        // Find the closest point on the curve network to p, via dCN's polyline BVH, then
        // evaluated exactly on the spline. Not refined past that; caller does the Newton step.
        int closestPoint(const Eigen::Vector3d& p, const Polynet::dcurvenet* dCN, cnBindData& bind, bool snap = true, double snapTol = 1e-6) const;
        // Newton iterations to refine an initial guess t-value to get true closest point
        double optimizeT(double t, int s, const Eigen::Vector3d& p, int max_iter = 10) const;

        // Get bounding box diagonal length
        double getBBoxDiag() const;

        friend class Polynet::dcurvenet;
    protected:
        // No class inheritance
    private:
        // --------- INITIALIZATION -----------
        // Create new control
        int addControl(Eigen::Vector3d pos);

        // Add a spline given the start, end, and two tangent endpoints
        std::pair<int, int> addSpline(int start, int end, Eigen::Vector3d t0, Eigen::Vector3d t1, const Mesh::mesh* M);

        // --------- ITERATORS -----------
        // Get the adjacent tangent vectors to a control vertex
        std::vector<Eigen::Vector3d> ctrlAdjTans(int c) const;

        // Returns a list of the splines adjacent to a control 
        std::vector<int> ctrlAdjSplines(int c) const;

        // Returns vertices adjacent to a spline
        std::vector<int> splineAdjCtrls(int s) const;

        //  --------- OTHER -----------
        // Sort the halfedges of a control to be CCW
        // Should only be performed AFTER assigning normals to all verts
        int sortAdjHE(int c);
        // Compute what kind of vertex each control is using the valence of splines
        int assignCtrlType(int c);
        // Trace out curves
        int traceCurves();
        // Helper to find next traced spline
        int nextHEFromControl(int curr_he, int curr_end);

        // Get the 4 control points of a spline
        void splineCtrlPts(int s, Eigen::Vector3d& c0, Eigen::Vector3d& c1, Eigen::Vector3d& c2, Eigen::Vector3d& c3) const;
        // Compute the length of the diagonal of the bounding box
        void computeBBoxDiag();

        // Store attributes as lists
        std::vector<Control> C;
        std::vector<HalfEdge> HE;
        std::vector<CubicSpline> S;
        std::vector<Curve> Crv;

        // Mesh-projection data, parallel to C and HE respectively (empty if no mesh is used)
        std::vector<projData> vertData;
        std::vector<projData> tanData;

        // Map from input control index to output control index
        std::map<int, int> inputCtoC;
        // Map from input tangent index to output halfedge index
        std::map<int, int> inputTtoHE;

        // User sampling parameter
        int alpha = 5;
        // True once a mesh has been bound (at construction or later)
        bool setMesh = false;
        // Naive fallback sampling density (arc length per sample) when no mesh is available
        // TODO: Replace with a curvature-based heuristic
        double defaultSample = 0.1;

        // AABB Diagonal length
        double bboxDiag = 0.0;
};

}   // namespace Curvenet