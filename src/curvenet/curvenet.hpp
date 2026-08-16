// curvenet.hpp
#pragma once

#include "curvenet_types.hpp"
#include "curvenet_types.hpp"
#include "mesh/mesh.hpp"
#include <Eigen/Core>
#include <vector>
#include <array>
#include <map>

namespace DCurvenet {
    class dcurvenet;
}

namespace Curvenet {

// Curvenet class consisting of cubic bezier splines
// NOTE: This class is NOT meant to be an editable curve network. It is purely static and has the purpose of supporting a profilemover object.
class curvenet {
    public:
        // Constructor takes four points [start, tangent 1, tangent 2, end], and associated normals
        // NOTE: Constructor assumes you already have no duplicates in your inputs
        curvenet(const std::vector<Eigen::Vector3d>& Controls, const std::vector<Eigen::Vector3d>& Tangents, const std::vector<std::array<int, 4>>& Splines, const Mesh::mesh& M, int alpha = 5);
        // Empty constructor
        curvenet();

        // --------- UPDATE CURVENET -----------
        void updateCurveNet(std::vector<Eigen::Vector3d> Controls, std::vector<Eigen::Vector3d> Tangents);
        // Assign weight to a control
        int assignWeight(int c, bool fixed_w = true, double w = 1.0);
        void resetWeights();
        // Check if all weights are free or not
        // If all weights are free, then return them all to 1
        int validWeights();

        // --------- GETTERS -----------
        const int numControls() const { return C.size(); }
        const int numSplines() const { return S.size(); }
        const int numCurves() const { return Crv.size(); }
        // TODO: REMOVE the below getters and use friend classes instead
        const std::vector<Control>& controls() const { return C; }
        const std::vector<HalfEdge>& halfedges() const { return HE; }
        const std::vector<CubicSpline>& splines() const { return S; }
        const std::vector<Curve>& curves() const { return Crv; }
        std::vector<int> controlLocalSplineIdx(int c, int s) const;

        // --------- SAMPLING -----------
        // Sample a bezier curve at time t
        Eigen::Vector3d tSampleBezier(const Eigen::Vector3d& c0, const Eigen::Vector3d& c1, const Eigen::Vector3d& c2, const Eigen::Vector3d& c3, double t) const;
        Eigen::Vector3d tSampleBezier(int s, double t) const;
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
        // Compute number of samples to take on a spline given a user parameter alpha
        int computeNumSamples(double arclen);

        // --------- EDITING -----------
        // Exposed position edits
        int editControlPos(int c, Eigen::Vector3d pos);
        // Exposed normal augmentation
        // NOTE: Normals should only applied to vertices once, by projection onto the mesh
        int editControlN(int c, Eigen::Vector3d normal);

        // --------- OTHER -----------
        int ctrlProjDataFromMesh(const Mesh::mesh& m);
        int tanProjDataFromMesh(const Mesh::mesh& m);
        int sortAdjHEAll();
        int assignCtrlTypeAll();

        friend class DCurvenet::dcurvenet;
    protected:
        // No class inheritance
    private:
        // --------- INITIALIZATION -----------
        // Create new control
        int addControl(Eigen::Vector3d pos);

        // Add a spline given the start, end, and two tangent endpoints
        std::pair<int, int> addSpline(int start, int end, Eigen::Vector3d t0, Eigen::Vector3d t1);

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

        // Store attributes as lists
        std::vector<Control> C;
        std::vector<HalfEdge> HE;
        std::vector<CubicSpline> S;
        std::vector<Curve> Crv;

        // Map from input control index to output control index
        std::map<int, int> inputCtoC;
        // Map from input tangent index to output halfedge index
        std::map<int, int> inputTtoHE;

        // User sampling parameter
        int alpha = 5;
        double meanE = 0.0;
};

}   // namespace Curvenet