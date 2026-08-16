// dcurvenet.hpp
#pragma once

#include "dcurvenet_types.hpp"
#include "curvenet/curvenet.hpp"
#include "mesh/mesh.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <map>
#include <glm/glm.hpp>
#include <glm/vec3.hpp>

namespace Mesh {
    class cutmesh;
}

namespace Curvenet {
    class curvenet;
}

namespace DCurvenet {

// Discrete curvenet (i.e., polylines)
class dcurvenet {
    public:
        // NOTE: Vertices are copied directly from the Control of the curve network (CN),
        //       meaning they are in the same order and have the same corresponding indices.
        //       Curves are similar.
        // Takes the original curvenet and discretizes it
        // Alpha is the user-inputted sampling parameter
        dcurvenet(Curvenet::curvenet* CN, Mesh::mesh* M);
        // Initialize with empty constructor
        dcurvenet();

        // --------- GETTERS -----------
        const int numVerts() const { return V.size(); }
        const int numHalfedges() const { return HE.size(); }
        const int numCurves() const { return C.size(); }
        const std::vector<Vert>& verts() const { return V; }
        const std::vector<HalfEdge>& halfedges() const { return HE; }
        const std::vector<Curve>& curves() const { return C; }

        // --------- POLYSCOPE VIZ -----------
        int polyscopeFormat(Eigen::MatrixXd& Verts, 
                            std::vector<std::array<int, 2>>& Edges, 
                            std::vector<glm::vec3>& posEdgeTangents,
                            std::vector<glm::vec3>& posEdgeBinormals,
                            std::vector<glm::vec3>& posEdgeNormals,
                            std::vector<glm::vec3>& negEdgeTangents,
                            std::vector<glm::vec3>& negEdgeBinormals,
                            std::vector<glm::vec3>& negEdgeNormals,
                            std::vector<double>& weights) const;

        
        // --------- RUNTIME COMPUTATION -----------
        // Update with new curvenet positions and local frames
        void updateDiscCurveNet();
        // Move a single vertex to a new position
        int moveVert(int v, Eigen::Vector3d new_pos); 
        // Compute the deformation gradient on an edge given a new scaled frame
        Eigen::Matrix3d computeHEDefGrad(int he);
        // Compute deformation gradients on all halfedges
        int computeDefGradAll();
        // Pre-compute maps
        int computedCNMats(Eigen::MatrixXd& f_dCN_flat, Eigen::MatrixXd& x_dCN);
        // Compute dCN positions matrix
        int computedCNVerts(Eigen::MatrixXd& x_dCN);

        // Propagate weights along curvenet
        int propagateWeights();

        friend class Mesh::cutmesh;    // Friend class to access curvenet variables
    protected:
        // No inherited classes
    private:
        // --------- INITIALIZATION -----------
        // Add a vertex
        int addVert(Curvenet::Control ctrl, int ctrl_idx = -1);
        int addVert(Eigen::Vector3d new_pos, Eigen::Vector3d new_n = Eigen::Vector3d::Zero(), int ctrl_idx = -1, int ctrl_type = -1, int adjSize = 0);
        // Add an edge and return the index of the new edge
        int addEdge(int origin, int dest, int prev_he0 = -1, int next_he1 = -1, int c = -1);
        // Add a curve
        int addCurve(int crv);
        // Rewire incoming/outgoing halfedges of an intersection vertex such that topology is correct
        int rewireVertAdjHE(int v);
        // Compute projection data
        int computeProjData();

        // --------- SCALED FRAME COMPUTATION -----------
        // For controls, computes their corner normals and widths. For non-intersections, this method does nothing (return -1)
        int vertCornerNormalsWidths(int v, std::vector<curveDeformData>& curveData);
        // Corner normals on all vertices
        int allCornerNormalsAndWidths(std::vector<curveDeformData>& curveData);
        // Transport corner normals and widths from the two end corners of a curve
        int transportNWOnCurve(int c, const curveDeformData& cData);
        // Transport normals and widths for all curves
        int transportNormalsAndWidths(const std::vector<curveDeformData>& curveData);

        // Compute local frame on a curve
        int computeScaledFrameOnCurve(int c);
        // Compute scaled frames on all curves
        int computeScaledFrames();
        // Initialize new frames as copies of old
        int copyFrameToRest(int he);
        int copyFramesToRest();
        // Validate that frames are not zero or NaN
        int validateFrames();
        
        // --------- SCALED FRAME HELPERS -----------
        
        // Accumulates rotation matrices and lengths by tracing from a starting halfedge to an end vertex
        double accumulateRotations(int start_he, int end_v, std::vector<Eigen::Matrix3d>& rots, std::vector<double>& lens);
        // Compute torsion
        double computeTorsion(Eigen::Vector3d n_1, Eigen::Vector3d n_k, Eigen::Matrix3d Om_k, Eigen::Vector3d t_k);
        
        // --------- OTHER -----------
        // Compute number of samples to take for a given
        int computeNumSamples(double arclen);

        // Check if a halfedge is the positive or negative side
        bool isPositiveHalfedge(int he) const;

        // Attributes as lists
        std::vector<Vert> V;
        std::vector<HalfEdge> HE;
        std::vector<Edge> E;
        std::vector<Curve> C;
        
        // Map input vertex index to local vertex index
        std::map<int, int> inputCtoV;
        // Map input curve index to local curve index
        std::map<int, int> inputCrvToC;

        // Pointer to parent curvenet
        Curvenet::curvenet* CN;
        // Pointer to parent mesh
        Mesh::mesh* M;
};

}   // namespace DCurvenet