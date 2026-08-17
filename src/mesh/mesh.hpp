// mesh.hpp
#pragma once

#include "mesh_types.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/StdVector>
#include <vector>
#include <map>


namespace Mesh {

// Template mesh class that augments our standard mesh setup with some quantities we need
/* 
File Descriptors:
    - mesh_types.hpp: Structs for primal objects used by the mesh class (Verts, Edges, Faces, HalfEdges, etc.)
    - mesh.cpp: Mesh instance initialization functions and getters
    - mesh_utils.cpp: Standard mesh operations (projection, editing, etc.)
    - mesh_iter.cpp: Standard mesh iterators and queries (vertex-edge mapping, vertex umbrellas, edge-face mapping, etc.)
*/
class mesh {
    public:
        // Constructor
        mesh(const std::vector<Eigen::Vector3d>& V_List, const std::vector<std::vector<int>>& F_List);
        mesh();

        // Project a vertex onto the mesh
        // mesh_utils.cpp
        // A version which returns a vertProjData object
        vertProjData computeVProjection(const Eigen::Vector3d& v, Eigen::Vector3d& proj, bool snap = true, bool fast = true) const;
        // A generalization of computeVProjection that also computes other data if the caller wants
        int computeVBinding(const Eigen::Vector3d& p, meshBindData& bind, bool snap = true, bool fast = true) const;

        // Recover a point location from stored coordinates
        int recoverCoords(int elType, int elIdx, const Eigen::VectorXd& coords, Eigen::Vector3d& p, bool fast = true);

        // Evaluate the basis weights of the verts spanning the mesh element a projData sits on
        int evaluateBasis(const vertProjData& proj, const Eigen::VectorXd& coords, std::vector<std::pair<int, double>>& basis) const;
        // Compute the mesh laplacian on verts
        int computeLaplacian(Eigen::SparseMatrix<double>& L);
        // Compute the lumped mass matrix (as a vector) on verts
        int computeMass(Eigen::VectorXd& A);
        // Get verts as a matrix
        int vertsAsMatrix(Eigen::MatrixXd& Verts);
        // Directly overwrite a vertex position (no BVH/normal recompute)
        void setVertPos(int v, const Eigen::Vector3d& pos);

        // Getters
        Eigen::Vector3d getNormal(vertProjData projData) const;
        Eigen::Vector3d getNormal(int elType, int elIdx) const;
        Eigen::Vector3d getVNormal(int v) const;
        Eigen::Vector3d getENormal(int e) const;
        Eigen::Vector3d getFNormal(int f) const;
        int getNumActiveV() const;

        // Get mean edge length
        double getMeanE() const;
        double getSquaredMeanE() const;
        // Get bounding box diagonal length
        double getBBoxDiag() const;

        friend class cutmesh;   // Let cutmesh read its internals :)

    protected:
        // INITIALIZATION HELPERS
        // Internal function to precompute normals and areas on mesh structures
        double computeFVectorArea(int f, Eigen::Vector3d& fN);  // 1 face
        void computeFNormalsAreas();        // All faces
        // Compute edge normals
        int computeENormal(int e, Eigen::Vector3d& eN, bool weight_fN = true);   // 1 edge
        void computeENormals(bool weight_fN = true);       // All edges
        // weight_fN weights by adjacent face areas
        void computeVNormalsAreas(bool weight_fN = true);       // All vertices
        // Clear all mesh attributes
        bool clearMesh();
        // Create a new vertex but do NOT insert it
        Vert createVertex(Eigen::Vector3d pos, Eigen::Vector3d n);
        HalfEdge createHalfEdge(bool boundary = false, 
                                    int twin = -1, 
                                    int dest = -1, 
                                    int edge = -1, 
                                    int face = -1, 
                                    int next = -1, 
                                    int prev = -1, 
                                    int dCN_idx = -1);
        Edge createEdge(int he = -1, Eigen::Vector3d n = Eigen::Vector3d::Zero());
        bool copyVertex(int v, Vert& new_vert);
        bool copyHalfEdge(int he, HalfEdge& new_he);
        bool copyEdge(int e, Edge& new_e);

        void countNumActive();


        // ------------- ITERATORS + QUERYING (mesh_iter.cpp) -----------------

        // Returns a CCW list of a vertex's OUTGOING halfedge indices
        std::vector<int> vertAdjHEs(int v) const;
        // Returns a CCW list of all incoming AND outgoing halfedge indices
        std::vector<int> vertAllHEs(int v) const;
        // Returns a CCW list of a vertex's adjacent vertices
        std::vector<int> vertAdjVerts(int v) const;
        // Returns a list of vertices in a loop from a given halfedge
        std::vector<int> vertLoop(int he) const;
        // Returns a CCW list of a vertex's adjacent faces
        std::vector<int> vertAdjFaces(int v) const;

        // Returns the endpoints of an edge in an arbitrary order.
        std::vector<int> edgeAdjVerts(int e) const;
        // Returns the adjacent face(s) of an edge (-1 indicates boundary)
        std::vector<int> edgeAdjFaces(int e) const;

        // Returns the halfedge index given the face index and edge index
        int halfedgeAtFaceEdge(int f, int e) const;
        // Get an arbitrary halfedge loop
        std::vector<int> halfedgeLoop(int he) const;

        // Returns a CCW list of a face's vertices
        std::vector<Eigen::Vector3d> faceAdjVerts(int f, bool origin = false) const;
        std::vector<Eigen::Vector3d> adjVerts(std::vector<int> vertIdxs) const;
        // Returns a CCW list of a face's vertex indices
        std::vector<int> faceAdjVertIdxs(int f, bool origin = false) const;
        // Returns a CCW list of a face's half edges
        std::vector<int> faceAdjHalfEdges(int f) const;

        // Adjacent faces given the element type and its index
        std::vector<int> adjFaces(int elType, int elIdx) const;

        // Returns the outgoing boundary HE if a vertex is a boundary vertex, else returns -1
        int vertIsBoundary(int v) const;

        // ------------- ATTRIBUTES -----------------
        // List of primal mesh elements
        std::vector<Vert> V;
        std::vector<HalfEdge> HE;
        std::vector<Edge> E;
        std::vector<Face> F;
        // Counters for convenience
        int active_v = 0;
        int active_e = 0;
        int active_f = 0;
        // Vertex-to-Edge Map for easy indexing
        std::map<std::pair<int, int>, int> vertPairToHE;

        // Mean edge length on mesh
        double meanE;
        // AABB Diagonal length
        double bboxDiag;

    private:
        // ------------- INITIALIZATION (mesh.cpp)  -----------------
        // Main init function
        bool initHalfEdgeMesh(const std::vector<Eigen::Vector3d>& V_List, const std::vector<std::vector<int>>& F_List);
        // Internal function to precompute height functions on both planar/nonplanar faces
        Eigen::VectorXd computeFaceHeight(int f) const;

        // Computes mean edge length on the mesh
        void computeMeanE();
        // Compute the length of the diagonal of the bounding box.
        void computeBBoxDiag();
        // Compute bounding volume hierarchy
        int computeBVH();
        int buildBVHNode(const std::vector<int>& faces, int depth);


        // ------------- UTILITIES (mesh_utils.cpp) -----------------
        // Projection helpers
        // Squared distance of point to AABB
        double pointAABBDist2(const Eigen::Vector3d& p, int box) const;
        // Closest point of a point to a face
        int closestPointOnFace(int f,
                            const Eigen::Vector3d& p,
                            Eigen::Vector3d& proj,
                            vertProjData& projData,
                            bool snap) const;
        // Find closest face given the BVH
        int closestFaceBVH(const Eigen::Vector3d& p,
                        Eigen::Vector3d& proj,
                        vertProjData& projData,
                        bool snap) const;
        int computeBindCoords(int elType, int elIdx, const Eigen::Vector3d& proj, Eigen::VectorXd& coords) const;
        int computeBindFrame(int elType, int elIdx, Eigen::Matrix3d& frame) const;

        // GEODESICS
        // Optional default input parameters for traced intersection vertices
        int traceGeodesic(const Vert& start,
                        const vertProjData& startData,
                        const Vert& end,
                        const vertProjData& endData,
                        Eigen::Vector3d prevDirec,
                        vertProjData prevData,
                        std::vector<Vert>& tracedVerts,
                        std::vector<vertProjData>& tracedProjData,
                        int depth = 0,
                        int max_depth = 5,
                        bool recompute = false,
                        bool fast = true);
        // Slow Termination check for traceGeodesic: Check if the end is visible from the start
        // on a shared face
        bool testVisibility(int f, Eigen::Vector3d start, Eigen::Vector3d end, double eps = 1e-4);

        int rayCastOnFace(int f, 
                        Eigen::Vector3d start, 
                        Eigen::Vector3d direc, 
                        Eigen::Vector3d& hit, 
                        vertProjData& hitData,
                        double eps = 1e-12);
        // Compute the next walk element given that we intersected with an edge
        int nextEl_Edge(int e, const vertProjData& originData, const Eigen::Vector3d& prev_direc, 
                    Eigen::Vector3d& next_direc, vertProjData& nextData, bool bdy_snap = true);
        // Figures out which next attribute to walk on given the very start is on an edge
        int nextEl_EdgeStart(int e, const vertProjData& originData, const Eigen::Vector3d& start_direc,
                           Eigen::Vector3d& next_direc, vertProjData& nextData, bool bdy_snap, double eps = 1e-6);
        // Helper for next edge start
        int snapWalkToEdge(int e, const Eigen::Vector3d& direc, Eigen::Vector3d& next_direc,
                         vertProjData& nextData, double eps = 1e-6) const;
        // Compute the next walk element given that we intersected with a vertex
        int nextEl_Vert(int v, const vertProjData& originData, 
                    const Eigen::Vector3d& prev_direc, Eigen::Vector3d& next_direc, 
                    vertProjData& nextData, bool bdy_snap = true, double eps = 1e-6);
        // Figures out which next attribute to walk on given the start is on a vert
        int nextEl_VertStart(int v, const vertProjData& originData, const Eigen::Vector3d& start_direc,
                           Eigen::Vector3d& next_direc, vertProjData& nextData, bool bdy_snap, double eps = 1e-6);

        // BVH
        std::vector<AABB> BVH;
        int max_depth = 5;
};

}   // namespace Mesh