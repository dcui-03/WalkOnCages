#include "mesh.hpp"

#include "../utils/decUtils.hpp"
#include "../utils/utils.hpp"
#include <Eigen/Core>
#include <vector>
#include <map>
#include <algorithm>
#include <utility>
#include <iostream>

// Mesh class functions for initialization

namespace Mesh {

// Constructor takes the projected curvenet and the mesh, and produces a cut-mesh
mesh::mesh(const std::vector<Eigen::Vector3d>& V_List, const std::vector<std::vector<int>>& F_List) {
    if (!initHalfEdgeMesh(V_List, F_List)) {
        throw std::runtime_error("Failed to initialize halfedge mesh.");
    }
    computeFNormalsAreas();
    computeENormals();
    computeVNormalsAreas();
    computeMeanE();
    computeBBoxDiag();
    if (computeBVH() != 1) {
        throw std::runtime_error("Failed to build BVH.");
    }
    
    return;
}

mesh::mesh() {

}

// Initializes the half edge mesh (Verts, Edges, Faces, Halfedges) from a vertex and face list
bool mesh::initHalfEdgeMesh(const std::vector<Eigen::Vector3d>& V_List, const std::vector<std::vector<int>>& F_List) {
    clearMesh();
    // Guard against empty meshes
    if (V_List.empty()) {
        std::cout << "EMPTY INPUT!" << std::endl;
        return false;
    }

    // 1. Initialize All Vertices
    V.resize(V_List.size());
    for (int v = 0; v < V_List.size(); v++) {
        V[v].pos = V_List[v];
    }
    active_v = V.size();

    std::vector<std::vector<int>> orientedF = F_List;

    if (!Utils::orientFacesConsistently(orientedF)) {
        std::cout << "Unable to orient faces." << std::endl;
        return false;
    }

    // 2. Initialize faces, edges, and interior halfedges
    F.resize(orientedF.size());
    // The number of face corners is a maximum on the number of edges needed
    int numFaceCorners = 0;
    for (const std::vector<int>& faceVerts : orientedF) {
        if (faceVerts.size() < 3) { // Not a valid face
            return false;
        }
        for (int vi : faceVerts) {
            if (vi < 0 || vi >= V_List.size()) {
                std::cout << "Invalid vertex index" << std::endl;
                return false;
            }
        }

        numFaceCorners += faceVerts.size();
    }
    HE.reserve(2 * numFaceCorners);
    E.reserve(numFaceCorners);

    for (int f = 0; f < orientedF.size(); f++) {
        const std::vector<int>& fVerts = orientedF[f];
        const int fSize = static_cast<int>(fVerts.size());

        // Temporary list of face HE's
        std::vector<int> faceHEs(fSize, -1);
        // Create interior halfedges for this face
        for (int i = 0; i < fSize; i++) {
            int vi = fVerts[i];
            int vj = fVerts[(i + 1) % fSize];
            // Keys for the vertex pair to HE map
            std::pair<int, int> dirKey = std::make_pair(vi, vj);
            std::pair<int, int> oppKey = std::make_pair(vj, vi);
            // Check that we don't already have this edge. If so, then there's a duplicate
            if (vertPairToHE.find(dirKey) != vertPairToHE.end()) {
                std::cout << "Duplicate edge found" << std::endl;
                return false;
            }
            // Insert the new halfedge into the list and set its attributes
            int heIdx = HE.size();
            HE.emplace_back();

            HE[heIdx].dest = vj;
            HE[heIdx].face = f;

            faceHEs[i] = heIdx;
            vertPairToHE[dirKey] = heIdx;

            // Store one outgoing halfedge for vi
            if (V[vi].he == -1) {
                V[vi].he = heIdx;
            }

            // Check whether the opposite direction already exists
            auto oppIt = vertPairToHE.find(oppKey);
            if (oppIt != vertPairToHE.end()) {  // If so, then hook the two halfedges up
                int oppHE = oppIt->second;
                // If the opposite halfedge already has a twin, then more than 2 face touch the edge (non-manifold)
                if (HE[oppHE].twin != -1) {
                    return false;
                }
                // If we run into an invalid or non-existent edge, then something is also wrong
                int eIdx = HE[oppHE].edge;
                if (eIdx < 0 || eIdx >= E.size()) {
                    return false;
                }
                HE[heIdx].twin = oppHE;
                HE[oppHE].twin = heIdx;
                HE[heIdx].edge = eIdx;
            } else {    // If not, add in this new edge.
                int eIdx = E.size();
                E.emplace_back();

                E[eIdx].he = heIdx;
                HE[heIdx].edge = eIdx;
            }
        }

        // Wire next/prev HE's around the face.
        for (int i = 0; i < fSize; ++i) {
            int he = faceHEs[i];

            HE[he].next = faceHEs[(i + 1) % fSize];
            HE[he].prev = faceHEs[(i + fSize - 1) % fSize];
        }

        F[f].he = faceHEs[0];
    }

    // 3. Create boundary halfedges
    std::vector<int> boundaryHEs;
    const int numInteriorHEs = HE.size();
    // Iterate over interior halfedges and find any that have no twin (i.e., twin = -1)
    for (int he = 0; he < numInteriorHEs; ++he) {
        if (HE[he].twin != -1) {
            continue;
        }
        // Get the starting and ending vertices
        int u = HE[HE[he].prev].dest;
        int v = HE[he].dest;
        // Get an index and insert the new boundary half edge in + attributes
        int bhe = HE.size();
        HE.emplace_back();

        HE[bhe].dest = u;
        HE[bhe].twin = he;
        HE[bhe].edge = HE[he].edge;
        HE[bhe].face = -1;
        HE[bhe].boundary = true;

        HE[he].twin = bhe;

        // Prefer boundary outgoing halfedge for boundary vertices (for easy querying)
        V[v].he = bhe;

        std::pair<int, int> bKey = std::make_pair(v, u);
        // If the boundary halfedge already exists somehow, then something is wrong
        if (vertPairToHE.find(bKey) != vertPairToHE.end()) {
            std::cout << "boundary halfedge already exists" << std::endl;
            return false;
        }
        vertPairToHE[bKey] = bhe;
        boundaryHEs.push_back(bhe);
    }

    // 4. Connect next/prev for boundary halfedges
    std::map<int, int> boundaryOutgoingFromVertex;  // Store a map with the outgoing HE from each bdy vertex
    // Iterate over boundary halfedges
    for (int bhe : boundaryHEs) {
        // Boundary halfedge origin is the dest of its twin
        int origin = HE[HE[bhe].twin].dest;
        // If a vertex has more than one outgoing boundary halfedge, then it must be nonmanifold
        if (boundaryOutgoingFromVertex.find(origin) != boundaryOutgoingFromVertex.end()) {
            std::cout << "Nonmanifold edge found" << std::endl;
            return false;
        }
        boundaryOutgoingFromVertex[origin] = bhe;
    }
    // Iterate over boundary vertices and find the associated next/prev
    for (int bhe : boundaryHEs) {
        int dest = HE[bhe].dest;
        // make sure that there exists an associated next vertex (i.e., valid boundary configuration)
        auto nextIt = boundaryOutgoingFromVertex.find(dest);
        if (nextIt == boundaryOutgoingFromVertex.end()) {
            return false;
        }
        int nextBHE = nextIt->second;

        HE[bhe].next = nextBHE;
        HE[nextBHE].prev = bhe;
    }

    countNumActive();
    std::cout << "Num. active vertices: " << active_v << std::endl;
    std::cout << "Num. active edges: " << active_e << std::endl;
    std::cout << "Num. active faces: " << active_f << std::endl;
    return true;    // success!
}

// Clear all mesh attributes
bool mesh::clearMesh() {
    // Clear all lists
    V.clear();
    HE.clear();
    E.clear();
    F.clear();
    vertPairToHE.clear();
    // Reset number of vertices
    active_v = 0;
    active_e = 0;
    active_f = 0;
    // Reset mesh qualities
    meanE = 0.0;
    bboxDiag = 0.0;
    return true;
}

// Getters
Eigen::Vector3d mesh::getNormal(Utils::projData projData) const {
    if (projData.elType == 0) {
        return getVNormal(projData.elIdx);
    } else if (projData.elType == 1) {
        return getENormal(projData.elIdx);
    } else if (projData.elType == 2) {
        return getFNormal(projData.elIdx);
    } else {
        return Eigen::Vector3d::Zero();
    }
}
Eigen::Vector3d mesh::getNormal(int elType, int elIdx) const {
    if (elType == 0) {
        return getVNormal(elIdx);
    } else if (elType == 1) {
        return getENormal(elIdx);
    } else if (elType == 2) {
        return getFNormal(elIdx);
    } else {
        return Eigen::Vector3d::Zero();
    }
}
Eigen::Vector3d mesh::getVNormal(int v) const {
    if (v < 0 || v >= V.size()) {
        return Eigen::Vector3d::Zero();
    }
    return V[v].n;
}
Eigen::Vector3d mesh::getENormal(int e) const {
    if (e < 0 || e >= E.size()) {
        return Eigen::Vector3d::Zero();
    }
    return E[e].n;
}
Eigen::Vector3d mesh::getFNormal(int f) const {
    if (f < 0 || f >= F.size()) {
        return Eigen::Vector3d::Zero();
    }
    Eigen::Vector3d n = F[f].n;
    return n;
}

// Get mean edge length
double mesh::getMeanE() const {
    return meanE;
}

int mesh::getNumActiveV() const {
    return active_v;
}

Eigen::VectorXd mesh::computeFaceHeight(int f) const {
    const std::vector<Eigen::Vector3d> fVertsPos = faceAdjVerts(f);
    int fSize = fVertsPos.size();
    Eigen::VectorXd faceH = Eigen::VectorXd::Zero(fSize);
    // Special handling for triangles (must be planar)
    if (fSize == 3) {
        faceH.setZero();
        return faceH;
    }
    // 1. compute barycenter and face normal
    Eigen::Vector3d faceCenter = DECUtils::computeBarycenter(fVertsPos);
    Eigen::Vector3d faceN = F[f].n;
    // 2. Project face vertices onto the Newell plane and grab height
    std::vector<Eigen::Vector2d> proj_v(fSize);
    for (int v = 0; v < fSize; v++) {
        Eigen::Vector3d proj3d = Utils::projectPointOntoPlane(faceN, faceCenter, fVertsPos[v]);
        // vector from old point to plane point
        Eigen::Vector3d heightVec = fVertsPos[v] - proj3d;
        // Get 2D version
        double height = heightVec.norm();   // How far we are from the plane
        if (height <= 1e-6) {   // If we are on/close to the surface, just snap to the plane
            faceH(v) = 0.0;
        } else if ((heightVec.normalized()).dot(faceN) > 0.0) { // We are above the plane
            faceH(v) = height;
        } else {    // We are below the plane
            faceH(v) = -1.0 * height;
        }
    }
    return faceH;
}

// Function which computes a single face's normal/area
double mesh::computeFVectorArea(int f, Eigen::Vector3d& fN) {
    std::vector<Eigen::Vector3d> fVertsPos = faceAdjVerts(f);
    double area = DECUtils::vectorArea(fVertsPos, fN); // TODO: Be careful about degenerate normals!
    return area;
}

// Internal function to precompute normals on all mesh structures
void mesh::computeFNormalsAreas() {
    for (int f = 0; f < F.size(); f++) {
        if(!F[f].active) {
            continue;
        }
        Eigen::Vector3d fN = Eigen::Vector3d::Zero();
        double fArea = computeFVectorArea(f, fN);
        if (fArea <= 1e-8) {    // Throw an error if we have a degenerate face
            std::cout << "mesh::computeFNormalsAreas(): One face is degnerate." << std::endl;
            //throw std::runtime_error("mesh::computeFNormalsAreas(): One face is degnerate.");
            //return;
        }
        F[f].fArea = fArea;
        F[f].n = fN;
    }
    return;
}

// Compute edge normals
int mesh::computeENormal(int e, Eigen::Vector3d& eN, bool weight_fN) {
    std::vector<int> eFaces = edgeAdjFaces(e);
    eN = Eigen::Vector3d::Zero();
    if (eFaces[0] != -1) {
        eN += F[eFaces[0]].n;
    }
    if (eFaces[1] != -1) {
        eN += F[eFaces[1]].n;
    }
    eN.normalize();  // TODO: needs safe normalization
    return 1;
}
void mesh::computeENormals(bool weight_fN) {
    for (int e = 0; e < E.size(); e++) {
        if (!E[e].active) {
            continue;
        }
        Eigen::Vector3d eN = Eigen::Vector3d::Zero();
        computeENormal(e, eN, weight_fN);
        E[e].n = eN;
    }
    return;
}

// Function which computes a single vertex's normal/area
// Function which computes ALL vertex normals and areas
void mesh::computeVNormalsAreas(bool weight_fN) {
    // Reset vertex normals and areas
    for (int v = 0; v < V.size(); v++) {
        V[v].n = Eigen::Vector3d::Zero();
        V[v].vArea = 0.0;
    }

    // Accumulate face normals onto their adjacent vertices
    for (int f = 0; f < F.size(); f++) {
        if (!F[f].active) {
            continue;
        }

        std::vector<int> fVerts = faceAdjVertIdxs(f);
        const int n = fVerts.size();
        if (n < 3) {
            continue;
        }

        const double baryArea = F[f].fArea / static_cast<double>(n);

        for (int vi : fVerts) {
            if (vi < 0 || vi >= V.size() || !V[vi].active) {
                continue;
            }

            if (weight_fN) {
                V[vi].n += baryArea * F[f].n;
            } else {
                V[vi].n += F[f].n;
            }

            V[vi].vArea += baryArea;
        }
    }

    // Normalize accumulated normals
    for (int v = 0; v < V.size(); v++) {
        if (!V[v].active) {
            continue;
        }
        const double len = V[v].n.norm();
        V[v].n /= len;
    }
    return;
}

// Computes mean edge length on the mesh
void mesh::computeMeanE() {
    meanE = 0.0;
    if (active_e == 0) {
        return;
    }
    for (int e = 0; e < E.size(); ++e) {
        if (!E[e].active) {
            continue;
        }
        std::vector<int> eVerts = edgeAdjVerts(e);
        meanE += (V[eVerts[0]].pos - V[eVerts[1]].pos).norm();
    }
    meanE /= active_e;
    return;
}

double mesh::getSquaredMeanE() const {
    double sqMeanE = 0.0;
    if (active_e == 0) {
        return 0.0;
    }
    for (int e = 0; e < E.size(); ++e) {
        if (!E[e].active) {
            continue;
        }
        std::vector<int> eVerts = edgeAdjVerts(e);
        sqMeanE += (V[eVerts[0]].pos - V[eVerts[1]].pos).squaredNorm();
    }
    sqMeanE /= active_e;
    return sqMeanE;
}

// Get bounding box diagonal length
double mesh::getBBoxDiag() const {
    return bboxDiag;
}

// Computes the diagonal length of the mesh's AABB
void mesh::computeBBoxDiag() {
    bool found = false;     // Safety
    Eigen::Vector3d minV, maxV;

    for (int v = 0; v < V.size(); v++) {
        if (!V[v].active) {
            continue;
        }
        if (!found) {
            minV = V[v].pos;
            maxV = V[v].pos;
            found = true;
        } else {
            minV = minV.cwiseMin(V[v].pos);
            maxV = maxV.cwiseMax(V[v].pos);
        }
    }
    bboxDiag = found ? (maxV - minV).norm() : 0.0;
    return;
}

// Compute bounding volume hierarchy
int mesh::computeBVH() {
    BVH.clear();

    // Get active face list
    std::vector<int> activeFaces;
    activeFaces.reserve(F.size());
    for (int f = 0; f < F.size(); f++) {
        if (F[f].active) {
            activeFaces.push_back(f);
        }
    }
    if (activeFaces.empty()) {
        return -1;
    }

    int root = buildBVHNode(activeFaces, 0);
    if (root < 0) {
        return -1;
    }
    // Get the true max depth
    max_depth = 0;
    for (int i = 0; i < BVH.size(); i++) {
        max_depth = std::max(BVH[i].depth, max_depth);
    }
    return 1;
}

int mesh::buildBVHNode(const std::vector<int>& faces, int depth) {
    if (faces.empty()) {
        return -1;
    }
    int max_depth = 16;     // Hard stop so we don't recurse too much
    int leafFaceCount = 50; // Maximum number of faces a leaf node can have
    int nodeIdx = static_cast<int>(BVH.size());
    BVH.emplace_back();

    BVH[nodeIdx].depth = depth;
    BVH[nodeIdx].leaf = false;
    BVH[nodeIdx].faces.clear();
    BVH[nodeIdx].children.clear();

    // 1. Compute tight AABB around all faces in this node
    bool found = false;
    Eigen::Vector3d minV;
    Eigen::Vector3d maxV;
    for (int f_idx = 0; f_idx < faces.size(); f_idx++) {
        int f = faces[f_idx];
        if (f < 0 || f >= F.size() || !F[f].active) {
            continue;
        }
        std::vector<Eigen::Vector3d> fVerts = faceAdjVerts(f);
        for (int v_idx = 0; v_idx < fVerts.size(); v_idx++) {
            const Eigen::Vector3d& p = fVerts[v_idx];
            if (!found) {
                minV = p;
                maxV = p;
                found = true;
            } else {
                minV = minV.cwiseMin(p);
                maxV = maxV.cwiseMax(p);
            }
        }
    }
    if (!found) {
        return -1;
    }
    BVH[nodeIdx].bdyVerts = {minV, maxV};

    // 2. Leaf stopping criteria
    if (depth >= max_depth || faces.size() <= leafFaceCount) {
        BVH[nodeIdx].leaf = true;
        BVH[nodeIdx].faces = faces;
        return nodeIdx;
    }

    // 3. Choose split axis using longest AABB extent
    Eigen::Vector3d extent = maxV - minV;
    int axis = 0;
    extent.maxCoeff(&axis);
    // Degenerate box: cannot split meaningfully
    if (extent[axis] <= 1e-12) {
        BVH[nodeIdx].leaf = true;
        BVH[nodeIdx].faces = faces;
        return nodeIdx;
    }

    // 4. Compute face centroids along selected axis
    // Note since we only split along one axis, we only need to take the average along that axis
    std::vector<std::pair<double, int>> centroidFacePairs;
    centroidFacePairs.reserve(faces.size());
    for (int f_idx = 0; f_idx < faces.size(); f_idx++) {
        int f = faces[f_idx];
        if (f < 0 || f >= F.size() || !F[f].active) {
            continue;
        }
        std::vector<Eigen::Vector3d> fVerts = faceAdjVerts(f);
        if (fVerts.empty()) {
            continue;
        }
        Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
        for (int v_idx = 0; v_idx < fVerts.size(); v_idx++) {
            centroid += fVerts[v_idx];
        }
        centroid /= static_cast<double>(fVerts.size());
        centroidFacePairs.push_back(std::make_pair(centroid[axis], f));
    }
    // We've hit the leaf face count. Terminate at this node
    if (centroidFacePairs.size() <= leafFaceCount) {
        BVH[nodeIdx].leaf = true;
        for (int i = 0; i < centroidFacePairs.size(); i++) {
            BVH[nodeIdx].faces.push_back(centroidFacePairs[i].second);
        }
        return nodeIdx;
    }

    // 5. Sort faces by centroid coordinate
    std::sort(centroidFacePairs.begin(), centroidFacePairs.end());
    // 6. Split face list in half
    int mid = centroidFacePairs.size() / 2;
    // Safety: Degenerate case
    if (mid <= 0 || mid >= centroidFacePairs.size()) {
        BVH[nodeIdx].leaf = true;
        for (int i = 0; i < centroidFacePairs.size(); i++) {
            BVH[nodeIdx].faces.push_back(centroidFacePairs[i].second);
        }
        return nodeIdx;
    }

    // Distinguish the left and right side faces
    std::vector<int> leftFaces;
    std::vector<int> rightFaces;
    leftFaces.reserve(mid);
    rightFaces.reserve(centroidFacePairs.size() - mid);
    for (int i = 0; i < centroidFacePairs.size(); i++) {
        if (i < mid) {
            leftFaces.push_back(centroidFacePairs[i].second);
        } else {
            rightFaces.push_back(centroidFacePairs[i].second);
        }
    }

    // 7. Recurse to get children
    int leftChild = buildBVHNode(leftFaces, depth + 1);
    int rightChild = buildBVHNode(rightFaces, depth + 1);
    if (leftChild < 0 || rightChild < 0) {  // If either child encounters an error, then return as leaf
        BVH[nodeIdx].leaf = true;
        BVH[nodeIdx].faces = faces;
        BVH[nodeIdx].children.clear();
        return nodeIdx;
    }
    BVH[nodeIdx].children.push_back(leftChild);
    BVH[nodeIdx].children.push_back(rightChild);

    // Clear any non-leaf nodes' faces for storage (we will never need them)
    BVH[nodeIdx].faces.clear();
    BVH[nodeIdx].leaf = false;
    return nodeIdx;
}


// Create a new vertex but do NOT insert it
Vert mesh::createVertex(Eigen::Vector3d pos, Eigen::Vector3d n) {
    Vert v;
    v.pos = pos;
    v.n = n;
    return v;
}

HalfEdge mesh::createHalfEdge(bool boundary, int twin, int dest, int edge, int face, int next, int prev, int dCN_idx) {
    HalfEdge he;
    he.boundary = boundary;
    he.twin = twin;
    he.dest = dest;
    he.edge = edge;
    he.face = face;
    he.next = next;
    he.prev = prev;
    he.dCN_idx = dCN_idx;
    return he;
}

Edge mesh::createEdge(int he, Eigen::Vector3d n) {
    Edge e;
    e.he = he;
    e.n = n;
    return e;
}
bool mesh::copyVertex(int v, Vert& new_vert) {
    if (v < 0 || v > V.size()) {
        return false;
    }
    const Vert& vert = V[v];
    new_vert = createVertex(vert.pos, vert.n);
    return true;
}
bool mesh::copyHalfEdge(int he, HalfEdge& new_he) {
    if (he < 0 || he > HE.size()) {
        return false;
    }
    const HalfEdge& halfedge = HE[he];
    new_he = createHalfEdge(halfedge.boundary, halfedge.twin, halfedge.dest, halfedge.edge, halfedge.face, halfedge.next, halfedge.prev, halfedge.dCN_idx);
    return true;
}
bool mesh::copyEdge(int e, Edge& new_e) {
    if (e < 0 || e > E.size()) {
        return false;
    }
    const Edge& edge = E[e];
    new_e = createEdge(edge.he, edge.n);
    return true;
}

// Count number of active elements
void mesh::countNumActive() {
    active_v = 0;
    active_e = 0;
    active_f = 0;
    for (int f = 0; f < F.size(); f++) {
        if (F[f].active) {
            active_f++;
        }
    }

    for (int e = 0; e < E.size(); e++) {
        if (E[e].active) {
            active_e++;
        }
    }
    
    for (int v = 0; v < V.size(); v++) {
        if (V[v].active) {
            active_v++;
        }
    }
    return;
}

}   // namespace Mesh