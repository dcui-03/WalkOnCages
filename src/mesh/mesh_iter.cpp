// File for mesh iterators and querying
#include "mesh.hpp"

#include <Eigen/Core>
#include <vector>


namespace Mesh {

// Returns a CCW list of a vertex's OUTGOING halfedge indices
std::vector<int> mesh::vertAdjHEs(int v) const {
    if (v < 0 || v >= V.size() || V[v].he < 0) {
        return {};
    }
    std::vector<int> outgoingHEs;
    const int he0 = V[v].he;
    int he_curr = he0;
    do {
        outgoingHEs.push_back(he_curr);
        he_curr = HE[HE[he_curr].prev].twin;
    } while (he_curr != he0);
    return outgoingHEs;
}

// Returns a CCW list of ALL of a vertex's incoming and outgoin halfedge indices
std::vector<int> mesh::vertAllHEs(int v) const {
    if (v < 0 || v >= V.size() || V[v].he < 0) {
        return {};
    }
    std::vector<int> adjHEs;
    const int he0 = V[v].he;
    int he_curr = he0;
    do {
        adjHEs.push_back(he_curr);
        adjHEs.push_back(HE[he_curr].prev);
        he_curr = HE[HE[he_curr].prev].twin;
    } while (he_curr != he0);
    return adjHEs;
}

// Returns a CCW list of a vertex's adjacent vertices
std::vector<int> mesh::vertAdjVerts(int v) const {
    if (v < 0 || v >= V.size() || V[v].he < 0) {
        return {};
    }
    std::vector<int> adjHE = vertAdjHEs(v);
    std::vector<int> adjVerts(adjHE.size());
    for (int he = 0; he < adjHE.size(); he++) {
        adjVerts[he] = HE[adjHE[he]].dest;
    }
    return adjVerts;
}

// Get an arbitrary vertex loop starting from a certain halfedge index
std::vector<int> mesh::vertLoop(int he) const {
    int max_search = active_e/2;
    int iter = 0;
    std::vector<int> vLoop;
    int he_curr = he;
    do {
        vLoop.push_back(HE[he_curr].dest);
        he_curr = HE[he_curr].next;
        iter++;
    } while((he_curr != he) && (iter <= max_search));
    return vLoop;
}

// Returns a CCW list of a vertex's adjacent faces
// INCLUDES BOUNDARY if a boundary is adjacent
std::vector<int> mesh::vertAdjFaces(int v) const {
    if (v < 0 || v >= V.size() || V[v].he < 0) {
        return {};
    }
    std::vector<int> adjHE = vertAdjHEs(v);
    std::vector<int> adjFaces;
    for (int he = 0; he < adjHE.size(); he++) {
        int f = HE[adjHE[he]].face;
        bool faceExists = false;
        // Check that we have not already recorded this face
        for (int fi = 0; fi < adjFaces.size(); fi++) {
            if (f == adjFaces[fi]) {
                faceExists = true;
                break;
            }
        }
        if (faceExists) {   // If so, skip this case.
            continue;
        }
        adjFaces.push_back(f);
    }
    return adjFaces;
}

// Returns the endpoints of an edge in an arbitrary order.
std::vector<int> mesh::edgeAdjVerts(int e) const {
    int he = E[e].he;
    std::vector<int> v_pair = {HE[HE[he].prev].dest, HE[he].dest};
    return v_pair;
}

// Returns the adjacent face(s) of an edge (-1 indicates boundary)
std::vector<int> mesh::edgeAdjFaces(int e) const {
    int he0 = E[e].he;
    int he1 = HE[he0].twin;
    std::vector<int> v_pair = {HE[he0].face, HE[he1].face};
    return v_pair;
}

// Returns the halfedge index given the face index and edge index
// If both halfedges point to the same face, returns an arbitrary one
int mesh::halfedgeAtFaceEdge(int f, int e) const {
    int he0 = E[e].he;
    int he1 = HE[he0].twin;
    if (HE[he0].face == f) {
        return he0;
    }
    return he1;
}

// Get an arbitrary halfedge loop starting from a certain index
std::vector<int> mesh::halfedgeLoop(int he) const {
    int max_search = active_e/2;
    int iter = 0;
    std::vector<int> heLoop;
    int he_curr = he;
    do {
        heLoop.push_back(he_curr);
        he_curr = HE[he_curr].next;
        iter++;
    } while((he_curr != he) && (iter <= max_search));

    return heLoop;
}

// Returns a CCW list of a face's vertices
// NOTE: In case of scrambled vertex ordering (ex. interior loops), it's safest to do this by halfedge
// origin flag: order this based on the halfedge origins vs. as halfedge dests
std::vector<Eigen::Vector3d> mesh::faceAdjVerts(int f, bool origin) const {
    std::vector<int> fVerts = faceAdjVertIdxs(f, origin);
    std::vector<Eigen::Vector3d> fVertsPos(fVerts.size());
    for (int v = 0; v < fVerts.size(); v++) {
        fVertsPos[v] = V[fVerts[v]].pos;
    }
    return fVertsPos;
}
// If given adjacent indices
std::vector<Eigen::Vector3d> mesh::adjVerts(std::vector<int> vertIdxs) const {
    std::vector<Eigen::Vector3d> vertsPos(vertIdxs.size());
    for (int v = 0; v < vertIdxs.size(); v++) {
        vertsPos[v] = V[vertIdxs[v]].pos;
    }
    return vertsPos;
}

// Returns a CCW list of face vertex indices
std::vector<int> mesh::faceAdjVertIdxs(int f, bool origin) const {
    std::vector<int> fHalfEdges = faceAdjHalfEdges(f);
    std::vector<int> fVerts(fHalfEdges.size());
    for (int he = 0; he < fHalfEdges.size(); he++) {
        int he_idx = he;
        if (origin) {
            he_idx = (he + fHalfEdges.size() - 1)%fHalfEdges.size();
        }
        fVerts[he] = HE[fHalfEdges[he_idx]].dest;
    }
    return fVerts;
}

// Returns a CCW list of a face's halfedges
std::vector<int> mesh::faceAdjHalfEdges(int f) const {
    int max_search = active_e / 2;
    std::vector<int> fHalfEdges;
    const int he0 = F[f].he;
    int he_curr = he0;
    int iter = 0;
    do {
        fHalfEdges.push_back(he_curr);
        he_curr = HE[he_curr].next;
        iter++;
    } while(he_curr != he0 && iter < max_search);
    return fHalfEdges;
}

std::vector<int> mesh::adjFaces(int elType, int elIdx) const {
    std::vector<int> adjF;
    if (elType == 0) {
        adjF = vertAdjFaces(elIdx);
    } else if (elType == 1) {
        adjF = edgeAdjFaces(elIdx);
    } else if (elType == 2) {
        adjF.push_back(elIdx);
    }
    return adjF;
}

// Returns the outgoing boundary HE if a vertex is a boundary vertex, else returns -1
int mesh::vertIsBoundary(int v) const {
    if (v < 0 || v >= V.size() || V[v].he < 0) {
        return {};
    }
    // Hard check: Check all outgoing halfedges
    // Safety in case the soft boundary halfedge rule is accidentally violated
    std::vector<int> adjHE = vertAdjHEs(v);
    for (int he = 0; he < adjHE.size(); he++) {
        if (HE[adjHE[he]].boundary) {
            return 1;
        }
    }
    return -1;
}

}   // namespace Mesh