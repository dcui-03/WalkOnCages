#include "mesh.hpp"

#include "../utils/decUtils.hpp"
#include "../utils/utils.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include <limits>
#include <queue>
#include <utility>
#include <algorithm>
#include <iostream>

// Utility functions for tracing mesh geodesics

namespace Mesh {

// Trace a "straightest" geodesic (ish) from the start vert to the end
int mesh::traceGeodesic(const Vert& start, 
                  const Vert& end, 
                  Eigen::Vector3d prevDirec,
                  vertProjData prevData, 
                  std::vector<Vert>& tracedVerts,
                  int depth,
                  const int max_depth,
                  bool recompute,
                  bool fast) {
    double eps = 1e-6;
    // First, do some simple tests for termination
    // It's good to have these to catch tiny directional drift
    if (prevData.elType == end.projData.elType && prevData.elIdx == end.projData.elIdx) {   // The next mesh element is exactly the goal
        return 1;
    } else if (start.projData.elType == 0 && end.projData.elType == 0) {    // Both are vertices
        if (vertPairToHE.find(std::make_pair(start.projData.elIdx, end.projData.elIdx)) != vertPairToHE.end()) {
            return 1;
        }
    } else if (start.projData.elType == 0 && end.projData.elType == 1) {    // Start is vertex, end is edge
        // Check if either end of the edge is the vertex
        if (HE[E[end.projData.elIdx].he].dest == start.projData.elIdx || HE[HE[E[end.projData.elIdx].he].twin].dest == start.projData.elIdx) {
            return 1;
        }
    } else if (start.projData.elType == 1 && end.projData.elType == 0) {    // Start is edge, end is vertex
        // Check if either end of the edge is the vertex
        if (HE[E[start.projData.elIdx].he].dest == end.projData.elIdx || HE[HE[E[start.projData.elIdx].he].twin].dest == end.projData.elIdx) {
            return 1;
        }
    } else if (start.projData.elType == 1 && end.projData.elType == 1) {    // Both are on edges
        // Check if they share an edge
        if (start.projData.elIdx == end.projData.elIdx) {
            return 1;
        }
    }

    // Otherwise need to do a face-wise check
    // Find shared faces between start and end, if any
    std::vector<int> startAdjF = adjFaces(start.projData.elType, start.projData.elIdx);
    std::vector<int> endAdjF = adjFaces(end.projData.elType, end.projData.elIdx);
    std::vector<int> sharedAdjF;
    for (int i = 0; i < startAdjF.size(); i++) {
        if (startAdjF[i] == -1) {
            continue;
        }
        for (int j = 0; j < endAdjF.size(); j++) {
            if (endAdjF[j] == -1) {
                continue;
            }
            if (startAdjF[i] == endAdjF[j]) {
                sharedAdjF.push_back(startAdjF[i]);
            }
        }
    }
    // Check if there are any shared faces (termination condition)
    if (sharedAdjF.size() >= 1) {   // Share at least one face
        if (fast) { // Fast version is a simple face check
            return 1;
        } else {    // Slow version does a visibility check
            for (int f_idx : sharedAdjF) {
                if (testVisibility(f_idx, start.pos, end.pos)) {
                    return 1;
                }
            }
        }
    }
    // If we still haven't found the vert after searching the max depth, assume that we are going in the wrong direc
    if (++depth > max_depth) {
        return -1;
    }
    std::cout << "Reached depth " << depth << std::endl;
    /*
    std::cout << "Start vertex: " << start.pos[0] << ", "
                                  << start.pos[1] << ", "
                                  << start.pos[2] << "; End vertex: "
                                  << end.pos[0] << ", "
                                  << end.pos[1] << ", "
                                  << end.pos[2] << std::endl;
    std::cout << "Previous direction: " << prevDirec[0] << ", "
                                        << prevDirec[1] << ", "
                                        << prevDirec[2] << std::endl;
    */
    // We have to walk, so we now compute the next walk direction
    vertProjData nextData;
    Eigen::Vector3d nextDirec;
    prevDirec.normalize();
    if (recompute && start.projData.elType == 2) {
        Eigen::Vector3d refDirec = (end.pos - start.pos);
        nextData.elType = prevData.elType;
        nextData.elIdx = prevData.elIdx;
        // If projection is degenerate, defer to valid previous version
        if (Utils::projectVectorOntoTangentPlane(F[nextData.elIdx].n, refDirec, nextDirec) <= eps) {
            nextDirec = prevDirec;
        } else if (nextDirec.dot(prevDirec) < 0.0) {    // Else do a soft check to make sure we're going the right way
            nextDirec *= -1;
        }
    } else {
        int success = -1;
        if (start.projData.elType == 0) {   // We are on a vert
            if (prevData.elType == 0) {
                success = nextEl_VertStart(start.projData.elIdx, prevData, prevDirec, nextDirec, nextData, true);
            } else {
                // std::cout << "Landed on vertex. Computing next direction" << std::endl;
                success = nextEl_Vert(start.projData.elIdx, prevData, prevDirec, nextDirec, nextData, true);
            }
        } else if (start.projData.elType == 1) {    // We are on an edge
            if (prevData.elType == 1) {             // Just started walking
                // std::cout << "Landed on edge. Computing next direction for prev is edge" << std::endl;
                success = nextEl_EdgeStart(start.projData.elIdx, prevData, prevDirec, nextDirec, nextData, true);
            } else if (prevData.elType == 2) {  // We landed on an edge from a face
                // std::cout << "Landed on edge. Computing next direction for prev is face" << std::endl;
                success = nextEl_Edge(start.projData.elIdx, prevData, prevDirec, nextDirec, nextData, true);
            } else {
                std::cout << "Landed on edge. Invalid next direction" << std::endl;
            }
        } else {        // We are on a face
            nextData = prevData;
            success = 1;
            // std::cout << "Landed on face. Computing next direction by projection" << std::endl;
            if (Utils::projectVectorOntoTangentPlane(F[nextData.elIdx].n, prevDirec, nextDirec) == -1) {
                return -1;
            }
        }
        if (success != 1) { // Could not find a new direction
            std::cout << "Could not find a suitable next direction" << std::endl;
            return -1;
        }
    }
    if (nextData.elType < 0) {  // No valid next direction
        std::cout << "No valid next direction" << std::endl;
        return -1;
    }
    /*
    std::cout << "Next direction: " << nextDirec[0] << ", "
                                        << nextDirec[1] << ", "
                                        << nextDirec[2] << std::endl;
    */
    // Case 1: Check if walk direction is on an edge, simply grab the other end vertex of the edge
    if (nextData.elType == 1) {
        // Opt to actually compute the matching direction instead of assuming the start is at a vertex
        // This way we can handle degenerate cases where the initial walk direction is on an edge
        int he = E[nextData.elIdx].he;
        Eigen::Vector3d heVec = (V[HE[he].dest].pos - V[HE[HE[he].twin].dest].pos).normalized();
        // Orient ourselves correctly
        int next;
        if (nextDirec.dot(heVec) >= 0.0) {
            next = HE[he].dest;
        } else {
            next = HE[HE[he].twin].dest;
        }
        // std::cout << "Walking along edge. Next vert found: " << next << std::endl;
        Vert nextVert = createVertex(V[next].pos, V[next].n, 2, -1, 0, next);
        // Recurse
        tracedVerts.push_back(nextVert);
        return traceGeodesic(nextVert, end, nextDirec, nextData, tracedVerts, depth, max_depth, true, fast);
    }

    // If we reached this point, we are definitely walking on a face
    // Project walk direction onto specified direction
    Eigen::Vector3d hit;
    vertProjData hit_Data;
    // std::cout << "Walking along face " << nextData.elIdx << std::endl;
    if (rayCastOnFace(nextData.elIdx, start.pos, nextDirec, hit, hit_Data) == -1) {
        std::cout << "Raycasting failed." << std::endl;
        return -1;
    }
    std::cout << "Raycast succeeded. Recursing." << std::endl;
    // Create a new vertex at intersection and append to list
    Vert nextVert = createVertex(hit, getNormal(hit_Data), 2, -1, hit_Data);

    tracedVerts.push_back(nextVert);
    return traceGeodesic(nextVert, end, nextDirec, nextData, tracedVerts, depth, max_depth, true, fast);
}

// Test whether the start and end are visible from each other on a particular face
// Returns true if so.
bool mesh::testVisibility(int f, Eigen::Vector3d start, Eigen::Vector3d end, double eps) {
    std::vector<Eigen::Vector3d> fVerts = faceAdjVerts(f);
    if (fVerts.size() == 3) {   // Triangles must be convex
        return true;
    }
    Eigen::Vector3d t1, t2;
    Utils::buildPlaneBasis(F[f].n, t1, t2);
    // For simplicity, assume tangent plane is centered on the start
    std::vector<Eigen::Vector2d> projFVerts(fVerts.size());
    Eigen::Vector2d start2D = Eigen::Vector2d::Zero();
    Eigen::Vector2d end2D = Utils::convertTo2D(end, start, t1, t2);
    for (int v = 0; v < fVerts.size(); v++) {
        Eigen::Vector3d proj3D = Utils::projectPointOntoPlane(F[f].n, start, fVerts[v]);
        projFVerts[v] = Utils::convertTo2D(proj3D, start, t1, t2);
    }
    Eigen::Vector2d direc = end2D - start2D;
    double t_end = direc.norm();
    direc.normalize();

    std::vector<double> intersections;
    // Raycast to find nearest segment, tracking the t-vals of intersections
    for (int v = 0; v < projFVerts.size(); v++) {
        int v_p1 = (v+projFVerts.size()+1) % projFVerts.size();
        double t, u;
        // Raycast to find nearest segment
        if (Utils::raycastToSegment2D(start2D, direc, projFVerts[v], projFVerts[v_p1], t, u)) {
            // Filter out any that are behind
            if (t >= eps) {
                intersections.push_back(t);
            }
        }
    }

    // If intersects none then we are outside the face and/or pointing the wrong way... not good, but check anyways
    if (intersections.size() == 0) {
        return false;
    }
    // Otherwise, check if we hit anything first
    for (double t : intersections) {
        // If we are w/in an epsilon, then ignore
        if (t <= t_end - eps) {
            return false;
        }
    }
    return true;
}

// Find the next intersection point while walking on a particular face
int mesh::rayCastOnFace(int f, Eigen::Vector3d start, Eigen::Vector3d direc, Eigen::Vector3d& hit, vertProjData& hitData, double eps) {
    double tol = 1e-5 * bboxDiag;
    Eigen::Vector3d projDirec;
    Utils::projectVectorOntoTangentPlane(F[f].n, direc, projDirec);
    std::vector<int> fVertIdxs = faceAdjVertIdxs(f);
    std::vector<Eigen::Vector3d> fVerts = adjVerts(fVertIdxs);

    Eigen::Vector3d t1, t2;
    Utils::buildPlaneBasis(F[f].n, t1, t2);
    // For simplicity, assume tangent plane is centered on the start
    std::vector<Eigen::Vector2d> projFVerts(fVerts.size());
    Eigen::Vector2d start2D = Eigen::Vector2d::Zero();
    Eigen::Vector2d direc2D = Utils::convertTo2D(start + projDirec, start, t1, t2);
    for (int v = 0; v < fVerts.size(); v++) {
        Eigen::Vector3d proj3D = Utils::projectPointOntoPlane(F[f].n, start, fVerts[v]);
        projFVerts[v] = Utils::convertTo2D(proj3D, start, t1, t2);
    }
    direc2D.normalize();

    // intersections are a pair of local segment index and the u parameter along that segment
    std::vector<std::pair<int, double>> intersections;
    std::vector<double> intersections_t;
    // Raycast to find nearest segment, tracking the t-vals of intersections
    for (int v = 0; v < projFVerts.size(); v++) {
        int v_p1 = (v+projFVerts.size()+1) % projFVerts.size();
        double t, u;
        // Raycast to find nearest segment
        if (Utils::raycastToSegment2D(start2D, direc2D, projFVerts[v], projFVerts[v_p1], t, u)) {
            // Filter out any that are behind
            if (t >= eps) {
                intersections.push_back(std::make_pair(v, u));
                intersections_t.push_back(t);
            }
        }
    }

    // If intersects none then we are outside the face and/or pointing the wrong way... not good, but check anyways
    if (intersections.size() == 0) {
        return -1;
    }
    // Otherwise, check what we hit first
    int nearest_idx = 0;
    double nearest = intersections_t[0];
    for (int i = 1; i < intersections_t.size(); i++) {
        if (intersections_t[i] < nearest) {
            nearest = intersections_t[i];
            nearest_idx = i;
        }
    }
    // Compute the intersection
    int v = intersections[nearest_idx].first;
    double u = intersections[nearest_idx].second;
    hit = fVerts[v] + u * (fVerts[(v+1)%fVerts.size()] - fVerts[v]);
    // Apply snapping as necessary
    hitData.elType = 1;
    hitData.elIdx = HE[vertPairToHE[std::make_pair(fVertIdxs[v], fVertIdxs[(v+1)%fVertIdxs.size()])]].edge;
    int v_next = (v+1) % fVerts.size();
    if ((hit - fVerts[v]).norm() <= tol) {
        hit = fVerts[v];
        hitData.elType = 0;
        hitData.elIdx = fVertIdxs[v];
    } else if ((hit - fVerts[v_next]).norm() <= tol) {
        hit = fVerts[v_next];
        hitData.elType = 0;
        hitData.elIdx = fVertIdxs[v_next];
    }
    return 1;
}

// Compute the next walk element given that we intersected with an edge
// Returns the next 
int mesh::nextEl_Edge(int e, const vertProjData& originData, const Eigen::Vector3d& prev_direc, 
                    Eigen::Vector3d& next_direc, vertProjData& nextData, bool bdy_snap) {
                        // Get adjacent faces
    nextData.elType = 2;
    if (originData.elType != 2 || originData.elIdx < 0 || originData.elIdx >= F.size()) {
        return -1;
    }
    std::vector<int> adjF = edgeAdjFaces(e);
    nextData.elIdx = adjF[0];
    if (originData.elIdx == adjF[0]) {
        nextData.elIdx = adjF[1];
    }
    // Projection step for safety
    Eigen::Vector3d proj_direc;
    double valid = Utils::projectVectorOntoTangentPlane(F[originData.elIdx].n, prev_direc, proj_direc);
    if (valid <= 0.0) {
        return -1;
    }
    // Handle boundary case first
    // Snap to the better boundary edge direction
    if (nextData.elIdx == -1) {
        if (bdy_snap) {
            nextData.elIdx = e;
            nextData.elType = 1;
            // Find the better fit direction 
            int v1 = HE[E[e].he].dest;
            int v0 = HE[HE[E[e].he].twin].dest;
            if (((V[v1].pos - V[v0].pos).normalized()).dot(proj_direc) >= 0) {
                next_direc = (V[v1].pos - V[v0].pos).normalized();
            } else {
                next_direc = (V[v0].pos - V[v1].pos).normalized();
            }
            return 1;
        } else {    // Hit a dead end
            return -1;
        }
    }
    // Otherwise, we just rotate onto the plane of the adajcent face
    // Get the halfedge of the incoming face
    int he0 = E[e].he;
    if (HE[he0].face == originData.elIdx) {
        he0 = HE[he0].twin;
    }
    // Compute the local tangent basis of the origin face
    Eigen::Vector3d tangent = ((V[HE[he0].dest].pos - V[HE[HE[he0].twin].dest].pos).normalized());
    // Normalized for numerical safety
    Eigen::Vector3d biN_origin = tangent.cross(F[originData.elIdx].n).normalized();
    Eigen::Vector3d biN_next = tangent.cross(F[nextData.elIdx].n).normalized();
    Eigen::Matrix3d rot = Utils::computeRotation(biN_origin, biN_next);
    next_direc = rot * proj_direc;
    return 1;
}

// Find the next walk direction given we are starting from an edge
int mesh::nextEl_EdgeStart(int e, const vertProjData& originData, const Eigen::Vector3d& start_direc,
                           Eigen::Vector3d& next_direc, vertProjData& nextData, bool bdy_snap, double eps) {
    nextData.elType = -1;
    nextData.elIdx = -1;
    next_direc.setZero();
    if (e < 0 || e >= E.size()) {
        return -1;
    }
    if (originData.elType != 1 || originData.elIdx != e) {
        return -1;
    }
    // Get adjacent attributes
    int he0 = E[e].he;
    int he1 = HE[he0].twin;
    int v0 = HE[he1].dest;
    int v1 = HE[he0].dest;
    Eigen::Vector3d tangent = (V[v1].pos - V[v0].pos).normalized();
    int f0 = HE[he0].face;
    int f1 = HE[he1].face;
    double dot0 = 0.0;
    double dot1 = 0.0;

    // Check the edge itself
    // First, check if the desired direction is essentially along the edge
    Eigen::Vector3d edgeProj;
    double edgeProjLen = Utils::projectVectorOntoTangentPlane(E[e].n, start_direc, edgeProj);
    Eigen::Vector3d binorm0 = (E[e].n.cross(tangent)).normalized();
    Eigen::Vector3d binorm1 = -1 * binorm0;
    if (edgeProjLen > eps) {
        double side = edgeProj.dot(binorm0);
        if (std::abs(side) <= 1e-5) {
            return snapWalkToEdge(e, start_direc, next_direc, nextData, eps);
        } 
    } else {
        return -1;
    }

    // The face whose outward binormal aligns more with the walk direction is the face we are coming from (prev)
    dot0 = edgeProj.dot(binorm0);
    dot1 = edgeProj.dot(binorm1);
    vertProjData tempOrigin;
    nextData.elType = 2;
    if (dot0 > dot1) {
        nextData.elIdx = f0;
    } else {
        nextData.elIdx = f1;
    }

    if (nextData.elIdx == -1) { // Snap to boundary edge if needed
        if (bdy_snap) {
            nextData.elType = 1;
            nextData.elIdx = e;
            if (tangent.dot(edgeProj) >= 0.0) {
                next_direc = tangent;
            } else {
                next_direc = -1 * tangent;
            }
        } else {
            return -1;
        }
    }

    // We must be on a face; get the new direction
    Utils::projectVectorOntoTangentPlane(F[nextData.elIdx].n, edgeProj, next_direc);
    return 1;
}

// Helper for next edge that performs edge snapping
int mesh::snapWalkToEdge(int e, const Eigen::Vector3d& direc, Eigen::Vector3d& next_direc,
                         vertProjData& nextData, double eps) const {
    nextData.elType = -1;
    nextData.elIdx = -1;
    next_direc.setZero();
    // Get adjacent attributes
    int he = E[e].he;
    int twin = HE[he].twin;
    int v0 = HE[twin].dest;
    int v1 = HE[he].dest;

    Eigen::Vector3d edgeVec = V[v1].pos - V[v0].pos;
    if (edgeVec.norm() <= eps) {
        return -1;
    }
    edgeVec.normalize();

    // Compute which side of the edge 
    if (direc.dot(edgeVec) >= 0.0) {
        next_direc = edgeVec;
    } else {
        next_direc = -1*edgeVec;
    }

    nextData.elType = 1;
    nextData.elIdx = e;
    return 1;
}

// Compute the next walk element given that we intersected with a vertex
int mesh::nextEl_Vert(int v, const vertProjData& originData, 
                    const Eigen::Vector3d& prev_direc, Eigen::Vector3d& next_direc, 
                    vertProjData& nextData, bool bdy_snap, double eps) {
    nextData.elType = 2;
    // Gather all of the adjacent halfedges and faces
    std::vector<int> adjHE = vertAdjHEs(v);
    // Figure out which halfedge "matches" the face origin and also find the local umbrella angle sum
    int localIdx = -1;
    double angleSum = 0.0;
    std::vector<double> angleBuckets(adjHE.size());
    std::vector<Eigen::Vector3d> cornerNormals(adjHE.size());
    for (int he_idx = 0; he_idx < adjHE.size(); he_idx++) {
        int he0 = adjHE[he_idx];
        int he1 = adjHE[(he_idx + 1) % adjHE.size()];
        // Compute angle between the two vectors by first computing an axis
        Eigen::Vector3d vec0 = V[HE[he0].dest].pos - V[v].pos;
        Eigen::Vector3d vec1 = V[HE[he1].dest].pos - V[v].pos;
        Eigen::Vector3d axis = vec0.cross(vec1);
        if (axis.norm() <= eps) {
            axis = F[HE[he0].face].n;
        }
        cornerNormals[he_idx] = axis.normalized();
        // Get the positive signed angle
        angleBuckets[he_idx] = Utils::signedAngle(vec0, vec1, axis, true);
        angleSum += angleBuckets[he_idx];
        // Get the matching face as a local halfedge index in adjHE
        if (originData.elType == 1 && HE[adjHE[he_idx]].edge == originData.elIdx) {
            localIdx = he_idx;
        } else if (originData.elType == 2 && HE[adjHE[he_idx]].face == originData.elIdx) {
            localIdx = he_idx;
        }
    }
    if (localIdx == -1) {   // Did not find a match
        return -1;
    }
    // Compute the projected walk direction onto the corner normal's plane
    Eigen::Vector3d proj_direc;
    double valid = Utils::projectVectorOntoTangentPlane(cornerNormals[localIdx], prev_direc, proj_direc);
    if (valid <= 0.0) {
        return -1;
    }
    double angleToFace = Utils::signedAngle(V[HE[adjHE[localIdx]].dest].pos - V[v].pos, -1 * prev_direc, cornerNormals[localIdx], true);

    // Get the halfway angle
    double targetAngle = angleSum / 2.0;
    int stopHE = -1;
    targetAngle += angleToFace; // Add this on so we can start accumulating angles from the halfedge
    // Iterate over faces to get the target angle
    for (int f = 0; f < adjHE.size(); f++) {
        int f_curr = (localIdx + f) % adjHE.size();
        if (targetAngle <= angleBuckets[f_curr]) {  // Found the face
            stopHE = f_curr;
            break;
        } else {    // Keep going
            targetAngle -= angleBuckets[f_curr];
        }
    }
    if (stopHE < 0) {   // Just in case
        return -1;
    }
    // Use the computed face and remnant angle to determine the vector
    Eigen::Vector3d heVec = (V[HE[adjHE[stopHE]].dest].pos - V[v].pos).normalized();
    // To make this well-posed, compute corner normal as the average of adjacent corner normals
    if (cornerNormals[stopHE].norm() <= eps) {
        cornerNormals[stopHE] = (cornerNormals[(stopHE + 1)%adjHE.size()] + cornerNormals[(stopHE + adjHE.size() - 1)%adjHE.size()]) / 2;
    }
    Eigen::Matrix3d rot = Utils::computeRotation(cornerNormals[stopHE], targetAngle);
    next_direc = rot * heVec;

    Eigen::Vector3d heVec1 = (V[HE[adjHE[(stopHE+1)%adjHE.size()]].dest].pos - V[v].pos).normalized();
    nextData.elIdx = HE[adjHE[stopHE]].face;
    // Process boundary if we hit one
    if (nextData.elIdx == -1) {
        if (bdy_snap) {
            nextData.elType = 1;
            // Compute best adjacent edge to snap to
            if (next_direc.dot(heVec) >= next_direc.dot(heVec1)) {
                next_direc = heVec;
                nextData.elIdx = HE[adjHE[stopHE]].edge;
            } else {
                next_direc = heVec1;
                nextData.elIdx = HE[adjHE[(stopHE+1)%adjHE.size()]].edge;
            }
            return 1;
        } else {    // Failure
            return -1;
        }
    }
    // Check for edge snaps
    if (next_direc.normalized().dot(heVec) >= 1.0 - eps) {
        next_direc = heVec;
        nextData.elIdx = HE[adjHE[stopHE]].edge;
        nextData.elType = 1;
        return 1;
    } else if (next_direc.normalized().dot(heVec1) >= 1.0 - eps) {
        next_direc = heVec1;
        nextData.elIdx = HE[adjHE[(stopHE+1)%adjHE.size()]].edge;
        nextData.elType = 1;
        return 1;
    }
    return 1;
}

int mesh::nextEl_VertStart(int v, const vertProjData& originData, const Eigen::Vector3d& start_direc,
                           Eigen::Vector3d& next_direc, vertProjData& nextData, bool bdy_snap, double eps) {
    nextData.elType = -1;
    nextData.elIdx = -1;
    next_direc.setZero();
    if (v < 0 || v >= V.size()) {
        return -1;
    }
    if (originData.elType != 0 || originData.elIdx != v) {
        return -1;
    }
    // Get adjacent attributes
    Eigen::Vector3d t0, t1;
    Utils::buildPlaneBasis(V[v].n, t0, t1);
    std::vector<int> adjHE = vertAdjHEs(v);
    std::vector<double> adjAngles(adjHE.size());
    for (int he = 0; he < adjHE.size(); he++) {
        double theta;
        if (!Utils::directionAngleInPlane(V[v].pos, V[HE[adjHE[he]].dest].pos, V[v].n, t0, t1, theta)) {
            return -1;
        }
        adjAngles[he] = theta;
    }
    double start_theta;
    if (!Utils::directionAngleInPlane(V[v].pos, V[v].pos + start_direc, V[v].n, t0, t1, start_theta)) {
        return -1;
    }
    int face_he = -1;
    for (int he = 0; he < adjHE.size(); he++) { // Find which "corner" we are in
        double curr_angle = adjAngles[he];
        double next_angle = adjAngles[(he+1) % adjHE.size()];

        double d = std::abs(adjAngles[he] - start_theta);
        d = std::min(d, 2.0 * M_PI - d);
        if (d <= eps) { // i.e., we are on the edge of the halfedge
            next_direc = (V[HE[adjHE[he]].dest].pos - V[v].pos).normalized();
            nextData.elType = 1;
            nextData.elIdx = HE[adjHE[he]].edge;
            return 1;
        }
        // Wrap-around check
        double temp_start = start_theta;
        if (adjAngles[he] > next_angle) {
            next_angle += 2.0*M_PI;
            if (temp_start < adjAngles[he]) {
                temp_start += 2.0*M_PI;
            }
        }
        // Once we find the right one, return
        if (temp_start < next_angle && temp_start > adjAngles[he]) {
            face_he = he;
        }
    }
    if (face_he == -1) {    // NOTE: face_he is an index local to adjHE
        return -1;
    }

    nextData.elType = 2;
    nextData.elIdx = HE[adjHE[face_he]].face;
    // To get the estimated new direction, project onto vertex tangent plane, and then project again onto the corner normal plane
    Eigen::Vector3d v_ProjDirec;
    Utils::projectVectorOntoTangentPlane(V[v].n, start_direc, v_ProjDirec);
    int next_he = (face_he + 1) % adjHE.size();
    Eigen::Vector3d tangent0 = (V[HE[adjHE[face_he]].dest].pos - V[v].pos).normalized();
    Eigen::Vector3d tangent1 = (V[HE[adjHE[next_he]].dest].pos - V[v].pos).normalized();
    Eigen::Vector3d corner_normal = tangent0.cross(tangent1);
    if (corner_normal.norm() <= eps) {
        if (nextData.elIdx != -1) {
            corner_normal = F[nextData.elIdx].n;
        } else {
            corner_normal = V[v].n;
        }
    } else {
        corner_normal.normalize();
    }
    if (Utils::projectVectorOntoTangentPlane(corner_normal, v_ProjDirec, next_direc) <= eps) {
        return -1;
    };
    next_direc.normalize();
    // Handle the boundary case
    if (HE[adjHE[face_he]].face == -1) {    // Handle boundary case
        if (bdy_snap) {
            nextData.elType = 1;
            if (tangent0.dot(next_direc) >= tangent1.dot(next_direc)) {
                nextData.elIdx = HE[adjHE[face_he]].edge;
            } else {
                nextData.elIdx = HE[adjHE[next_he]].edge;
            }
        } else {
            return -1;
        }
    }
    
    return 1;
    
}

}   // namespace Mesh