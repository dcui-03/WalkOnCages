#include "polynet.hpp"

#include <vector>
#include <algorithm>

namespace Polynet {
    // Compute what kind of vertex this is using the valence of adjacent edges
    int polynet::assignVertType(int v) {
        std::vector<int> adjHE = V[v].adjHE;
        // NOTE: Treat self-loops as not endpoints
        V[v].vType = std::min(static_cast<int>(adjHE.size()), 3);
        return V[v].vType;
    }
    int polynet::assignVertTypeAll() {
        for (int v = 0; v < V.size(); v++) {
            if (V[v].active) {
                assignVertType(v);
            }
        }
        return 1;
    }

    // Trace edge chains into curves
    int polynet::traceCurves() {
        C.clear();
        std::vector<bool> edgeFound(E.size(), false);
        // Search from endpoint and intersection vertices
        for (int v = 0; v < V.size(); v++) {
            int vType = V[v].vType;
            if (vType == 2) {  // Skip all loops for now
                continue;
            }
            const std::vector<int> adjHE = V[v].adjHE;
            // Start from each outgoing halfedge and trace until we hit a stop point
            for (int he = 0; he < adjHE.size(); he++) {
                int e0 = HE[adjHE[he]].edge;
                if (edgeFound[e0]) {   // Skip any edges we've already seen
                    continue;
                }

                // Otherwise, need to trace a curve
                int c = C.size();
                C.emplace_back();
                C[c].start = v;
                C[c].he_start = adjHE[he];
                bool curveEnd = false;
                int counter = 0;
                int curr_he = adjHE[he];
                // Trace edges until we hit an intersection or an endpoint
                do {
                    int e = HE[curr_he].edge;
                    E[e].curve = c;
                    edgeFound[e] = true;
                    counter++;
                    int curr_end = HE[curr_he].dest;
                    // Go to next edge
                    if (V[curr_end].vType == 1 || V[curr_end].vType == 3) {
                        curveEnd = true;
                        C[c].end = curr_end;
                        C[c].he_end = HE[curr_he].twin;
                    } else {    // Must be 2 outgoing HE's from this one. Pick the one we haven't gone to yet
                        int next_he = nextHEFromVert(curr_he, curr_end);
                        if (next_he < 0) {
                            return -1;
                        }
                        curr_he = next_he;
                    }
                } while (!curveEnd && counter < E.size());
                if (!curveEnd) {    // Error check
                    return -1;
                }
            }
        }

        // Start from remaining edges, which must form closed loops
        for (int v = 0; v < V.size(); v++) {
            int vType = V[v].vType;
            if (vType == 1 || vType == 3) {  // Skip endpoints and intersections
                continue;
            }
            const std::vector<int> adjHE = V[v].adjHE;
            // Start from each outgoing halfedge and trace until we hit the start point again
            for (int he = 0; he < adjHE.size(); he++) {
                int e0 = HE[adjHE[he]].edge;
                if (edgeFound[e0]) {   // Skip any edges we've already seen
                    continue;
                }

                // Otherwise, need to trace a curve
                int c = C.size();
                C.emplace_back();
                C[c].start = v;
                C[c].he_start = adjHE[he];
                bool curveEnd = false;
                int counter = 0;
                int curr_he = adjHE[he];
                // Trace edges until we hit the start again
                do {
                    int e = HE[curr_he].edge;
                    E[e].curve = c;
                    edgeFound[e] = true;
                    counter++;
                    int curr_end = HE[curr_he].dest;
                    // Go to next edge
                    if (curr_end == v) {    // Hit the start point again
                        curveEnd = true;
                        C[c].end = curr_end;
                        C[c].he_end = HE[curr_he].twin;
                    } else {    // Must be 2 outgoing HE's from this one. Pick the one we haven't gone to yet
                        int next_he = nextHEFromVert(curr_he, curr_end);
                        if (next_he < 0) {
                            return -1;
                        }
                        curr_he = next_he;
                    }
                } while (!curveEnd && counter < E.size());
                if (!curveEnd) {    // Error check
                    return -1;
                }
            }
        }
        return 1;
    }

    // Helper for traceCurves
    // Compute the next halfedge from a degree 2 vertex
    int polynet::nextHEFromVert(int curr_he, int curr_end) {
        if (curr_end < 0 || curr_end >= V.size()) {
            return -1;
        }
        if (V[curr_end].vType != 2) {
            return -1;
        }
        const std::vector<int>& adjHE = V[curr_end].adjHE;
        if (adjHE.size() != 2) {
            return -1;
        }
        int back_he = HE[curr_he].twin;
        if (adjHE[0] == back_he) {
            return adjHE[1];
        }
        if (adjHE[1] == back_he) {
            return adjHE[0];
        }
        // The halfedge we arrived on does not actually end at this vertex according to the vertex's adjacency list
        return -1;
    }
}   // namespace Polynet
