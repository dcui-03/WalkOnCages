#include "curvenet.hpp"

#include "utils/utils.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <algorithm>


namespace Curvenet {
    // Returns a list of the tangents to a specified control
    std::vector<Eigen::Vector3d> curvenet::ctrlAdjTans(int c) const {
        assert(c >= 0 && c < static_cast<int>(C.size()));
        std::vector<Eigen::Vector3d> adjT;
        for (int he = 0; he < C[c].adjHE.size(); he++) {
            adjT.push_back(HE[C[c].adjHE[he]].tan);
        }
        return adjT;
    }

    // Returns a list of the splines adjacent to a control 
    std::vector<int> curvenet::ctrlAdjSplines(int c) const {
        assert(c >= 0 && c < static_cast<int>(C.size()));
        std::vector<int> adjS;
        std::vector<int> adjHE = C[c].adjHE;
        for (int he = 0; he < adjHE.size(); he++) {
            if (std::find(adjS.begin(), adjS.end(), HE[adjHE[he]].s) == adjS.end()) {
                adjS.push_back(HE[adjHE[he]].s);
            }
        }
        return adjS;
    }

    // Returns vertices adjacent to a spline
    std::vector<int> curvenet::splineAdjCtrls(int s) const {
        assert(s >= 0 && s < static_cast<int>(S.size()));
        std::vector<int> adjC;
        int he0 = S[s].he;
        adjC.push_back(HE[he0].origin);   // First vertex
        if (HE[HE[he0].twin].origin != HE[he0].origin) {
            adjC.push_back(HE[HE[he0].twin].origin);
        }
        return adjC;
    }

    // Because we sort the outgoing halfedges in CCW order,
    // we can get the local indices of adjHE s.t. the adjacent halfedge is part of the specified spline
    std::vector<int> curvenet::controlLocalSplineIdx(int c, int s) const {
        assert(c >= 0 && c < static_cast<int>(C.size()));
        assert(s >= 0 && s < static_cast<int>(S.size()));
        std::vector<int> local_idxs;
        const std::vector<int> cAdjHE = C[c].adjHE;
        for (int he = 0; he < cAdjHE.size(); he++) {
            if (HE[cAdjHE[he]].s == s) {
                local_idxs.push_back(he);
            }
        }
        return local_idxs;
    }
}   // namespace Curvenet