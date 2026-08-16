#include "dcurvenet.hpp"

#include "utils/utils.hpp"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>
#include <utility>
#include <iostream>

namespace DCurvenet {

// Runtime operators
int dcurvenet::computedCNMats(Eigen::MatrixXd& f_dCN_flat, Eigen::MatrixXd& x_dCN) {
    int num_HE = numHalfedges();
    int num_V  = numVerts();

    if (f_dCN_flat.rows() != num_HE || f_dCN_flat.cols() != 9) {
        f_dCN_flat.resize(num_HE, 9);
    }

    if (x_dCN.rows() != num_V || x_dCN.cols() != 3) {
        x_dCN.resize(num_V, 3);
    }

    #pragma omp parallel for
    for (int he = 0; he < num_HE; he++) {
        f_dCN_flat.row(he) =
            Utils::flattenMatrix3d(HE[he].defData.defGrad).transpose();
    }

    #pragma omp parallel for
    for (int v = 0; v < num_V; v++) {
        x_dCN.row(v) = V[v].new_pos.transpose();
    }

    return 1;
}

}   // namespace DCurvenet