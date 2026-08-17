#include "dcurvenet.hpp"

namespace Polynet {

    bool dcurvenet::isPositiveHalfedge(int he) const {
        return E[HE[he].edge].he == he;
    }

}   // namespace Polynet