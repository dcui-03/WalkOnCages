// dcurvenet_types.hpp
#pragma once

// File with dcurvenet-only linkage data (parallel to polynet's V/C, see polynet_types.hpp)

namespace Polynet {

    // Per-vertex data linking a dcurvenet vert back to its source curvenet control, if any
    struct dCNVertData {
        int cn_idx = -1;    // curvenet index if coincident with a control vertex
        int cn_type = -1;   // curvenet vertex type if coincident with a control vertex
    };
}   // namespace Polynet
