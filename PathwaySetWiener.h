// PathwaySetWiener.h
#pragma once

#include "PathwaySet.h"
#include "PathwayWiener.h"
#include <string>

class Grid2D;

class PathwaySetWiener
{
public:
    enum class Release {
        LeftUniform,
        LeftFluxWeighted,
        CenterPoint
    };

    static PathwaySet initialize(
        size_t nPaths,
        Release release,
        Grid2D* grid
    );

    static void trackAll(
        PathwaySet& paths,
        Grid2D* grid,
        const WienerParams& wp,
        int max_steps = 2000000,
        const std::string& qx_name = "qx",
        const std::string& qy_name = "qy"
    );
};
