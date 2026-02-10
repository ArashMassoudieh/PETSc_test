// PathwayWiener.h
#pragma once

#include "Pathway.h"
#include <random>
#include <string>

class Grid2D;

enum class WienerMode {
    Off = 0,
    W1D_X,   // noise only in x
    W1D_Y,   // noise only in y
    W2D      // noise in x and y (ellipse if Dx != Dy)
};

struct WienerParams
{
    WienerMode mode = WienerMode::Off;

    double D = 0.01;        // diffusion coefficient
    double dt = 1e-3;       // fixed timestep
    unsigned long seed = 0; // RNG seed
    double rx = 1;
    double ry = 0.1;
};

// Add-on tracker: does NOT modify Particle/Pathway definitions.
class PathwayWiener
{
public:
    static void trackParticleWiener(
        Pathway& path,
        Grid2D* grid,
        const WienerParams& wp,
        int max_steps = 2000000,
        const std::string& qx_name = "qx",
        const std::string& qy_name = "qy"
    );
};
