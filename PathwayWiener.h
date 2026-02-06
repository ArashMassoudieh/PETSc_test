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

    double Dx = 0.0;        // diffusion in x  (used in W1D_X and W2D)
    double Dy = 0.0;        // diffusion in y  (used in W1D_Y and W2D)
    double dt = 1e-3;       // fixed timestep
    unsigned long seed = 0; // RNG seed
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
