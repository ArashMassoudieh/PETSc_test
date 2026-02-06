// PathwayWiener.cpp
#include "PathwayWiener.h"
#include "grid.h"
#include <cmath>
#include <stdexcept>

static inline double safe_sqrt(double x)
{
    return (x > 0.0) ? std::sqrt(x) : 0.0;
}

void PathwayWiener::trackParticleWiener(
    Pathway& path,
    Grid2D* grid,
    const WienerParams& wp,
    int max_steps,
    const std::string& qx_name,
    const std::string& qy_name
){
    if (!grid) throw std::runtime_error("trackParticleWiener: grid is null");
    if (path.empty()) throw std::runtime_error("trackParticleWiener: pathway is empty");
    if (wp.dt <= 0.0) throw std::runtime_error("trackParticleWiener: dt must be positive");
    if (wp.Dx < 0.0 || wp.Dy < 0.0) throw std::runtime_error("trackParticleWiener: Dx,Dy must be >= 0");
    if (wp.mode == WienerMode::Off) throw std::runtime_error("trackParticleWiener: mode=Off");

    const double Lx = grid->Lx();
    const double Ly = grid->Ly();

    Particle cur = path.last();
    double t = cur.t();

    std::mt19937_64 rng(wp.seed);
    std::normal_distribution<double> N01(0.0, 1.0);

    // Set stochastic amplitudes based on selected mode
    double sx = 0.0, sy = 0.0;
    if (wp.mode == WienerMode::W1D_X) {
        sx = safe_sqrt(2.0 * wp.Dx * wp.dt);
        sy = 0.0;
    } else if (wp.mode == WienerMode::W1D_Y) {
        sx = 0.0;
        sy = safe_sqrt(2.0 * wp.Dy * wp.dt);
    } else { // W2D
        sx = safe_sqrt(2.0 * wp.Dx * wp.dt);
        sy = safe_sqrt(2.0 * wp.Dy * wp.dt);
    }

    int step = 0;
    while (step < max_steps) {

        // already completed
        if (cur.x() >= Lx) {
            cur.setPosition(Lx, cur.y());
            cur.setT(t);
            cur.setActive(true);
            path.addParticle(cur);
            break;
        }

        // Drift (velocity at current position)
        double vx = 0.0, vy = 0.0;
        try {
            auto vel = grid->getVelocityAt(cur.x(), cur.y(), qx_name, qy_name);
            vx = vel.first;
            vy = vel.second;
        } catch (...) {
            cur.setActive(false);
            path.addParticle(cur);
            break;
        }

        const double dWx = N01(rng);
        const double dWy = N01(rng);

        const double new_x = cur.x() + vx * wp.dt + sx * dWx;
        const double new_y = cur.y() + vy * wp.dt + sy * dWy;

        t += wp.dt;

        // Right boundary => completed
        if (new_x >= Lx) {
            cur.setPosition(Lx, new_y);
            cur.setT(t);
            cur.setVelocity(vx, vy);
            cur.setActive(true);
            path.addParticle(cur);
            break;
        }

        // Top/bottom => failed
        if (new_y < 0.0 || new_y > Ly) {
            cur.setPosition(new_x, new_y);
            cur.setT(t);
            cur.setVelocity(vx, vy);
            cur.setActive(false);
            path.addParticle(cur);
            break;
        }

        // Normal step
        cur.setPosition(new_x, new_y);
        cur.setT(t);
        cur.setVelocity(vx, vy);
        cur.setActive(true);
        path.addParticle(cur);

        ++step;
    }

    if (step >= max_steps) {
        Particle last = path.last();
        last.setActive(false);
    }
}
