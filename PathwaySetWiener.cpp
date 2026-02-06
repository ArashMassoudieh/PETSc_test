// PathwaySetWiener.cpp
#include "PathwaySetWiener.h"
#include "grid.h"
#include <random>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <iostream>

PathwaySet PathwaySetWiener::initialize(
    size_t nPaths,
    Release release,
    Grid2D* grid
){
    if (!grid) throw std::runtime_error("PathwaySetWiener::initialize: grid is null");

    PathwaySet ps((int)nPaths);

    const double Lx = grid->Lx();
    const double Ly = grid->Ly();
    const double dy = grid->dy();
    const int ny = grid->ny();

    std::random_device rd;
    std::mt19937 gen(rd());

    if (release == Release::CenterPoint) {
        const double x0 = 0.5 * Lx;
        const double y0 = 0.5 * Ly;
        for (size_t i = 0; i < nPaths; ++i) {
            Pathway p((int)i);
            p.addParticle(x0, y0, 0.0);
            ps.addPathway(p);
        }
        return ps;
    }

    if (release == Release::LeftUniform) {
        std::uniform_real_distribution<double> dist_y(0.0, Ly);
        for (size_t i = 0; i < nPaths; ++i) {
            const double y0 = dist_y(gen);
            Pathway p((int)i);
            p.addParticle(0.0, y0, 0.0);
            ps.addPathway(p);
        }
        return ps;
    }

    if (release == Release::LeftFluxWeighted) {
        const std::vector<double>& qx = grid->flux("qx");

        std::vector<double> flux_values;
        std::vector<double> y_positions;

        for (int j = 0; j < ny; ++j) {
            const double y = (j + 0.5) * dy;
            const int idx = j * (grid->nx() + 1) + 0;  // left boundary
            const double flux_val = qx[idx];
            if (flux_val > 0.0) {
                flux_values.push_back(flux_val);
                y_positions.push_back(y);
            }
        }
        if (flux_values.empty()) throw std::runtime_error("LeftFluxWeighted: no positive flux at left boundary");

        std::vector<double> cdf(flux_values.size());
        std::partial_sum(flux_values.begin(), flux_values.end(), cdf.begin());
        const double total = cdf.back();
        for (double& v : cdf) v /= total;

        std::uniform_real_distribution<double> U01(0.0, 1.0);
        std::uniform_real_distribution<double> offset(-dy/2.0, dy/2.0);

        for (size_t i = 0; i < nPaths; ++i) {
            const double r = U01(gen);
            auto it = std::lower_bound(cdf.begin(), cdf.end(), r);
            size_t k = (size_t)std::distance(cdf.begin(), it);
            if (k >= y_positions.size()) k = y_positions.size() - 1;

            double y0 = y_positions[k] + offset(gen);
            y0 = std::max(0.0, std::min(Ly, y0));

            Pathway p((int)i);
            p.addParticle(0.0, y0, 0.0);
            ps.addPathway(p);
        }
        return ps;
    }

    throw std::runtime_error("Unknown Release mode");
}

void PathwaySetWiener::trackAll(
    PathwaySet& paths,
    Grid2D* grid,
    const WienerParams& wp,
    int max_steps,
    const std::string& qx_name,
    const std::string& qy_name
){
    if (!grid) throw std::runtime_error("PathwaySetWiener::trackAll: grid is null");
    if (wp.mode == WienerMode::Off) throw std::runtime_error("PathwaySetWiener::trackAll: mode=Off");

    int completed = 0, failed = 0;

    for (size_t i = 0; i < paths.size(); ++i) {
        try {
            WienerParams wpp = wp;
            wpp.seed = wp.seed + 1315423911UL * (unsigned long)(i + 1);

            PathwayWiener::trackParticleWiener(paths[i], grid, wpp, max_steps, qx_name, qy_name);

            if (!paths[i].empty() && paths[i].last().isActive()) completed++;
            else failed++;

        } catch (const std::exception& e) {
            failed++;
            std::cerr << "Wiener tracking error pathway " << i << ": " << e.what() << "\n";
        }
    }

    std::cout << "Wiener tracking complete: " << completed << " completed, "
              << failed << " failed\n";
}
