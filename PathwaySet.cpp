// PathwaySet.cpp
#include "PathwaySet.h"
#include "grid.h"
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <random>
#include <iostream>
#include <iomanip>
#include <sstream>

PathwaySet::PathwaySet()
{
}

PathwaySet::PathwaySet(int reserve_size)
{
    pathways_.reserve(reserve_size);
}

void PathwaySet::addPathway(const Pathway& p)
{
    pathways_.push_back(p);
}

void PathwaySet::addPathway(int id)
{
    pathways_.emplace_back(id);
}

void PathwaySet::Initialize(size_t number_of_pathways, Weighting weighting, Grid2D* grid)
{
    if (!grid) {
        throw std::runtime_error("Initialize: grid pointer is null");
    }

    pathways_.clear();
    pathways_.reserve(number_of_pathways);

    const double Ly = grid->Ly();
    const double dy = grid->dy();
    const int ny = grid->ny();

    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());

    if (weighting == Weighting::Uniform) {
        // Uniform random distribution in y
        std::uniform_real_distribution<double> dist_y(0.0, Ly);

        for (size_t i = 0; i < number_of_pathways; ++i) {
            double y = dist_y(gen);
            Pathway path(i);
            path.addParticle(0.0, y, 0.0);  // x=0, y=random, t=0
            pathways_.push_back(path);
        }
    }
    else if (weighting == Weighting::FluxWeighted) {
        // Get qx flux at left boundary (x = 0)
        const std::vector<double>& qx = grid->flux("qx");

        // Collect positive fluxes and their y-positions
        std::vector<double> flux_values;
        std::vector<double> y_positions;

        for (int j = 0; j < ny; ++j) {
            double y = (j + 0.5) * dy;  // Cell center y-coordinate
            int idx = j * (grid->nx() + 1) + 0;  // qx index at left boundary (i=0)
            double flux_val = qx[idx];

            if (flux_val > 0.0) {
                flux_values.push_back(flux_val);
                y_positions.push_back(y);
            }
        }

        if (flux_values.empty()) {
            throw std::runtime_error("Initialize: no positive flux at left boundary");
        }

        // Create cumulative distribution function (CDF)
        std::vector<double> cdf(flux_values.size());
        std::partial_sum(flux_values.begin(), flux_values.end(), cdf.begin());
        double total_flux = cdf.back();

        // Normalize CDF
        for (auto& val : cdf) {
            val /= total_flux;
        }

        // Sample pathways based on flux weighting
        std::uniform_real_distribution<double> dist_uniform(0.0, 1.0);

        for (size_t i = 0; i < number_of_pathways; ++i) {
            double r = dist_uniform(gen);

            // Find which bin this sample falls into
            auto it = std::lower_bound(cdf.begin(), cdf.end(), r);
            size_t idx = std::distance(cdf.begin(), it);

            if (idx >= y_positions.size()) {
                idx = y_positions.size() - 1;
            }

            // Get y position (with small random offset within the cell)
            double y_base = y_positions[idx];
            std::uniform_real_distribution<double> dist_offset(-dy/2.0, dy/2.0);
            double y = y_base + dist_offset(gen);

            // Clamp to domain
            y = std::max(0.0, std::min(Ly, y));

            Pathway path(i);
            path.addParticle(0.0, y, 0.0);  // x=0, y=flux-weighted, t=0
            pathways_.push_back(path);
        }
    }
    else {
        throw std::runtime_error("Initialize: unknown weighting type");
    }

    std::cout << "Initialized " << number_of_pathways << " pathways at left boundary" << std::endl;
}

void PathwaySet::trackAllPathways(Grid2D* grid, double dx_step,
                                  const std::string& qx_name,
                                  const std::string& qy_name)
{
    if (!grid) {
        throw std::runtime_error("trackAllPathways: grid pointer is null");
    }

    std::cout << "Tracking " << pathways_.size() << " pathways..." << std::endl;

    int completed = 0;
    int failed = 0;

    for (size_t i = 0; i < pathways_.size(); ++i) {
        try {
            pathways_[i].trackParticle(grid, dx_step, qx_name, qy_name);

            if (pathways_[i].last().isActive()) {
                completed++;
            } else {
                failed++;
            }

            // Progress reporting
            if ((i + 1) % 100 == 0 || i == pathways_.size() - 1) {
                std::cout << "  Tracked " << (i + 1) << "/" << pathways_.size()
                << " pathways (completed: " << completed
                << ", failed: " << failed << ")" << std::endl;
            }
        } catch (const std::exception& e) {
            std::cerr << "Error tracking pathway " << i << ": " << e.what() << std::endl;
            failed++;
        }
    }

    std::cout << "Tracking complete: " << completed << " completed, "
              << failed << " failed" << std::endl;
}

Pathway& PathwaySet::operator[](size_t i)
{
    return pathways_[i];
}

const Pathway& PathwaySet::operator[](size_t i) const
{
    return pathways_[i];
}

Pathway& PathwaySet::at(size_t i)
{
    return pathways_.at(i);
}

const Pathway& PathwaySet::at(size_t i) const
{
    return pathways_.at(i);
}

std::vector<Pathway>::iterator PathwaySet::begin()
{
    return pathways_.begin();
}

std::vector<Pathway>::iterator PathwaySet::end()
{
    return pathways_.end();
}

std::vector<Pathway>::const_iterator PathwaySet::begin() const
{
    return pathways_.begin();
}

std::vector<Pathway>::const_iterator PathwaySet::end() const
{
    return pathways_.end();
}

size_t PathwaySet::size() const
{
    return pathways_.size();
}

bool PathwaySet::empty() const
{
    return pathways_.empty();
}

void PathwaySet::reserve(size_t n)
{
    pathways_.reserve(n);
}

void PathwaySet::clear()
{
    pathways_.clear();
}

size_t PathwaySet::countActive() const
{
    size_t count = 0;
    for (const auto& p : pathways_) {
        if (!p.empty() && p.last().isActive()) {
            ++count;
        }
    }
    return count;
}

size_t PathwaySet::countCompleted() const
{
    size_t count = 0;
    for (const auto& p : pathways_) {
        if (!p.empty() && p.last().isActive() && p.last().x() >= p.last().x()) {
            ++count;
        }
    }
    return count;
}

void PathwaySet::removeIncomplete()
{
    pathways_.erase(
        std::remove_if(pathways_.begin(), pathways_.end(),
                       [](const Pathway& p) {
                           return p.empty() || !p.last().isActive();
                       }),
        pathways_.end()
        );
}

void PathwaySet::writeToFile(const std::string& filename) const
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    file << "# PathwaySet with " << pathways_.size() << " pathways\n";
    file << "# pathway_id num_points path_length travel_time active\n";

    for (const auto& p : pathways_) {
        if (!p.empty()) {
            file << p.id() << " " << p.size() << " "
                 << p.pathLength() << " " << p.duration() << " "
                 << (p.last().isActive() ? 1 : 0) << "\n";
        }
    }
    file.close();
}

void PathwaySet::readFromFile(const std::string& filename)
{
    // Note: This only reads the summary, not full pathway data
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    pathways_.clear();
    std::string line;

    // Skip headers
    while (std::getline(file, line) && line[0] == '#') {}

    // This is a simplified version - full implementation would need
    // to read individual pathway files
    file.close();
}

void PathwaySet::writeAllVTK(const std::string& prefix) const
{
    for (const auto& p : pathways_) {
        if (!p.empty()) {
            std::ostringstream filename;
            filename << prefix << "_pathway_"
                     << std::setw(6) << std::setfill('0') << p.id() << ".vtk";
            p.writeVTK(filename.str());
        }
    }
}

// In PathwaySet.cpp - update writeCombinedVTK:
void PathwaySet::writeCombinedVTK(const std::string& filename) const
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    // Count total points
    size_t total_points = 0;
    for (const auto& p : pathways_) {
        total_points += p.size();
    }

    // Write VTK header
    file << "# vtk DataFile Version 3.0\n";
    file << "All pathways\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";
    file << "POINTS " << total_points << " double\n";

    // Write all points
    for (const auto& p : pathways_) {
        for (const auto& particle : p) {
            file << particle.x() << " " << particle.y() << " 0.0\n";
        }
    }

    // Write line connectivity
    file << "\nLINES " << pathways_.size() << " "
         << (total_points + pathways_.size()) << "\n";

    size_t point_offset = 0;
    for (const auto& p : pathways_) {
        file << p.size();
        for (size_t i = 0; i < p.size(); ++i) {
            file << " " << (point_offset + i);
        }
        file << "\n";
        point_offset += p.size();
    }

    // Write pathway IDs as cell data
    file << "\nCELL_DATA " << pathways_.size() << "\n";
    file << "SCALARS pathway_id int\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& p : pathways_) {
        file << p.id() << "\n";
    }

    // Add point data for velocities and time
    file << "\nPOINT_DATA " << total_points << "\n";

    // Write time
    file << "SCALARS time double\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& p : pathways_) {
        for (const auto& particle : p) {
            file << particle.t() << "\n";
        }
    }

    // Write velocity vectors
    file << "\nVECTORS velocity double\n";
    for (const auto& p : pathways_) {
        for (const auto& particle : p) {
            file << particle.qx() << " " << particle.qy() << " 0.0\n";
        }
    }

    file.close();
}
double PathwaySet::meanPathLength() const
{
    if (pathways_.empty()) return 0.0;

    double sum = 0.0;
    size_t count = 0;

    for (const auto& p : pathways_) {
        if (!p.empty()) {
            sum += p.pathLength();
            ++count;
        }
    }

    return count > 0 ? sum / count : 0.0;
}

double PathwaySet::meanTravelTime() const
{
    if (pathways_.empty()) return 0.0;

    double sum = 0.0;
    size_t count = 0;

    for (const auto& p : pathways_) {
        if (!p.empty()) {
            sum += p.duration();
            ++count;
        }
    }

    return count > 0 ? sum / count : 0.0;
}

std::pair<double, double> PathwaySet::pathLengthRange() const
{
    if (pathways_.empty()) return {0.0, 0.0};

    double min_len = std::numeric_limits<double>::max();
    double max_len = 0.0;

    for (const auto& p : pathways_) {
        if (!p.empty()) {
            double len = p.pathLength();
            min_len = std::min(min_len, len);
            max_len = std::max(max_len, len);
        }
    }

    return {min_len, max_len};
}

std::pair<double, double> PathwaySet::travelTimeRange() const
{
    if (pathways_.empty()) return {0.0, 0.0};

    double min_time = std::numeric_limits<double>::max();
    double max_time = 0.0;

    for (const auto& p : pathways_) {
        if (!p.empty()) {
            double time = p.duration();
            min_time = std::min(min_time, time);
            max_time = std::max(max_time, time);
        }
    }

    return {min_time, max_time};
}
