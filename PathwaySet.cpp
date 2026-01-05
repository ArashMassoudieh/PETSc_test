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

PathwaySet PathwaySet::sampleParticlePairs(double Delta_x, size_t num_samples) const
{
    if (pathways_.empty()) {
        throw std::runtime_error("sampleParticlePairs: PathwaySet is empty");
    }

    // Filter pathways that are long enough to accommodate Delta_x
    std::vector<size_t> valid_pathway_indices;
    for (size_t i = 0; i < pathways_.size(); ++i) {
        if (!pathways_[i].empty()) {
            double path_length = pathways_[i].last().x() - pathways_[i].first().x();
            if (path_length >= Delta_x) {
                valid_pathway_indices.push_back(i);
            }
        }
    }

    if (valid_pathway_indices.empty()) {
        throw std::runtime_error("sampleParticlePairs: no pathways long enough for Delta_x");
    }

    // Create result PathwaySet with 2 pathways (one for x, one for x+Delta_x)
    PathwaySet result(2);
    Pathway pathway_at_x(0);      // Particles at location x
    Pathway pathway_at_x_plus(1); // Particles at location x + Delta_x

    // Reserve space
    pathway_at_x.reserve(num_samples);
    pathway_at_x_plus.reserve(num_samples);

    // Random number generators
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> pathway_dist(0, valid_pathway_indices.size() - 1);

    // Sample particle pairs
    for (size_t sample = 0; sample < num_samples; ++sample) {
        // Randomly select a valid pathway
        size_t random_pathway_idx = valid_pathway_indices[pathway_dist(gen)];
        const Pathway& selected_pathway = pathways_[random_pathway_idx];

        try {
            // Extract a random pair from this pathway
            auto pair = selected_pathway.extractRandomPairWithSeparation(Delta_x);

            // Add to respective pathways
            pathway_at_x.addParticle(pair.first);
            pathway_at_x_plus.addParticle(pair.second);

        } catch (const std::exception& e) {
            std::cerr << "Warning: Failed to sample from pathway "
                      << random_pathway_idx << ": " << e.what() << std::endl;
            // Continue to next sample
        }
    }

    // Add the two pathways to result
    result.addPathway(pathway_at_x);
    result.addPathway(pathway_at_x_plus);

    std::cout << "Sampled " << pathway_at_x.size() << " particle pairs with separation Delta_x = "
              << Delta_x << std::endl;

    return result;
}

double PathwaySet::calculateCorrelation(size_t pathway1_idx, size_t pathway2_idx,
                                        const std::string& quantity) const
{
    if (pathway1_idx >= pathways_.size() || pathway2_idx >= pathways_.size()) {
        throw std::runtime_error("calculateCorrelation: pathway index out of range");
    }

    const Pathway& path1 = pathways_[pathway1_idx];
    const Pathway& path2 = pathways_[pathway2_idx];

    if (path1.size() != path2.size()) {
        throw std::runtime_error("calculateCorrelation: pathways must have equal number of particles");
    }

    if (path1.empty()) {
        throw std::runtime_error("calculateCorrelation: pathways are empty");
    }

    size_t n = path1.size();

    // Extract values based on quantity
    std::vector<double> values1(n), values2(n);
    for (size_t i = 0; i < n; ++i) {
        if (quantity == "qx") {
            values1[i] = path1[i].qx();
            values2[i] = path2[i].qx();
        } else if (quantity == "qy") {
            values1[i] = path1[i].qy();
            values2[i] = path2[i].qy();
        } else {
            throw std::runtime_error("calculateCorrelation: unknown quantity '" + quantity + "'");
        }
    }

    // Calculate means
    double mean1 = 0.0, mean2 = 0.0;
    for (size_t i = 0; i < n; ++i) {
        mean1 += values1[i];
        mean2 += values2[i];
    }
    mean1 /= n;
    mean2 /= n;

    // Calculate correlation coefficient
    double numerator = 0.0;
    double variance1 = 0.0;
    double variance2 = 0.0;

    for (size_t i = 0; i < n; ++i) {
        double dev1 = values1[i] - mean1;
        double dev2 = values2[i] - mean2;
        numerator += dev1 * dev2;
        variance1 += dev1 * dev1;
        variance2 += dev2 * dev2;
    }

    // Check for zero variance
    if (variance1 < 1e-15 || variance2 < 1e-15) {
        std::cerr << "Warning: calculateCorrelation has near-zero variance" << std::endl;
        return 0.0;
    }

    double correlation = numerator / std::sqrt(variance1 * variance2);

    return correlation;
}
