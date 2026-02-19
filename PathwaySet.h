// PathwaySet.h
#ifndef PATHWAYSET_H
#define PATHWAYSET_H

#include "Pathway.h"
#include <vector>
#include <string>
#include <utility>
#include "TimeSeries.h"
#include <gsl/gsl_rng.h>

class Grid2D;  // Forward declaration

class PathwaySet {
public:
    // Weighting enum for initialization
    enum class Weighting {
        Uniform,
        FluxWeighted
    };

    // Constructors
    PathwaySet();
    explicit PathwaySet(int reserve_size);

    // Add pathways
    void addPathway(const Pathway& p);
    void addPathway(int id);

    // Initialize pathways at left boundary
    void Initialize(size_t number_of_pathways, Weighting weighting, Grid2D* grid);
    void InitializeAtOrigin(size_t number_of_pathways);
    // Particle tracking
    void trackAllPathways(Grid2D* grid, double dx_step,
                          const std::string& qx_name = "qx",
                          const std::string& qy_name = "qy");

    void trackAllPathwaysWithDiffusion(Grid2D* grid, double dx_step,
                                       double D, unsigned long seed = 0,
                                       const std::string& qx_name = "qx",
                                       const std::string& qy_name = "qy");

    // Access pathways
    Pathway& operator[](size_t i);
    const Pathway& operator[](size_t i) const;
    Pathway& at(size_t i);
    const Pathway& at(size_t i) const;

    // Iterators
    std::vector<Pathway>::iterator begin();
    std::vector<Pathway>::iterator end();
    std::vector<Pathway>::const_iterator begin() const;
    std::vector<Pathway>::const_iterator end() const;

    // Size and capacity
    size_t size() const;
    bool empty() const;
    void reserve(size_t n);
    void clear();

    // Pathway management
    size_t countActive() const;
    size_t countCompleted() const;  // Pathways that exited right boundary
    void removeIncomplete();  // Remove pathways that didn't complete

    // File I/O
    void writeToFile(const std::string& filename) const;
    void readFromFile(const std::string& filename);
    void writeAllVTK(const std::string& prefix) const;  // Write each pathway as separate VTK
    void writeCombinedVTK(const std::string& filename) const;  // All pathways in one file

    // Statistics
    double meanPathLength() const;
    double meanTravelTime() const;
    std::pair<double, double> pathLengthRange() const;  // min, max
    std::pair<double, double> travelTimeRange() const;  // min, max

    // Sample particle pairs for correlation analysis
    PathwaySet sampleParticlePairs(double Delta_x, size_t num_samples) const;

    // Calculate correlation between corresponding particles in two pathways
    double calculateCorrelation(size_t pathway1_idx, size_t pathway2_idx,
                                const std::string& quantity = "qx") const;

    TimeSeries<double> trackDiffusion(double dt, const double &rx, const double &ry, const double &D);

private:
    std::vector<Pathway> pathways_;
};

#endif // PATHWAYSET_H
