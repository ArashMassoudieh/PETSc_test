// Pathway.h
#ifndef PATHWAY_H
#define PATHWAY_H

#include "Particle.h"
#include <vector>
#include <string>
#include <random>

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
    std::mt19937_64 rng;
};

class Grid2D;

class Pathway {
public:
    // Constructors
    Pathway();
    explicit Pathway(int id);

    // Factory method for reserving size
    static Pathway withReservedSize(int reserve_size);
    static Pathway withReservedSize(int reserve_size, int id);

    // Add particle to pathway
    void addParticle(const Particle& p);
    void addParticle(double x, double y, double t);

    // Access particles
    Particle& operator[](size_t i);
    const Particle& operator[](size_t i) const;
    Particle& at(size_t i);
    const Particle& at(size_t i) const;

    // Iterators
    std::vector<Particle>::iterator begin();
    std::vector<Particle>::iterator end();
    std::vector<Particle>::const_iterator begin() const;
    std::vector<Particle>::const_iterator end() const;

    // Size
    size_t size() const;
    bool empty() const;
    void reserve(size_t n);
    void clear();

    // Getters
    int id() const { return id_; }
    void setId(int id) { id_ = id; }

    // Get first and last particle
    const Particle& first() const;
    const Particle& last() const;
    Particle& first();
    Particle& last();

    // Get start and end positions
    std::pair<double, double> startPosition() const;
    std::pair<double, double> endPosition() const;

    // Get time range
    double startTime() const;
    double endTime() const;
    double duration() const;

    // Get total path length
    double pathLength() const;

    // File I/O
    void writeToFile(const std::string& filename) const;
    void readFromFile(const std::string& filename);
    void writeVTK(const std::string& filename) const;

    // Particle tracking
    void trackParticle(Grid2D* grid, double dx_step,
                       const std::string& qx_name = "qx",
                       const std::string& qy_name = "qy");

    double trackParticleDiffusion(const double &dt, const double &rx, const double &ry, const double &D);

    // Check if pathway has uniform x-spacing
    bool isUniformX(double tolerance = 1e-6) const;
    void cacheUniformXStatus(double tolerance = 1e-6);

    // Setters for uniform x status
    void setUniformX(bool is_uniform) { uniform_x_ = is_uniform; }
    void setUniformXChecked(bool checked) { uniform_x_checked_ = checked; }

    // Interpolate particle properties at given x-location
    Particle interpolateAtX(double target_x, double tolerance = 1e-10) const;

    // Extract two random particles separated by distance Delta_x
    std::pair<Particle, Particle> extractRandomPairWithSeparation(double Delta_x) const;

    void setWienerParams(WienerParams _wp) {wp= _wp;}
private:
    std::vector<Particle> particles_;
    int id_;
    bool uniform_x_;
    bool uniform_x_checked_;
    WienerParams wp;
};

#endif // PATHWAY_H
