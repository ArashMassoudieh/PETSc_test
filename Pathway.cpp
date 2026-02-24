// Pathway.cpp
#include "Pathway.h"
#include <fstream>
#include <stdexcept>
#include <cmath>
#include "grid.h"
#include <iostream>

Pathway::Pathway()
    : id_(-1), uniform_x_(false), uniform_x_checked_(false)
{
}

Pathway::Pathway(int id)
    : id_(id), uniform_x_(false), uniform_x_checked_(false)
{
}

Pathway Pathway::withReservedSize(int reserve_size)
{
    Pathway path;
    path.particles_.reserve(reserve_size);
    return path;
}

Pathway Pathway::withReservedSize(int reserve_size, int id)
{
    Pathway path(id);
    path.particles_.reserve(reserve_size);
    return path;
}

void Pathway::addParticle(const Particle& p)
{
    particles_.push_back(p);
}

void Pathway::addParticle(double x, double y, double t)
{
    Particle p(x, y, id_);
    p.setT(t);
    particles_.push_back(p);
}

Particle& Pathway::operator[](size_t i)
{
    return particles_[i];
}

const Particle& Pathway::operator[](size_t i) const
{
    return particles_[i];
}

Particle& Pathway::at(size_t i)
{
    return particles_.at(i);
}

const Particle& Pathway::at(size_t i) const
{
    return particles_.at(i);
}

std::vector<Particle>::iterator Pathway::begin()
{
    return particles_.begin();
}

std::vector<Particle>::iterator Pathway::end()
{
    return particles_.end();
}

std::vector<Particle>::const_iterator Pathway::begin() const
{
    return particles_.begin();
}

std::vector<Particle>::const_iterator Pathway::end() const
{
    return particles_.end();
}

size_t Pathway::size() const
{
    return particles_.size();
}

bool Pathway::empty() const
{
    return particles_.empty();
}

void Pathway::reserve(size_t n)
{
    particles_.reserve(n);
}

void Pathway::clear()
{
    particles_.clear();
}

const Particle& Pathway::first() const
{
    if (particles_.empty()) {
        throw std::runtime_error("Pathway::first() - pathway is empty");
    }
    return particles_.front();
}

const Particle& Pathway::last() const
{
    if (particles_.empty()) {
        throw std::runtime_error("Pathway::last() - pathway is empty");
    }
    return particles_.back();
}

Particle& Pathway::first()
{
    if (particles_.empty()) {
        throw std::runtime_error("Pathway::first() - pathway is empty");
    }
    return particles_.front();
}

Particle& Pathway::last()
{
    if (particles_.empty()) {
        throw std::runtime_error("Pathway::last() - pathway is empty");
    }
    return particles_.back();
}

std::pair<double, double> Pathway::startPosition() const
{
    if (particles_.empty()) {
        throw std::runtime_error("Pathway::startPosition() - pathway is empty");
    }
    return {particles_.front().x(), particles_.front().y()};
}

std::pair<double, double> Pathway::endPosition() const
{
    if (particles_.empty()) {
        throw std::runtime_error("Pathway::endPosition() - pathway is empty");
    }
    return {particles_.back().x(), particles_.back().y()};
}

double Pathway::startTime() const
{
    if (particles_.empty()) {
        throw std::runtime_error("Pathway::startTime() - pathway is empty");
    }
    return particles_.front().t();
}

double Pathway::endTime() const
{
    if (particles_.empty()) {
        throw std::runtime_error("Pathway::endTime() - pathway is empty");
    }
    return particles_.back().t();
}

double Pathway::duration() const
{
    if (particles_.empty()) {
        return 0.0;
    }
    return endTime() - startTime();
}

double Pathway::pathLength() const
{
    if (particles_.size() < 2) {
        return 0.0;
    }

    double length = 0.0;
    for (size_t i = 1; i < particles_.size(); ++i) {
        double dx = particles_[i].x() - particles_[i-1].x();
        double dy = particles_[i].y() - particles_[i-1].y();
        length += std::sqrt(dx*dx + dy*dy);
    }

    return length;
}

void Pathway::writeToFile(const std::string& filename) const
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    file << "# pathway_id: " << id_ << "\n";
    file << "# t x y qx qy\n";  // Updated header
    for (const auto& p : particles_) {
        file << p.t() << " " << p.x() << " " << p.y()
        << " " << p.qx() << " " << p.qy() << "\n";  // Added qx, qy
    }
    file.close();
}

void Pathway::readFromFile(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    particles_.clear();
    std::string line;

    // Read header to get pathway ID
    std::getline(file, line);
    if (line.find("pathway_id:") != std::string::npos) {
        size_t pos = line.find(":");
        if (pos != std::string::npos) {
            id_ = std::stoi(line.substr(pos + 1));
        }
    }

    // Skip second header line if present
    std::getline(file, line);
    if (line[0] != '#') {
        file.seekg(0);
    }

    double t, x, y;
    while (file >> t >> x >> y) {
        addParticle(x, y, t);
    }

    file.close();
}

void Pathway::writeVTK(const std::string& filename) const
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    // Write VTK header for polyline
    file << "# vtk DataFile Version 3.0\n";
    file << "Particle pathway " << id_ << "\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";
    file << "POINTS " << particles_.size() << " double\n";

    for (const auto& p : particles_) {
        file << p.x() << " " << p.y() << " 0.0\n";
    }

    // Define the line connectivity
    file << "\nLINES 1 " << (particles_.size() + 1) << "\n";
    file << particles_.size();
    for (size_t i = 0; i < particles_.size(); ++i) {
        file << " " << i;
    }
    file << "\n";

    // Add time as scalar data
    file << "\nPOINT_DATA " << particles_.size() << "\n";
    file << "SCALARS time double\n";
    file << "LOOKUP_TABLE default\n";
    for (const auto& p : particles_) {
        file << p.t() << "\n";
    }

    file << "\nVECTORS velocity double\n";
    for (const auto& p : particles_) {
        file << p.qx() << " " << p.qy() << " 0.0\n";
    }

    file.close();
}


double Pathway::trackParticleDiffusion(const double &dt, const double &rx, const double &ry, const double &D, std::mt19937_64 &rng)
{
    if (particles_.empty()) {
        throw std::runtime_error("trackParticleDiffusion: pathway has no starting particle");
    }

    Particle myParticle = particles_[0];
    std::normal_distribution<double> N01(0.0, 1.0);
    double t = 0;

    while (myParticle.x() * myParticle.x() / (rx * rx) +
               myParticle.y() * myParticle.y() / (ry * ry) < 1.0)
    {
        t += dt;
        const double dWx = N01(rng);
        const double dWy = N01(rng);
        myParticle.setX(myParticle.x() + dWx * std::sqrt(2.0 * D * dt));
        myParticle.setY(myParticle.y() + dWy * std::sqrt(2.0 * D * dt));
        myParticle.setT(t);
        this->addParticle(myParticle);
    }
    return myParticle.t();
}

void Pathway::trackParticle(Grid2D* grid, double dx_step,
                            const std::string& qx_name,
                            const std::string& qy_name)
{
    if (!grid) {
        throw std::runtime_error("trackParticle: grid pointer is null");
    }

    if (particles_.empty()) {
        throw std::runtime_error("trackParticle: pathway has no starting particle");
    }

    if (dx_step <= 0.0) {
        throw std::runtime_error("trackParticle: dx_step must be positive");
    }

    // Get domain bounds
    const double Lx = grid->Lx();
    const double Ly = grid->Ly();

    // Start from the last particle in the pathway
    Particle current = particles_.back();
    double t = current.t();

    const int max_steps = 1000000;  // Safety limit
    int step_count = 0;

    // Track until particle exits right boundary
    while (current.x() < Lx && step_count < max_steps) {

        // Get velocity at current position
        double vx, vy;
        try {
            auto vel = grid->getVelocityAt(current.x(), current.y(), qx_name, qy_name);
            vx = vel.first;
            vy = vel.second;
        } catch (const std::exception& e) {
            // Particle left domain (top, bottom boundaries)
            current.setActive(false);
            particles_.push_back(current);
            break;
        }

        // Calculate velocity magnitude
        double v_mag = std::sqrt(vx*vx + vy*vy);

        if (v_mag < 1e-12) {
            // Velocity is essentially zero - particle is stuck
            current.setActive(false);
            particles_.push_back(current);
            break;
        }

        // Calculate dt needed to move distance dx_step
        double dt = dx_step / v_mag;

        // Update position
        double new_x = current.x() + vx * dt;
        double new_y = current.y() + vy * dt;

        // Check if particle exits right boundary
        if (new_x >= Lx) {
            // Particle exits - interpolate to exact exit point
            double t_exit = (Lx - current.x()) / vx;
            new_x = Lx;
            new_y = current.y() + vy * t_exit;
            t += t_exit;

            current.setPosition(new_x, new_y);
            current.setT(t);
            current.setVelocity(vx, vy);  // Set velocity BEFORE pushing
            current.setActive(true);  // Successfully exited
            particles_.push_back(current);
            break;
        }

        // Check if particle exits top or bottom boundary
        if (new_y < 0.0 || new_y > Ly) {
            current.setActive(false);
            particles_.push_back(current);
            break;
        }

        // Update particle - do ALL updates before pushing
        t += dt;
        current.setPosition(new_x, new_y);
        current.setT(t);
        current.setVelocity(vx, vy);  // MOVED: Set velocity after position update
        particles_.push_back(current);

        step_count++;
    }

    if (step_count >= max_steps) {
        std::cerr << "Warning: Particle tracking reached maximum steps ("
                  << max_steps << ") for pathway " << id_ << std::endl;
        current.setActive(false);
    }
    setUniformX(true);
}

void Pathway::trackParticleWithDiffusion(Grid2D* grid, double dx_step,
                                         double D, gsl_rng* rng,
                                         const std::string& qx_name,
                                         const std::string& qy_name)
{
    if (!grid)
        throw std::runtime_error("trackParticleWithDiffusion: grid pointer is null");
    if (particles_.empty())
        throw std::runtime_error("trackParticleWithDiffusion: pathway has no starting particle");
    if (dx_step <= 0.0)
        throw std::runtime_error("trackParticleWithDiffusion: dx_step must be positive");

    const double Lx = grid->Lx();
    const double Ly = grid->Ly();
    const int nx = grid->nx();
    const int ny = grid->ny();
    const double gdx = grid->dx();
    const double gdy = grid->dy();

    const bool has_Dy = grid->hasFlux("D_y");

    // Cap: diffusive jump <= cap_factor * dy
    const double cap_factor = 2.0;

    Particle current = particles_.back();
    double t = current.t();

    const int max_steps = 1000000;
    int step_count = 0;

    while (current.x() < Lx && step_count < max_steps) {

        // --- 1. Get velocity at current position ---
        double vx, vy;
        try {
            auto vel = grid->getVelocityAt(current.x(), current.y(), qx_name, qy_name);
            vx = vel.first;
            vy = vel.second;
        } catch (const std::exception& e) {
            current.setActive(false);
            particles_.push_back(current);
            break;
        }

        double v_mag = std::sqrt(vx*vx + vy*vy);

        if (v_mag < 1e-12) {
            current.setActive(false);
            particles_.push_back(current);
            break;
        }

        // --- 2. Get local D_y AND its y-gradient (BEFORE dt) ---
        double D_y_local = D;
        double dDy_dy = 0.0;   // Ito drift correction: d(D_y)/dy

        if (has_Dy) {
            const std::vector<double>& D_y = grid->flux("D_y");

            double px = current.x();
            double py = current.y();

            int i = static_cast<int>(px / gdx);
            int j = static_cast<int>(py / gdy);
            i = std::min(i, nx - 1);
            j = std::min(j, ny - 1);

            double xi  = (px - i * gdx) / gdx;
            double eta = (py - j * gdy) / gdy;
            xi  = std::max(0.0, std::min(1.0, xi));
            eta = std::max(0.0, std::min(1.0, eta));

            int idx_bottom = j * nx + i;
            int idx_top    = (j + 1) * nx + i;

            // --- Interpolate D_y at particle position ---
            if (i < nx - 1) {
                int idx_bottom_right = j * nx + (i + 1);
                int idx_top_right    = (j + 1) * nx + (i + 1);

                double Dy_left  = (1.0 - eta) * D_y[idx_bottom]       + eta * D_y[idx_top];
                double Dy_right = (1.0 - eta) * D_y[idx_bottom_right] + eta * D_y[idx_top_right];
                D_y_local = (1.0 - xi) * Dy_left + xi * Dy_right;
            } else {
                D_y_local = (1.0 - eta) * D_y[idx_bottom] + eta * D_y[idx_top];
            }
            D_y_local = std::max(0.0, D_y_local);

            // --- Compute dD_y/dy via finite difference on face values ---
            // D_y at face j and face j+1 are the values bounding this cell.
            // Gradient within cell (i,j): (D_y[top] - D_y[bottom]) / dy
            // Interpolate in x for consistency.
            double grad_left = (D_y[idx_top] - D_y[idx_bottom]) / gdy;
            if (i < nx - 1) {
                int idx_bottom_right = j * nx + (i + 1);
                int idx_top_right    = (j + 1) * nx + (i + 1);
                double grad_right = (D_y[idx_top_right] - D_y[idx_bottom_right]) / gdy;
                dDy_dy = (1.0 - xi) * grad_left + xi * grad_right;
            } else {
                dDy_dy = grad_left;
            }
        }

        // --- 3. Adaptive dt: cap based on local diffusion ---
        double dt_adv = dx_step / v_mag;

        double D_eff = std::max(D, D_y_local);
        double dt_cap = (D_eff > 1e-15)
                            ? (cap_factor * cap_factor * gdy * gdy) / (2.0 * D_eff)
                            : 1e10;

        double dt = std::min(dt_adv, dt_cap);

        // --- 4. Compute displacements ---
        double adv_dx = vx * dt;
        double adv_dy = vy * dt;

        // Longitudinal diffusion (constant D)
        double diff_dx = std::sqrt(2.0 * D * dt) * gsl_ran_gaussian_ziggurat(rng, 1.0);

        // Transverse: Ito drift correction + diffusive noise
        // SDE: dy = (dD_y/dy) dt + sqrt(2 D_y) dW
        double drift_dy = dDy_dy * dt;
        double noise_dy = std::sqrt(2.0 * D_y_local * dt) * gsl_ran_gaussian_ziggurat(rng, 1.0);
        double diff_dy = drift_dy + noise_dy;

        double new_x = current.x() + adv_dx + diff_dx;
        double new_y = current.y() + adv_dy + diff_dy;

        // --- 5. Boundary conditions ---

        // Right boundary exit
        if (new_x >= Lx) {
            double t_exit = (Lx - current.x()) / vx;
            new_x = Lx;
            new_y = current.y() + vy * t_exit;
            t += t_exit;

            current.setPosition(new_x, new_y);
            current.setT(t);
            current.setVelocity(vx, vy);
            current.setActive(true);
            particles_.push_back(current);
            break;
        }

        // Left boundary: reflect
        if (new_x < 0.0) {
            new_x = -new_x;
        }

        // Top/bottom boundaries: reflect
        if (new_y < 0.0) {
            new_y = -new_y;
        } else if (new_y > Ly) {
            new_y = 2.0 * Ly - new_y;
        }

        t += dt;
        current.setPosition(new_x, new_y);
        current.setT(t);
        current.setVelocity(vx, vy);
        particles_.push_back(current);

        step_count++;
    }

    if (step_count >= max_steps) {
        std::cerr << "Warning: Particle tracking reached maximum steps ("
                  << max_steps << ") for pathway " << id_ << std::endl;
        current.setActive(false);
    }
    setUniformX(false);
    setUniformXChecked(true);
}



bool Pathway::isUniformX(double tolerance) const
{
    if (particles_.size() < 3) {
        return true;  // Too few points to be non-uniform
    }

    // Check if we've already computed this
    if (uniform_x_checked_) {
        return uniform_x_;
    }

    // Sample only first 10 increments (or fewer if pathway is short)
    size_t n_samples = std::min(size_t(10), particles_.size() - 1);

    // Get first increment as reference
    double dx_ref = particles_[1].x() - particles_[0].x();

    // Check if all sampled increments match the reference
    for (size_t i = 1; i < n_samples; ++i) {
        double dx = particles_[i+1].x() - particles_[i].x();
        if (std::abs(dx - dx_ref) > tolerance) {
            return false;
        }
    }

    return true;
}

void Pathway::cacheUniformXStatus(double tolerance)
{
    uniform_x_ = isUniformX(tolerance);
    uniform_x_checked_ = true;
}

Particle Pathway::interpolateAtX(double target_x, double tolerance) const
{
    if (particles_.empty()) {
        throw std::runtime_error("interpolateAtX: pathway is empty");
    }

    if (particles_.size() == 1) {
        return particles_[0];
    }

    // Check if target_x is within pathway bounds
    double x_min = particles_.front().x();
    double x_max = particles_.back().x();

    if (target_x < x_min - tolerance || target_x > x_max + tolerance) {
        throw std::runtime_error("interpolateAtX: target_x outside pathway bounds");
    }

    // Clamp to bounds if very close
    if (target_x < x_min) target_x = x_min;
    if (target_x > x_max) target_x = x_max;

    size_t i1, i2;

    if (uniform_x_ && uniform_x_checked_) {
        // Fast path: use direct indexing for uniform spacing
        double dx = particles_[1].x() - particles_[0].x();

        if (std::abs(dx) < tolerance) {
            // All particles at same x (degenerate case)
            return particles_[0];
        }

        // Calculate which segment contains target_x
        double frac_index = (target_x - x_min) / dx;
        i1 = static_cast<size_t>(std::floor(frac_index));

        // Clamp to valid range
        if (i1 >= particles_.size() - 1) {
            i1 = particles_.size() - 2;
        }
        i2 = i1 + 1;
    }
    else {
        // Slow path: binary search for non-uniform spacing
        // Find two consecutive particles bracketing target_x
        bool found = false;
        for (size_t i = 0; i < particles_.size() - 1; ++i) {
            if (particles_[i].x() <= target_x && target_x <= particles_[i+1].x()) {
                i1 = i;
                i2 = i + 1;
                found = true;
                break;
            }
        }

        if (!found) {
            throw std::runtime_error("interpolateAtX: failed to find bracketing particles");
        }
    }

    const Particle& p1 = particles_[i1];
    const Particle& p2 = particles_[i2];

    // Calculate interpolation weight
    double dx = p2.x() - p1.x();

    if (std::abs(dx) < tolerance) {
        // Particles at essentially same x-location
        return p1;
    }

    double weight = (target_x - p1.x()) / dx;

    // Interpolate using particle operators
    return p1 + (p2 - p1) * weight;
}

std::pair<Particle, Particle> Pathway::extractRandomPairWithSeparation(double Delta_x) const
{
    if (particles_.empty()) {
        throw std::runtime_error("extractRandomPairWithSeparation: pathway is empty");
    }

    double x_min = particles_.front().x();
    double x_max = particles_.back().x();
    double path_length = x_max - x_min;

    if (Delta_x > path_length) {
        throw std::runtime_error("extractRandomPairWithSeparation: Delta_x exceeds pathway length");
    }

    if (Delta_x < 0.0) {
        throw std::runtime_error("extractRandomPairWithSeparation: Delta_x must be positive");
    }

    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(x_min, x_max - Delta_x);

    // Generate random starting x position
    double x1 = dist(gen);
    double x2 = x1 + Delta_x;

    // Interpolate particles at these locations
    Particle p1 = interpolateAtX(x1);
    Particle p2 = interpolateAtX(x2);

    return {p1, p2};
}
