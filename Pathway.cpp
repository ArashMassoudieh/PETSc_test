// Pathway.cpp
#include "Pathway.h"
#include <fstream>
#include <stdexcept>
#include <cmath>
#include "grid.h"
#include <iostream>

Pathway::Pathway()
    : id_(-1)
{
}

Pathway::Pathway(int id)
    : id_(id)
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

// Add to Pathway.cpp

// In Pathway.cpp - update trackParticle function:

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
}
