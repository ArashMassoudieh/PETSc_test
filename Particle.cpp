// Particle.cpp
#include "Particle.h"

Particle::Particle()
    : x_(0.0), y_(0.0), active_(true), id_(-1), t_(0.0), qx_(0.0), qy_(0.0)
{
}

Particle::Particle(double x, double y, int id)
    : x_(x), y_(y), active_(true), id_(id), t_(0.0), qx_(0.0), qy_(0.0)
{
}

void Particle::setPosition(double x, double y) {
    x_ = x;
    y_ = y;
}

void Particle::move(double dx, double dy) {
    x_ += dx;
    y_ += dy;
}

// Operator overloads
Particle Particle::operator+(const Particle& other) const {
    Particle result;
    result.x_ = x_ + other.x_;
    result.y_ = y_ + other.y_;
    result.t_ = t_ + other.t_;
    result.qx_ = qx_ + other.qx_;
    result.qy_ = qy_ + other.qy_;
    result.active_ = active_;  // Keep the state of first particle
    result.id_ = id_;          // Keep the ID of first particle
    return result;
}

Particle Particle::operator-(const Particle& other) const {
    Particle result;
    result.x_ = x_ - other.x_;
    result.y_ = y_ - other.y_;
    result.t_ = t_ - other.t_;
    result.qx_ = qx_ - other.qx_;
    result.qy_ = qy_ - other.qy_;
    result.active_ = active_;  // Keep the state of first particle
    result.id_ = id_;          // Keep the ID of first particle
    return result;
}

Particle Particle::operator*(double scalar) const {
    Particle result;
    result.x_ = x_;              // Position unchanged
    result.y_ = y_;              // Position unchanged
    result.t_ = t_ * scalar;     // Time scaled
    result.qx_ = qx_ * scalar;   // Velocity scaled
    result.qy_ = qy_ * scalar;   // Velocity scaled
    result.active_ = active_;    // Keep the state
    result.id_ = id_;            // Keep the ID
    return result;
}

// Friend function for scalar * Particle
Particle operator*(double scalar, const Particle& p) {
    return p * scalar;
}

