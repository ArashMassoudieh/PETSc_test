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
