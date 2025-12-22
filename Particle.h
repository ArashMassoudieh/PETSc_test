// Particle.h
#ifndef PARTICLE_H
#define PARTICLE_H

class Particle {
public:
    // Constructors
    Particle();
    Particle(double x, double y, int id = -1);

    // Getters
    double x() const { return x_; }
    double y() const { return y_; }
    bool isActive() const { return active_; }
    int id() const { return id_; }
    double qx() const { return qx_; }
    double qy() const { return qy_; }

    // Setters
    void setPosition(double x, double y);
    void setX(double x) { x_ = x; }
    void setY(double y) { y_ = y; }
    void setActive(bool active) { active_ = active; }
    void setId(int id) { id_ = id; }
    void setQx(double qx) { qx_ = qx; }
    void setQy(double qy) { qy_ = qy; }
    void setVelocity(double qx, double qy) { qx_ = qx; qy_ = qy; }

    // Movement
    void move(double dx, double dy);

    double t() const { return t_; }

    // Existing setters...
    void setT(double t) { t_ = t; }

private:
    double x_;          // x-coordinate
    double y_;          // y-coordinate
    double t_;          // time
    bool active_;       // Whether particle is still being tracked
    int id_;            // Particle identifier
    double qx_;         // v_x
    double qy_;         // v_y
};

#endif // PARTICLE_H
