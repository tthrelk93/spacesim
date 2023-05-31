#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include <cmath>
#include <QTimer>
#include <QElapsedTimer>


// Gravitational constant
const double G = 6.67430e-11;

// Define a 3D point
struct Point3D {
    double x, y, z;
};

// Define a 3D object with position and radius
struct Object3D {
    Point3D position;
    double radius;
};

// Define an Octree node
class OctreeNode {
public:
    Point3D center;
    double size;
    std::vector<Object3D> objects;
    OctreeNode* children[8] = {nullptr};

    OctreeNode(const Point3D& center, double size);
    ~OctreeNode();
    bool fits(const Object3D& object) const;
    void insert(const Object3D& object);
    void subdivide();
    bool checkCollisions() const;
};

class Spacecraft {
public:
    double mass;
    double radius;  // Add this line
    double x, y, z;  // Position
    double vx, vy, vz;  // Velocity

    Spacecraft(double mass, double radius, double x, double y, double z, double vx, double vy, double vz);
    void updatePosition(double dt);
    void updateVelocity(double fx, double fy, double fz, double dt, QElapsedTimer& debugTimer);
};

class CelestialBody {
public:
    double mass;
    double radius;
    double x, y, z;
    double vx, vy, vz;

    CelestialBody(double mass, double radius, double x, double y, double z, double vx, double vy, double vz);
    void updatePosition(double dt);
    void updateVelocity(double fx, double fy, double fz, double dt);
};


class Universe {
public:
    std::vector<Spacecraft> spacecrafts;
    std::vector<CelestialBody> bodies;
    OctreeNode octree;  // The universe as an Octree

    Universe();
    void addSpacecraft(const Spacecraft &s);
    void removeSpacecraft(int index);
    void addBody(const CelestialBody &b);
    void removeBody(int index);
    bool checkCollisions();
    void update(double dt, QElapsedTimer& debugTimer);
    void fireEngine(int spacecraftIndex, double fx, double fy, double fz, double dt);
};

std::vector<Universe> runSimulation(int numSteps);


#endif // SIMULATION_H
