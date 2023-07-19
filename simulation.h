#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include <cmath>
#include <QTimer>
#include <QElapsedTimer>
#include <unordered_map>
#include <string>
#include <QQueue>
#include <QVector3D>

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
    double renderSize;

    Spacecraft(double mass, double radius, double x, double y, double z, double vx, double vy, double vz);
    void updatePosition(double dt);
    void updateVelocity(double fx, double fy, double fz, double dt, QElapsedTimer& debugTimer);
};

class CelestialBody {
public:
    std::string name;
    double mass;
    double radius;
    QVector3D position;
    QVector3D velocity;   
    double objectType = 0.0;
    double renderSize;
    bool isMoon;
    QQueue<QVector3D> pastPositions;
    QVector3D force;
    //std::vector<CelestialBody> moons;
    // Default constructor
    CelestialBody() : name(""), mass(0), radius(0), position(QVector3D()), velocity(QVector3D()), objectType(0), renderSize(0), isMoon(false)  {}
    CelestialBody(std::string name, double mass, double radius, QVector3D position, QVector3D velocity, double objectType, double renderSize, bool isMoon);
      
    void updatePosition(double dt);

    void updateVelocity(QVector3D force, double dt);
};



class Universe {
public:
    std::vector<Spacecraft> spacecrafts;
    std::vector<CelestialBody> bodies;
    OctreeNode octree;  // The universe as an Octree
     std::unordered_map<std::string, CelestialBody> bodyMap;  // Hash map for celestial bodies

    Universe();
     QVector3D calculateForce(const CelestialBody &body, float dt = 0.0, QVector3D dv = QVector3D(0, 0, 0));
    void calculateAndApplyForces();
    void addSpacecraft(const Spacecraft &s);
    void removeSpacecraft(int index);
    void addBody(const CelestialBody &b);
    void removeBody(int index);
    bool checkCollisions();
    void computeVelocity(float dt);
    void updateLocation(float dt);
    QVector3D partial_step(QVector3D &f, QVector3D &df, double scale,float dt);
    QVector3D calculate_single_body_acceleration(float dt, CelestialBody &target, bool isRungeKutta);
    void update(float dt, QElapsedTimer& debugTimer, bool isRungeKutta);
    void fireEngine(int spacecraftIndex, double fx, double fy, double fz, double dt);
};

std::vector<Universe> runSimulation(int numSteps);


#endif // SIMULATION_H
