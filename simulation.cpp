#include "simulation.h"
#include <QApplication>
#include <QOpenGLWindow>
#include <QOpenGLFunctions>
#include <iostream>
#include <QTimer>
#include <QElapsedTimer>

const bool debug = true;
OctreeNode::OctreeNode(const Point3D& center, double size)
    : center(center), size(size) {
    if (debug) {
        std::cout << "Created OctreeNode at (" << center.x << ", " << center.y << ", " << center.z << ") with size " << size << std::endl;
    }
}

OctreeNode::~OctreeNode() {
    for (int i = 0; i < 8; ++i) {
        delete children[i];
    }
}


bool OctreeNode::fits(const Object3D& object) const {
    double halfSize = size / 2;

    // Check if object fits within node boundaries
    bool fitsInNode = object.position.x - object.radius >= center.x - halfSize &&
                      object.position.x + object.radius <= center.x + halfSize &&
                      object.position.y - object.radius >= center.y - halfSize &&
                      object.position.y + object.radius <= center.y + halfSize &&
                      object.position.z - object.radius >= center.z - halfSize &&
                      object.position.z + object.radius <= center.z + halfSize;
    // Check if object is small enough relative to node size
    bool isSmallEnough = object.radius <= size;  // Adjust this fraction as needed

    if (debug) {
        std::cout << "Checking if object at (" << object.position.x << ", " << object.position.y << ", " << object.position.z << ") with radius " << object.radius
                  << " fits in node at (" << center.x << ", " << center.y << ", " << center.z << ") with size " << size << " fitsInNode: " << fitsInNode << ", " << " isSmallEnough: " << isSmallEnough << std::endl;
    }

    return fitsInNode && isSmallEnough;
}


// Insert an object into the Octree
void OctreeNode::insert(const Object3D& object) {

    // If this node doesn't have children, and it's either empty or the object fits, store the object here
    if (children[0] == nullptr && (objects.empty() || fits(object))) {
        objects.push_back(object);
        return;
    }

    // If this node doesn't have children, but the object doesn't fit, subdivide
    if (children[0] == nullptr) {
        subdivide();
    }

    for (int i = 0; i < 8; ++i) {
        //if (children[i]->fits(object)) {
            std::cout << "Inserting object at (" << object.position.x << ", " << object.position.y << ", " << object.position.z << ") "
                      << "into child " << i << " of node at (" << center.x << ", " << center.y << ", " << center.z << ")" << std::endl;
            children[i]->insert(object);
           // break;
       // }
    }
}

// Subdivide this node into 8 children
void OctreeNode::subdivide() {
    std::cout << "here  ";
    for (int i = 0; i < 8; ++i) {
        double x = center.x + size / 4 * ((i & 1) ? 1 : -1);
        double y = center.y + size / 4 * ((i & 2) ? 1 : -1);
        double z = center.z + size / 4 * ((i & 4) ? 1 : -1);
        children[i] = new OctreeNode({x, y, z}, size / 2);
        if (debug) {
            std::cout << "Created child " << i << " of node at (" << center.x << ", " << center.y << ", " << center.z << ")"
                      << " with center at (" << x << ", " << y << ", " << z << ") and size " << size / 2 << std::endl;
        }
    }
    if (debug) {
        std::cout << "Subdividing node at (" << center.x << ", " << center.y << ", " << center.z << ")" << std::endl;
    }
    // Move any objects that fit into children into those children
    for (const Object3D& object : objects) {
        for (int i = 0; i < 8; ++i) {
            if (children[i]->fits(object)) {
                children[i]->insert(object);
            }

        }
    }
    objects.clear();
}

// Check for collisions
bool OctreeNode::checkCollisions() const {
    // If there are multiple objects in this node, check for collision
    if (objects.size() > 1) {
        std::cout << "Checking collisions in node at (" << center.x << ", " << center.y << ", " << center.z << ") with " << objects.size() << " objects" << std::endl;
        for (int i = 0; i < objects.size(); i++) {
            for (int j = i + 1; j < objects.size(); j++) {
                double dx = objects[i].position.x - objects[j].position.x;
                double dy = objects[i].position.y - objects[j].position.y;
                double dz = objects[i].position.z - objects[j].position.z;
                double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (distance < objects[i].radius + objects[j].radius) {
                    // The objects are closer than the sum of their radii - collision!
                    std::cout << "Collision detected between objects at (" << objects[i].position.x << ", " << objects[i].position.y << ", " << objects[i].position.z << ") and (" << objects[j].position.x << ", " << objects[j].position.y << ", " << objects[j].position.z << ")" << std::endl;
                    return true;
                }
            }
        }
    }

    // If there's a single object in this node, but it doesn't fit, there's a collision
    if (objects.size() == 1 && !fits(objects[0])) {
        std::cout << "Single object in node at (" << center.x << ", " << center.y << ", " << center.z << ") does not fit" << std::endl;
        return true;
    }

    // If there are no objects in this node, but there are children, check them for collisions
    if (children[0] != nullptr) {
        std::cout << "Checking collisions in children of node at (" << center.x << ", " << center.y << ", " << center.z << ")" << std::endl;
        for (int i = 0; i < 8; ++i) {
            if (children[i]->checkCollisions()) {
                return true;
            }
        }
    }
    // No collisions
    return false;
}




Spacecraft::Spacecraft(double mass, double radius, double x, double y, double z, double vx, double vy, double vz)
    : mass(mass), radius(radius), x(x), y(y), z(z), vx(vx), vy(vy), vz(vz) {}

// ... (implement the rest of the Spacecraft methods here) ...
// Update position based on velocity
void Spacecraft::updatePosition(double dt) {
    x += vx * dt;
    y += vy * dt;
    z += vz * dt;
}

// Update velocity based on force
void Spacecraft::updateVelocity(double fx, double fy, double fz, double dt, QElapsedTimer& debugTimer) {
    vx += fx / mass * dt;
    vy += fy / mass * dt;
    vz += fz / mass * dt;
    if(debug && debugTimer.elapsed() % 500 < 16){
        qDebug() << "craft velocity (vx, vy, vz):" << vx << vy << vz;
    }
}

CelestialBody::CelestialBody(double mass, double radius, double x, double y, double z, double vx, double vy, double vz)
    : mass(mass), radius(radius), x(x), y(y), z(z), vx(vx), vy(vy), vz(vz) {

}

void CelestialBody::updatePosition(double dt) {
    x += vx * dt;
    y += vy * dt;
    z += vz * dt;
}

void CelestialBody::updateVelocity(double fx, double fy, double fz, double dt) {
    vx += fx / mass * dt;
    vy += fy / mass * dt;
    vz += fz / mass * dt;
}


Universe::Universe()
    : octree({0, 0, 0}, 800) {  // The size of the Octree is now 1.0
    if (debug) {
        std::cout << "Created Universe with Octree root node at (0, 0, 0) and size 800 " << std::endl;
    }
}


// ... (implement the rest of the Universe methods here) ...
void Universe::addSpacecraft(const Spacecraft &s) {
    spacecrafts.push_back(s);
    octree.insert({{s.x, s.y, s.z}, s.radius});
}

// Remove a spacecraft
void Universe::removeSpacecraft(int index) {
    // You'll need to implement OctreeNode::remove to remove objects from the Octree
    //octree.remove({{spacecrafts[index].x, spacecrafts[index].y, spacecrafts[index].z}, spacecrafts[index].radius});
    spacecrafts.erase(spacecrafts.begin() + index);
}

// Add a celestial body
void Universe::addBody(const CelestialBody &b) {
    bodies.push_back(b);
    octree.insert({{b.x, b.y, b.z}, b.radius});
}

// Remove a celestial body
void Universe::removeBody(int index) {
    // You'll need to implement OctreeNode::remove to remove objects from the Octree
    //octree.remove({{bodies[index].x, bodies[index].y, bodies[index].z}, bodies[index].radius});
    bodies.erase(bodies.begin() + index);
}

// Check for collisions
bool Universe::checkCollisions() {
    return octree.checkCollisions();
}

void Universe::update(double dt, QElapsedTimer& debugTimer) {
    double scaleFactor = 1.0e6;

    if(debug && debugTimer.elapsed() % 500 < 16){
        qDebug() << "";
        qDebug() << "------values at each time step dt------" << dt;
    }

    // Calculate forces for spacecrafts
    for (Spacecraft &s : spacecrafts) {
        for (const CelestialBody &b : bodies) {
            double dx = b.x - s.x;
            double dy = b.y - s.y;
            double dz = b.z - s.z;
            double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            double force = G * s.mass * b.mass / (dist * dist * dist);

            // Scale the force to the scaled units
            double force_scaled = force / (scaleFactor * scaleFactor);

            s.updateVelocity(force_scaled * dx, force_scaled * dy, force_scaled * dz, dt, debugTimer);
        }
    }

    // Calculate forces for celestial bodies
    for (CelestialBody &body1 : bodies) {
        for (const CelestialBody &body2 : bodies) {
            if (&body1 != &body2) {  // Don't calculate force of body on itself
                double dx = body2.x - body1.x;
                double dy = body2.y - body1.y;
                double dz = body2.z - body1.z;
                double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                double force = G * body1.mass * body2.mass / (dist * dist * dist);

                // Scale the force to the scaled units
                double force_scaled = force / (scaleFactor * scaleFactor);

                body1.updateVelocity(force_scaled * dx, force_scaled * dy, force_scaled * dz, dt);
            }
        }
    }

    // Check for collisions
    if (checkCollisions() && debugTimer.elapsed() % 500 < 16){
        // If a collision is detected, print a message and exit the program
        std::cout << "Collision detected. Exiting program." << std::endl;
        exit(0);
    }

    // Update positions
    for (Spacecraft &s : spacecrafts) {
        s.updatePosition(dt);
    }
    for (CelestialBody &b : bodies) {
        b.updatePosition(dt);
    }

    // Update Octree
    // ...

    if(debug && debugTimer.elapsed() % 500 < 16){
        qDebug() << "";
    }
}


// Simulate engine firing
void Universe::fireEngine(int spacecraftIndex, double fx, double fy, double fz, double dt) {
    Spacecraft &s = spacecrafts[spacecraftIndex];
   // s.updateVelocity(fx, fy, fz, dt); add in qtimer
}

std::vector<Universe> runSimulation(int numSteps) {
    // ... (implement the runSimulation function here) ...
}
