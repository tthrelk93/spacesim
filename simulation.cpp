#include "simulation.h"
#include <QApplication>
#include <QOpenGLWindow>
#include <QOpenGLFunctions>
#include <iostream>
#include <QTimer>
#include <QElapsedTimer>
#include <unordered_map>
#include <string>
#include <QVector3D>

const bool debug = true;
OctreeNode::OctreeNode(const Point3D& center, double size)
    : center(center), size(size) {
    if (debug) {
       // std::cout << "Created OctreeNode at (" << center.x << ", " << center.y << ", " << center.z << ") with size " << size << std::endl;
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
        //std::cout << "Checking if object at (" << object.position.x << ", " << object.position.y << ", " << object.position.z << ") with radius " << object.radius
                  //<< " fits in node at (" << center.x << ", " << center.y << ", " << center.z << ") with size " << size << " fitsInNode: " << fitsInNode << ", " << " isSmallEnough: " << isSmallEnough << std::endl;
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
           // std::cout << "Inserting object at (" << object.position.x << ", " << object.position.y << ", " << object.position.z << ") "
                   //   << "into child " << i << " of node at (" << center.x << ", " << center.y << ", " << center.z << ")" << std::endl;
            children[i]->insert(object);
           // break;
       // }
    }
}

// Subdivide this node into 8 children
void OctreeNode::subdivide() {

    for (int i = 0; i < 8; ++i) {
        double x = center.x + size / 4 * ((i & 1) ? 1 : -1);
        double y = center.y + size / 4 * ((i & 2) ? 1 : -1);
        double z = center.z + size / 4 * ((i & 4) ? 1 : -1);
        children[i] = new OctreeNode({x, y, z}, size / 2);
        if (debug) {
          //  std::cout << "Created child " << i << " of node at (" << center.x << ", " << center.y << ", " << center.z << ")"
                     // << " with center at (" << x << ", " << y << ", " << z << ") and size " << size / 2 << std::endl;
        }
    }
    if (debug) {
      //  std::cout << "Subdividing node at (" << center.x << ", " << center.y << ", " << center.z << ")" << std::endl;
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
      //  std::cout << "Checking collisions in node at (" << center.x << ", " << center.y << ", " << center.z << ") with " << objects.size() << " objects" << std::endl;
        for (int i = 0; i < objects.size(); i++) {
            for (int j = i + 1; j < objects.size(); j++) {
                double dx = objects[i].position.x - objects[j].position.x;
                double dy = objects[i].position.y - objects[j].position.y;
                double dz = objects[i].position.z - objects[j].position.z;
                double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
                if (distance < objects[i].radius + objects[j].radius) {
                    // The objects are closer than the sum of their radii - collision!
                  //  std::cout << "Collision detected between objects at (" << objects[i].position.x << ", " << objects[i].position.y << ", " << objects[i].position.z << ") and (" << objects[j].position.x << ", " << objects[j].position.y << ", " << objects[j].position.z << ")" << std::endl;
                    return true;
                }
            }
        }
    }

    // If there's a single object in this node, but it doesn't fit, there's a collision
    if (objects.size() == 1 && !fits(objects[0])) {
       // std::cout << "Single object in node at (" << center.x << ", " << center.y << ", " << center.z << ") does not fit" << std::endl;
        return true;
    }

    // If there are no objects in this node, but there are children, check them for collisions
    if (children[0] != nullptr) {
       // std::cout << "Checking collisions in children of node at (" << center.x << ", " << center.y << ", " << center.z << ")" << std::endl;
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
       // qDebug() << "craft velocity (vx, vy, vz):" << vx << vy << vz;
    }
}

CelestialBody::CelestialBody(std::string name, double mass, double radius, QVector3D position, QVector3D velocity, double objectType, double renderSize, bool isMoon)
    : name(name), mass(mass), radius(radius), position(position), velocity(velocity), objectType(objectType), renderSize(renderSize), isMoon(isMoon) {
        qDebug() << "in constructor for: " << name; 
    }

void CelestialBody::updatePosition(double dt) {
    position += velocity * dt;
}

void CelestialBody::updateVelocity(QVector3D force, double dt) {
    velocity += force * dt / mass;
}


Universe::Universe()
    : octree({0, 0, 0}, 4000) {  // The size of the Octree is now 1.0
    if (debug) {
       // std::cout << "Created Universe with Octree root node at (0, 0, 0) and size 4000 " << std::endl;
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
    octree.insert({{b.position.x(), b.position.y(), b.position.z()}, b.radius});
    qDebug() << "inserting body into bodyMap: " << b.name;
    bodyMap.insert({b.name, b});
}

// Remove a celestial body
void Universe::removeBody(int index) {
    // You'll need to implement OctreeNode::remove to remove objects from the Octree
    //octree.remove({{bodies[index].position.x(), bodies[index].position.y(), bodies[index].position.z()}, bodies[index].radius});
    bodyMap.erase(bodies[index].name);
    bodies.erase(bodies.begin() + index);
}


// Check for collisions
bool Universe::checkCollisions() {
    return octree.checkCollisions();
}

bool fireDebugs = false;
void Universe::update(float dt, QElapsedTimer& debugTimer, bool isRungeKutta) {
    
    if(isRungeKutta){
        computeVelocity(dt);
        updateLocation(dt);
    } else {
        // Half step for velocity
        for (CelestialBody &body : bodies) {
            QVector3D acceleration = calculate_single_body_acceleration(dt, body, false);
            body.velocity += 0.5 * dt * acceleration;
        }

        // Full step for position
        for (CelestialBody &body : bodies) {
            body.position += dt * body.velocity;
        }

        // Half step for velocity
        for (CelestialBody &body : bodies) {
            QVector3D acceleration = calculate_single_body_acceleration(dt, body, false);
            body.velocity += 0.5 * dt * acceleration;
        }
    }
}


QVector3D Universe::calculate_single_body_acceleration(float dt, CelestialBody& target, bool isRungeKutta)
{
    double G_const = 6.67408e-11; //m3 kg - 1 s - 2
    QVector3D acceleration{ 0, 0, 0 };
    QVector3D velocity_update{ 0, 0, 0 };
    QVector3D location_update{ 0, 0, 0 };
    
    if(isRungeKutta){
       for (CelestialBody &externalBody : bodies) {
           
                if (&externalBody != &target) {

                    QVector3D k1{ 0, 0, 0 };
                    QVector3D k2{ 0, 0, 0 };
                    QVector3D k3{ 0, 0, 0 };
                    QVector3D k4{ 0, 0, 0 };

                    double r = (pow((target.position.x() - externalBody.position.x()), 2) +
                                pow((target.position.y() - externalBody.position.y()), 2) +
                                pow((target.position.z() - externalBody.position.z()), 2));

                    r = sqrt(r);

                    auto tmp = G_const * externalBody.mass / (r*r*r);

                    //k1 - acceleration at current location
                    k1 = (externalBody.position - target.position) * tmp;

                    //k2 - acceleration 0.5 timesteps in the future based on k1 acceleration value
                    velocity_update = partial_step(target.velocity, k1, 0.5, dt);
                    location_update = partial_step(target.position, velocity_update, 0.5, dt);
                    k2 = (externalBody.position - location_update) * tmp;

                    //k3 acceleration 0.5 timesteps in the future using k2 acceleration
                    velocity_update = partial_step(target.velocity, k2, 0.5, dt);
                    location_update = partial_step(target.position, velocity_update, 0.5, dt);
                    k3 = (externalBody.position - location_update) * tmp;

                    //k4 - location 1 timestep in the future using k3 acceleration
                    velocity_update = partial_step(target.velocity, k3, 1, dt);
                    location_update = partial_step(target.position, velocity_update, 1, dt);
                    k4 = (externalBody.position - location_update) * tmp;

                    acceleration += (k1 + k2 * 2 + k3 * 2 + k4) / 6;
                    
                 
                }
       }
           
        return acceleration;
} else {
    for (CelestialBody &externalBody : bodies) {
       
        if (&externalBody != &target) {
            QVector3D r = target.position - externalBody.position;
            double dist = r.length();
            acceleration -= G_const * externalBody.mass / (dist*dist*dist) * r;
        }
    }

    return acceleration;
}
       

}

QVector3D Universe::partial_step(QVector3D &f, QVector3D &df, double scale,float dt)
{
    return QVector3D{
                        static_cast<float>(f.x() + df.x() * dt * scale),
                        static_cast<float>(f.y() + df.y() * dt * scale),
                        static_cast<float>(f.z() + df.z() * dt * scale)
    };
}


void Universe::computeVelocity(float dt)
{
 for (CelestialBody &body1 : bodies) {
		QVector3D acceleration = Universe::calculate_single_body_acceleration(dt, body1, true);
		body1.velocity += acceleration * dt;
	}

}

void Universe::updateLocation(float dt)
{
   for (CelestialBody &body1 : bodies) {

        body1.position += body1.velocity * dt;
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
