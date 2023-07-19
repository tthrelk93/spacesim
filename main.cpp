#include <GL/glew.h>
#include </Users/thomasthrelkeld/Qt/6.6.0/macos/lib/QtWidgets.framework/Versions/A/Headers/QApplication>
#include <QOpenGLWindow>
#include <QOpenGLFunctions>
#include <QTimer>
#include <QElapsedTimer>
#include <glm/glm.hpp>  
#include <SDL.h>
#include <QKeyEvent>
#include <QWheelEvent>
#include <QSlider>
#include <QVBoxLayout>
#include <QOpenGLWidget>
#include <QPainter>
#include <unordered_map>
#include <string>
#include <QPainter>
#include <QString>
#include <QLabel>
#include <QCheckBox>
#include <QListWidget>
#include <QListWidgetItem>
#include "simulation.h"
#include <QVector3D>
#include <iostream>
#include <QQuaternion>


class MyOpenGLWindow : public QOpenGLWidget, protected QOpenGLFunctions {

private:
    Universe universe;  // Your Universe object
    float cameraRadius = 0.0f;
    float changeFocusCameraRadiusScale = 0.2f;
    //These renderRadii do not effects the physics at all. The are simply visual. The value for the sun was arbitrarily set to a value that seemed visually large enough
    //and then the other sizes are based on the realistic percentage of the sun's size that each celestial body is in real life
    float sunRenderRadius = 8550000000.0f;
    float earthRenderRadius = 900000000.0f;
    float moonRenderRadius = 24000000.0f;
    float venusRenderRadius = 855000000.0f;
    float mercuryRenderRadius = 300000000.0f;
    float marsRenderRadius = 400000000.0f;
    float jupiterRenderRadius = 9873000000.0f;
    float saturnRenderRadius = 8226000000.0f;
    float uranusRenderRadius = 3582000000.0f;
    float neptuneRenderRadius = 3474000000.0f;

    bool topDownView = false;
    bool mousePress = false;

    float angleX = 0.0f;
    float angleY = 0.0f;
    float angleZ = 0.0f;
    QPoint lastMousePos;
    float scaleFactor = 1.0;  // The scale factor you used to scale down the universe. Currently set to 1 because it gave the most realistic physics
    const bool debug = true;
    const bool SUN_RELATIVE_LSR = false;
    const bool USE_RUNGE_KUTTA = false;
    int MAX_PAST_POSITIONS = 100;
    QElapsedTimer debugTimer;
    float timeScale = 1.0f;
    const double GRAVITATIONAL_CONSTANT = 6.67430e-11;  // in m^3 kg^-1 s^-2
    CelestialBody* selectedBody;
     double SUN_MASS = 1.989e30;  // in kilograms


public:
     QTimer simulationTimer;
     QTimer renderTimer;
     double visRadiusScale = 1/30000.0f; // use smallest celestial body
    float miniMapSize = 400.0f; // Adjust the size of the mini map as needed
    float miniMapMargin = 10.0f; // Adjust the margin of the mini map as needed
     QVector3D cameraPos;
    QListWidget* listWidget;

    MyOpenGLWindow(QWidget *parent = nullptr)
        : QOpenGLWidget(parent) {
        connect(&simulationTimer, &QTimer::timeout, this, &MyOpenGLWindow::updateUniverse);
        connect(&renderTimer, &QTimer::timeout, this, &MyOpenGLWindow::updateWrapper);
        listWidget = new QListWidget;
        debugTimer.start();
        addSun();
        addMercury();
        addVenus();
        addEarth();
        addMars();
       // addJupiter();
        //addSaturn();
       // addNeptune();
       // addUranus();

    }


  

    void addSun(){
        const double SUN_RADIUS = 6.9634e8;  // in meters
        QVector3D SUN_VELOCITY(0, 0, 0);
        if(SUN_RELATIVE_LSR){
            SUN_VELOCITY = QVector3D(-9700.0, 20000.0, 7000.0);  // in meters per second (approximate value for solar apex motion)
        }
        // Scale the velocity components based on the scaleFactor
        QVector3D scaledSunVelocity = SUN_VELOCITY / scaleFactor;
        // Set the initial position of the Sun
        QVector3D initialSunPosition(0.0, 0.0, 0.0);  // Set the initial position of the Sun
        // Set the Sun's velocity and position
        CelestialBody sun("Sun", SUN_MASS, SUN_RADIUS / scaleFactor, initialSunPosition, scaledSunVelocity, 0.0, sunRenderRadius*visRadiusScale, false);
        universe.addBody(sun);
        selectedBody = &sun;
        qDebug() << "Adding Sun with position:" << initialSunPosition << "and velocity:" << scaledSunVelocity;

    }

    void addMercury(){
        // Constants for Mercury
        const double MERCURY_MASS = 3.3011e23;  // in kilograms
        const double MERCURY_RADIUS = 2.4397e6;  // in meters
        const double SUN_MERCURY_DISTANCE = 5.791e10;  // in meters

        // Calculate the velocity of Mercury
        double mercury_velocity_magnitude = sqrt(GRAVITATIONAL_CONSTANT * SUN_MASS / SUN_MERCURY_DISTANCE) / scaleFactor;
        double mercury_inclination = 7.0 * M_PI / 180.0;

        // Set the Mercury's velocity and position
        QVector3D mercury_position(SUN_MERCURY_DISTANCE / scaleFactor * cos(mercury_inclination), 0, SUN_MERCURY_DISTANCE / scaleFactor * sin(mercury_inclination));
        QVector3D mercury_velocity(0, -mercury_velocity_magnitude, 0);
        CelestialBody mercury("Mercury", MERCURY_MASS, MERCURY_RADIUS / scaleFactor, mercury_position, mercury_velocity, 1.0, mercuryRenderRadius*visRadiusScale, false);
        universe.addBody(mercury);
        qDebug() << "Adding Mercury with position:" << mercury_position << "and velocity:" << mercury_velocity;
    }

       void addVenus(){
       // Constants for Venus
       const double VENUS_MASS = 4.8675e24;  // in kilograms
       const double VENUS_RADIUS = 6.0518e6;  // in meters
       const double SUN_VENUS_DISTANCE = 1.0821e11;  // in meters
       
       double venus_velocity_magnitude = sqrt(GRAVITATIONAL_CONSTANT * SUN_MASS / SUN_VENUS_DISTANCE) / scaleFactor;
       double venus_inclination = 3.4 * M_PI / 180.0;
       QVector3D venus_position(SUN_VENUS_DISTANCE / scaleFactor * cos(venus_inclination), 0, SUN_VENUS_DISTANCE / scaleFactor * sin(venus_inclination));
        QVector3D venus_velocity(0, -venus_velocity_magnitude, 0);
       // Set the Venus's velocity and position
       CelestialBody venus("Venus", VENUS_MASS, VENUS_RADIUS / scaleFactor, venus_position, venus_velocity, 2.0, venusRenderRadius*visRadiusScale, false);
       universe.addBody(venus);
       qDebug() << "Adding Venus with position:" << venus_position << "and velocity:" << venus_velocity;
   }

       void addEarth(){
       // Constants for the Earth
       const double EARTH_MASS = 5.972e24;  // in kilograms
       const double EARTH_RADIUS = 6.371e6;  // in meters
       const double EARTH_SUN_DISTANCE = 1.496e11;  // in meters
       // Constants for the Moon
       const double MOON_MASS = 7.342e22;  // in kilograms
       const double MOON_RADIUS = 1.737e6;  // in meters
       // Distance from the Earth to the Moon
       const double EARTH_MOON_DISTANCE = 3.844e8;  // in meters
       // Calculate the velocity of Earth
        double earth_velocity_magnitude = sqrt(GRAVITATIONAL_CONSTANT * SUN_MASS / EARTH_SUN_DISTANCE) / scaleFactor;
        double earth_inclination = 0.0 * M_PI / 180.0;

        QVector3D earth_position(EARTH_SUN_DISTANCE / scaleFactor * cos(earth_inclination), 0, EARTH_SUN_DISTANCE / scaleFactor * sin(earth_inclination));
        QVector3D earth_velocity(0, -earth_velocity_magnitude * cos(earth_inclination), -earth_velocity_magnitude * sin(earth_inclination));
       // QVector3D earth_velocity(0, -earth_velocity_magnitude, 0);
        CelestialBody earth("Earth", EARTH_MASS, EARTH_RADIUS / scaleFactor, earth_position, earth_velocity, 3, earthRenderRadius*visRadiusScale, false);

        
        double moon_velocity_magnitude = sqrt(GRAVITATIONAL_CONSTANT * EARTH_MASS / EARTH_MOON_DISTANCE) / scaleFactor;
        double moon_inclination = 5.14 * M_PI / 180.0;
        QVector3D moon_position = earth_position + QVector3D(EARTH_MOON_DISTANCE / scaleFactor * cos(moon_inclination), 0, EARTH_MOON_DISTANCE / scaleFactor * sin(moon_inclination));
        QVector3D moon_velocity = QQuaternion::fromAxisAndAngle(QVector3D(0, 1, 0), earth_inclination * 180.0 / M_PI) * QVector3D(0, moon_velocity_magnitude * cos(moon_inclination), moon_velocity_magnitude * sin(moon_inclination));

        CelestialBody moon("Moon", MOON_MASS, MOON_RADIUS / scaleFactor, moon_position, earth_velocity + moon_velocity, 3.1, moonRenderRadius*visRadiusScale, true);   
        //earth.moons.push_back(moon);

        universe.addBody(moon);
        universe.addBody(earth);
        
        qDebug() << "Adding Earth with position:" << earth_position << "and velocity:" << earth_velocity;
         qDebug() << "Adding Moon with position:" << moon_position << "and velocity:" << moon_velocity;
   }

    void addMars(){
       // Constants for Mars
    const double MARS_MASS = 0.64171e24;  // in kilograms
    const double MARS_RADIUS = 3.3895e6;  // in meters
    const double SUN_MARS_DISTANCE = 2.2792e11;  // in meters

    // Constants for Phobos
    const double PHOBOS_MASS = 1.0659e16;  // in kilograms
    const double PHOBOS_RADIUS = 1.1e4;  // in meters
    const double MARS_PHOBOS_DISTANCE = 9.377e6;  // in meters

    // Constants for Deimos
    const double DEIMOS_MASS = 1.4762e15;  // in kilograms
    const double DEIMOS_RADIUS = 6.2e3;  // in meters
    const double MARS_DEIMOS_DISTANCE = 2.346e7;  // in meters

    // Calculate the velocity of Mars
    double mars_velocity_magnitude = sqrt(GRAVITATIONAL_CONSTANT * SUN_MASS / SUN_MARS_DISTANCE) / scaleFactor;
    double mars_inclination = 1.85 * M_PI / 180.0;

    QVector3D mars_position(SUN_MARS_DISTANCE / scaleFactor * cos(mars_inclination), 0, SUN_MARS_DISTANCE / scaleFactor * sin(mars_inclination));
    QVector3D mars_velocity(0, -mars_velocity_magnitude, 0);
    CelestialBody mars("Mars", MARS_MASS, MARS_RADIUS / scaleFactor, mars_position, mars_velocity, 4.0, marsRenderRadius*visRadiusScale, false);

   // Calculate the velocity of Phobos
double phobos_velocity_magnitude = sqrt(GRAVITATIONAL_CONSTANT * MARS_MASS / MARS_PHOBOS_DISTANCE) / scaleFactor;
double phobos_inclination = 1.08 * M_PI / 180.0;  // Inclination of Phobos's orbit
QVector3D phobos_position = mars_position + QVector3D(MARS_PHOBOS_DISTANCE / scaleFactor * cos(phobos_inclination), 0, MARS_PHOBOS_DISTANCE / scaleFactor * sin(phobos_inclination));
QVector3D phobos_velocity = QQuaternion::fromAxisAndAngle(QVector3D(0, 1, 0), mars_inclination * 180.0 / M_PI) * QVector3D(0, phobos_velocity_magnitude * cos(phobos_inclination), phobos_velocity_magnitude * sin(phobos_inclination));
phobos_velocity = phobos_velocity - mars_velocity;  // Transform to the reference frame of the Sun
CelestialBody phobos("Phobos", PHOBOS_MASS, PHOBOS_RADIUS / scaleFactor, phobos_position, phobos_velocity, 4.1, moonRenderRadius*visRadiusScale, true);

// Calculate the velocity of Deimos
double deimos_velocity_magnitude = sqrt(GRAVITATIONAL_CONSTANT * MARS_MASS / MARS_DEIMOS_DISTANCE) / scaleFactor;
double deimos_inclination = 1.79 * M_PI / 180.0;  // Inclination of Deimos's orbit
QVector3D deimos_position = mars_position + QVector3D(MARS_DEIMOS_DISTANCE / scaleFactor * cos(deimos_inclination), 0, MARS_DEIMOS_DISTANCE / scaleFactor * sin(deimos_inclination));
QVector3D deimos_velocity = QQuaternion::fromAxisAndAngle(QVector3D(0, 1, 0), mars_inclination * 180.0 / M_PI) * QVector3D(0, deimos_velocity_magnitude * cos(deimos_inclination), deimos_velocity_magnitude * sin(deimos_inclination));
deimos_velocity = deimos_velocity - mars_velocity;  // Transform to the reference frame of the Sun
CelestialBody deimos("Deimos", DEIMOS_MASS, DEIMOS_RADIUS / scaleFactor, deimos_position, deimos_velocity, 4.2, moonRenderRadius*visRadiusScale, true);


    universe.addBody(phobos);
    universe.addBody(deimos);
    universe.addBody(mars);
    qDebug() << "Adding Mars with position:" << mars_position << "and velocity:" << mars_velocity;
    qDebug() << "Adding Phobos with position:" << phobos_position << "and velocity:" << phobos_velocity;
    qDebug() << "Adding Deimos with position:" << deimos_position << "and velocity:" << deimos_velocity;
    }
      

   void addJupiter(){
       // Constants for Jupiter
       const double JUPITER_MASS = 1898.19e24;  // in kilograms
       const double JUPITER_RADIUS = 69.911e6;  // in meters
       const double SUN_JUPITER_DISTANCE = 7.7857e11;  // in meters

        double jupiter_velocity_magnitude = sqrt(GRAVITATIONAL_CONSTANT * SUN_MASS / SUN_JUPITER_DISTANCE) / scaleFactor;
       double jupiter_inclination = 1.304 * M_PI / 180.0;
       QVector3D jupiter_position(SUN_JUPITER_DISTANCE / scaleFactor * cos(jupiter_inclination), 0, SUN_JUPITER_DISTANCE / scaleFactor * sin(jupiter_inclination));
        QVector3D jupiter_velocity(0, -jupiter_velocity_magnitude, 0);
       // Set the Jupiter's velocity and position
       CelestialBody jupiter("Jupiter", JUPITER_MASS, JUPITER_RADIUS / scaleFactor, jupiter_position, jupiter_velocity, 5.0, jupiterRenderRadius*visRadiusScale, false);
       universe.addBody(jupiter);
       qDebug() << "Adding Jupiter with position:" << jupiter_position << "and velocity:" << jupiter_velocity;
   }

   void addSaturn(){
       // Constants for Saturn
       const double SATURN_MASS = 568.34e24;  // in kilograms
       const double SATURN_RADIUS = 58.232e6;  // in meters
       const double SUN_SATURN_DISTANCE = 1.4335e12;  // in meters
       double saturn_velocity_magnitude = sqrt(GRAVITATIONAL_CONSTANT * SUN_MASS / SUN_SATURN_DISTANCE) / scaleFactor;
       double saturn_inclination = 2.485 * M_PI / 180.0;
        QVector3D saturn_position(SUN_SATURN_DISTANCE / scaleFactor * cos(saturn_inclination), 0, SUN_SATURN_DISTANCE / scaleFactor * sin(saturn_inclination));
        QVector3D saturn_velocity(0, -saturn_velocity_magnitude, 0);
      
       // Set the Saturn's velocity and position
       CelestialBody saturn("Saturn", SATURN_MASS, SATURN_RADIUS / scaleFactor, saturn_position, saturn_velocity, 6.0, saturnRenderRadius*visRadiusScale, false);
       universe.addBody(saturn);
       qDebug() << "Adding Saturn with position:" << saturn_position << "and velocity:" << saturn_velocity;
   }

   void addNeptune(){
       // Constants for Neptune
       const double NEPTUNE_MASS = 102.413e24;  // in kilograms
       const double NEPTUNE_RADIUS = 24.622e6;  // in meters
       const double SUN_NEPTUNE_DISTANCE = 4.4951e12;  // in meters
       double neptune_velocity_magnitude = sqrt(GRAVITATIONAL_CONSTANT * SUN_MASS / SUN_NEPTUNE_DISTANCE) / scaleFactor;
       double neptune_inclination = 1.767 * M_PI / 180.0;
        QVector3D neptune_position(SUN_NEPTUNE_DISTANCE / scaleFactor * cos(neptune_inclination), 0, SUN_NEPTUNE_DISTANCE / scaleFactor * sin(neptune_inclination));
        QVector3D neptune_velocity(0, -neptune_velocity_magnitude, 0);
       // Set the Neptune's velocity and position
       CelestialBody neptune("Neptune", NEPTUNE_MASS, NEPTUNE_RADIUS / scaleFactor, neptune_position, neptune_velocity, 7.0, neptuneRenderRadius*visRadiusScale, false);
       universe.addBody(neptune);
       qDebug() << "Adding Neptune with position:" << neptune_position << "and velocity:" << neptune_velocity;
   }

   void addUranus(){
       // Constants for Uranus
       const double URANUS_MASS = 86.813e24;  // in kilograms
       const double URANUS_RADIUS = 25.362e6;  // in meters
       const double SUN_URANUS_DISTANCE = 2.8725e12;  // in meters
        double uranus_velocity_magnitude = sqrt(GRAVITATIONAL_CONSTANT * SUN_MASS / SUN_URANUS_DISTANCE) / scaleFactor;
       double uranus_inclination = 0.772 * M_PI / 180.0;
       QVector3D uranus_position(SUN_URANUS_DISTANCE / scaleFactor * cos(uranus_inclination), 0, SUN_URANUS_DISTANCE / scaleFactor * sin(uranus_inclination));
        QVector3D uranus_velocity(0, -uranus_velocity_magnitude, 0);
       // Set the Uranus's velocity and position
       CelestialBody uranus("Uranus", URANUS_MASS, URANUS_RADIUS / scaleFactor, uranus_position, uranus_velocity, 8.0, uranusRenderRadius*visRadiusScale, false);
       universe.addBody(uranus);
       qDebug() << "Adding Uranus with position:" << uranus_position << "and velocity:" << uranus_velocity;
   }




protected:

    void initializeGL() override {
        initializeOpenGLFunctions();
        // Enable depth testing
        glEnable(GL_DEPTH_TEST);
        // Set the clear color to black
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);  // Change to a gray color
        // Generate sphere vertices
        std::vector<GLfloat> sphereVertices = generateSphereVertices(1.0f, 64, 32);
        // Create a vertex buffer object (VBO)
        GLuint vbo;
        glGenBuffers(1, &vbo);
        // Bind the VBO and upload the vertex data
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, sphereVertices.size() * sizeof(GLfloat), &sphereVertices[0], GL_STATIC_DRAW);
        // Set up the projection matrix
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        // Set up the model-view matrix
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        updateCameraPosition();
    }

    void wheelEvent(QWheelEvent *event) override {
        // Get the number of degrees the wheel was rotated
        float degrees = event->angleDelta().y() / 8.0f;
        // Convert degrees to steps (a full step is 15 degrees)
        float steps = degrees / 15.0f;
        // Adjust the camera radius based on the number of steps
        cameraRadius -= steps * 15000.0f;  // Adjust the factor as needed


    }

    void mousePressEvent(QMouseEvent *event) override {
        if (event->button() == Qt::LeftButton) {
            mousePress = true;
            lastMousePos = event->pos();
        }
    }

    void mouseMoveEvent(QMouseEvent* event) override {
        if (mousePress) {
            QPoint delta = event->pos() - lastMousePos;
            angleX += delta.x() * 0.05f;
            angleY += delta.y() * 0.05f;
            //angleY = std::max(-359.9f, std::min(359.9f, angleY));


            //updateCameraPosition();

            lastMousePos = event->pos();
        }
    }

    void mouseReleaseEvent(QMouseEvent *event) override {
        if (event->button() == Qt::LeftButton) {
            mousePress = false;
        }
    }

    void updateCameraPosition() {

        if(!topDownView){

            float x = cameraRadius * std::sin(angleY * M_PI / 180.0f) * std::cos(angleX * M_PI / 180.0f);
            float y = cameraRadius * std::sin(angleY * M_PI / 180.0f) * std::sin(angleX * M_PI / 180.0f);
            float z = cameraRadius * std::cos(angleY * M_PI / 180.0f);

            // Add the position of the object to the camera position
            cameraPos = QVector3D(selectedBody->position.x() * visRadiusScale + x, selectedBody->position.y() * visRadiusScale + y, selectedBody->position.z() * visRadiusScale + z);
        } else {
            double zVal = 50000;
            // For top-down view, set the camera's position to match the selected body's position and add an offset in the z direction
            cameraPos = QVector3D(selectedBody->position.x() * visRadiusScale, selectedBody->position.y() * visRadiusScale, zVal);

        }

        update(); // Request a redraw of the scene
    }


    void paintGL() override {
        
        // Clear the screen
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        // Enable depth testing for correct rendering order
       glEnable(GL_DEPTH_TEST);

        if(selectedBody->name == ""){
            for (CelestialBody& body : universe.bodies) {
                if(body.name == "Sun"){
                    selectedBody = &body;
                    break;
                }
            }
       }

        // Draw the celestial bodies
        for (CelestialBody& body : universe.bodies) {

            glPushMatrix();

            if (body.objectType == 0.0) {
                 glColor3f(1.0f, 1.0f, 0.0f);  // Set the color to yellow for the sun
            } else if (body.objectType == 1.0) {
                glColor3f(0.0f, 1.0f, 0.0f);  // Set the color to green for Mercury           
            } else if (body.objectType == 2.0) {
                 glColor3f(1.0f, 0.0f, 1.0f);  // Set the color to magenta for Venus
            } else if (body.objectType == 3.0){
                  glColor3f(0.0f, 0.0f, 1.0f);  // Set the color to blue for Earth
            } else if (body.objectType == 4.0){
                 glColor3f(0.5f, 0.0f, 0.0f);  // Set the color to orange for Mars
            } else if (body.objectType == 5.0){
                 glColor3f(1.0f, 0.5f, 0.0f);  // Set the color to dark red for Jupiter
            } else if (body.objectType == 6.0){
                 glColor3f(0.5f, 0.5f, 0.0f);  // Set the color to olive for Saturn
            } else if (body.objectType == 7.0){
                 glColor3f(0.0f, 0.0f, 0.5f);  // Set the color to navy for Neptune               
            } else if (body.objectType == 8.0){
                glColor3f(0.0f, 0.5f, 0.5f);  // Set the color to teal for Uranus
            } else {
                glColor3f(1.0f, 1.0f, 1.0f);  // Set the color to white for Moons
            }

            QVector3D scalePosition(body.position.x() * visRadiusScale, body.position.y() * visRadiusScale, body.position.z() * visRadiusScale);

            glTranslatef(scalePosition.x(), scalePosition.y(), scalePosition.z());
            drawSphere(body.renderSize, 20, 20, body.objectType, false);  // Pass the radius of the celestial body
            glPopMatrix();

            glBegin(GL_LINE_STRIP);
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            float alpha = 1.0f;  // Start with full opacity
            int count = 0;  // Count of past positions
            for (const QVector3D& pastPosition : body.pastPositions) {
                alpha -= 0.005;  // Increase the rate of decrease
                alpha = std::max(alpha, 0.0f);  // Ensure alpha never goes below 0
                
                if (body.objectType == 0.0) {
                 glColor4f(1.0f, 1.0f, 0.0f, alpha);  // Set the color to yellow for the sun
            } else if (body.objectType == 1.0) {
                glColor4f(0.0f, 1.0f, 0.0f, alpha);  // Set the color to green for Mercury           
            } else if (body.objectType == 2.0) {
                 glColor4f(1.0f, 0.0f, 1.0f, alpha);  // Set the color to magenta for Venus
            } else if (body.objectType == 3.0){
                  glColor4f(0.0f, 0.0f, 1.0f, alpha);  // Set the color to blue for Earth
            } else if (body.objectType == 4.0){
                 glColor4f(0.5f, 0.0f, 0.0f, alpha);  // Set the color to orange for Mars
            } else if (body.objectType == 5.0){
                 glColor4f(1.0f, 0.5f, 0.0f, alpha);  // Set the color to dark red for Jupiter
            } else if (body.objectType == 6.0){
                 glColor4f(0.5f, 0.5f, 0.0f, alpha);  // Set the color to olive for Saturn
            } else if (body.objectType == 7.0){
                 glColor4f(0.0f, 0.0f, 0.5f, alpha);  // Set the color to navy for Neptune               
            } else if (body.objectType == 8.0){
                glColor4f(0.0f, 0.5f, 0.5f, alpha);  // Set the color to teal for Uranus
            } else {
                glColor4f(1.0f, 1.0f, 1.0f, alpha);  // Set the color to white for Moons
            }
                glVertex3f(pastPosition.x()  * visRadiusScale, pastPosition.y()  * visRadiusScale, pastPosition.z() * visRadiusScale);
                 count++;
                  // Limit the number of past positions
                if (count >= MAX_PAST_POSITIONS) {
                    break;
                }
            }
             glEnd();
            if(body.name == selectedBody->name){
                 selectedBody = &body;
            }
        }

        QVector3D scalePos(selectedBody->position.x() * visRadiusScale, selectedBody->position.y() * visRadiusScale, selectedBody->position.z() * visRadiusScale);
            // Set up the projection matrix for the top-down view or perspective view
            if (topDownView) {
                 // Set up the projection matrix for the top-down view
                 glMatrixMode(GL_PROJECTION);
                 glLoadIdentity();

               double zoomFactor = 0.005;  // Adjust this value to zoom in or out

                double left = (double)((this->width()*2)*-1) / zoomFactor;
                double right = (double)(this->width()*2) / zoomFactor;
                double bottom = (double)((this->width()*2)*-1) / zoomFactor;
                double top = (double)(this->width()*2) / zoomFactor;
                 double zNear = 1.0;  // Set near clipping plane based on camera distance
                 double zFar = 10000000000000000000.0;  // Set far clipping plane based on camera distance


                 glOrtho(left, right, bottom, top, zNear, zFar);
                 glMatrixMode(GL_MODELVIEW);// Set up the model-view matrix for top-down view
                 glLoadIdentity();

                 gluLookAt(cameraPos.x(), cameraPos.y(), cameraPos.z(),
                           scalePos.x() , scalePos.y(), scalePos.z(),
                           0.0f, 0.0f, 1.0f);

            } else {

                 glMatrixMode(GL_PROJECTION);// Set up the projection matrix for the perspective view
                 glLoadIdentity();
                 double fovy = 90;  // Field of view angle, in degrees, in the y direction
                 double aspect = (double)this->width() / (double)this->height();  // Aspect ratio (width to height)
                 double zNear = 1.0;  // Distance to near clipping plane;  // Distance to near clipping plane
                 double zFar = 10000000000000000000;  // Distance to far clipping plane
                 gluPerspective(fovy, aspect, zNear, zFar);

                 glMatrixMode(GL_MODELVIEW);
                 glLoadIdentity();

                 gluLookAt(cameraPos.x(), cameraPos.y(), cameraPos.z(),
                           scalePos.x() , scalePos.y(), scalePos.z(),
                           0.0f, 0.0f, 1.0f);

            }

    }


    void drawSphere(double r, int lats, int longs, int objectType, bool isMiniMapObj) {
         glPushMatrix();  // Push the current matrix onto the stack

        double visibleRadius = r;
        //qDebug() << "vis radius: " << visibleRadius;
        glScalef(visibleRadius, visibleRadius, visibleRadius);  // Scale

        int i, j;
        for(i = 0; i <= lats; i++) {
            double lat0 = M_PI * (-0.5 + (double) (i - 1) / lats);
            double z0  = sin(lat0);
            double zr0 =  cos(lat0);

            double lat1 = M_PI * (-0.5 + (double) i / lats);
            double z1 = sin(lat1);
            double zr1 = cos(lat1);

            glBegin(GL_QUAD_STRIP);
            for(j = 0; j <= longs; j++) {
                double lng = 2 * M_PI * (double) (j - 1) / longs;
                double x = cos(lng);
                double y = sin(lng);

                glNormal3f(x * zr0, y * zr0, z0);
                glVertex3f(x * zr0, y * zr0, z0);
                glNormal3f(x * zr1, y * zr1, z1);
                glVertex3f(x * zr1, y * zr1, z1);
            }
            glEnd();
        }
        glPopMatrix();  // Pop the current matrix off the stack
    }
    std::vector<GLfloat> generateSphereVertices(float radius, unsigned int rings, unsigned int sectors) {
        std::vector<GLfloat> vertices;
        float const R = 1.0f / (float)(rings - 1);
        float const S = 1.0f / (float)(sectors - 1);
        for (unsigned int r = 0; r < rings; ++r) {
            for (unsigned int s = 0; s < sectors; ++s) {
                float const y = sin(-M_PI_2 + M_PI * r * R);
                float const x = cos(2 * M_PI * s * S) * sin(M_PI * r * R);
                float const z = sin(2 * M_PI * s * S) * sin(M_PI * r * R);

                vertices.push_back(x * radius);
                vertices.push_back(y * radius);
                vertices.push_back(z * radius);
            }
        }
        return vertices;
    }

    void resizeGL(int w, int h) override {
        // Set the viewport to cover the entire widget
        glViewport(0, 0, w, h);
        // Set up the projection matrix
        glMatrixMode(GL_PROJECTION);// Set up the projection matrix for the perspective view
        glLoadIdentity();
        double fovy = 90;  // Field of view angle, in degrees, in the y direction
        double aspect = (double)this->width() / (double)this->height();  // Aspect ratio (width to height)
        double zNear = 1.0;  // Distance to near clipping plane;  // Distance to near clipping plane
        double zFar = 10000000000000000000;  // Distance to far clipping plane
        gluPerspective(fovy, aspect, zNear, zFar);               
    }

public slots:
    void updateWrapper() {
         update();
    }
    void toggleTopDownView(int state) {
         topDownView = (state == Qt::Checked);
         
         updateCameraPosition();

    }

    void changeFocus(QListWidgetItem* current, QListWidgetItem* previous) {
         // Extract the body name from the item text
         if (current) {
            QString bodyName = current->text().split(":").first();
            // Find the body in the universe
            if (universe.bodyMap.count(bodyName.toStdString()) > 0) {
                selectedBody = &universe.bodyMap[bodyName.toStdString()];
               
                updateCameraPosition();
            }
         }
    }


    void adjustTimeScale(float value) {
       timeScale = static_cast<float>(value);  // Adjust the scaling factor as needed
    }


    void updateUniverse() {
       
       universe.update(timeScale/1,debugTimer, USE_RUNGE_KUTTA); // Assuming 60 frames per second

       listWidget->clear();
       // Add the text for each celestial body
       for (const CelestialBody& body : universe.bodies) {
            //todo do i need to apply the scale here?
            QString text = QString("%1: Position (%2, %3, %4), Velocity (%5, %6, %7)")
                               .arg(QString::fromStdString(body.name))
                               .arg(body.position.x() * visRadiusScale).arg(body.position.y()* visRadiusScale ).arg(body.position.z() * visRadiusScale)
                               .arg(body.velocity.x()).arg(body.velocity.y()).arg(body.velocity.z());
            listWidget->addItem(text);
       }
       // Update the list of past positions for each celestial body
       for (CelestialBody& body : universe.bodies) {
            body.pastPositions.push_front(body.position);  // Add the current position to the front of the list
            if(body.name == selectedBody->name){
                selectedBody = &body;
            }
       }

       // Update the selected body's position
       if (universe.bodyMap.count(selectedBody->name) > 0) {
            selectedBody = &universe.bodyMap[selectedBody->name];
       }
       updateCameraPosition();
    }

};

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);

    // Create a QWidget to serve as the main window
    QWidget window;
    window.resize(1000, 1000);

    // Create a QGridLayout to arrange the widgets
    QGridLayout* layout = new QGridLayout(&window);

    // Create the MyOpenGLWindow
    MyOpenGLWindow* openglWindow = new MyOpenGLWindow;

    openglWindow->simulationTimer.start(100);
    openglWindow->renderTimer.start(16);
    // Create the QSlider
    QSlider* timeSlider = new QSlider(Qt::Horizontal);
    timeSlider->setRange(1, 100000);  // Adjust the range as needed
    timeSlider->setValue(1);

    // Create the QCheckBox for toggle switch
    QCheckBox* topDownCheckBox = new QCheckBox("Top-Down View");
    layout->addWidget(topDownCheckBox, 3, 1, 1, 1);
    //topDownCheckBox->setChecked(true);
    //topDownCheckBox->setCheckState(Qt::CheckState(true));

    // Connect the stateChanged signal of the QCheckBox to the toggleTopDownView slot
    QObject::connect(topDownCheckBox, &QCheckBox::stateChanged, openglWindow, &MyOpenGLWindow::toggleTopDownView);

    openglWindow->listWidget->setFixedSize(600, 230);  // Adjust the width and height as needed

    // Add the MyOpenGLWindow to the layout, spanning 2 rows and 2 columns
    layout->addWidget(openglWindow, 0, 0, 2, 2);

    // Add the QSlider to the layout in the third row
    layout->addWidget(timeSlider, 2, 0, 1, 2);

    // Add the QListWidget to the layout in the top right corner
    layout->addWidget(openglWindow->listWidget, 0, 1, 1, 1);
    openglWindow->listWidget->setFixedSize(620, 70);  // Adjust the width and height as needed
    // Check if the signal is being emitted
    QObject::connect(openglWindow->listWidget, &QListWidget::currentItemChanged, openglWindow, &MyOpenGLWindow::changeFocus);

    // Create a QLabel to display the camera coordinates
    QLabel* cameraLabel = new QLabel(&window);
    layout->addWidget(cameraLabel, 0, 0, 1, 1);
    cameraLabel->setAlignment(Qt::AlignTop | Qt::AlignLeft);

    QObject::connect(timeSlider, &QSlider::valueChanged, [openglWindow](int value) {
        float floatValue = static_cast<float>(value);
        openglWindow->adjustTimeScale(floatValue);
    });
    // Show the main window
    window.show();

    // Set up a timer to update the camera coordinates label every 200 milliseconds
    QTimer cameraTimer;
    QObject::connect(&cameraTimer, &QTimer::timeout, [&openglWindow, &cameraLabel]() {
        QString cameraText = QString("Camera Position: (%1, %2, %3)")
                                 .arg(openglWindow->cameraPos.x() )
                                 .arg(openglWindow->cameraPos.y() )
                                 .arg(openglWindow->cameraPos.z() );
        cameraLabel->setText(cameraText);
    });
    cameraTimer.start(200);  // Update every 200 milliseconds

    return app.exec();
}



