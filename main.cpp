#include <GL/glew.h>
#include <QApplication>
#include <QOpenGLWindow>
#include <QOpenGLFunctions>
#include <QTimer>
#include <QElapsedTimer>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
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
#include <QListWidget>
#include <QListWidgetItem>
#include "simulation.h"

class MyOpenGLWindow : public QOpenGLWidget, protected QOpenGLFunctions {

private:
    Universe universe;  // Your Universe object
    float cameraX = 0;
    float cameraY = 0;
    float cameraZ = 5000.0f;  // Adjust this value as needed

    bool mousePress = false;
    float angleX = 0.0f;
    float angleY = 0.0f;
    float angleZ = 0.0f;
    QPoint lastMousePos;

    float scaleFactor = 1.0f;
    // Constants for the Earth
    const double EARTH_MASS = 5.972e24;  // in kilograms
    const double EARTH_RADIUS = 6.371e6;  // in meters
    // Constants for the Moon
    const double MOON_MASS = 7.342e22;  // in kilograms
    const double MOON_RADIUS = 1.737e6;  // in meters
    // Distance from the Earth to the Moon
    const double EARTH_MOON_DISTANCE = 3.844e8;  // in meters
    // Low Earth orbit altitude
    const double LOW_EARTH_ORBIT_ALTITUDE = 2.0e6;  // in meters
    const bool debug = true;
    QElapsedTimer debugTimer;
    int timeScale = 1;
    const double GRAVITATIONAL_CONSTANT = 6.67430e-11;  // in m^3 kg^-1 s^-2
    CelestialBody selectedBody;

public:
    QListWidget* listWidget;
    MyOpenGLWindow(QWidget *parent = nullptr)
        : QOpenGLWidget(parent) {
listWidget = new QListWidget;
        debugTimer.start();
        // Add some celestial bodies
        scaleFactor = 1.0e15;  // The scale factor you used to scale down the universe

        // Constants for the Sun
        const double SUN_MASS = 1.989e30;  // in kilograms
        const double SUN_RADIUS = 6.9634e8;  // in meters
        const double EARTH_SUN_DISTANCE = 1.496e11;  // in meters

        // Create the Sun
        CelestialBody sun("Sun", SUN_MASS, SUN_RADIUS / scaleFactor, 0, 0, 0, 0, 0, 0, 0);
       // universe.bodyMap["Sun"] = sun;
        universe.addBody(sun);
        selectedBody = sun;

        // Calculate the Earth's initial velocity around the Sun

        double earth_position_scaled = EARTH_SUN_DISTANCE / scaleFactor;
         double moon_position_scaled = (EARTH_MOON_DISTANCE + EARTH_SUN_DISTANCE) / scaleFactor;
        // Calculate the Earth's initial velocity around the Sun
         // Calculate the Earth's initial velocity around the Sun
         double earth_velocity = std::sqrt(GRAVITATIONAL_CONSTANT * SUN_MASS / EARTH_SUN_DISTANCE);

         // Set the Earth's velocity and position
         CelestialBody earth("Earth", EARTH_MASS, EARTH_RADIUS / scaleFactor, earth_position_scaled, 0, 0, 0, earth_velocity / scaleFactor, 0, 1);
         universe.addBody(earth);

         // Calculate the Moon's initial velocity around the Earth
         double moon_velocity = std::sqrt(GRAVITATIONAL_CONSTANT * EARTH_MASS / EARTH_MOON_DISTANCE);

         // The Moon's total velocity is the sum of its velocity around the Earth and the Earth's velocity around the Sun
         double moon_total_velocity = moon_velocity /*+ earth_velocity*/;

         // Set the Moon's velocity and position
         CelestialBody moon("Moon", MOON_MASS, MOON_RADIUS / scaleFactor, moon_position_scaled, 0, 0, 0, moon_total_velocity / scaleFactor, 0, 2);
         universe.addBody(moon);


        // Calculate the distance from the center of the Earth
        double distance = EARTH_RADIUS + LOW_EARTH_ORBIT_ALTITUDE;

        // Scale down the distance
        double distance_scaled = distance / scaleFactor;

        // Calculate the initial position of the spacecraft (in low Earth orbit)
        double initial_position = EARTH_RADIUS + LOW_EARTH_ORBIT_ALTITUDE;

        // Add an offset to the initial position
        double offset = 2.0e8;  // Adjust this value as needed
        initial_position += offset;

        double initial_position_scaled = initial_position / scaleFactor;

        // Calculate the gravitational force from the Earth
        double force_earth = GRAVITATIONAL_CONSTANT * EARTH_MASS / std::pow(initial_position_scaled, 2);

        double velocity = std::sqrt(force_earth * initial_position_scaled); //todo use init pos or scaled pos

        double velocity_scaled = velocity / scaleFactor;

        double dx = initial_position_scaled - earth.x;
        double dy = 0 - earth.y;  // Assuming the Earth's y-coordinate is 0


        // Get a vector that's orthogonal to the displacement vector
        double vx = dy;
        double vy = -dx;

        // Normalize the velocity vector (make it have a length of 1)
        double length = std::sqrt(vx*vx + vy*vy);
        vx /= length;
        vy /= length;

        // Scale the velocity vector by the calculated velocity
        vx *= velocity_scaled;
        vy *= velocity_scaled;

        // Set the spacecraft's velocity
        //Spacecraft spacecraft(1000, 1, initial_position_scaled, 0, 0, vx, vy, 0);
        //universe.addSpacecraft(Spacecraft(1000, 0.5, initial_position_scaled, 0, 0, vx, vy, 0));

        if(debug){
            qDebug() << "spacecraft starting Distance (scaled):" << distance_scaled;
            qDebug() << "spacecraft starting vx, vy: " << vx << ", " << vy;
            qDebug() << "spacecraft starting Initial position:" << initial_position;
            qDebug() << "all objects Initial positions:" << sun.x << earth.x << moon.x;
        }
        // Add the spacecraft to the universe
        //universe.addSpacecraft(spacecraft);

    }


protected:

    void initializeGL() override {
        initializeOpenGLFunctions();

        // Enable depth testing
        glEnable(GL_DEPTH_TEST);

        // Set the clear color to black
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);  // Change to a gray color

        // Generate sphere vertices
        std::vector<GLfloat> sphereVertices = generateSphereVertices(1.0f, 16, 32);

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

        CelestialBody& sun = universe.bodyMap["Sun"];

        gluLookAt(cameraX, cameraY, cameraZ,  // Camera position
                  sun.x, sun.y, sun.z,  // Look-at point (average position)
                  0.0, 0.0, 1.0); // Up vector


    }

    void wheelEvent(QWheelEvent *event) override {
        // Get the number of degrees the wheel was rotated
        float degrees = event->angleDelta().y() / 8.0f;
        // Convert degrees to steps (a full step is 15 degrees)
        float steps = degrees / 15.0f;
        // Adjust the camera radius based on the number of steps
        cameraZ -= steps * 2000.0f;  // Adjust the factor as needed
        // Ensure the radius stays within a certain range
        cameraZ = std::max(1000.0f, std::min(1000000.0f, cameraZ));  // Adjust the range as needed
    }


    void mousePressEvent(QMouseEvent *event) override {
        if (event->button() == Qt::LeftButton) {
            mousePress = true;
            lastMousePos = event->pos();
        }
    }

    void mouseReleaseEvent(QMouseEvent *event) override {
        if (event->button() == Qt::LeftButton) {
            mousePress = false;
        }
    }

    void mouseMoveEvent(QMouseEvent *event) override {
        if (mousePress) {
            QPoint delta = event->pos() - lastMousePos;
            angleX += delta.x();
            angleY += delta.y();

            lastMousePos = event->pos();
        }
    }


    double calculateMaxDistance() {
        double maxDistance = 0.0;
        for (const CelestialBody& body1 : universe.bodies) {
            for (const CelestialBody& body2 : universe.bodies) {
                if (&body1 != &body2) {
                    double dx = body1.x - body2.x;
                    double dy = body1.y - body2.y;
                    double dz = body1.z - body2.z;
                    double distance = sqrt(dx*dx + dy*dy + dz*dz);
                    if (distance > maxDistance) {
                        maxDistance = distance;
                    }
                }
            }
        }
        return maxDistance;
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


    void paintGL() override {
        QPainter painter(this);

        // Clear the screen
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Set up the projection matrix
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        double fovy = 90.0;  // Field of view angle, in degrees, in the y direction
        double aspect = (double)this->width() / (double)this->height();  // Aspect ratio (width to height)
        double zNear = 0.01;  // Distance to near clipping plane
        double zFar = 10000000000000.0;  // Distance to far clipping plane

        gluPerspective(fovy, aspect, zNear, zFar);

        // Set up the model-view matrix
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();


        float dx = selectedBody.x - cameraX;
        float dy = selectedBody.y - cameraY;
        float dz = selectedBody.z - cameraZ;

        // The up vector will be a rotation of 90 degrees around the y-axis
        float upVectorX = dz;
        float upVectorY = dy;
        float upVectorZ = dx;

        gluLookAt(cameraX, cameraY, cameraZ,  // Camera position
                  selectedBody.x, selectedBody.y, selectedBody.z,  // Look-at point
                  upVectorX, upVectorY, upVectorZ); // Up vector


         glRotatef(angleX, 0.0f, 1.0f, 0.0f);
         glRotatef(angleY, 1.0f, 0.0f, 0.0f);

        // Draw celestial bodies and spacecrafts

        for (const CelestialBody& body : universe.bodies) {
            glPushMatrix();
            glTranslatef(body.x * scaleFactor, body.y * scaleFactor, body.z * scaleFactor);
            if(body.objectType == 0){
                glColor3f(1.0f, 1.0f, 0.0f);  // Set the color to red
            } else if(body.objectType == 1) {
                glColor3f(0.0f, 0.0f, 1.0f);  // Set the color to red
            } else if(body.objectType == 2){
                 glColor3f(1.0f, 1.0f, 1.0f);  // Set the color to red
            }
            drawSphere(body.radius * scaleFactor, 20, 20, body.objectType);  // Pass the radius of the celestial body
            glPopMatrix();
            glBegin(GL_LINE_STRIP);
            float alpha = 1.0f;  // Start with full opacity
            for (const QVector3D& pastPosition : body.pastPositions) {
                 glColor4f(1.0f, 1.0f, 1.0f, alpha);  // White color with alpha opacity
                 glVertex3f(pastPosition.x() * scaleFactor, pastPosition.y() * scaleFactor, pastPosition.z() * scaleFactor);
                 //alpha *= 0.95f;  // Decrease the opacity for the next segment
            }
            glEnd();
        }

        for (const Spacecraft& craft : universe.spacecrafts) {
            glPushMatrix();
            glTranslatef(craft.x * scaleFactor, craft.y * scaleFactor, craft.z * scaleFactor);
            glColor3f(0.0f, 1.0f, 0.0f);  // Set the color to green
            drawSphere(0.01 * scaleFactor, 20, 20, 3);  // Use a small constant value for the spacecraft size
            glPopMatrix();
        }




        // End using QPainter
        painter.end();
    }

    void drawSphere(double r, int lats, int longs, int objectType) {
        double visibleRadius = 0;
        if(objectType == 0){
            visibleRadius = 2000; //std::max(r, 0.1 * scaleFactor);  // Use a minimum radius for visibility
        } else if(objectType == 1) {
           visibleRadius = 18.3; //std::max(r, 0.1 * scaleFactor);  // Use a minimum radius for visibility
        } else if (objectType == 2){
           visibleRadius = 4.98;
        } else {
           visibleRadius = 1; //change this to a realistic value for a craft
        }
        glPushMatrix();  // Push the current matrix onto the stack
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

    void resizeGL(int w, int h) override {
        // Set the viewport to cover the entire widget
        glViewport(0, 0, w, h);

        // Set up the projection matrix
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        double fovy = 90.0;  // Field of view angle, in degrees, in the y direction
        double aspect = (double)w / (double)h;  // Aspect ratio (width to height)
        double zNear = 0.01;  // Distance to near clipping plane
        double zFar = 10000000000000.0;  // Distance to far clipping plane

        gluPerspective(fovy, aspect, zNear, zFar);

        // Calculate average position of all celestial bodies
        float dx = selectedBody.x - cameraX;
        float dy = selectedBody.y - cameraY;
        float dz = selectedBody.z - cameraZ;

        // The up vector will be a rotation of 90 degrees around the y-axis
        float upVectorX = -dz;
        float upVectorY = 0.0f;
        float upVectorZ = dx;

        gluLookAt(cameraX, cameraY, cameraZ,  // Camera position
                  selectedBody.x, selectedBody.y, selectedBody.z,  // Look-at point
                  upVectorX, upVectorY, upVectorZ); // Up vector
    }



public slots:
    void changeFocus(QListWidgetItem* current, QListWidgetItem* previous) {
        // Extract the body name from the item text
        if (current) {
            QString bodyName = current->text().split(":").first();
            // Find the body in the universe

            if (universe.bodyMap.count(bodyName.toStdString()) > 0) {

                selectedBody = universe.bodyMap[bodyName.toStdString()];
            }
        }
    }

    void adjustTimeScale(int value) {
        // Adjust the time scale based on the slider value
        // This is just an example; you'll need to replace this with your own logic
        timeScale = static_cast<double>(value) / 0.00005;
    }

    void updateUniverse() {
        // Update your universe here
        universe.update(timeScale/1200.0,debugTimer); // Assuming 60 frames per second

        // Calculate the scale factor
        double maxDistance = 0.0;
        for (const CelestialBody& body : universe.bodies) {
            for (const CelestialBody& otherBody : universe.bodies) {
                double distance = std::sqrt(std::pow(body.x - otherBody.x, 2) +
                                            std::pow(body.y - otherBody.y, 2) +
                                            std::pow(body.z - otherBody.z, 2));
                maxDistance = std::max(maxDistance, distance);
            }
        }
        double minScaleFactor = 0.01;  // Adjust as needed
        double maxScaleFactor = 1.0;  // Adjust as needed
        scaleFactor = this->width() / maxDistance;
        scaleFactor = std::max(minScaleFactor, std::min(maxScaleFactor, static_cast<double>(scaleFactor)));
        if (scaleFactor < 0.1) scaleFactor = 0.1;  // Set a lower limit for the scale factor

        if(debug && debugTimer.elapsed() % 500 < 16){
            for (const CelestialBody& body : universe.bodies) {
                qDebug() << "body pos:" << body.x << body.y << body.z;
            }
            for (const Spacecraft& craft : universe.spacecrafts) {
                qDebug() << "craft pos:" << craft.x << craft.y << craft.z;
            }
        }
        listWidget->clear();
        // Add the text for each celestial body
        for (const CelestialBody& body : universe.bodies) {
            QString text = QString("%1: Position (%2, %3, %4), Velocity (%5, %6, %7)")
                               .arg(QString::fromStdString(body.name)).arg(body.x).arg(body.y).arg(body.z)
                               .arg(body.vx).arg(body.vy).arg(body.vz);
            listWidget->addItem(text);
        }
        // Update the list of past positions for each celestial body
        for (CelestialBody& body : universe.bodies) {
            body.pastPositions.push_front(QVector3D(body.x, body.y, body.z));  // Add the current position to the front of the list
            while (body.pastPositions.size() > 1000) {  // Limit the size of the list
                body.pastPositions.pop_back();  // Remove the oldest position
            }
        }
        // Redraw the scene
        update();
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

    // Create the QSlider
    QSlider* timeSlider = new QSlider(Qt::Horizontal);
    timeSlider->setRange(1, 500);  // Adjust the range as needed
    timeSlider->setValue(1);


    openglWindow->listWidget->setFixedSize(600, 200);  // Adjust the width and height as needed

    // Add the MyOpenGLWindow to the layout, spanning 2 rows and 2 columns
    layout->addWidget(openglWindow, 0, 0, 2, 2);

    // Add the QSlider to the layout in the third row
    layout->addWidget(timeSlider, 2, 0, 1, 2);

    // Add the QListWidget to the layout in the top right corner
    layout->addWidget(openglWindow->listWidget, 0, 1, 1, 1);
    openglWindow->listWidget->setFixedSize(600, 200);  // Adjust the width and height as needed
    // Check if the signal is being emitted
    QObject::connect(openglWindow->listWidget, &QListWidget::currentItemChanged, openglWindow, &MyOpenGLWindow::changeFocus);


    // Connect the QSlider's valueChanged signal to a slot in MyOpenGLWindow
    // You'll need to add a slot named adjustTimeScale to MyOpenGLWindow for this to work
    QObject::connect(timeSlider, &QSlider::valueChanged, openglWindow, &MyOpenGLWindow::adjustTimeScale);

    // Show the main window
    window.show();

    // Set up the timer as before
    QTimer timer;
    QObject::connect(&timer, &QTimer::timeout, openglWindow, &MyOpenGLWindow::updateUniverse);
    timer.start(16);  // 60 frames per second

    return app.exec();
}



