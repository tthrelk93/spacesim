#ifndef VECTOR3_H
#define VECTOR3_H

class Vector3 {
public:
    double x, y, z;

    Vector3();
    Vector3(double x, double y, double z);

    Vector3 operator+(const Vector3& v) const;
    Vector3 operator-(const Vector3& v) const;
    Vector3 operator*(double scalar) const;
    Vector3 operator/(double scalar) const;

    double magnitude() const;
};

#endif // VECTOR3_H
