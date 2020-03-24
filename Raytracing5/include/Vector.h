#pragma once
#include <math.h>
#include <vector>

class Vector {
    public:
        Vector(double x=0, double y=0, double z=0) {
            coord[0] =x;
            coord[1] =y;
            coord[2] =z;
        }
    //creation d'un accesseur constant et non constant
    const double& operator[](int i) const { return coord[i];}
    double& operator[](int i) { return coord[i];}

    //fonction qui renvoi la norme au carré d'un vecteur
    double getNorm2() {
        return coord[0]*coord[0] + coord[1]*coord[1] + coord[2]*coord[2];
    }
    //fonction qui calcule la norme dun vecteur
    void normalize() {
        double norm = sqrt(getNorm2());
        coord[0] /= norm;
        coord[1] /= norm;
        coord[2] /= norm;
    }
    Vector getNormalized() {
        Vector result(*this);
        result.normalize();
        return result;
    }

    Vector& operator+=(const Vector& b) {
        coord[0] += b[0];
        coord[1] += b[1];
        coord[2] += b[2];
        return *this;
    }

    private:
        double coord[3];
};

//addition de deux vecteurs
Vector operator+(const Vector& a, const Vector& b);
Vector operator-(const Vector& a, const Vector& b);
Vector operator*(double a, const Vector& b) ;
Vector operator*(const Vector& a, const Vector& b);
Vector operator*(const Vector& a,double b) ;
Vector operator/(const Vector& a, double b);
double dot( const Vector& a, const Vector& b);
Vector cross( const Vector& a, const Vector& b);

Vector random_cos(const Vector& N);


//classe qui va gerer les rayons
class Ray {
public:
    //un rayon a un origine et une direction
    Ray(const Vector& o, const Vector& d) : origin(o), direction(d) {};
    Vector origin, direction;
};



