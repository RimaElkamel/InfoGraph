#include "Vector.h"
#include <random>

#define M_PI 3.1415926535897932

static std::default_random_engine engine;
static std::uniform_real_distribution<double> unifom(0,1);

Vector operator+(const Vector& a, const Vector& b) {
    return Vector(a[0]+b[0], a[1]+b[1], a[2]+b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
    return Vector(a[0]-b[0], a[1]-b[1], a[2]-b[2]);
}
Vector operator*(double a, const Vector& b) {
    return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const Vector& b) {
    return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}
Vector operator*(const Vector& a,double b) {
    return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator/(const Vector& a, double b) {
    return Vector(a[0]/b, a[1]/b, a[2]/b);
}
double dot( const Vector& a, const Vector& b) {  //produit scalaire entre deux nombres
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
Vector cross( const Vector& a, const Vector& b) {
    return Vector(a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0] );
}
Vector random_cos(const Vector& N) {
    double r1 = unifom(engine);
    double r2 = unifom(engine);
    Vector direction_aleatoire_repere_local(cos(2*M_PI*r1)*sqrt(1 - r2),sin(2*M_PI*r1)*sqrt(1 - r2), sqrt(r2)) ;
    Vector aleatoire(unifom(engine) - 0.5, unifom(engine) - 0.5, unifom(engine) - 0.5);
    Vector tangeante1 = cross(N, aleatoire); tangeante1.normalize();
    Vector tangeante2 = cross(tangeante1, N);

    return direction_aleatoire_repere_local[2]*N + direction_aleatoire_repere_local[0]*tangeante1 + direction_aleatoire_repere_local[1]*tangeante2;

}


