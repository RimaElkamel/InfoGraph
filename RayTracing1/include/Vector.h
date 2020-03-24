#ifndef VECTOR_H
#define VECTOR_H
#include <math.h>
#include <vector>
class Vector {
 public :
     Vector(double x=0, double y=0, double z=0){
         coords[0] = x;
         coords[1] = y;
         coords[2] = z;
     }
      double operator[](int i) const {
        return coords[i];
     }
     double operator[](int i){
        return coords[i];
     }

     //La norme au carré d'un vecteur (pour faire des comparaisons ultérieures)
     double getNormSquare() const{
         return coords[0] * coords[0] + coords[1] * coords[1] + coords[2] * coords[2] ;
     }
     void operator*=(double a){
        coords[0]*=a;
        coords[1]*=a;
        coords[2]*=a;
     }
     void operator/=(double a){
        coords[0]/=a;
        coords[1]/=a;
        coords[2]/=a;
     }

     //Calcul la norme d'un vecteur
     void normalize(){
         double aux = sqrt(getNormSquare());
         coords[0] /= aux;
         coords[1] /= aux;
         coords[2] /= aux;
     }

     Vector getNormalized(){
        Vector result(*this);
        result.normalize();
        return result;
     }

     Vector& operator+=(const Vector& b){
         coords[0]+=b[0];
         coords[1]+=b[1];
         coords[2]+=b[2];
         return *this;
     }
private :
     double coords[3];
 };

 Vector operator+(const Vector& A, const Vector& B );

 Vector operator-(const Vector& A, const Vector& B );

 Vector operator*(double k, const Vector& B );
 Vector operator*(const Vector& B , double k);
 Vector operator*(const Vector& B , const Vector& A);

 Vector operator/(const Vector& B , double k);

 //Produit scalaire entre deux vecteurs
  double dot(const Vector& A, const Vector& B);


//Classe pour les rayons
class Ray{
public :
     Ray(const Vector& C, Vector u):C(C), u(u){};
     //C: origine du vecteur et u est la direction du vecteur
     Vector C, u;
 };

class Sphere{
public:
    Sphere(const Vector& O, double R,const Vector &c): O(O),R(R),colorS(c){};

    //Tester s'il y a une intersection entre un rayon et une sphère
bool if_intersect(const Ray& r, Vector& P, Vector& N, double& t){
    // resout l esuqation du second degré
    double a=1;
    double b= 2*dot(r.u, r.C - O);
    double Rcarre = R*R;
    double c= (r.C - O).getNormSquare()-Rcarre;

    double delta = b*b - 4*a*c;
    if (delta<0) return false;

    double racineDelta = sqrt(delta);
    double t1 = (-b-racineDelta)/(2*a);
    double t2 = (-b+racineDelta)/(2*a);

    if (t2 < 0) return false;

    if(t1>0) t=t1;
    else t=t2;

    P = r.C + t* r.u;
    N=(P-O).getNormalized();
    return true;
}

    //O: origine de la sphère et R est le rayon de la sphère
    // colorS : la couleur de la sphère
    Vector O;
    double R;
    Vector colorS;

};

class Scene{
public:
    Scene(){};
    void ajouterSphere(const Sphere& s) { spheres.push_back(s);}

    //sphere_indice pour renvoyer l indice de la sphère interceptée
    bool if_intersect(const Ray& r, Vector& P, Vector& N, int &sphere_indice){
        bool has_inter = false;
        double min_t= 1E99;
        for(int i=0; i<spheres.size(); i++)
        {
            Vector P1,N1;
            double t;
            bool has_inter1 = spheres[i].if_intersect(r,P1,N1,t);
            if(has_inter1){
                has_inter=true;
                if(t<min_t){
                    min_t=t;
                    P = P1;
                    N = N1;
                    sphere_indice = i;
                }
            }
        }
        return has_inter;
    }

    std::vector<Sphere> spheres;
};
#endif // VECTOR_H
