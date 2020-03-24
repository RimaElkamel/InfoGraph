#include <omp.h>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <limits>
#include <math.h>
#include <random>

#define MATH_PI 3.1415922653589

std::default_random_engine engine;
std::uniform_real_distribution<double> distrib(0,1);

inline double sqr(double x){
    return x*x;
}

void save_image(const char* filename, const unsigned char* tableau, int w, int h) { // (0,0) is top-left corner

    FILE *f;

    int filesize = 54 + 3 * w*h;

    unsigned char bmpfileheader[14] = { 'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0 };
    unsigned char bmpinfoheader[40] = { 40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0 };
    unsigned char bmppad[3] = { 0,0,0 };

    bmpfileheader[2] = (unsigned char)(filesize);
    bmpfileheader[3] = (unsigned char)(filesize >> 8);
    bmpfileheader[4] = (unsigned char)(filesize >> 16);
    bmpfileheader[5] = (unsigned char)(filesize >> 24);

    bmpinfoheader[4] = (unsigned char)(w);
    bmpinfoheader[5] = (unsigned char)(w >> 8);
    bmpinfoheader[6] = (unsigned char)(w >> 16);
    bmpinfoheader[7] = (unsigned char)(w >> 24);
    bmpinfoheader[8] = (unsigned char)(h);
    bmpinfoheader[9] = (unsigned char)(h >> 8);
    bmpinfoheader[10] = (unsigned char)(h >> 16);
    bmpinfoheader[11] = (unsigned char)(h >> 24);

    f = fopen(filename, "wb");
    fwrite(bmpfileheader, 1, 14, f);
    fwrite(bmpinfoheader, 1, 40, f);
    unsigned char *row = new unsigned char[w * 3];
    for (int i = 0; i < h; i++)
    {
        for (int j = 0; j < w; j++) {
            row[j * 3] = tableau[(w*(h - i - 1) * 3) + j * 3+2];
            row[j * 3+1] = tableau[(w*(h - i - 1) * 3) + j * 3+1];
            row[j * 3+2] = tableau[(w*(h - i - 1) * 3) + j * 3];
        }
        fwrite(row, 3, w, f);
        fwrite(bmppad, 1, (4 - (w * 3) % 4) % 4, f);
    }
    fclose(f);
    delete[] row;
}

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
     void normalize(){
         double aux = sqrt(getNormSquare());
         coords[0] /= aux;
         coords[1] /= aux;
         coords[2] /= aux;
     }
     Vector& operator+=(const Vector& b){
         coords[0]+=b[0];
         coords[1]+=b[1];
         coords[2]+=b[2];
         return *this;
     }
     double coords[3];
 };
 Vector operator+(const Vector& A, const Vector& B ){
     return Vector(A[0]+B[0], A[1]+B[1], A[2]+B[2]);
 }
 Vector operator-(const Vector& A, const Vector& B ){
     return Vector(A[0]-B[0], A[1]-B[1], A[2]-B[2]);
 }

 Vector operator*(double k, const Vector& B ){
     return Vector(k*B[0], k*B[1], k*B[2]);
 }
 Vector operator*(const Vector& B , double k){
     return Vector(k*B[0], k*B[1], k*B[2]);
 }
 Vector operator*(const Vector& B , const Vector& A){
     return Vector(A[0]*B[0], A[1]*B[1], A[2]*B[2]);
 }
 Vector operator/(const Vector& B , double k){
     return Vector(B[0]/k, B[1]/k, B[2]/k);
 }

 double dot(const Vector& A, const Vector& B){
     return A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
 }

 Vector cr(const Vector& A, const Vector& B){
     return Vector(A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2], A[0] * B[1] - A[1] * B[0]);
 }

Vector random_cos(const Vector& N){
     //double alea1 = distrib(engine[omp_get_thread_num()]);
     double alea1 = distrib(engine);
     //double alea2 = distrib(engine[omp_get_thread_num()]);
     double alea2 = distrib(engine);
     Vector alea_dir_L(cos(2 * MATH_PI * alea1) * sqrt(1 - alea2), sin(2 * MATH_PI * alea1) * sqrt(1-alea2),sqrt(alea2));
     Vector aleatoir(distrib(engine)+0.5,distrib(engine)+0.5,distrib(engine)+0.5);
     Vector tang1 = cr(N,aleatoir);
     tang1.normalize();
     Vector tang2 = cr(tang1, N);
     return alea_dir_L[2] * N + alea_dir_L[0] * tang1 + alea_dir_L[1] * tang2;
}


class Ray{
public :
     Ray(const Vector& C, Vector u):C(C), u(u){};
     Vector C, u;
 };

class Object{
 public :
        Object(){};
 	virtual bool intersect(const Ray& d, Vector& P, Vector& N, double& t) const =0;

        Vector albedo;
        bool mirro,transp;
};


class Sphere : public Object{
public:
    Sphere(const Vector& O, double R, const Vector& albedo, bool mirro = false, bool transp = false): O(O),R(R){this->albedo = albedo;
this->mirro = mirro;
this->transp = transp;};


    bool intersect(const Ray& r, Vector& P, Vector& N, double& t) const {
        //2*t^2+b*t+c
        double a = 1;
        double b = 2 * dot(r.u, r.C - O);
        double c =(r.C-O).getNormSquare() -R*R;

        double delta = b*b - 4*a*c;
        if(delta<0) return false;
        double sqrtDelta = sqrt(delta);
        double t1 = (-b + sqrtDelta)/(2*a);
        if(t1<0) return false;

        double t0 = (-b - sqrtDelta)/(2*a);
        if(t0 > 0){
            t = t0;
        }else{
            t = t1;
        }
        P = r.C + t*r.u;
        N = (P - O);
	N.normalize();
        return true;
    }
    Vector O;
    double R;

};

class Triangle : public Object{
public :
Vector A, B, C;
Triangle(const Vector& A, const Vector& B, const Vector& C, const Vector& albedo, bool mirro = false, bool transp = false) :A(A), B(B), C(C){
	this->albedo = albedo;
	this->mirro = mirro;
	this->transp = transp;
};

 bool intersect(const Ray& d, Vector& P, Vector& N, double& t) const {
	N = cr(B - A, C - A);
        N.normalize();
        t = dot(C-d.C,N) / dot(d.u,N);
        if(t<0) return false;
        P = d.C + t*d.u;
        Vector u = B-A;
        Vector v = C-A;
        Vector w = P-A;
        double m11 = u.getNormSquare();
        double m12 = dot(u,v);
        double m22 = v.getNormSquare();
        double detm = m11*m22 - m12*m12;
        double b11 = dot(w,u);
        double b21 = dot(w, v);
        double detb=b11*m22 - b21*m12;
        double beta = detb/detm;
        double g12 = b11;
        double g22 = b21;
        double detg = m11*g22 - m12*g12;
        double gamma = detg/detm;
        double alpha = 1-beta-gamma;
        if(alpha<0 || alpha>1) return false;
        if(beta<0 || beta>1) return false;
        if(gamma<0 || gamma>1) return false;

        return true;
}

};

class Scene{
public :
    Scene(){};

    void ajouterSphere(const Sphere& s){
        objects.push_back(&s);
    }
    void ajouterTriangle(const Triangle& s){
        objects.push_back(&s);
    }

    bool intersect(const Ray& r, Vector& P, Vector& N, int& indice, double& t) const {
        bool has_inter = false;
        t = std::numeric_limits<double>::max();

        for(int i=0; i<objects.size();i++){
           Vector Pp, Np;
           double tp;
           bool interlocal = objects[i]->intersect(r,Pp, Np, tp);
           if(interlocal){
                has_inter = true;
                if(tp < t){
                    t = tp;
                    P = Pp;
                    N = Np;
                    indice = i;
                }
           }
        }
       return has_inter;
    }

    Vector renduCouleur(const Ray& r,int numrebond){
        if(numrebond<0) return Vector(0.,0.,0.);
        Vector P, N;
        int indice_sphere;
	double t;
        bool has_intersection = intersect(r, P, N, indice_sphere,t);

	Vector intensite_pixel=Vector(0.,0.,0.);
        if(has_intersection){
            if(indice_sphere == 0){
                   return L->albedo * intensiteL / (4* MATH_PI * L->R * L->R);
            }

            if(objects[indice_sphere]->mirro){
                Vector direction_miroir = r.u - 2*dot(r.u,N)*N;
                Ray rray_miroir(P+0.0001*N, direction_miroir);
                intensite_pixel = renduCouleur(rray_miroir, numrebond-1);
            }

            else {
		  if(objects[indice_sphere]->transp){
		        double n1=1;
		        double n2=1.4;
		        Vector NforTransp(N);
		        if(dot(r.u,N)>0){
		            std::swap(n1,n2);
		            NforTransp = -1*N;
		        }
		        double delta = 1-sqr(n1/n2)*(1-sqr(dot(r.u,NforTransp)));
                         if(delta>0){
                               Vector direction_refracte = (n1/n2)*(r.u - dot(r.u, NforTransp)*NforTransp) - NforTransp*sqrt(delta);
                               Ray rayon_refracte(P - 0.0001*NforTransp,direction_refracte);
                               intensite_pixel = renduCouleur(rayon_refracte, numrebond-1);
                          }

		    }
	            else{
                             //eclairage indirect
                             Vector axePO = (P - L->O);
                             axePO.normalize();
                             Vector dir_aleatoire = random_cos(axePO);
                             Vector point_aleatoire = dir_aleatoire * L->R + L->O;
                             Vector wi = (point_aleatoire - P);
                             wi.normalize();
                             double d_lightp = (point_aleatoire - P).getNormSquare();
                             Vector Np = dir_aleatoire;


			    Ray rayon_lum(P+0.01*N, wi);
			    Vector Pprime, Nprime;
			    int indiceprime;
			    double tprime;
                            bool has_inter = intersect(rayon_lum,Pprime, Nprime, indiceprime,tprime);

                              if(has_inter && tprime*tprime < d_lightp*0.99){
				 intensite_pixel = Vector(0.,0.,0.);
			    }else{
                             intensite_pixel = (intensiteL / (4*MATH_PI*d_lightp)*std::max(0.,dot(N,wi))*dot(Np,-1*wi) /dot(axePO, dir_aleatoire))* objects[indice_sphere]->albedo ;               }

                             Vector direction_aleatoir = random_cos(N);
                             Ray rray_aleatoir(P+0.0001*N, direction_aleatoir);
                             intensite_pixel += renduCouleur(rray_aleatoir, numrebond-1) * objects[indice_sphere]->albedo;

		       }
        	}
	}
        return intensite_pixel;
    }
    std::vector<const Object*> objects;
    Sphere* L ;
    double intensiteL;
};

int main()
{
    int W = 512;
    int H = 512;

    Scene s;
    Sphere s_lumiere(Vector(15, 70, -30),15,Vector(1.,1.,1.));
    Sphere s1(Vector(0., 0., -55.), 10, Vector(1.,1.,1.));
    //Sphere s2(Vector(-15., 0., -35.), 10, Vector(1.,1.,1.),false,true);
    //Sphere s3(Vector(15., 0., -75.), 10, Vector(1.,1.,1.),true);
    Sphere ssol(Vector(0., -2000-20, 0.), 2000, Vector(1.,1.,1.));
    Sphere splafond(Vector(0., 2000+100, 0.), 2000, Vector(1.,1.,1.));
    Sphere smurgauche(Vector(-2000-50, 0., 0.), 2000, Vector(0.,1.,0.));
    Sphere smurdroit(Vector(2000+50, 0., 0.), 2000, Vector(0.,0.,1.));
    Sphere smurfond(Vector(0., 0., -2000-100), 2000, Vector(0.,1.,1.));

    Triangle tri(Vector(-10,-10,-20),Vector(10,-10,-20),Vector(10,10,-20),Vector(1,0,1));
    s.ajouterSphere(s_lumiere);
    s.ajouterSphere(s1);
    //s.ajouterSphere(s2);
    //s.ajouterSphere(s3);
    s.ajouterSphere(ssol);
    s.ajouterSphere(splafond);
    s.ajouterSphere(smurgauche);
    s.ajouterSphere(smurdroit);
    s.ajouterSphere(smurfond);
    s.ajouterTriangle(tri);

    s.L = &s_lumiere;
    s.intensiteL = 10000000000;

    double alpha = 60 * MATH_PI /180;
    double d = W / (2 * tan(alpha/2.));
    const int nbr_ray = 50;
    Vector pos_lum(0.,0.,0.);
    double focus_distance = 55;
    double aperture=5.;
    std::vector<unsigned char> img(W*H * 3, 0);
#pragma cmp parallel for schedule(dynamic,1)
    for(int i = 0 ; i<H ; i++){
        for(int j = 0; j < W; j++){

                Vector I = Vector(0.,0.,0.);

		for(int k =0;k<nbr_ray;k++){
                        //Box Muller
                        double alea1 = distrib(engine);
                        double alea2 = distrib(engine);
                        double R = sqrt(-2 * log(alea1));
                        Vector u(j-W / 2 + 0.5 + R * cos(2 * MATH_PI * alea2), -i+H /2 + 0.5 + R * sin(2 * MATH_PI * alea2), -d);
                        u.normalize();
                        Vector destination = pos_lum + focus_distance * u;
                        Vector new_origin = pos_lum + Vector((distrib(engine)-0.5)*aperture,(distrib(engine)-0.5)*aperture,0);
                        Vector new_direction = (destination - new_origin);
                        new_direction.normalize();
                        Ray r(new_origin, new_direction);
                 	I += s.renduCouleur(r, 2)/nbr_ray;
                }

                img[(i*W + j)*3 + 0] = std::min(255., pow(I[0], 0.45));
                img[(i*W + j)*3 + 1] = std::min(255., pow(I[1], 0.45));
                img[(i*W + j)*3 + 2] = std::min(255., pow(I[2], 0.45));
           }
        }
    save_image("imagetr.png", &img[0], W, H);

    return 0;
}
