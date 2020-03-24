#include <vector>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <random>
#include "Geometry.h"
#define M_PI 3.1415926535897932
#include "omp.h"
std::default_random_engine engine;
std::uniform_real_distribution<double> unifom(0,1);



void save_image(const char* filename, const unsigned char* pixels, int W, int H) { // stored as RGB

    unsigned char bmpfileheader[14] = { 'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0 };
    unsigned char bmpinfoheader[40] = { 40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0 };
    unsigned char bmppad[3] = { 0,0,0 };

    int filesize = 54 + 3* W*H;
    bmpfileheader[2] = (unsigned char)(filesize);
    bmpfileheader[3] = (unsigned char)(filesize >> 8);
    bmpfileheader[4] = (unsigned char)(filesize >> 16);
    bmpfileheader[5] = (unsigned char)(filesize >> 24);

    bmpinfoheader[4] = (unsigned char)(W);
    bmpinfoheader[5] = (unsigned char)(W >> 8);
    bmpinfoheader[6] = (unsigned char)(W >> 16);
    bmpinfoheader[7] = (unsigned char)(W >> 24);
    bmpinfoheader[8] = (unsigned char)(H);
    bmpinfoheader[9] = (unsigned char)(H >> 8);
    bmpinfoheader[10] = (unsigned char)(H >> 16);
    bmpinfoheader[11] = (unsigned char)(H >> 24);

    FILE* f;
    f = fopen(filename, "wb");
    fwrite(bmpfileheader, 1, 14, f);
    fwrite(bmpinfoheader, 1, 40, f);
    std::vector<unsigned char> bgr_pixels(W*H * 3);
    for (int i = 0; i < W*H; i++) {
        bgr_pixels[i * 3] = pixels[i*3 + 2];
        bgr_pixels[i * 3+1] = pixels[i*3 + 1];
        bgr_pixels[i * 3+2] = pixels[i*3];
    }
    for(int i =0; i<H; i++) {
        fwrite(&bgr_pixels[0] + (W*(H - i -1) * 3), 3, W, f);
        fwrite(bmppad, 1, (4 - (W * 3) % 4) % 4, f);
    }
    fclose(f);
}

Vector renduCouleur(const Ray& r, const Scene& s, int nbrebonds) {
    if(nbrebonds == 0) return Vector(0, 0, 0);
    Vector P, N;
    int sphere_id;
    double t;
    bool has_inter = s.intersection(r, P, N, sphere_id, t);

    Vector intensite_pixel(0, 0, 0);
    if(has_inter) {
        if(s.objects[sphere_id]->miroir) {
            Vector direction_miroir = r.direction - 2*dot(N, r.direction)*N;
            Ray rayon_miroir(P+0.001*N, direction_miroir);
            intensite_pixel = renduCouleur(rayon_miroir, s, nbrebonds - 1);
        }  else if(s.objects[sphere_id]->transparent) {
            double n1 = 1;
            double n2 = 1.3;
            Vector normale_pour_transporence(N);
            if(dot(r.direction, N)>0) {//on sort de la sphere
                n1 = 1.3;
                n2 = 1;
                normale_pour_transporence = -1*N;
            }
            double radical = 1 - pow(n1/n2,2)* (1 - pow(dot(normale_pour_transporence, r.direction), 2));
            if(radical>0) {
                Vector direction_refracte = (n1/n2)*(r.direction - dot(r.direction, normale_pour_transporence)*normale_pour_transporence) - normale_pour_transporence* sqrt(radical) ;
                Ray rayon_refracte(P - 0.01*normale_pour_transporence, direction_refracte);
                intensite_pixel = renduCouleur(rayon_refracte, s, nbrebonds - 1);
            }
        }
        else {
            Vector axeOP = (P - s.lumiere->O).getNormalized();
            Vector dir_aleatoire = random_cos(axeOP);
            Vector point_aleatoire = dir_aleatoire * s.lumiere->R + s.lumiere->O ;
            Vector wi = (point_aleatoire -P ).getNormalized();
            double d_light2 = (point_aleatoire - P).getNorm2();
            Vector Np = dir_aleatoire;

            Ray ray_light(P + 0.01*N, wi);
            Vector P_light, N_light;
            int sphere_id_light;
            double t_light;
            bool has_inter_light = s.intersection(ray_light, P_light, N_light, sphere_id_light, t_light);

            if(has_inter_light && t_light*t_light<d_light2*0.99) {
                intensite_pixel = Vector(0,0,0);
            } else {
                intensite_pixel = (s.intensite_lumiere / (4*M_PI*d_light2)* std::max(0., dot(N, wi))* dot(Np, -1*wi) / dot(axeOP, dir_aleatoire)) * s.objects[sphere_id]->albedo;
            }
            Vector dir_alea = random_cos(N);
            Ray rayon_aleatoire(P+0.001*N, dir_alea);
            intensite_pixel += renduCouleur(rayon_aleatoire, s, nbrebonds - 1) * s.objects[sphere_id]->albedo;
        }
    }

    return intensite_pixel;
}

int main()
{
    int W = 64;
    int H = 64;
    const int nrays = 80;
    double fov = 60 * M_PI / 180;

    //creation d'une source lumineuse
    Sphere slum(Vector(15, 70, -30), 15, Vector(1., 1., 1.));
    Sphere s2(Vector(0, -2000 - 20, 0), 2000, Vector(1, 1, 1)); //sol de couleur blanc
    Sphere s3(Vector(0, 2000 + 100, 0), 2000, Vector(1, 1, 1)); //plafond de couleur blanc
    Sphere s4(Vector(-2000 - 50, 0 , 0), 2000, Vector(0, 1, 0)); //mur gauche de couleur blanc
    Sphere s5(Vector(2000 + 50, 0,  0), 2000, Vector(1, 0, 0)); //mur droit de couleur rouge
    Sphere s6(Vector(0, 0, -2000 - 100), 2000, Vector(0, 1, 1)); //mur au fond de couleur cyant

    Geometry g("img/BeautifulGirl.obj", 10, Vector(0, 0, -55), Vector(1., 1., 1.));

    Scene s;
    s.addSphere(slum);

    s.addGeometry(g);
    s.addSphere(s2);
    s.addSphere(s3);
    s.addSphere(s4);
    s.addSphere(s5);
    s.addSphere(s6);
    s.lumiere = &slum;
    s.intensite_lumiere = 1000000000000;
    Vector pos_lum(0.,0.,0.);
    double focus_distance = 55;
    double aperture = 0.5;


    std::vector<unsigned char> image(W*H * 3);

    #pragma cmp parallel for
    for(int i = 0; i<H; i++) {
        for(int j =0; j<W; j++) {
            printf("(%d, %d )", i,j);
            Vector color(0., 0., 0.);
            for(int k=0; k<nrays; k++) {

                double alea1 = unifom(engine);
                double alea2 = unifom(engine);
                double R = sqrt(-2*log(alea1));

                Vector direction(j - W/2 + 0.5 + R*cos(2*M_PI*alea2), i - H/2 + 0.5 + R*sin(2*M_PI*alea2), -W/(2*tan(fov/2)));
                direction.normalize();

                Vector destination = pos_lum + focus_distance * direction;
                Vector new_origin = pos_lum + Vector((unifom(engine) - 0.5) * aperture, (unifom(engine) - 0.5) * aperture, 0);

                Ray r(new_origin, (destination - new_origin).getNormalized());

                color += renduCouleur(r, s, 1)/nrays;

            }

           image[((H - i - 1)*W + j) * 3 + 0] = std::min(255., std::max(0., std::pow(color[0], 1/2.2)));
           image[((H - i - 1)*W + j) * 3 + 1] = std::min(255., std::max(0., std::pow(color[1], 1/2.2)));
           image[((H - i - 1)*W + j) * 3 + 2] = std::min(255., std::max(0., std::pow(color[2], 1/2.2)));
    }}
    save_image("img1.png", &image[0], W, H);

    return 0;
}
