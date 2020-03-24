using namespace std;
#include <iostream>
#include <stdio.h>
#include <vector>
#include "Vector.h"

//Pour l'utilisation de pow
#include <cmath>
#define MATH_PI 3.1415922653589

//Pour la génération des nombres aléatoires
#include <random>
std::default_random_engine engine;
std::uniform_real_distribution<double> distrib(0,1);
//Fonction pour sauvegarder une image déjà donnée
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

Vector renduCouleur(const Ray &r, const Scene& s, int nombreRebond){
        if (nombreRebond == 0) return Vector(0,0,0);
        Vector pos_lum(30,60, -50);
        double intensite_lum = 1000000000;
        Vector P, N;
        int sphere_indice;
        double t;
        bool has_inter = s.if_intersect(r, P, N, sphere_indice,t);
        Vector inten_pix(0.,0.,0.);
        if(has_inter){
                // Si la sphère a un caractére miroir
                if(s.spheres[sphere_indice].if_miroir){
                    Vector miroir_direction = r.u - 2*dot(N, r.u)*N;
                    Ray miroir_ray(P + 0.001*N, miroir_direction);
                    inten_pix =  renduCouleur(miroir_ray ,s, nombreRebond - 1);
                }
                else if (s.spheres[sphere_indice].if_transparent){
                double n1 = 1;
                double n2=1.3;
                Vector normal_transparent(N);
                //tester si le rayon est en train d entrer ou sortir de la sphère
                if(dot(r.u,N)>0){
                    // on sort
                    n1 = 1.3;
                    n2 = 1;
                    normal_transparent = -1*N;
                }
                // Utilisation de la formule de réfraction
                double carreRef = (n1/n2)*(n1/n2);
                double carreDot = dot(normal_transparent,r.u)*dot(normal_transparent,r.u);
                double res = 1-carreRef*(1-carreDot);
                if(res>0){
                    Vector transparent_direction = (n1/n2)*(r.u-dot(r.u,normal_transparent)*normal_transparent)- normal_transparent * (sqrt(res));
                //Rayon dans la direction réfractée pour avoir la propriété de transparence
                    Ray transparent_ray(P - 0.001*normal_transparent, transparent_direction);
                    inten_pix =  renduCouleur(transparent_ray ,s, nombreRebond - 1);
                }
                }

        else{
                // eclairage direct
                Ray rayon_ombre(P+0.01*N, (pos_lum-P).getNormalized());
                Vector P_ombre, N_ombre;
                int sphere_indice_ombre;
                double t_ombre;
                bool has_inter_lum = s.if_intersect(rayon_ombre,P_ombre,N_ombre,sphere_indice_ombre, t_ombre);
                double distance_lumiere = (pos_lum-P).getNormSquare();
                if(has_inter_lum &&(t_ombre*t_ombre)< distance_lumiere){
                inten_pix = Vector(0,0,0);
                }
                else{
                Vector color_of_sphere= s.spheres[sphere_indice].colorS/MATH_PI;
                Vector num = color_of_sphere * intensite_lum*std::max(0.,dot((pos_lum-P).getNormalized(),N));
                double denom = (pos_lum-P).getNormSquare();
                inten_pix = num/ denom;

                 //img[(i*W + j)*3 + 0] = 255;
                //composante verte du pixel
                 //img[(i*W + j)*3 + 1] = 255;
                //composante bleue du pixel
                // img[(i*W + j)*3 + 2] = 255;
                }
                //Partie eclairage indirect :: repère local et global
                double alea1 = distrib(engine);
                double alea2 = distrib(engine);

                Vector alea_direction_L(cos(2*MATH_PI*alea1)*sqrt(1-alea2),sin(2*MATH_PI*alea1)*sqrt(1-alea2),sqrt(alea2));
                Vector alea_vector(distrib(engine)-0.5,distrib(engine)-0.5,distrib(engine)-0.5);
                Vector tang1 = cr(N, alea_vector);
                tang1.normalize();
                Vector tang2 = cr(tang1, N);

                Vector alea_direction_G = alea_direction_L[0]*tang1 + alea_direction_L[1]*tang2;


                Ray alea_ray(P + 0.001*N, alea_direction_G);
                inten_pix +=  renduCouleur(alea_ray ,s, nombreRebond - 1) * s.spheres[sphere_indice].colorS;

                }
        }
        return inten_pix;
}
int main()
{
    int W = 1024;
    int H = 1024;
    const int nbr_ray = 80;
    //Sphere sphere1(Vector(0,0,-55),7, Vector(1,0,1),true);
    Sphere sphere1(Vector(0,0,-55),5, Vector(0,1,0));
    Sphere s2(Vector(-15., 0., -55.), 10, Vector(1.,1.,1.));
    Sphere sphere2(Vector(0,22,-55),5, Vector(0,0,1));
    Sphere sphere3(Vector(0,12,-55),5, Vector(1,0,1));
    Sphere sphere4(Vector(15,0,-55),10, Vector(1,0,0));
    Sphere sphereAuSol(Vector(0,-2000-20,0),2000, Vector(1,1,1));
    Sphere sphereAuplafond(Vector(0,2000+100,0),2000, Vector(1.,1.,1.));
    Sphere sphereGauche(Vector(-2000-50,0,0),2000, Vector(0.,1.,0.));
    Sphere sphereDroite(Vector(2000+50,0,0),2000, Vector(0.,0.,1.));
    Sphere sphereAuFond(Vector(0,0,-2000-100),2000, Vector(0.,1.,1.));

    Scene scene;
    //scene.ajouterSphere(sphere1);
    scene.ajouterSphere(s2);
    //scene.ajouterSphere(sphere2);
    //scene.ajouterSphere(sphere3);
    scene.ajouterSphere(sphere4);
    scene.ajouterSphere(sphereAuSol);
    scene.ajouterSphere(sphereAuplafond);
    scene.ajouterSphere(sphereGauche);
    scene.ajouterSphere(sphereDroite);
    scene.ajouterSphere(sphereAuFond);

    Vector pos_lum(15,30, -40);
    double intensite_lum = 1000000;

    std::vector<unsigned char> img(W*H * 3);

    //parcourir les pixels de l'image
#pragma omp parallel for
    for(int i = 0 ; i<H ; i++){
        for(int j = 0; j < W; j++){

        Vector u(j-W/2, i-H/2, -W/(2*tan((60* MATH_PI/180)/2)));
        u.normalize();

        Ray r(Vector(0,0,0),u);


        Vector couleur(0,0,0);
        //pour avoir une image mojns bruitée
        for(int k =0;k<nbr_ray;k++)
            couleur+=renduCouleur(r, scene,5)/nbr_ray;

        double int_pixel_1= std::max(0.,std::pow(couleur[0],1/2.2));
        double int_pixel_2= std::max(0.,std::pow(couleur[1],1/2.2));
        double int_pixel_3= std::max(0.,std::pow(couleur[2],1/2.2));
        //composante rouge du pixel
        img[((H-i-1)*W + j)*3 + 0] = std::min(255.,int_pixel_1);
        //composante verte du pixel
         img[((H-i-1)*W + j)*3 + 1] = std::min(255.,int_pixel_2);
        //composante bleue du pixel
         img[((H-i-1)*W + j)*3 + 2]= std::min(255.,int_pixel_3);
        }
        }
save_image("image2.png", &img[0], W, H);
    return 0;
}
