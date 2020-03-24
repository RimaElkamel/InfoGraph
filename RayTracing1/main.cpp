using namespace std;
#include <iostream>
#include <stdio.h>
#include <vector>
#include "Vector.h"

#define MATH_PI 3.1415922653589

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


int main()
{
    int W = 1024;
    int H = 1024;

    Sphere sphere1(Vector(0,0,-55),5, Vector(0,1,0));
    Sphere sphere2(Vector(0,20,-55),5, Vector(0,0,1));
    Sphere sphere3(Vector(0,10,-55),5, Vector(1,0,1));
    Sphere sphere4(Vector(0,-10,-55),5, Vector(1,1,0));
    Sphere sphereAuSol(Vector(0,-2000-20,0),2000, Vector(1,1,1));
    Sphere sphereAuplafond(Vector(0,2000+100,0),2000, Vector(1,1,1));
    Sphere sphereGauche(Vector(-2000-50,0,0),2000, Vector(0,1,0));
    Sphere sphereDroite(Vector(2000+50,0,0),2000, Vector(0,0,1));
    Sphere sphereAuFond(Vector(0,0,-2000-100),2000, Vector(0,1,1));


    Scene scene;
    scene.ajouterSphere(sphere1);
    scene.ajouterSphere(sphere2);
    scene.ajouterSphere(sphere3);
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
    for(int i = 0 ; i<H ; i++){
        for(int j = 0; j < W; j++){

        Vector u(j-W/2, i-H/2, -W/(2*tan((60* MATH_PI/180)/2)));
        u.normalize();

        Ray r(Vector(0,0,0),u);
        Vector P, N;
        int sphere_indice;
        bool has_inter = scene.if_intersect(r, P, N, sphere_indice);

        Vector inten_pix(0,0,0);
        if(has_inter){
        Vector color_of_sphere= scene.spheres[sphere_indice].colorS;
        inten_pix = color_of_sphere * intensite_lum*std::max(0.,dot((pos_lum-P).getNormalized(),N))/ (pos_lum-P).getNormSquare();

         //img[(i*W + j)*3 + 0] = 255;
        //composante verte du pixel
         //img[(i*W + j)*3 + 1] = 255;
        //composante bleue du pixel
        // img[(i*W + j)*3 + 2] = 255;
        }
        double int_pixel_1= std::max(0.,inten_pix[0]);
        double int_pixel_2= std::max(0.,inten_pix[1]);
        double int_pixel_3= std::max(0.,inten_pix[2]);
        //composante rouge du pixel
        img[((H-i-1)*W + j)*3 + 0] = std::min(255.,int_pixel_1);
        //composante verte du pixel
         img[((H-i-1)*W + j)*3 + 1] = std::min(255.,int_pixel_2);
        //composante bleue du pixel
         img[((H-i-1)*W + j)*3 + 2]= std::min(255.,int_pixel_3);
        }
        }
save_image("image4p.png", &img[0], W, H);
    return 0;
}
