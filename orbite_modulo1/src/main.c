#include "utils.h"
#include "physics.h"
#include "io.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int stamp(char *arg, char *mode) {
    return strcmp(arg, mode) == 0;
}

int main(int argc, char *argv[]) {

    /*Lettura variabili da terminale*/
    double M_s= atof(argv[1]);
    int Np= atoi(argv[2]);
    double delta_t=atof(argv[3]);           //perchè delta_t è blu?

    /*Funzione di lettura dei parametri*/
    GeneralData all;
    read_params(&all);

    all.Ms = M_s;
    all.Np = Np;
    all.delta_t = delta_t;
    /*Allocazione di memoria e lettura dei raggi*/
    all.r = malloc(Np * sizeof(double));
    if (!all.r) return 1;
    for (int i = 0; i < Np; i++) {
        all.r[i] = atof(argv[4 + i]);
    }

    PlanetState *planets = malloc(Np * sizeof(PlanetState));
    if (!planets) return 1;

    init_planets(&all, planets);

    /*APRIMAO I FILE*/
    FILE **files = open_binary_files(Np);                             //!!!!!!!

    /*INTEGRAZIONE DI EULERO*/
    //Calcolo del periodo a partire dal pianeta più distante
    double r_max = all.r[Np-1];
    double T = 2 * M_PI * sqrt(pow(r_max,3) / (all.G_scaled * M_s));
    int nsteps = (int)(T / delta_t);            
    for (int n = 0; n < nsteps; n++) {          // era nstep-1 ma chat dice di toglierlo
        for (int i = 0; i < Np; i++) {
            //qua si scrive il file binario
            OutputRecord rec;
            rec.t = n * all.delta_t;
            rec.x = planets[i].x;
            rec.y = planets[i].y;
            rec.vx = planets[i].vx;
            rec.vy = planets[i].vy;
            rec.ax = planets[i].ax;
            rec.ay = planets[i].ay;
            write_record(files[i], rec);
        }
        Eulero(&all, planets);
    }


    /*CHIUDIAMO I FILE*/
    close_binary_files(Np, files);

    free(all.r);
    free(planets);

    return 0;
}
