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
    /*Controllo sulla lettura da terminale*/
    if (argc < 2) {
        printf("Uso: %s SUN oppure TRAPPIST1 skip\n", argv[0]);
        return 1;
    }

    /*Lettura dei parametri*/
    GeneralData all = {0};

    if (stamp(argv[1], "SUN")) {
        read_params("params_sun.txt", &all);
    } 
    else if (stamp(argv[1], "TRAPPIST1")) {
        read_params("params_trappist1.txt", &all);
    } 
    else {
        printf("Errore: scrivere SUN oppure TRAPPIST1\n");
        return 1;
    }

    PlanetState *planets = malloc(all.Np * sizeof(PlanetState));
    if (!planets) return 1;

    init_planets(&all, planets);

    /*APRIMAO I FILE*/
    FILE **files = open_binary_files(all.Np);                             

    /*INTEGRAZIONE DI EULERO*/
    //Calcolo del periodo a partire dal pianeta più distante
    double r_max = all.r[all.Np-1];
    double T = 2 * M_PI * sqrt(pow(r_max,3) / (all.G_scaled * all.Ms));
    int nsteps = (int)(T / all.delta_t);            
    for (int n = 0; n < nsteps; n++) {          
        for (int i = 0; i < all.Np; i++) {
            //Scrittura del file binario
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
    close_binary_files(all.Np, files);

    free(all.r);
    free(planets);

    return 0;
}
