#include "io.h"
#include "physics.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>


void read_params(GeneralData *gd) {
    FILE *fp = fopen("params.txt", "r");        
    if (!fp) {                                  
        printf("Errore: impossibile aprire params.txt\n");
        exit(1);
    }
    char name[64];                  
    double value;                  

    while (fscanf(fp, "%s %lf", name, &value) == 2) {           //scanf cerca di leggere una stringa e un numero, se riesce a leggere entrambi il ciclo continua

        if (strcmp(name, "N_particles") == 0) {                         
            gd->Np = value;
        }
        else if (strcmp(name, "L_grid") == 0) {
            gd->L_grid = value;
        }
        else if (strcmp(name, "dt") == 0) {
            gd->dt = value;
        }
        else if (strcmp(name, "t_end") == 0) {
            gd->t_end = value;
        }
        
    }
    fclose(fp);                                 
}

void write_file(const char *filename, Particella *p, int N) {
    FILE *fp = fopen(filename, "w");

    if (fp == NULL) {
        printf("Errore apertura file %s\n", filename);
        return;
    }

    for (int i = 0; i < N; i++) {
        fprintf(fp, "%f %f %f %f %f %f\n", p[i].Pos, p[i].Vel, p[i].Acc, p[i].Density, p[i].Pressure, p[i].U);
    }         

    fclose(fp);
}
