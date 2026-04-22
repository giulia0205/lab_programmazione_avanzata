#ifndef IO_H
#define IO_H

#include <stdio.h>

typedef struct {
    double Mnorm, Lnorm, Tnorm;

    double Ms;        // massa della stella
    double delta_t;   // passo temporale
    double G_scaled;  // costante gravitazionale nelle unità interne

    int Np;
    double *r;
} GeneralData;

//dichiarazione della struct che rappresenta un record da salvare nel file binario
typedef struct {
    double t, x, y, vx, vy, ax, ay;   
} OutputRecord;

//dichiarazione delle funzioni di io.c
void read_params(const char *filename, GeneralData *gd);
FILE **open_binary_files(int Np);
void write_record(FILE *file, OutputRecord rec);
void close_binary_files(int Np, FILE **files);

#endif