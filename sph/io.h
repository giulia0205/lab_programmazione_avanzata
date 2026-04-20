#ifndef IO_H
#define IO_H

#include "physics.h"

typedef struct {

    int Np;        
    double L_grid;       
    double dt;
    double t_end;

} GeneralData;

void read_params(GeneralData *gd);

void write_file(const char *filename, Particella *p, int N);

#endif
