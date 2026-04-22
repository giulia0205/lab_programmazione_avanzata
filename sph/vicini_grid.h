#ifndef VICINI_GRID_H
#define VICINI_GRID_H

#include "physics.h"

typedef struct {
    Particella **head;   /* array di teste */
    int ncell;
    double xmin;
    double xmax;
    double delta;
    double L;
} Grid1D;


int cell_index(double x, const Grid1D *g);
void build_grid(Grid1D *g, Particella *p, int N);

void compute_h_grid(Particella *p, int N, const Grid1D *g);
void compute_density_grid(Particella *p, int N, const Grid1D *g);
void compute_acc_dU_grid(Particella *p, int N, const Grid1D *g);

#endif