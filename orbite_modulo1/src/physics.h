#ifndef PHYSICS_H
#define PHYSICS_H

#include "utils.h"

typedef struct {
    double x, y;     
    double vx, vy;   
    double ax, ay;   
} PlanetState;

void compute_acceleration(const GeneralData *gd, PlanetState *p);
void init_planets(const GeneralData *gd, PlanetState *planets);
void Eulero(const GeneralData *gd, PlanetState *planets);

#endif