#include "physics.h"            
#include <math.h>


//ATTENZIONE: PERCHÈ A VOLTE *P E ALTRE * PLANETS?????? non importa è solo il puntatore

/*Calcolo dell'accelerazione*/
void compute_acceleration(const GeneralData *gd, PlanetState *p) {
    double r2 = p->x * p->x + p->y * p->y;
    double r  = pow(r2,0.5);
    p->ax = - gd->G_scaled * gd->Ms * p->x / (r*r*r);
    p->ay = - gd->G_scaled * gd->Ms * p->y / (r*r*r);
}

/*CONDIZIONI INIZIALI: Trasformazione dei dati di input in stato fisico adatto all'integrazione*/
void init_planets(const GeneralData *gd, PlanetState *planets) {
    for (int i = 0; i < gd->Np; i++) {
        double r = gd->r[i];
        planets[i].x  = -r;
        planets[i].y  = 0.0;
        planets[i].vx = 0.0;
        planets[i].vy = pow(gd->G_scaled * gd->Ms / r, 0.5);
        compute_acceleration(gd, &planets[i]);             
    }
}

/*Step per l'integrazione di Eulero*/
void Eulero(const GeneralData *gd, PlanetState *planets) {
    for (int i = 0; i < gd->Np; i++) {
        PlanetState *p = &planets[i];          
        double vy_old = p->vy;
        double vx_old = p->vx;
        compute_acceleration(gd, p);
        p->vx = vx_old + p->ax * gd->delta_t;
        p->vy = vy_old + p->ay * gd->delta_t;
        p->x  = p->x + vx_old * gd->delta_t;
        p->y  = p->y + vy_old * gd->delta_t;

    }
}



