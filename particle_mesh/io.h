#ifndef IO_H                 
#define IO_H

#include "physics.h"

typedef struct {
    double Vnorm, Mnorm, Lnorm, Tnorm;          //unità di misura

    double delta_t;     //passo temporale iniziale
    double t_final;     //tempo finale della simulazione 

    double G_scaled;  

    int N_p;             //numero particelle 
    int N_g;             //numero celle Fourier
    int N_tsteps_est;    //stima iniziale del numero di timestep

    double L_Box;                               //grandezza della box
    double delta_rho_max;                       //contrasto di densità massimo
} GeneralData;


void read_params(GeneralData *gd);  
Particle *allocate_particles(int N); 
void write_ic_text(const char *filename, Particle *part, int N);
void write_ic_binary(const char *filename, Particle *part, int N);
void write_histogram(const char *filename, Particle *part, int N, double L, int nbins);
Particle *read_ic_text(const char *filename, int *N);
Particle *read_ic_bin(const char *filename, int *N);
void write_snapshot(int step, Particle *part, int Np, double *rho, int Ng, double L);
void write_final_data(Particle *part, int Np, double *rho, double *pot, double *force, int Ng, double L);

#endif