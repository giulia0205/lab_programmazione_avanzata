#ifndef PHYSICS_H
#define PHYSICS_H

typedef struct {
    double x;     
    double v;   
    double mass;   
} Particle;

/*Distribuzione di particelle*/
double uniform_rand();
double pdf(double x, double L, double A);
void generate_distribution(Particle *part, int N, double L, double A);

/*Densità*/
void density(double *rho, Particle *part, int Np, int Ng, double L);

/*Potenziale*/
void compute_potential_fft(double *rho, double *pot, int Ng, double L);

/*Forza*/
void compute_force_grid(double *force, double *pot, int Ng, double L);
void force_to_particles(double *acc, double *force, Particle *part, int Np, int Ng, double L);

/*Leapfrog*/
double compute_timestep(Particle *part, int Np, int Ng, double L);
void leapfrog_start(Particle *part, double *acc, int Np, double dt);
void leapfrog_drift(Particle *part, int Np, double dt, double L);
void leapfrog_kick(Particle *part, double *acc, int Np, double dt);

#endif