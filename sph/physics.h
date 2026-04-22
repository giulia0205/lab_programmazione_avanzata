#ifndef PHYSICS_H
#define PHYSICS_H

typedef struct Particella {
    double Mass;
    double Pos;
    double Vel;
    double Acc;
    double Pressure;
    double Density;
    double U;
    double dU;
    double h;
    int ghost;   /* 0 = particella fisica, 1 = ghost */
    struct Particella *next;
} Particella;

void init_particles_shock_tube(Particella *p, int N, int Ng_sx, int Ng_dx, double xmin, double xmax);

double compute_h_mean(Particella *p, int N);
double kernel_cubic_1d(double r, double h);
double grad_kernel_cubic_1d(double xi, double xj, double h);
double hij_mean(double hi, double hj);
double equation_of_state(double rho, double u);
double artificial_viscosity(const Particella *pi, const Particella *pj, double hij);

void compute_pressure(Particella *p, int N) ;

void leapfrog_half_kick(Particella *p, int N, double dt);
void leapfrog_drift(Particella *p, int N, double dt);

#endif