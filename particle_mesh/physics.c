#include "physics.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <fftw3.h>
#define PI 3.14159265358979323846

/*DISTRIBUZIONE INIZIALE DELLE PARTICELLE*/
double uniform_rand() {
    return rand() / (RAND_MAX + 1.0);
}

double pdf(double x, double L, double A) {
    return 1.0 - A * cos(2.0 * PI * x / L);
}

void generate_distribution(Particle *part, int N, double L, double A) {
    /*Controllo che la pdf resti sempre non negativa*/
    if (fabs(A) > 1.0) {
        printf("Errore: A_deltaPar deve soddisfare |A| <= 1\n");
        exit(1);
    }

    int count = 0;
    double pmax = 1.0 + fabs(A);

    while (count < N) {
        double x_try = L * uniform_rand();
        double y_try = pmax * uniform_rand();

        if (y_try < pdf(x_try, L, A)) {
            part[count].x = x_try;
            part[count].v = 0.0;
            part[count].mass = L / (double)N;
            count++;
        }
    }
}

/*DENSITÀ*/
/*Densità NGP*/
#ifdef NGP
void density(double *rho, Particle *part, int Np, int Ng, double L) {
    double dx = L / Ng;
    /*Condizione iniziale di densità nulla*/
    for (int i = 0; i < Ng; i++) {
        rho[i] = 0.0;
    }

    for (int p = 0; p < Np; p++) {
        double x = part[p].x;
        /*Periodocità della box*/
        while (x < 0.0) x += L;
        while (x >= L) x -= L;
        /*Ricerca della cella più vicina*/
        int cell = (int)(x / dx);
        /*Calcolo della massa totale in ciascuna cella*/
        rho[cell] += part[p].mass;
    }
    /*Calcolo della densità a partire dalla massa totale in ciascuna cella*/
    for (int i = 0; i < Ng; i++) {
        rho[i] /= dx;
    }
}
#endif

/*Densità CIC*/
#ifdef CIC
void density(double *rho, Particle *part, int Np, int Ng, double L) {
    double dx = L / Ng;
    /*Condizione iniziale di densità nulla*/
    for (int i = 0; i < Ng; i++) {
        rho[i] = 0.0;
    }

    for (int p = 0; p < Np; p++) {
        double x = part[p].x;
        /*Periodicità della box*/
        while (x < 0.0) x += L;
        while (x >= L) x -= L;
        /*Costruzione della griglia centrata sui centri delle celle*/
        double u = x / dx - 0.5;            //centro della cella i-esima
        int i0 = (int)floor(u);             //floor() prende la parte intera inferiore, quindi i0 è la cella di sinistra
        double frac = u - i0;               //parte frazionaria che dice quanto siamo vicini al centro della cella
        /*Indici delle celle coinvolte*/
        int i_left = i0;
        int i_right = i0 + 1;
        /*Periodicità sugli indici*/
        if (i_left < 0) i_left += Ng;
        if (i_left >= Ng) i_left -= Ng;

        if (i_right < 0) i_right += Ng;
        if (i_right >= Ng) i_right -= Ng;
        /*Pesi sulla massa*/
        double w_left = 1.0 - frac;
        double w_right = frac;
        /*Calcolo delle masse totali sulle due celle*/
        rho[i_left] += part[p].mass * w_left;
        rho[i_right] += part[p].mass * w_right;
    }
    /*Calcolo della densità a partire dalla massa totale in ciascuna cella*/
    for (int i = 0; i < Ng; i++) {
        rho[i] /= dx;
    }
}
#endif

/*Densità TSC*/
#ifdef TSC
void density(double *rho, Particle *part, int Np, int Ng, double L) {
    double dx = L / Ng;
    /*Condizione iniziale di densità nulla*/
    for (int i = 0; i < Ng; i++) {
        rho[i] = 0.0;
    }

    for (int p = 0; p < Np; p++) {
        double x = part[p].x;
        /*Periodicità della box*/
        while (x < 0.0) x += L;
        while (x >= L) x -= L;
        /*Costruzione della griglia centrata sui centri delle celle*/
        double u = x / dx - 0.5;
        int j = (int)floor(u + 0.5);        //ricerca dell'indice della cella centrale
        double s = u - j;                   //distanza della particella dal centro della cella espressa in unità di cella
        /*Indici delle celle coinvolte*/
        int i_left = j - 1;
        int i_center = j;
        int i_right = j + 1;
        /*Periodicità sugli indici*/
        if (i_left < 0) i_left += Ng;
        if (i_left >= Ng) i_left -= Ng;

        if (i_center < 0) i_center += Ng;
        if (i_center >= Ng) i_center -= Ng;

        if (i_right < 0) i_right += Ng;
        if (i_right >= Ng) i_right -= Ng;
        /*Calcolo dei pesi per la massa*/
        double w_left = 0.5 * (0.5 - s) * (0.5 - s);
        double w_center = 0.75 - s * s;
        double w_right = 0.5 * (0.5 + s) * (0.5 + s);
        /*Calcolo della massa totale in ciascuna cella*/
        rho[i_left] += part[p].mass * w_left;
        rho[i_center] += part[p].mass * w_center;
        rho[i_right] += part[p].mass * w_right;
    }
    /*Calcolo della densità a partire dalla massa totale in ciascuna cella*/
    for (int i = 0; i < Ng; i++) {
        rho[i] /= dx;
    }
}
#endif

/*CALCOLO DEL POTENZIALE*/
void compute_potential_fft(double *rho, double *pot, int Ng, double L) {
    fftw_complex *kDensity, *kPot;
    fftw_plan fft_fwd, fft_bck;

    kDensity = fftw_alloc_complex(Ng);
    kPot = fftw_alloc_complex(Ng);

    if (kDensity == NULL || kPot == NULL) {
        printf("Errore allocazione FFTW\n");
        exit(1);
    }

    fft_fwd = fftw_plan_dft_r2c_1d(Ng, rho, kDensity, FFTW_ESTIMATE);
    fft_bck = fftw_plan_dft_c2r_1d(Ng, kPot, pot, FFTW_ESTIMATE);

    /*FFT diretta della densità*/
    fftw_execute(fft_fwd);

    /*Modo k=0*/
    kPot[0][0] = 0.0;
    kPot[0][1] = 0.0;

    /*Fattore fondamentale dei numeri d'onda*/
    double knorm = 2.0 * PI / L;

    /*Potenziale in Fourier*/
    for (int i = 1; i < Ng/2 + 1; i++) {
        double k = i * knorm;

        kPot[i][0] = -kDensity[i][0] / (k * k);
        kPot[i][1] = -kDensity[i][1] / (k * k);
    }

    /*FFT inversa*/
    fftw_execute(fft_bck);

    /*Normalizzazione*/
    for (int i = 0; i < Ng; i++) {
        pot[i] /= Ng;
    }

    fftw_destroy_plan(fft_fwd);
    fftw_destroy_plan(fft_bck);
    fftw_free(kDensity);
    fftw_free(kPot);
}


/*CALCOLO DELLA FORZA*/
/*Calcolo della forza sulla griglia a partire dal potenziale*/
void compute_force_grid(double *force, double *pot, int Ng, double L) {
    double dx = L / Ng;

    for (int i = 0; i < Ng; i++) {
        int ip = i + 1;
        int im = i - 1;

        /*Periodicità degli indici*/
        if (ip >= Ng) ip -= Ng;
        if (im < 0) im += Ng;

        force[i] = -(pot[ip] - pot[im]) / (2.0 * dx);
    }
}

/*Calcolo della forza su ciascuna particella attraverso NGP*/
#ifdef NGP
void force_to_particles(double *acc, double *force, Particle *part, int Np, int Ng, double L) {
    double dx = L / Ng;

    for (int p = 0; p < Np; p++) {
        double x = part[p].x;

        while (x < 0.0) x += L;
        while (x >= L) x -= L;

        int cell = (int)(x / dx);

        acc[p] = force[cell];
    }
}
#endif

/*Calcolo della forza su ciascuna particella attraverso CIC*/
#ifdef CIC
void force_to_particles(double *acc, double *force, Particle *part, int Np, int Ng, double L) {
    double dx = L / Ng;

    for (int p = 0; p < Np; p++) {
        double x = part[p].x;

        while (x < 0.0) x += L;
        while (x >= L) x -= L;

        double u = x / dx - 0.5;
        int i0 = (int)floor(u);
        double frac = u - i0;

        int i_left = i0;
        int i_right = i0 + 1;

        if (i_left < 0) i_left += Ng;
        if (i_left >= Ng) i_left -= Ng;

        if (i_right < 0) i_right += Ng;
        if (i_right >= Ng) i_right -= Ng;

        double w_left = 1.0 - frac;
        double w_right = frac;

        acc[p] = w_left * force[i_left] + w_right * force[i_right];
    }
}
#endif

/*Calcolo della forza su ciascuna particella attraverso TSC*/
#ifdef TSC
void force_to_particles(double *acc, double *force, Particle *part, int Np, int Ng, double L) {
    double dx = L / Ng;

    for (int p = 0; p < Np; p++) {
        double x = part[p].x;

        while (x < 0.0) x += L;
        while (x >= L) x -= L;

        double u = x / dx - 0.5;
        int j = (int)floor(u + 0.5);
        double s = u - j;

        int i_left = j - 1;
        int i_center = j;
        int i_right = j + 1;

        if (i_left < 0) i_left += Ng;
        if (i_left >= Ng) i_left -= Ng;

        if (i_center < 0) i_center += Ng;
        if (i_center >= Ng) i_center -= Ng;

        if (i_right < 0) i_right += Ng;
        if (i_right >= Ng) i_right -= Ng;

        double w_left = 0.5 * (0.5 - s) * (0.5 - s);
        double w_center = 0.75 - s * s;
        double w_right = 0.5 * (0.5 + s) * (0.5 + s);

        acc[p] = w_left * force[i_left] + w_center * force[i_center] + w_right * force[i_right];
    }
}
#endif

/*DELTA_T*/
double compute_timestep(Particle *part, int Np, int Ng, double L) {
    double dx = L / Ng;
    double vmax = 0.0;
    double w = 0.5;

    for (int p = 0; p < Np; p++) {
        double vabs = fabs(part[p].v);
        if (vabs > vmax) vmax = vabs;
    }

    /*Caso iniziale: velocità nulle*/
    if (vmax == 0.0) {
        return w * dx;
    }

    return w * dx / vmax;
}

/*LEAPFROG*/
void leapfrog_drift(Particle *part, int Np, double dt, double L) {
    for (int p = 0; p < Np; p++) {
        part[p].x += part[p].v * dt;

        /*Periodicità*/
        while (part[p].x < 0.0) part[p].x += L;
        while (part[p].x >= L) part[p].x -= L;
    }
}

void leapfrog_kick(Particle *part, double *acc, int Np, double dt) {
    for (int p = 0; p < Np; p++) {
        part[p].v += acc[p] * dt * 0.5;
    }
}