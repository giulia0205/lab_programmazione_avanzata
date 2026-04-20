#include "physics.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stddef.h>

#define GAMMA 1.4
#define SIGMA_1D (2.0/3.0)
#define N_NGB     64
#define HALF_NGB  32
#define EPS       1.0e-12

/*CONDIZIONI INIZIALI SHOCK TUBE CON GHOST PARTICLES*/
void init_particles_shock_tube(Particella *p, int N, int Ng_sx, int Ng_dx, double xmin, double xmax) {
    const double x0 = 0.0;

    /*Stati iniziali richiesti*/
    const double rho_sx = 1.0;
    const double p_sx   = 1.0;
    const double v_sx   = 0.0;

    const double rho_dx = 0.125;
    const double p_dx   = 0.1;
    const double v_dx   = 0.0;

    /*Energia interna dall'equazione di stato*/
    const double u_sx = p_sx / ((GAMMA - 1.0) * rho_sx);
    const double u_dx = p_dx / ((GAMMA - 1.0) * rho_dx);

    /*Divisione del dominio fisico in parte sinistra e destra*/
    int n_sx = N / 2;
    int n_dx = N - n_sx;                 //definizione che tiene conto del caso in cui N è dispari

    double L_sx = x0 - xmin;
    double L_dx = xmax - x0;

    /*Spaziature uniformi nelle due metà fisiche*/
    double delta_x_sx = L_sx / (double)n_sx;
    double delta_x_dx = L_dx / (double)n_dx;

    /*Masse scelte per riprodurre le densità iniziali*/
    double m_sx = rho_sx * delta_x_sx;
    double m_dx = rho_dx * delta_x_dx;

    /*GHOST A SINISTRA*/
    for (int g = 0; g < Ng_sx; g++) {
        p[g].Mass     = m_sx;
        p[g].Pos      = xmin - (Ng_sx - g - 0.5) * delta_x_sx;
        p[g].Vel      = v_sx;
        p[g].Acc      = 0.0;
        p[g].Pressure = p_sx;
        p[g].Density  = rho_sx;
        p[g].U        = u_sx;
        p[g].dU       = 0.0;
        p[g].h        = EPS;
        p[g].ghost    = 1;
        p[g].next     = NULL;
    }

    /*ZONA FISICA SINISTRA*/
    for (int i = 0; i < n_sx; i++) {
        int k = Ng_sx + i;

        p[k].Mass     = m_sx;
        p[k].Pos      = xmin + (i + 0.5) * delta_x_sx;
        p[k].Vel      = v_sx;
        p[k].Acc      = 0.0;
        p[k].Pressure = p_sx;
        p[k].Density  = rho_sx;
        p[k].U        = u_sx;
        p[k].dU       = 0.0;
        p[k].h        = EPS;
        p[k].ghost    = 0;
        p[k].next     = NULL;
    }

    /*ZONA FISICA DESTRA*/
    for (int j = 0; j < n_dx; j++) {
        int k = Ng_sx + n_sx + j;

        p[k].Mass     = m_dx;
        p[k].Pos      = x0 + (j + 0.5) * delta_x_dx;
        p[k].Vel      = v_dx;
        p[k].Acc      = 0.0;
        p[k].Pressure = p_dx;
        p[k].Density  = rho_dx;
        p[k].U        = u_dx;
        p[k].dU       = 0.0;
        p[k].h        = EPS;
        p[k].ghost    = 0;
        p[k].next     = NULL;
    }

    /*GHOST A DESTRA*/
    for (int g = 0; g < Ng_dx; g++) {
        int k = Ng_sx + N + g;

        p[k].Mass     = m_dx;
        p[k].Pos      = xmax + (g + 0.5) * delta_x_dx;
        p[k].Vel      = v_dx;
        p[k].Acc      = 0.0;
        p[k].Pressure = p_dx;
        p[k].Density  = rho_dx;
        p[k].U        = u_dx;
        p[k].dU       = 0.0;
        p[k].h        = EPS;
        p[k].ghost    = 1;
        p[k].next     = NULL;
    }
}

double compute_h_mean(Particella *p, int N) {
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        sum += p[i].h;
    }
    return sum / (double)N;
}

/*Distanza assoluta r_ij = |x_i - x_j|*/
/*
double rij(double xi, double xj) {
    return fabs(xi - xj);
}

/*Differenza con segno x_i - x_j (serve per il gradiente)*/
/*
double dxij(double xi, double xj) {
    return xi - xj;
}
*/

/*Kernel cubico B-spline 1D*/
double kernel_cubic_1d(double r, double h) {
    /*Controllo sul valore di h*/
    if (h <= 0.0) {
        fprintf(stderr, "Errore: h <= 0 nel kernel\n");
        exit(EXIT_FAILURE);
    }

    double q = r / h;

    if (q >= 0.0 && q <= 1.0) {
        return (SIGMA_1D / h) * (1.0 - 1.5*q*q + 0.75*q*q*q);
    }
    else if (q > 1.0 && q <= 2.0) {
        double term = 2.0 - q;
        return (SIGMA_1D / h) * (0.25 * term * term * term);
    }
    else {
        return 0.0;
    }
}

/*Gradiente del kernel rispetto a x_i: dW/dx_i*/
double grad_kernel_cubic_1d(double xi, double xj, double h) {
    if (h <= 0.0) {
        fprintf(stderr, "Errore: h <= 0 nel gradiente del kernel\n");
        exit(EXIT_FAILURE);
    }

    double dx = xi - xj;
    double r  = fabs(dx);

    if (r == 0.0) {
        return 0.0;
    }

    double q = r / h;
    double sign;
    if (dx > 0.0) {
        sign = 1.0;
    } else {
        sign = -1.0;
    }
    double df_dq;

    if (q >= 0.0 && q <= 1.0) {
        df_dq = -3.0*q + 2.25*q*q;
    }
    else if (q > 1.0 && q <= 2.0) {
        df_dq = -0.75 * (2.0 - q) * (2.0 - q);
    }
    else {
        return 0.0;
    }

    return (SIGMA_1D / (h * h)) * df_dq * sign;
}

/*Media aritmetica delle smoothing lengths*/
double hij_mean(double hi, double hj) {
    return 0.5 * (hi + hj);
}

/*Equazione di stato*/
double equation_of_state(double rho, double u) {
    return (GAMMA - 1.0) * rho * u;
}

/*Calcolo della viscosità artificiale di Monaghan*/
double artificial_viscosity(const Particella *pi, const Particella *pj, double hij) {
    const double alpha = 1.0;               //controlla il termine lineare
    const double beta  = 2.0;               //controlla il termine quadratico
    const double eta   = 0.1;               //piccolo parametro di sicurezza per evitare divisioni problematiche quando due particelle sono molto vicine

    double xij = pi->Pos - pj->Pos;         //posizione relativa
    double vij = pi->Vel - pj->Vel;         //velocità relativa

    /*La viscosità artificiale si attiva solo in avvicinamento*/
    if (vij * xij >= 0.0) {
        return 0.0;
    }

    double rhoij = 0.5 * (pi->Density + pj->Density);       //densità media della coppia

    /*Velocità del suono relativa alle particelle*/
    double ci = sqrt(GAMMA * (GAMMA - 1.0) * pi->U);
    double cj = sqrt(GAMMA * (GAMMA - 1.0) * pj->U);
    double cij = 0.5 * (ci + cj);

    double muij = (hij * vij * xij) / (xij * xij + eta * eta * hij * hij);      //quantità ausiliaria

    return (-alpha * cij * muij + beta * muij * muij) / rhoij;
}

/*Calcolo della pressione attraverso l'equazione di stato*/
void compute_pressure(Particella *p, int N) {
    for (int i = 0; i < N; i++) {
        p[i].Pressure = equation_of_state(p[i].Density, p[i].U);
    }
}

/*LEAPFROG*/
void leapfrog_half_kick(Particella *p, int N, double dt) {
    for (int i = 0; i < N; i++) {
        p[i].Vel += 0.5 * dt * p[i].Acc;
        p[i].U   += 0.5 * dt * p[i].dU;
    }
}

void leapfrog_drift(Particella *p, int N, double dt) {
    for (int i = 0; i < N; i++) {
        p[i].Pos += dt * p[i].Vel;
    }
}
