/*
Integrate the evolution of the cosmological scale factor a(t).
Let's remember that, given a cosmological model, the evolution of the scale factor is given by the Friedmann-Lemaitre-Robertson-Walker equations.
Generalizing the equation for a Dark Energy component: ((da/dt)*(1/a))^2 = H_0^2 * ((1+z)^3 * Omega_m + (1+z)^2 * Omega_k + (1+z)^(3*(1+w)) * Omega_DE)
        z is the cosmological redshift defined by: 1+z = a_0/a, where a_0=1 when z=0;   
        Omega_k = 1 - Omega_m - Omega_DE
        Omega_m and Omega_DE are the matter and Dark Energy density parameters at z=0.
        H_0 is the Hubble constant at z=0 (use 70 km/s/Mpc; note that this determines the UdM of the equation);
        w gives the equation of state of the Dark Energy, p = w*rho: w=-1 for a LambdaCDM cosmological model.
Your code should accept all parameters from the command line.
To integrate a(t), use both your midpoint code, and the above adaptive-timestep routine. 
Try a number of different cosmological model; integrate a(t) between t=0 and t=14 Gyr. 
Use both codes and try to guess the difference in accuracy/computing time.
Make plots.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h" 

/*Definizione delle costanti*/
double t_max = 14.0;          // Gyr
double H_0_kmsMpc = 70.0;     // km/s/Mpc
int nsteps = 500000;

/*Definizione di variabili per adaptive step*/
int kmax, kount;
float dxsav, *xp, **yp;

/*Variabili globali per il modello cosmologico*/
double H_0, Omega_m, Omega_DE, w;

/*Conversione H0 = 70 km/s/Mpc in 1/Gyr */
double H_0_Gyr(double H_0_kmsMpc) {
    double Mpc_km = 3.085677581e19;   // km in un Mpc
    double sec_in_Gyr = 3.1536e16;    // secondi in un Gyr
    double H0 = H_0_kmsMpc / Mpc_km;  // s^-1
    return H0 * sec_in_Gyr;           // 1/Gyr
}

double Omega_k(double Omega_m, double Omega_DE) {
    return 1.0 - Omega_m - Omega_DE;
}

double da_dt(double H_0, double Omega_m, double Omega_DE, double w, double a) {
    double rad;

    if (a <= 0.0) {
        fprintf(stderr, "Errore: a <= 0 non fisico.\n");
        exit(EXIT_FAILURE);
    }

    rad = Omega_m * pow(a, -3.0)
        + Omega_k(Omega_m, Omega_DE) * pow(a, -2.0)
        + Omega_DE * pow(a, -3.0 * (1.0 + w));

    if (rad < 0.0) {
        fprintf(stderr, "Errore: radicando negativo per a = %g\n", a);
        exit(EXIT_FAILURE);
    }

    return a * H_0 * sqrt(rad);
}

/*Prototipi richiesti da odeint*/
void odeint(float ystart[], int nvar, float x1, float x2, float eps,
            float h1, float hmin, int *nok, int *nbad,
            void (*derivs)(float, float [], float []),
            void (*rkqs)(float [], float [], int, float *, float, float,
                         float [], float *, float *,
                         void (*)(float, float [], float [])));

void rkqs(float y[], float dydx[], int n, float *x, float htry, float eps,
          float yscal[], float *hdid, float *hnext,
          void (*derivs)(float, float [], float []));

void rkck(float y[], float dydx[], int n, float x, float h,
          float yout[], float yerr[],
          void (*derivs)(float, float [], float []));

/*RK4 standard per passo adattivo semplificato*/
void rk4_step(double *y, double *t, double h, double H0, double Omega_m, double Omega_DE, double w) {
    double k1 = da_dt(H0, Omega_m, Omega_DE, w, *y);
    double k2 = da_dt(H0, Omega_m, Omega_DE, w, *y + 0.5*h*k1);
    double k3 = da_dt(H0, Omega_m, Omega_DE, w, *y + 0.5*h*k2);
    double k4 = da_dt(H0, Omega_m, Omega_DE, w, *y + h*k3);
    *y += h/6.0 * (k1 + 2.0*k2 + 2.0*k3 + k4);
    *t += h;
}

void derivs(float t, float y[], float dydt[]) {
    double a, rad;
    (void)t;

    a = y[1];

    if (a <= 0.0 || !isfinite(a)) {
        fprintf(stderr, "Errore in derivs: a non fisico = %.15e\n", a);
        exit(EXIT_FAILURE);
    }

    rad = Omega_m * pow(a, -3.0)
        + (1.0 - Omega_m - Omega_DE) * pow(a, -2.0)
        + Omega_DE * pow(a, -3.0 * (1.0 + w));

    if (rad <= 0.0 || !isfinite(rad)) {
        fprintf(stderr, "Errore in derivs: rad non fisico = %.15e, a = %.15e\n", rad, a);
        exit(EXIT_FAILURE);
    }

    dydt[1] = (float)(a * H_0 * sqrt(rad));
}

int stamp(const char *arg, const char *mode) {
    return strcmp(arg, mode) == 0;
}

int main(int argc, char *argv[]) {

    if (argc < 5) {
        fprintf(stderr, "Uso: %s Omega_m Omega_DE w METODO\n", argv[0]);
        fprintf(stderr, "METODO = PUNTO_MEDIO oppure PASSO_ADATTIVO\n");
        return 1;
    }

    /*i parametri vengono letti da riga di comando*/
    Omega_m = atof(argv[1]);
    Omega_DE = atof(argv[2]);
    w = atof(argv[3]); 
    H_0 = H_0_Gyr(H_0_kmsMpc);
    double a_init = 1.0e-4;

    /*da riga di comando leggiamo anche quale metodo utilizzare*/

    if (stamp(argv[4], "PUNTO_MEDIO")) {
        /* con nsteps punti, il passo corretto è t_max/(nsteps-1) */
        double h = t_max / (double)(nsteps - 1); 
        double t = 0.0;
        double *a = malloc(nsteps * sizeof(double));

        if (a == NULL) {
            fprintf(stderr, "Errore allocazione memoria.\n");
            return 1;
        }

        a[0] = a_init;

        for (int n = 0; n < nsteps - 1; n++) {              
            double k1 = da_dt(H_0, Omega_m, Omega_DE, w, a[n]);
            double a_mid = a[n] + 0.5 * h * k1;             
            double k_mid = da_dt(H_0, Omega_m, Omega_DE, w, a_mid);
            a[n+1] = a[n] + h * k_mid;
        }

        /* scrittura del file di output */
        FILE *f = fopen("a_puntomedio.dat", "w");
        if (f == NULL) {
            fprintf(stderr, "Errore apertura file a_puntomedio.dat\n");
            free(a);
            return 1;
        }

        t = 0.0;
        for (int n = 0; n < nsteps; n++) {
            fprintf(f, "%.15e %.15e\n", t, a[n]);
            t += h;
        }

        fclose(f);
        free(a);
        return 0;
    }

    if (stamp(argv[4], "PASSO_ADATTIVO")) {
        int nvar = 1;         // numero di variabili
        float ystart[2];      // indice da 1
        ystart[1] = (float)a_init;  

        float t1 = 0.0f;
        float t2 = (float)t_max;
        float eps  = 1.0e-8f;
        float h1   = 1.0e-6f;
        float hmin = 1.0e-12f;

        int nok = 0, nbad = 0;
        int Nplot = 1000;     // numero di punti per il plot

        /* Variabili globali per odeint
           dxsav deve essere compatibile con kmax, altrimenti non salvi
           punti fino a t2 e l'interpolazione finale va fuori range.
        */
        kmax = Nplot + 10;
        kount = 0;
        dxsav = (t2 - t1) / (Nplot - 1);

        xp = vector(1, kmax);       
        yp = matrix(1, nvar, 1, kmax);

        /* Chiamata al driver ODE */
        odeint(ystart, nvar, t1, t2, eps, h1, hmin, &nok, &nbad, derivs, rkqs);

        FILE *fplot = fopen("a_adaptive_full.dat", "w");
        if (fplot == NULL) {
            fprintf(stderr, "Errore apertura file a_adaptive_full.dat\n");
            free_vector(xp, 1, kmax);
            free_matrix(yp, 1, nvar, 1, kmax);
            return 1;
        }

        for (int i = 0; i < Nplot; i++) {
            double tplot = i * (t2 - t1) / (double)(Nplot - 1);
            double aplot;

            if (i == 0 || tplot <= 0.0) {
                aplot = a_init;
            }
            else if (kount < 1) {
                aplot = a_init;
            }
            else if (tplot <= xp[1]) {
                double t0 = 0.0;
                double t1_ = xp[1];
                double a0 = a_init;
                double a1_ = yp[1][1];
                aplot = a0 + (a1_ - a0) * (tplot - t0) / (t1_ - t0);
            }
            else if (tplot >= xp[kount]) {
                double t0 = xp[kount];
                double t1_ = t2;
                double a0 = yp[1][kount];
                double a1_ = ystart[1];

                if (fabs(t1_ - t0) > 1.0e-12) {
                    aplot = a0 + (a1_ - a0) * (tplot - t0) / (t1_ - t0);
                } else {
                    aplot = a1_;
                }
            }
            else {
                int k;
                for (k = 1; k < kount; k++) {
                    if (xp[k] <= tplot && tplot <= xp[k+1]) break;
                }

                {
                    double t0 = xp[k];
                    double t1_ = xp[k+1];
                    double a0 = yp[1][k];
                    double a1_ = yp[1][k+1];
                    aplot = a0 + (a1_ - a0) * (tplot - t0) / (t1_ - t0);
                }
            }

            fprintf(fplot, "%.15e %.15e\n", tplot, aplot);
        }

        fclose(fplot);
        free_vector(xp, 1, kmax);
        free_matrix(yp, 1, nvar, 1, kmax);
        return 0;
    }

    fprintf(stderr, "Metodo non riconosciuto: %s\n", argv[4]);
    fprintf(stderr, "Usa PUNTO_MEDIO oppure PASSO_ADATTIVO\n");
    return 1;
}