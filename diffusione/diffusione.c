/*WALL BOUNDARY CONDITIONS*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 100

int main(int argc, char *argv[]) {
    double dx = 1.0 / NX;   /*NX-1*/
    double K= atof(argv[1]);                /*K viene letta da terminale*/
    //double dt = 0.5 * pow(dx,2)/ K;         /*condizione di stabilità per Eulero-------- lui ha dissato a 0,0001*/
    double dt = 0.4 * dx * dx / K;
    double tmax = 1.0;
    int nsteps = (int)(tmax / dt);

    /*creazione vettore u con le condizioni iniziali richieste*/
    double u[NX];
    double x[NX];
    for (int i = 0; i < NX; i++) {
        x[i] = (i + 0.5) * dx;   
        u[i] = 0.0;
    }

    for (int i = 45; i <= 54; i++) {
        u[i] = 1.0;
    }

    /*wall boundary conditions*/
    u[0]   = u[1];
    u[NX-1] = u[NX-2];

    /*discretizzazione temporale*/
    double u_new[NX];

    for (int n = 0; n < nsteps; n++) {     
        for (int i = 1; i < NX-1; i++) {
            u_new[i] = u[i] + dt * K * (u[i+1]-2*u[i]+u[i-1])/(dx*dx);
        }

        /* copio u_new in u */
        for (int i = 1; i < NX-1; i++) {
            u[i] = u_new[i];
        }

        /* wall boundary conditions */
        u[0] = u[1];
        u[NX-1] = u[NX-2];
    }


    /*scrittura del file di output*/
    FILE *f = fopen("diffusione_output.dat", "w");
    for (int i = 0; i < NX; i++) {
        double x = (i + 0.5) * dx;
        fprintf(f, "%f %f\n", x, u[i]);
    }
    fclose(f);

    return 0;
}

