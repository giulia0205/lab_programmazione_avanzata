/*WALL BOUNDARY CONDITIONS*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 100

int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("Uso: %s K\n", argv[0]);
        return 1;
    }

    double K = atof(argv[1]);                /*K viene letta da terminale*/

    double dx = 1.0 / NX;   
    //double dt = 0.4 * dx * dx / K;           /*condizione di stabilità per Eulero*/
    double dt = 0.0001;
    double tmax = 1.0;
    int nsteps = (int)(tmax / dt);

    /*Creazione vettore u con le condizioni iniziali richieste*/
    double u[NX];
    double x[NX];
    for (int i = 0; i < NX; i++) {
        x[i] = (i + 0.5) * dx;   
        u[i] = 0.0;
    }

    for (int i = 45; i <= 54; i++) {
        u[i] = 1.0;
    }

    /*Wall boundary conditions*/
    u[0]   = u[1];
    u[NX-1] = u[NX-2];

    /*Discretizzazione temporale*/
    double u_new[NX];

    for (int n = 0; n < nsteps; n++) {     
        for (int i = 1; i < NX-1; i++) {
            u_new[i] = u[i] + dt * K * (u[i+1]-2*u[i]+u[i-1])/(dx*dx);
        }

        /*Copio u_new in u */
        for (int i = 1; i < NX-1; i++) {
            u[i] = u_new[i];
        }

        /*Wall boundary conditions */
        u[0] = u[1];
        u[NX-1] = u[NX-2];
    }


    /*Scrittura del file di output*/
    FILE *f = fopen("diffusione_output.dat", "w");
    for (int i = 0; i < NX; i++) {
        double x = (i + 0.5) * dx;
        fprintf(f, "%f %f\n", x, u[i]);
    }
    fclose(f);

    return 0;
}

