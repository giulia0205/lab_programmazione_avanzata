#include "io.h"
#include "physics.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main() {
    GeneralData all;
    double *rho = NULL;
    double *pot = NULL;
    FILE *fp = NULL;

    read_params(&all);

    rho = malloc(all.N_g * sizeof(double));
    pot = malloc(all.N_g * sizeof(double));

    if (rho == NULL || pot == NULL) {
        printf("Errore allocazione memoria\n");
        free(rho);
        free(pot);
        return 1;
    }

    double dx = all.L_Box / all.N_g;
    //double k = 2.0 * M_PI / all.L_Box;
    double k = 2.0 * acos(-1.0) / all.L_Box;

    /* Costruzione di una densità sinusoidale di test:
       rho(x) = sin(2*pi*x/L)
    */
    for (int i = 0; i < all.N_g; i++) {
        double x = (i + 0.5) * dx;
        rho[i] = sin(k * x);
    }

    /* Calcolo del potenziale numerico */
    compute_potential_fft(rho, pot, all.N_g, all.L_Box);

    /* Scrittura su file del confronto */
    fp = fopen("pot_test.txt", "w");
    if (fp == NULL) {
        printf("Errore apertura file pot_test.txt\n");
        free(rho);
        free(pot);
        return 1;
    }

    for (int i = 0; i < all.N_g; i++) {
        double x = (i + 0.5) * dx;
        double phi_exact = -sin(k * x) / (k * k);

        fprintf(fp, "%.15e %.15e %.15e %.15e\n",
                x, rho[i], pot[i], phi_exact);
    }

    fclose(fp);

    printf("Test completato. File scritto: pot_test.txt\n");

    free(rho);
    free(pot);
    return 0;
}