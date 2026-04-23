#include "io.h"
#include "physics.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int stamp(char *arg, char *mode) {
    return strcmp(arg, mode) == 0;
}

int main(int argc, char *argv[]) {
    GeneralData all;            /* struttura dei parametri */
    Particle *part = NULL;      /* array di particelle */
    double *rho = NULL;         /* array della densità */

    /*Lettura parametri*/
    read_params(&all);
    double dt_init = all.t_final / all.N_tsteps_est;                //passo iniziale, usato poi come tetto massimo del timestep adattivo: 0.1
    all.delta_t = dt_init;

    /*Controllo sul passaggio della tipologia di file da leggere*/
    if (argc < 2) {
        printf("Uso: %s TESTO oppure %s BINARIO\n", argv[0], argv[0]);
        return 1;
    }

    /*LETTURA FILE CONDIZIONI INIZIALI*/
    if (stamp(argv[1],"TESTO")) {
        part = read_ic_text("i_c.txt", &all.N_p);
    }
    else if (stamp(argv[1],"BINARIO")) {
        part = read_ic_bin("i_c.bin", &all.N_p);
    }
    else {
        printf("Argomento non valido: usa TESTO oppure BINARIO\n");
        return 1;
    }

    if (part == NULL) {
        printf("Errore nella lettura delle condizioni iniziali\n");
        return 1;
    }


    /*CALCOLO DELLA DENSITÀ DALLA DISTRIBUZIONE DI PARTICELLE*/
    /*Allocazione della densità*/
    rho = malloc(all.N_g * sizeof(double));
    /*Controllo sull'allocazione*/
    if (rho == NULL) {
        printf("Errore allocazione densita'\n");
        free(part);
        return 1;
    }
    /*Calcolo della densità*/
    density(rho, part, all.N_p, all.N_g, all.L_Box);


    /*CALCOLO DEL POTENZIALE PASSANDO IN FOURIER*/
    /*Allocazione di materia*/
    double *pot = malloc(all.N_g * sizeof(double));
    /*Controllo sull'allocazione*/
    if (pot == NULL) {
        printf("Errore allocazione potenziale\n");
        free(rho);
        free(part);
        return 1;
    }
    /*Calcolo del potenziale*/
    compute_potential_fft(rho, pot, all.N_g, all.L_Box);


    /*CALCOLO DELLA FORZA*/
    /*Allocazione di materia*/
    double *force = malloc(all.N_g * sizeof(double));
    double *acc = malloc(all.N_p * sizeof(double));
    /*Controllo sull'allocazione*/
    if (force == NULL || acc == NULL) {
        printf("Errore allocazione forza/accelerazione\n");
        free(force);
        free(acc);
        free(pot);
        free(rho);
        free(part);
        return 1;
    }
    /*Calcolo della forza sulla griglia*/
    compute_force_grid(force, pot, all.N_g, all.L_Box);
    force_to_particles(acc, force, part, all.N_p, all.N_g, all.L_Box);

    write_snapshot(0, part, all.N_p, rho, all.N_g, all.L_Box);

    /*Loop temporale*/
    double t = 0.0;
    int step = 1;
    double dt = all.delta_t;

    while (t < all.t_final) {
        if (t + dt > all.t_final) {
            dt = all.t_final - t;
        }

        /*Kick*/
        leapfrog_kick(part, acc, all.N_p, dt);
        /*Drift*/
        leapfrog_drift(part, all.N_p, dt, all.L_Box);

        /*Calcolo delle nuove forze sulle nuove posizioni*/
        density(rho, part, all.N_p, all.N_g, all.L_Box);
        compute_potential_fft(rho, pot, all.N_g, all.L_Box);
        compute_force_grid(force, pot, all.N_g, all.L_Box);
        force_to_particles(acc, force, part, all.N_p, all.N_g, all.L_Box);

        /*Kick*/
        leapfrog_kick(part, acc, all.N_p, dt);

        t += dt;

        /*Salvataggio snapshot ogni 200 step*/
        if (step % 200 == 0) {
            write_snapshot(step, part, all.N_p, rho, all.N_g, all.L_Box);
        }

        all.delta_t = compute_timestep(part, all.N_p, all.N_g, all.L_Box);
        if (all.delta_t > dt_init) {
            all.delta_t = dt_init;
        }

        dt = all.delta_t;
        step++;
    }

    /*SCRITTURA OUTPUT FINALI*/
    write_final_data(part, all.N_p, rho, pot, force, all.N_g, all.L_Box);
    
    free(acc);
    free(force);
    free(pot);
    free(rho);
    free(part);
    return 0;
}


