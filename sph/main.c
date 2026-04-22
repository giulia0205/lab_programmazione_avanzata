#include "physics.h"
#include "vicini_qsort.h"
#include "vicini_grid.h"
#include "io.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

double elapsed_seconds(struct timespec t0, struct timespec t1) {
    return (double)(t1.tv_sec - t0.tv_sec) + 1.0e-9 * (double)(t1.tv_nsec - t0.tv_nsec);
}

int stamp(char *arg, char *mode) {
    return strcmp(arg, mode) == 0;
}

int main(int argc, char *argv[]) {

    Particella *p = NULL;

    /*Funzione di lettura del parameter file*/
    GeneralData all;
    read_params(&all);
    fprintf(stderr, "PARAMS: dt=%g t_end=%g Np=%d\n", all.dt, all.t_end, all.Np);

    const double xmin = -0.5 * all.L_grid;           
    const double xmax =  0.5 * all.L_grid;

    double t = 0.0;
    /*Variabili necessarie al calcolo del wall-clock time*/
    struct timespec tic, toc;
    double total_force_time = 0.0;
    double step_force_time = 0.0;

    /*Controllo sul passaggio della tipologia di sorting*/
    if (argc < 2) {
        printf("Uso: %s QSORT oppure %s GRID\n", argv[0], argv[0]);
        return 1;
    }

    /*Allocazione memoria particelle ( con aggiunta di particelle fantasma)*/
    int Ng_sx = 40;
    int Ng_dx = 40;
    int Ntot = Ng_sx + all.Np + Ng_dx;
    p = malloc(Ntot * sizeof(Particella));
    if (p == NULL) {
        fprintf(stderr, "Errore: allocazione memoria per le particelle fallita.\n");
        return 1;
    }

    /*CONDIZIONI INIZIALI SHOCK TUBE*/
    init_particles_shock_tube(p, all.Np, Ng_sx, Ng_dx, xmin, xmax);

    /*METODO QSORT*/
    if (stamp(argv[1],"QSORT")) {

        qsort(p, Ntot, sizeof(Particella), cmp_posizione);
        compute_h_qsort(p, Ntot);
        compute_density_qsort(p, Ntot);
        compute_pressure(p, Ntot);
        compute_acc_dU_qsort(p, Ntot);

        /* Salvataggio condizioni iniziali */
        write_file("output_qsort_00000.txt", p + Ng_sx, all.Np);

        int step = 1;

        while(t < all.t_end) {                 

            clock_gettime(CLOCK_MONOTONIC, &tic);

            /*RICERCA DEI VICINI*/

            leapfrog_half_kick(p + Ng_sx, all.Np, all.dt);
            leapfrog_drift(p + Ng_sx, all.Np, all.dt);

            for (int i = Ng_sx; i < Ng_sx + all.Np; i++) {
                if (p[i].Pos < xmin) {
                    p[i].Pos = xmin + (xmin - p[i].Pos);
                    p[i].Vel = -p[i].Vel;
                }
                else if (p[i].Pos > xmax) {
                    p[i].Pos = xmax - (p[i].Pos - xmax);
                    p[i].Vel = -p[i].Vel;
                }
            }

            qsort(p, Ntot, sizeof(Particella), cmp_posizione);
            /*Controllo sull'ordinamento*/
            if (!check_sorted_by_pos(p, Ntot)) {
                fprintf(stderr, "Errore: array non ordinato correttamente.\n");
                free(p);
                return 1;
            }

            /*SPH ADATTIVO*/
            compute_h_qsort(p, Ntot);
            compute_density_qsort(p, Ntot);
            compute_pressure(p, Ntot);
            compute_acc_dU_qsort(p, Ntot);

            leapfrog_half_kick(p + Ng_sx, all.Np, all.dt);

            /*Calcolo del wall-clock time*/
            clock_gettime(CLOCK_MONOTONIC, &toc);
            step_force_time = elapsed_seconds(tic, toc);
            total_force_time += step_force_time;

            /*Scrittura file ASCII*/ 
            if (step % 20 == 0) {
                char filename[64];
                sprintf(filename, "output_qsort_%05d.txt", step);
                write_file(filename, p + Ng_sx, all.Np);
            }

            t += all.dt;
            step++;
        }

        write_file("output_final_qsort.txt", p + Ng_sx, all.Np);
        printf("Wall-clock time totale QSORT = %.10f s\n", total_force_time);
    }

    /*METODO GRID*/
    else if (stamp(argv[1],"GRID")) {

        /*Calcolo di h_ref iniziale */
        qsort(p, Ntot, sizeof(Particella), cmp_posizione);
        compute_h_qsort(p, Ntot);
        double h_ref = compute_h_mean(p + Ng_sx, all.Np);

        /*Costruzione della griglia iniziale */
        Grid1D g0;
        g0.xmin  = p[0].Pos - 1.0e-12;
        g0.xmax  = p[Ntot - 1].Pos + 1.0e-12;
        g0.L     = g0.xmax - g0.xmin;
        g0.delta = h_ref;
        g0.ncell = (int)ceil(g0.L / g0.delta);
        g0.head = calloc(g0.ncell, sizeof(Particella *));
        if (g0.head == NULL) {
            fprintf(stderr, "Errore: allocazione griglia iniziale fallita.\n");
            free(p);
            return 1;
        }

        build_grid(&g0, p, Ntot);
        compute_h_grid(p, Ntot, &g0);
        compute_density_grid(p, Ntot, &g0);
        compute_pressure(p, Ntot);
        compute_acc_dU_grid(p, Ntot, &g0);

        /*Salvataggio condizioni iniziali*/
        write_file("output_grid_00000.txt", p + Ng_sx, all.Np);

        h_ref = compute_h_mean(p + Ng_sx, all.Np);

        free(g0.head);

        int step = 1;
                
        while(t < all.t_end) {    
            /*RICERCA DEI VICINI*/              
            Grid1D g;

            /*Definizione delle celle di ampiezza pari a <h> del passo precedente */
            g.xmin  = p[0].Pos - 1.0e-12;
            g.xmax  = p[Ntot - 1].Pos + 1.0e-12;
            g.L = g.xmax - g.xmin;
            g.delta = h_ref;                            //h_ref ci serve ad aggiornare delta ad ogni passo temporale

            /*Controllo sull'ampiezza della griglia*/
            if (g.delta <= 0.0) {
                fprintf(stderr, "Errore: delta_cell non valida.\n");
                free(p);
                return 1;
            }

            g.ncell = (int)ceil(g.L / g.delta);         //ceil arrotonda per eccesso

            /*Allocazione memoria arrray di teste*/
            g.head = malloc(g.ncell * sizeof(Particella *));
            if (g.head == NULL) {
                fprintf(stderr, "Errore: allocazione griglia fallita.\n");
                free(p);
                return 1;
            }

            clock_gettime(CLOCK_MONOTONIC, &tic);

            leapfrog_half_kick(p + Ng_sx, all.Np, all.dt);
            leapfrog_drift(p + Ng_sx, all.Np, all.dt);

            for (int i = Ng_sx; i < all.Np + Ng_sx; i++) {
                if (p[i].Pos < xmin) {
                    p[i].Pos = xmin + (xmin - p[i].Pos);   // riflessione geometrica
                    p[i].Vel = -p[i].Vel;                  // inversione velocità
                }
                else if (p[i].Pos > xmax) {
                    p[i].Pos = xmax - (p[i].Pos - xmax);
                    p[i].Vel = -p[i].Vel;
                }
            }

            qsort(p, Ntot, sizeof(Particella), cmp_posizione);

            /*Definizione delle celle di ampiezza pari a <h> del passo precedente */
            g.xmin  = p[0].Pos - 1.0e-12;
            g.xmax  = p[Ntot - 1].Pos + 1.0e-12;
            g.L = g.xmax - g.xmin;
            g.delta = h_ref;                            

            /*Controllo sull'ampiezza della griglia*/
            if (g.delta <= 0.0) {
                fprintf(stderr, "Errore: delta_cell non valida.\n");
                free(p);
                return 1;
            }

            g.ncell = (int)ceil(g.L / g.delta);         

            /*Allocazione memoria arrray di teste*/
            g.head = malloc(g.ncell * sizeof(Particella *));
            if (g.head == NULL) {
                fprintf(stderr, "Errore: allocazione griglia fallita.\n");
                free(p);
                return 1;
            }

            /*Costruzione griglia*/
            build_grid(&g, p, Ntot);

            /*SPH ADATTIVO*/
            compute_h_grid(p, Ntot, &g);
            compute_density_grid(p, Ntot, &g);
            compute_pressure(p, Ntot);
            compute_acc_dU_grid(p, Ntot, &g);

            leapfrog_half_kick(p + Ng_sx, all.Np, all.dt);

            /*Aggiorna della scala della griglia per il passo successivo*/
            h_ref = compute_h_mean(p + Ng_sx, all.Np);

            /*Calcolo del wall-clock time*/
            clock_gettime(CLOCK_MONOTONIC, &toc);
            step_force_time = elapsed_seconds(tic, toc);
            total_force_time += step_force_time;

            /*Scrittura file ASCII*/ 
            if (step % 20 == 0) {
                char filename[64];
                sprintf(filename, "output_grid_%05d.txt", step);
                write_file(filename, p + Ng_sx, all.Np);
            }

            free(g.head);
            t += all.dt;
            step++;
        }

        write_file("output_final_grid.txt", p + Ng_sx, all.Np);
        printf("Wall-clock time totale GRID = %.10f s\n", total_force_time);
    }

    else {
        printf("Argomento non valido: usa QSORT oppure GRID\n");
        free(p);
        return 1;
    }

    free(p);
    return  0;
}