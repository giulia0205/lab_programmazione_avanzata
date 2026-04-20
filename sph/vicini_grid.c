#include "vicini_grid.h"
#include "physics.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vicini_qsort.h"

#define GAMMA 1.4
#define SIGMA_1D (2.0/3.0)
#define N_NGB     64
#define HALF_NGB  32
#define EPS       1.0e-12

/*CREAZIONE DELLA GRIGLIA*/

int cell_index(double x, const Grid1D *g) {
    int k = (int)floor((x - g->xmin) / g->delta);

    /*Controllo sul fatto che l'indice trovato stia nel range di indici della griglia*/
    if (k < 0 || k >= g->ncell) {
        fprintf(stderr, "ERRORE: particella fuori griglia: x=%g xmin=%g xmax=%g delta=%g ncell=%d k=%d\n");
        exit(1);
    }

    return k;
}

void build_grid(Grid1D *g, Particella *p, int N) {
    /*Azzaramento dell'intera griglia*/
    for (int k = 0; k < g->ncell; k++) {
        g->head[k] = NULL;
    }

    for (int i = 0; i < N; i++) {
        /*Calcolo dell'indice della particella*/
        int k = cell_index(p[i].Pos, g);
        /*Inserimento in testa*/
        p[i].next = g->head[k];
        g->head[k] = &p[i];
    }
}

/*CALCOLO DELLE QUANTITÀ FISICHE*/

int cmp_double(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;

    if (da < db) return -1;
    if (da > db) return 1;
    return 0;
}

/*Calcolo della smoothing lenght*/
void compute_h_grid(Particella *p, int N, const Grid1D *g) {                     
    for (int i = 0; i < N; i++) {
    
        double left_dist[256];      //array contenente i vicini di sinistra
        double right_dist[256];     //array contenente i vicini di destra
        int nleft = 0;              //conteggio dei vicini trovati a sinistra
        int nright = 0;             //conteggio dei vicini trovati a destra

        /*Calcolo della cella di riferimento per ciascuna particella*/
        int k0 = cell_index(p[i].Pos, g);

        /*Espansione progressiva di celle fino a trovare abbastanza vicini*/
        for (int shell = 0; shell < g->ncell; shell++) {
            int found_enough = (nleft >= HALF_NGB && nright >= HALF_NGB);
            if (found_enough) break;

            if (shell == 0) {
                int kc = k0;

                for (Particella *pj = g->head[kc]; pj != NULL; pj = pj->next) {
                    if (pj == &p[i]) continue;

                    double dx = pj->Pos - p[i].Pos;

                    /* Controllo anti-wrap spurio */
                    if (fabs(dx) > 0.5 * g->L) continue;

                    if (dx < 0.0 && nleft < 256) {
                        left_dist[nleft++] = -dx;
                    } else if (dx > 0.0 && nright < 256) {
                        right_dist[nright++] = dx;
                    }
                }

            } else {
                int k_left  = (k0 - shell + g->ncell) % g->ncell;
                int k_right = (k0 + shell) % g->ncell;

                /* Cella a sinistra */
                for (Particella *pj = g->head[k_left]; pj != NULL; pj = pj->next) {
                    if (pj == &p[i]) continue;

                    double dx = pj->Pos - p[i].Pos;

                    /* Controllo anti-wrap spurio */
                    if (fabs(dx) > 0.5 * g->L) continue;

                    if (dx < 0.0 && nleft < 256) {
                        left_dist[nleft++] = -dx;
                    } else if (dx > 0.0 && nright < 256) {
                        right_dist[nright++] = dx;
                    }
                }

                /* Cella a destra: attenzione a non ricontarla se coincide con la sinistra */
                if (k_right != k_left) {
                    for (Particella *pj = g->head[k_right]; pj != NULL; pj = pj->next) {
                        if (pj == &p[i]) continue;

                        double dx = pj->Pos - p[i].Pos;

                        /* Controllo anti-wrap spurio */
                        if (fabs(dx) > 0.5 * g->L) continue;

                        if (dx < 0.0 && nleft < 256) {
                            left_dist[nleft++] = -dx;
                        } else if (dx > 0.0 && nright < 256) {
                            right_dist[nright++] = dx;
                        }
                    }
                }
            }
        }

        /*Ordinamento crescente delle distanze*/
        qsort(left_dist,  nleft,  sizeof(double), cmp_double);
        qsort(right_dist, nright, sizeof(double), cmp_double);

        /*Estrazione delle distanze rilevanti*/
        double d_left;
        double d_right;

        if (nleft >= HALF_NGB) {
            d_left = left_dist[HALF_NGB - 1];
        } else if (nleft > 0) {
            d_left = left_dist[nleft - 1];
        } else {
            d_left = EPS;
        }

        if (nright >= HALF_NGB) {
            d_right = right_dist[HALF_NGB - 1];
        } else if (nright > 0) {
            d_right = right_dist[nright - 1];
        } else {
            d_right = EPS;
        }

        /*Costruzione di h*/
        p[i].h = 0.5 * fmax(d_left, d_right);

        /*Protezione numerica*/
        if (p[i].h < EPS) p[i].h = EPS;
    }
}

/*Calcolo della densità*/
void compute_density_grid(Particella *p, int N, const Grid1D *g) {             
    for (int i = 0; i < N; i++) {
        double rho = 0.0;
        int k0 = cell_index(p[i].Pos, g);

        for (int dk = -1; dk <= 1; dk++) {
            int kc = (k0 + dk + g->ncell) % g->ncell;

            for (Particella *pj = g->head[kc]; pj != NULL; pj = pj->next) {
                double dx = p[i].Pos - pj->Pos;
                double r  = fabs(dx);

                /*Controllo anti-wrap spurio*/
                if (r > 0.5 * g->L) continue;

                if (r <= 2.0 * p[i].h) {
                    rho += pj->Mass * kernel_cubic_1d(r, p[i].h);
                }
            }
        }

        p[i].Density = fmax(rho, EPS);
    }
}

/*Calcolo acc e dU*/
void compute_acc_dU_grid(Particella *p, int N, const Grid1D *g) {          
    for (int i = 0; i < N; i++) {
        double acc = 0.0;
        double du  = 0.0;

        int k0 = cell_index(p[i].Pos, g);

        for (int dk = -1; dk <= 1; dk++) {
            int kc = (k0 + dk + g->ncell) % g->ncell;

            for (Particella *pj = g->head[kc]; pj != NULL; pj = pj->next) {
                if (pj == &p[i]) continue;

                double dx = p[i].Pos - pj->Pos;
                double r  = fabs(dx);

                /*Controllo anti-wrap spurio*/
                if (r > 0.5 * g->L) continue;

                double hij = hij_mean(p[i].h, pj->h);               // si è usata h media
                /*Controllo sul suppporto del kernel*/
                if (r > 2.0 * hij) continue;

                double gradW = grad_kernel_cubic_1d(p[i].Pos, pj->Pos, hij);
                double Pi_ij = artificial_viscosity(&p[i], pj, hij);

                double rhoi2 = p[i].Density * p[i].Density + EPS;
                double rhoj2 = pj->Density * pj->Density + EPS;

                double common = p[i].Pressure / rhoi2
                              + pj->Pressure / rhoj2
                              + Pi_ij;

                double vij = p[i].Vel - pj->Vel;

                acc -= pj->Mass * common * gradW;
                du  += 0.5 * pj->Mass * common * vij * gradW;
            }
        }

        p[i].Acc = acc;
        p[i].dU  = du;
    }
}
