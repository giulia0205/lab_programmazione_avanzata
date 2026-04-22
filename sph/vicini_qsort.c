#include "vicini_qsort.h"
#include "physics.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define GAMMA 1.4
#define SIGMA_1D (2.0/3.0)
#define N_NGB     64
#define HALF_NGB  32
#define EPS       1.0e-12


/*massimo e minimo fra due valori*/
int imax(int a, int b) { 
    if (a > b)
        return a;
    else
        return b;
}
int imin(int a, int b) {
    if (a < b)
        return a;
    else
        return b;
}

int cmp_posizione(const void *a, const void *b) {
    const Particella *pa = (const Particella *)a;
    const Particella *pb = (const Particella *)b;

    if (pa->Pos < pb->Pos) return -1;
    if (pa->Pos > pb->Pos) return 1;
    return 0;
}

/*Controllo sull'ordinamento*/
int check_sorted_by_pos(Particella *p, int N) {
    for (int i = 0; i < N - 1; i++) {
        if (p[i].Pos > p[i + 1].Pos) {
            return 0;
        }
    }
    return 1;
}

/*Calcolo della smoothing lenght per ciascuna particella*/
void compute_h_qsort(Particella *p, int N) {
    for (int i = 0; i < N; i++) {               
        double d_left  = 0.0;
        double d_right = 0.0;

        if (i >= HALF_NGB) {                                //controllo sull'esistenza della particella i-32
            d_left = p[i].Pos - p[i - HALF_NGB].Pos;
        }
        if (i + HALF_NGB < N) {                             //controllo sull'esistenza della particella i+32
            d_right = p[i + HALF_NGB].Pos - p[i].Pos;
        }

        if (i < HALF_NGB) {                                 //caso particelle [0,32]
            p[i].h = 0.5 * fmax(d_right, EPS);
        } 
        else if (i + HALF_NGB >= N) {                       //caso particelle [N-32,N]
            p[i].h = 0.5 * fmax(d_left, EPS);
        } 
        else {                                              //caso particelle [32,N+32]
            p[i].h = 0.5 * fmax(d_left, d_right);
        }

        /*Controllo di sicurezza*/
        if (p[i].h < EPS) {
            p[i].h = EPS;
        }
    }
}

/*Calcolo della densità*/
void compute_density_qsort(Particella *p, int N) {
    for (int i = 0; i < N; i++) {

        /*Le ghost particles non vengono aggiornate: mantengono i valori iniziali assegnati come condizioni al contorno*/
        if (p[i].ghost) continue;

        double rho = 0.0;

        /*Range delle particelle adiacenti*/
        int jmin = i;
        while (jmin > 0 && (p[i].Pos - p[jmin - 1].Pos) <= 2.0 * p[i].h) {
            jmin--;
        }

        int jmax = i;
        while (jmax < N - 1 && (p[jmax + 1].Pos - p[i].Pos) <= 2.0 * p[i].h) {
            jmax++;
        }

        /*Ciclo sulle particelle vicine*/
        for (int j = jmin; j <= jmax; j++) {
            double r = fabs(p[i].Pos - p[j].Pos);

            if (r <= 2.0 * p[i].h) {
                rho += p[j].Mass * kernel_cubic_1d(r, p[i].h);
            }
        }

        p[i].Density = fmax(rho, EPS);
    }
}

/*Calcolo dell'accelerazione e della derivata dell'energia interna (vengono calcolare assieme in quanto usano la stessa somma sui vicini)*/
void compute_acc_dU_qsort(Particella *p, int N) {
    for (int i = 0; i < N; i++) {

        /*Le ghost particles non evolvono: non aggiorniamo né accelerazione né energia interna*/
        if (p[i].ghost) {
            p[i].Acc = 0.0;
            p[i].dU  = 0.0;
            continue;
        }

        double acc = 0.0;
        double dU  = 0.0;

        int jmin = i;
        while (jmin > 0 && (p[i].Pos - p[jmin - 1].Pos) <= 2.0 * p[i].h) {
            jmin--;
        }

        int jmax = i;
        while (jmax < N - 1 && (p[jmax + 1].Pos - p[i].Pos) <= 2.0 * p[i].h) {
            jmax++;
        }

        for (int j = jmin; j <= jmax; j++) {
            if (j == i) continue;                                       //la particella non interagisce con sè stessa

            double hij = hij_mean(p[i].h, p[j].h);
            double rij = fabs(p[i].Pos - p[j].Pos);

            if (rij > 2.0 * hij) continue;                              //oltre al supporto compatto del kernel le particelle non contribuiscono

            /*Calcolo dei gradienti*/
            double gradW = grad_kernel_cubic_1d(p[i].Pos, p[j].Pos, hij);
            double Pi_ij = artificial_viscosity(&p[i], &p[j], hij);

            /*Calcolo dei quadrati delle densità*/
            double rhoi2 = p[i].Density * p[i].Density + EPS;          //+EPS non è fisica, è protezione numerica: serve a evitare divisioni per zero o numeri troppo piccoli.
            double rhoj2 = p[j].Density * p[j].Density + EPS;

            /*Calcolo del termine comune alle due equazioni*/
            double common = (p[i].Pressure / rhoi2) + (p[j].Pressure / rhoj2) + Pi_ij;

            /*Calcolo della velocità relativa*/
            double vij = p[i].Vel - p[j].Vel;

            /*Aggiornamento di acc e dU*/
            acc -= p[j].Mass * common * gradW;
            dU  += 0.5 * p[j].Mass * common * vij * gradW;
        }

        p[i].Acc = acc;
        p[i].dU  = dU;
    }
}