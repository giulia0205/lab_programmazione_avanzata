#ifndef VICINI_QSORT_H
#define VICINI_QSORT_H

#include "physics.h"

int imax(int a, int b);
int imin(int a, int b);

int cmp_posizione(const void *a, const void *b);
int check_sorted_by_pos(Particella *p, int N);

void compute_h_qsort(Particella *p, int N);
void compute_density_qsort(Particella *p, int N);
void compute_acc_dU_qsort(Particella *p, int N);


#endif