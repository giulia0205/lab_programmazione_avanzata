#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define SIZE 10000

void recursive_free(int *arr[], int index, int total) {                 /*Funzione ricorsiva per deallocare ciascun indice dell'array*/
    if (index == total) return;
    free(arr[index]);
    recursive_free(arr, index + 1, total);
}

int main(int argc, char *argv[]) {

    if (argc != 2) {
        printf("Uso: %s N\n", argv[0]);
        return 1;
    }

    int N = atoi(argv[1]);

    srand(time(NULL));

    /*Allocazione dei vettori */
    int *x = malloc(SIZE * sizeof(int));        
    int *y = malloc(SIZE * sizeof(int));
    int *z = malloc(SIZE * sizeof(int));

    /*I vettori vengono inizializzati a 0 */
    for (int i = 0; i < SIZE; i++) {
        x[i] = 0;
        y[i] = 0;
        z[i] = 0;
    }

    /*Riempimento vellore x*/
    for (int i = 0; i < N; i++) {
        int index = rand() % SIZE;                              /*Il numero casuale compreso fra 0 e 9999 viene generato*/
        x[index] = 1;
    }

    /*Riempimento vettore y*/
    for (int i = 0; i < N; i++) {
        int index = rand() % SIZE;
        y[index] = 1;
    }

    /*Riempimento vettore z*/
    for (int i = 0; i < SIZE; i++) {
        z[i] = x[i] + y[i];  
    }

    int count_equal_two = 0;

    for (int i = 0; i < SIZE; i++) {                            /*print degli elementi di z pari a 2*/
        if (z[i] == 2) {
            printf("z[%d] = 2\n", i);
            count_equal_two=count_equal_two + 1;
        }
    }

    printf("\nTotale elementi di z uguali a 2: %d\n", count_equal_two);

    /*Deallocazione ricorsiva*/
    int *arrays_to_free[3] = {x, y, z};
    recursive_free(arrays_to_free, 0, 3);

    return 0;
}
