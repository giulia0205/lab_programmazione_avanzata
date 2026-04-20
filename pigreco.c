#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/* Funzione per stimare pi greco*/
double stima_pi(int num_punti, double raggio) {
    int dentro_cerchio = 0;
    for (int i = 0; i < num_punti; i++) {
        /*Genera punto casuale tra -r e +r per x e y*/
        double x = ((double)rand() / RAND_MAX) * 2 * raggio - raggio;
        double y = ((double)rand() / RAND_MAX) * 2 * raggio - raggio;

        /*Controlla se il punto è dentro al cerchio*/
        if (x * x + y * y <= raggio * raggio) {
            dentro_cerchio++;
        }
    }
    double rapporto = (double)dentro_cerchio / num_punti;
    return rapporto * 4;  /*area cerchio / area quadrato = πr² / (2r)² = π / 4 */
}

int main() {
    int num_punti = 1000000; 
    double raggio = 1.0;
    /*Inizializza seme per numeri casuali: così ogni volta si generano numeri casuali differenti*/
    srand(time(NULL));
    double pi_stimato = stima_pi(num_punti, raggio);
    printf("Stima di pi greco con %d punti: %f\n", num_punti, pi_stimato);
    return 0;
}










