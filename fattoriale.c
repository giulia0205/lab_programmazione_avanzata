/*Calcolo del fattoriale attraverso ricorsione*/
#include <stdio.h>

double fattoriale(int n) {
    if (n == 0 || n == 1)
        return 1;
    else
        return n * fattoriale(n - 1);
}

int main() {
    int numero;

    printf("Numero intero non negativo: ");
    scanf("%d", &numero);

    printf("Il fattoriale di %d è %0f\n", numero, fattoriale(numero));

    return 0;
}

/*numero max= 170*/
