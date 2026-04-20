/*Equazione di secondo grado: dichiarare i coefficienti nel main e leggerli da una funzione
    ESERCIZIO SUI PUNTATORI*/
#include <stdio.h>
#include <math.h>


void mycoefficient(double *a, double *b, double *c) {
    printf("Inserisci il coefficiente a: ");
    scanf("%lf", a);
    printf("Inserisci il coefficiente b: ");
    scanf("%lf", b);
    printf("Inserisci il coefficiente c: ");
    scanf("%lf", c);
}


int main() {
    double a, b, c;
    mycoefficient(&a, &b, &c); 

    double Delta = b*b - 4*a*c;  
    
    if (Delta > 0) {
        double x1 = (-b + sqrt(Delta)) / (2*a);
        double x2 = (-b - sqrt(Delta)) / (2*a);
        printf("Due soluzioni reali: x1 = %.2f, x2 = %.2f\n", x1, x2);
        
    } else if (Delta == 0) {
        double x = -b / (2*a);
        printf("Soluzione reale doppia: x = %.2f\n", x);
    } else {
        
      printf("Non esistono soluzioni reali");
    }

    /*controllo sulla lettura dei parametri*/
    printf("%.5f\n", a);
    printf("%.5f\n", b);
    printf("%.5f\n", c);

    return 0;
}