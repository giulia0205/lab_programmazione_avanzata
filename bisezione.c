/*risoluzione dell'equazione x^x-e^-x+1=0 con il metodo della bisezione*/

#include <stdio.h>
#include <math.h>


/*Definizione della funzione nei due casi TEST e TRASCEND*/
double y(double x);                     /*comando da terminale: gcc -DTEST*/

#ifdef TEST
double y(double x) {
    double acoeff = 1.0;
    double bcoeff = 2.0;
    double ccoeff = -1.0;
    return acoeff * pow(x, 2) + bcoeff * x + ccoeff;
}
#endif

#ifdef TRASCEND
double y(double x) {
    return pow(x, x) - exp(-x) - 1.0;
}
#endif


double bisezione(double x1, double x2, double epsilon) {
    double c;
    double a = y(x1);
    double b = y(x2);

    if (a * b > 0) {                /*Controllo che nell'intervallo considerato ci sia una radice*/
        printf ("la radice non è presente in questo intervallo");  
        return NAN;
    }

    while (x2-x1>epsilon) {
        double a = y(x1);
        double b = y(x2);
        c = (x1 + x2) / 2;
        double fc=y(c);

        if (fc==0.) {
            return c;
        }

        if (fc * y(x1) < 0) {
            x2 = c;
        } else {
            x1 = c;
        }        
    }
return c;
}


/*metodo della bisezione*/
int main(){
    double radice=bisezione(0., 2., 0.01);
    printf ("%.5f\n", radice);
    return 0;
}



