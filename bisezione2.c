/*risoluzione di varie equazioni con il metodo della bisezione
    ATTRAVERSO PUNTATORI E INPUT-OUTPUT
    per compilare: gcc -o bisezione2.x bisezione2.c -lm
                    ./bisezione2.x TEST 1.0 2.0 1.0 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 


double parab(double x, double acoeff, double bcoeff, double ccoeff) {
    return acoeff * pow(x, 2) + bcoeff * x + ccoeff;
}


double trasc1(double x, double acoeff, double bcoeff, double ccoeff) {
    return pow(x, x) - exp(-x) - 1.0;
}


double trasc2(double x, double acoeff, double bcoeff, double ccoeff) {
    return pow(x,3) + x + 1.0;
}

double trasc3(double x, double acoeff, double bcoeff, double ccoeff) {
    return pow(x,3)+x+1;
}

/*Funzione di confronto di stringhe*/
int stamp(const char *s1, const char *s2) {
    return strcmp(s1, s2) == 0;
}



double bisezione(double x1, double x2, double epsilon, double (*myfunc)(double, double, double, double),double acoeff, double bcoeff, double ccoeff) {
    double c;
    double a = myfunc(x1, acoeff, bcoeff, ccoeff);
    double b = myfunc(x2,acoeff, bcoeff, ccoeff);

    if (a * b > 0) {                                        /* Controllo che nell'intervallo considerato ci sia una radice*/
        printf ("La radice non è presente in questo intervallo");  
        return NAN;
    }

    while (x2-x1>epsilon) {
        a = myfunc(x1, acoeff, bcoeff, ccoeff);
        b = myfunc(x2,acoeff, bcoeff, ccoeff);
        c = (x1 + x2) / 2;
        double fc=myfunc(c,acoeff, bcoeff, ccoeff);;
        
        if (fc==0.) {
            return c;
        }

        if (fc * myfunc(x1,acoeff, bcoeff, ccoeff) < 0) {
            x2 = c;
        } else {
            x1 = c;
        }        
    }
    
return c;
}



/*metodo della bisezione*/
int main(int argc, char *argv[]){  
    double (*myfunc)(double, double, double, double);   /*argc è una variabile che contiene il numero di argomenti passati al programma da linea di comando;*/ 
    double acoeff, bcoeff, ccoeff;                      /*argv è un array di stringhe (char*), che contiene tutti questi argomenti.*/
    if (stamp(argv[1],"TEST")) {
        myfunc=&parab;
        acoeff=atof(argv[2]);
        bcoeff=atof(argv[3]);
        ccoeff=atof(argv[4]);
    }
    if (stamp(argv[1],"TRASCEND1")) {
        myfunc=&trasc1;
        acoeff=0.;
        bcoeff=0.;
        ccoeff=0.;
    }
    if (stamp(argv[1],"TRASCEND2")) {
        myfunc=&trasc2;
        acoeff=0.;
        bcoeff=0.;
        ccoeff=0.;
    }
    if (stamp(argv[1],"TRASCEND3")) {
        myfunc=&trasc3;
        acoeff=0.;
        bcoeff=0.;
        ccoeff=0.;
    }

    double radice=bisezione(0., 2., 0.01, myfunc, acoeff, bcoeff, ccoeff);
    printf ("%.5f\n", radice);
    return 0;
}




