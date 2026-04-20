/*risoluzione di varie equazioni con il metodo della bisezione
    ATTRAVERSO FUNZIONE CON NUMERO DI ARGOMENTI VARIABILE
    per compilare: gcc -o bisezione3.x bisezione3.c -lm
                    ./bisezione3.x TEST 1.0 2.0 1.0 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <stdarg.h>     /*per scrivere la funzione con numero di argomenti variabile*/

double parab_var(double x, int count, ...) {
    va_list args;
    va_start(args, count);
    double a = va_arg(args, double);
    double b = va_arg(args, double);
    double c = va_arg(args, double);
    va_end(args);
    return a*x*x + b*x + c;
}

double trasc1_var(double x, int count, ...) {
    return pow(x, x) - exp(-x) - 1.0;
}

double trasc2_var(double x, int count, ...) {
    return pow(x,3) + x + 1.0;
}

double trasc3_var(double x, int count, ...) {
    return pow(x,3)+x+1;
}

/*Funzione di confronto di stringhe*/
int stamp(const char *s1, const char *s2) {
    return strcmp(s1, s2) == 0;
}

typedef double (*func_var_ptr)(double, int, ...);

double bisezione(double x1, double x2, double epsilon, func_var_ptr f, int count, ...) {
    va_list args;                               
    va_start(args, count);

    double params[10];                  /*massimo 10 argomenti*/
    for(int i=0; i<count; i++) {
        params[i] = va_arg(args, double);            /*va_arg serve per leggere count argomenti di tipo double e salvarli nell’array params[];*/
    va_end(args);                                    /*dopo aver finito di leggere gli argomenti, si chiude la lista con va_end*/
    }

    double c;
    double a = f(x1, count, params[0], params[1], params[2]);
    double b = f(x2, count, params[0], params[1], params[2]);

    if (a * b > 0) {                                        /*Controllo che nell'intervallo considerato ci sia una radice*/
        printf ("La radice non è presente in questo intervallo");  
        return NAN;
    }

    while ((x2 - x1) > epsilon) {
        a = f(x1, count, params[0], params[1], params[2]);
        b = f(x2, count, params[0], params[1], params[2]);
        c = (x1 + x2) / 2;
        double fc = f(c, count, params[0], params[1], params[2]);

        if (fc==0.) {
            return c;
        }

        if (fc * a < 0) {
            x2 = c;
        } else {
            x1 = c;
        }
    }
    return c;
}

int main(int argc, char *argv[]) {  
    if (argc < 2) {
        printf("Uso: %s TIPO [parametri]\n", argv[0]);
        return 1;
    }

    func_var_ptr myfunc;
    double acoeff;
    double bcoeff;
    double ccoeff;
    int count = 0;

    if (stamp(argv[1], "TEST")) {
        if (argc < 5) {
            printf("Mancano coefficienti per parabola\n");
            return 1;
        }
        myfunc = &parab_var;
        acoeff = atof(argv[2]);
        bcoeff = atof(argv[3]);
        ccoeff = atof(argv[4]);
        count=3;
    } 
    
    if (stamp(argv[1],"TRASCEND1")) {
        myfunc=&trasc1_var;
        count=0;
    }
    if (stamp(argv[1],"TRASCEND2")) {
        myfunc=&trasc2_var;
        count=0;
    }
    if (stamp(argv[1],"TRASCEND3")) {
        myfunc=&trasc3_var;
        count=0;
    }

    double radice;
    if (count == 3) {
        radice = bisezione(0., 2., 0.01, myfunc, count, acoeff, bcoeff, ccoeff);
    } else {
        radice = bisezione(0., 2., 0.01, myfunc, count);
    }

    printf("Radice: %.5f\n", radice);
    return 0;
}
