#ifndef UTILS_H                 /*mettere tutto negli if serve nel caso in cui ci siano diversi header nel progetto per cui alcune strutture*/
#define UTILS_H                 /*potrebbero essere chiamate più volte e quindi ridefinite*/


typedef struct {
    double Mnorm, Lnorm, Tnorm;

    double Ms;        // massa della stella
    double delta_t;   // passo temporale
    double G_scaled;  // costante gravitazionale nelle unità interne

    int Np;
    double *r;
} GeneralData;


void read_params(const char *filename, GeneralData *gd); //dichiarazione della funzione di lettura dei parametri

#endif