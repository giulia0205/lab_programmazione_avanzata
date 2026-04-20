#include "utils.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

/*Funzione che prende i valori del parameter file, li memorizza in una struttura GeneralData e restituisce tali valori nel main:
    1. aprire il file
    2. leggere ogni riga
    3. capire quale valore stai leggendo
    4. assegnarlo alla struttura
    5. chiudere il file

    la faccio anche calcolare G nelle unità corrette?
*/
void read_params(GeneralData *gd) {
    FILE *fp = fopen("params.txt", "r");        //Apertura file
    if (!fp) {                                  //Controllo di corretta apertura del file
        printf("Errore: impossibile aprire params.txt\n");
        exit(1);
    }
    char name[20];                  //name conterrà: “UdL”, “UdM”, “UdT”
    double value;                   //value conterrà i valori numerici

    while (fscanf(fp, "%s %lf", name, &value) == 2) {           //scanf cerca di leggere una stringa e un numero, se riesce a leggere entrambi il ciclo continua

        if (strcmp(name, "UdL") == 0) {                         //controlla quale parametro sta leggendo
            gd->Lnorm = value;
        }
        else if (strcmp(name, "UdM") == 0) {
            gd->Mnorm = value;
        }
        else if (strcmp(name, "UdT") == 0) {
            gd->Tnorm = value;
        }
        /*
        else {                                                    
            printf("Attenzione: parametro sconosciuto %s\n", name);
        }*/
        
    }

    /*Riscaliamo G*/
    const double G = 6.67430e-11;
    gd->G_scaled = G * pow(gd->Tnorm,2) * gd->Mnorm / pow(gd->Lnorm,3);

    fclose(fp);                                  //Chiusura file
}
