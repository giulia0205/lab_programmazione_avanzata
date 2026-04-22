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

void read_params(const char *filename, GeneralData *gd) {
    FILE *fp = fopen(filename, "r");        //Apertura file
    if (!fp) {                              //Controllo di corretta apertura del file
        printf("Errore: impossibile aprire %s\n", filename);
        exit(1);
    }

    char name[20];                  //name conterrà il nome del parametro letto dal file

    gd->r = NULL;                   //inizializziamo il puntatore a NULL
    gd->Np = 0;                     //inizializziamo il numero di pianeti a 0

    while (fscanf(fp, "%19s", name) == 1) {     //scanf cerca di leggere il nome del parametro; se ci riesce il ciclo continua

        if (strcmp(name, "UdL") == 0) {                         //controlla quale parametro sta leggendo
            fscanf(fp, "%lf", &gd->Lnorm);
        }
        else if (strcmp(name, "UdM") == 0) {
            fscanf(fp, "%lf", &gd->Mnorm);
        }
        else if (strcmp(name, "UdT") == 0) {
            fscanf(fp, "%lf", &gd->Tnorm);
        }
        else if (strcmp(name, "Ms") == 0) {
            fscanf(fp, "%lf", &gd->Ms);
        }
        else if (strcmp(name, "delta_t") == 0) {
            fscanf(fp, "%lf", &gd->delta_t);
        }
        else if (strcmp(name, "Np") == 0) {
            fscanf(fp, "%d", &gd->Np);

            gd->r = malloc(gd->Np * sizeof(double));            //allocazione dinamica dell'array dei raggi
            if (!gd->r) {
                printf("Errore allocazione memoria per i raggi\n");
                fclose(fp);
                exit(1);
            }
        }
        else if (strcmp(name, "radii") == 0) {
            if (gd->Np <= 0) {                                  //controlliamo di aver già letto Np
                printf("Errore: leggere Np prima di radii\n");
                fclose(fp);
                exit(1);
            }

            for (int i = 0; i < gd->Np; i++) {                  //lettura dei raggi
                if (fscanf(fp, "%lf", &gd->r[i]) != 1) {
                    printf("Errore: numero di raggi insufficiente nel file\n");
                    fclose(fp);
                    exit(1);
                }
            }
        }
    }

    /*Riscaliamo G*/
    const double G = 6.67430e-11;
    gd->G_scaled = G * pow(gd->Tnorm, 2) * gd->Mnorm / pow(gd->Lnorm, 3);

    fclose(fp);                                  //Chiusura file
}
