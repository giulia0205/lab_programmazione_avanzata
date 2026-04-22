#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "io.h"             

/*CONVERTITORE DA FILE BINARIO A FINE DI TESTO*/

int stamp(char *arg, char *mode) {
    return strcmp(arg, mode) == 0;
}


int main(int argc, char *argv[]) {
    /*Controllo sulla lettura da terminale*/
    if (argc < 3) {
        printf("Uso: %s SUN oppure TRAPPIST1 skip\n", argv[0]);
        return 1;
    }

    /*Scelta del file dei parametri*/
    GeneralData all = {0};

    if (stamp(argv[1], "SUN")) {
        read_params("params_sun.txt", &all);
    }
    else if (stamp(argv[1], "TRAPPIST1")) {
        read_params("params_trappist1.txt", &all);
    }
    else {
        printf("Errore: scrivere SUN oppure TRAPPIST1\n");
        return 1;
    }

    /*Lettura da terminale dello skip*/
    int skip = atoi(argv[2]);

    char binname[256];
    char txtname[256];

    for (int i = 0; i < all.Np; i++) {
        /*Costruzione dei nomi dei file come stringhe, in modo da inserirli dentro gli array binname e txtname*/
        snprintf(binname, sizeof(binname), "planet_%d.bin", i);
        snprintf(txtname, sizeof(txtname), "planet_%d.txt", i);

        /*Apertura file binario*/
        FILE *bin = fopen(binname, "rb");
        if (!bin) {
            printf("Errore apertura file binario %s\n", binname);
            continue;
        }

        /*Creazione e apertura del file testo di output*/
        FILE *txt = fopen(txtname, "w");
        if (!txt) {
            printf("Errore apertura file testo %s\n", txtname);
            fclose(bin);
            continue;
        }

        OutputRecord rec;
        int step = 0;

        /*Lettura del file binario:     non vogliamo leggere solo un record ma tutti */
        while (fread(&rec, sizeof(OutputRecord), 1, bin) == 1) {
            if (step % skip == 0) {
                fprintf(txt, "%e %e %e %e %e %e %e\n",
                        rec.t,
                        rec.x,
                        rec.y,
                        rec.vx,
                        rec.vy,
                        rec.ax,
                        rec.ay);
            }
            step++;
        }

        /*Chiusura dei file aperti*/
        fclose(bin);
        fclose(txt);

        printf("Convertito: %s -> %s\n", binname, txtname);
    }

    free(all.r);

    return 0;
}