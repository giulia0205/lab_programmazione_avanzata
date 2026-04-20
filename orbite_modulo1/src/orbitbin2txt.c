#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "io.h"             //serve per la struct OutputRecord


/*CONVERTITORE DA FILE BINARIO A FINE DI TESTO:     passiamo da terminale Np e lo skip
    -apre il file binairo
    -legge i record
    -scrive il file .txt
*/

int main(int argc, char *argv[]) {
    /*Controllo sulla lettura da terminale*/
    if (argc < 3) {
        printf("Uso: %s Np skip\n", argv[0]);
        return 1;
    }

    /*Lettura da terminale del numero di pianeti e dello skip*/
    int Np = atoi(argv[1]);
    int skip = atoi(argv[2]);

    char binname[256];
    char txtname[256];

    for (int i = 0; i < Np; i++) {
        /*Costruzione dei nomi dei file come stringhe, in modo da inserirli dentro gli array binname e txtname*/
        /*vuol dire:
            -scrivi dentro binname
            -al massimo per sizeof(binname) caratteri
            -la stringa formattata "planet_%d.bin"
            -sostituendo %d con il valore di i
        */
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

    return 0;
}
