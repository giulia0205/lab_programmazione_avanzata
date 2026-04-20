/*
- legge un file di parametri;
- genera N particelle con distribuzione scelta;
- assegna a ciascuna particella posizione e velocità iniziale;
- salva tutto in un file che poi il PM leggerà;
- fa istogramma per controllo della distribuzione.
*/

#include "io.h"
#include "physics.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int stamp(char *arg, char *mode) {
    return strcmp(arg, mode) == 0;
}

int main(int argc, char *argv[]) {
    /*Controllo sul passaggio della tipologia di file da scrivere*/
    if (argc < 2) {
        printf("Uso: %s TESTO oppure %s BINARIO\n", argv[0], argv[0]);
        return 1;
    }

    /*Lettura dei parametri*/
    GeneralData all;
    read_params(&all);

    /*Allocazione di memoria*/
    Particle *part = allocate_particles(all.N_p);

    /*Distribuzione iniziale delle particelle*/
    srand(12345);           //random seed                                                               
    generate_distribution(part, all.N_p, all.L_Box, all.delta_rho_max);

    /*Scrittura del file delle condizioni iniziali*/
    if (stamp(argv[1],"TESTO")) {
        write_ic_text("i_c.txt", part, all.N_p);
    }

    if (stamp(argv[1],"BINARIO")) {
        write_ic_binary("i_c.bin", part, all.N_p);
    }

    /*Istogramma di controllo sulla distribuzione*/
    write_histogram("hist_ic.txt", part, all.N_p, all.L_Box, all.N_g);

    free(part);
    return 0;
}


