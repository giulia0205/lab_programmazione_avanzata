
#include "io.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


// funzione che nomina e apre i file
FILE **open_binary_files(int Np) {                                   //la funzione deve restituire un puntatore di file, quindi è dichiarata con FILE **
    FILE**files = malloc(Np*sizeof(FILE*));
    char name[256];  
    for (int i = 0; i < Np; i++) {
        snprintf(name, sizeof(name), "planet_%d.bin", i);           //scrive dentro a "name" la variabile formattata
        files[i] = fopen(name, "wb");
    }
    return files;                                                   // la funzione restituisce solo il puntatore
}

// funzione che scrive un record in un file 
void write_record(FILE *file, OutputRecord rec) {
    fwrite(&rec, sizeof(OutputRecord), 1, file);
}

// funzione che chiude i file
void close_binary_files(int Np, FILE **files) {                             //ora lavoriamo con l'array di files
    for (int i = 0; i < Np; i++) {
        fclose(files[i]);
    }
    free(files);
}