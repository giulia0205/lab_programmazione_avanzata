#ifndef IO_H
#define IO_H

#include <stdio.h>

//dichiarazione della struct che rappresenta un record da salvare nel file binario
typedef struct {
    double t, x, y, vx, vy, ax, ay;   
} OutputRecord;

//dichiarazione delle funzioni di io.c
FILE **open_binary_files(int Np);
void write_record(FILE *file, OutputRecord rec);
void close_binary_files(int Np, FILE **files);

#endif