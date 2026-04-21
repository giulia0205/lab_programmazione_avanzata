/* PER PLOTTARE:
gnuplot> set terminal pngcairo size 800,600

Terminal type is now 'pngcairo'
Options are ' background "#ffffff" enhanced fontscale 1.0 size 800, 600 '
gnuplot> set output "sfere7.png"
gnuplot> set xrange [-6:6]
gnuplot> set yrange [-6:6]
gnuplot> plot "sfere00002_3.00000.txt" using 1:2 with points pt 7
gnuplot> set output
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


/*definizione di una struct per i punti*/
typedef struct {
    float x;
    float y;
} point;

/*usiamo coordinate polari così generiamo punti casuali all'intereno della circonferenza e non solo sul bordo*/ 
point random_point_in_circle(float radius) {
    float r = rand()/(RAND_MAX+1.);             
    float theta = 2*M_PI*rand()/(RAND_MAX+1.);
    point p;
    p.x = r * cosf(theta);
    p.y = r * sinf(theta);
    return p;
}

void generatore_punti_file(int index, int count, float radius) {
    char filename[100];          /*Dichiara un array di caratteri (una stringa C) chiamata filename, lunga 100 caratteri, che servirà a contenere il nome del file.*/
    snprintf(filename, sizeof(filename), "circonferenze%05d_%.5f.txt", index, radius);       /*Scrive una stringa formattata dentro filename*/

    FILE *f = fopen(filename, "w");      /*Apre il file appena nominato per scrittura ("w"); f è un puntatore a FILE*/

    for (int i = 0; i < count; i++) {
        point p = random_point_in_circle(radius);
        fprintf(f, "%f %f\n", p.x, p.y);
    }

    fclose(f);
}

int main(int argc, char *argv[]) {
    if (argc != 5) {                //il programma si aspetta 5 argomenti totali, se non sono 5 argomenti, entra nell’if e restituisce un messaggio di errore.
        fprintf(stderr, "Uso: %s N M r1 r2\n", argv[0]);
        return 1;
    }

    /*Lettura degli argomenti da linea di comando*/
    int N = atoi(argv[1]);
    int M = atoi(argv[2]);
    float r1 = atof(argv[3]);
    float r2 = atof(argv[4]);

    /*Inizializza seme per numeri casuali: così ogni volta si generano numeri casuali differenti*/
    srand(time(NULL));

    /*Vengono generati i due file*/
    generatore_punti_file(1, N, r1);
    generatore_punti_file(2, M, r2);

    printf("File generati con successo.\n");

    return 0;
}


