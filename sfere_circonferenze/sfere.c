/*rifare esercizio delle circonferenze in 3D*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


/*definizione di una struct per i punti*/
typedef struct {
    float x;
    float y;
    float z;
} point;

/*usiamo coordinate sferiche così generiamo punti casuali all'intereno della sfera e non solo sul bordo*/ 
point random_point_in_sphere(float radius) {
    float r = ((float) rand()) / (RAND_MAX + 1.);    
    float costheta = 2. * ((float) rand()) / (RAND_MAX + 1.) - 1.;      /*Si usa direttamente costheta come variabile in modo che la distribuzione dei punti sia uniforme*/
    float sintheta = sqrtf(1. - costheta * costheta);                   /*se invece usassimo theta avremmo i punti concentrati ai poli*/
    float phi = 2*M_PI*rand()/(RAND_MAX+1.);

    point p;
    p.x = r * sintheta * cosf(phi);
    p.y = r * sintheta * sinf(phi);
    p.z = r * costheta;
    return p;
}


void generatore_punti_file(int index, int count, float radius) {
    char filename[100];          /*Dichiara un array di caratteri (una stringa C) chiamata filename, lunga 100 caratteri, che servirà a contenere il nome del file.*/
    snprintf(filename, sizeof(filename), "sfere%05d_%.5f.txt", index, radius);       /*Scrive una stringa formattata dentro filename*/

    FILE *f = fopen(filename, "w");      /*Apre il file appena nominato per scrittura ("w"); f è un puntatore a FILE*/
    /*if (!f) {
        perror("Errore apertura file");
        exit(1);
    }*/

    for (int i = 0; i < count; i++) {
        point p = random_point_in_sphere(radius);
        fprintf(f, "%f %f %f\n", p.x, p.y, p.z);
    }

    fclose(f);
}

int main(int argc, char *argv[]) {
    /*if (argc != 5) {                il programma si aspetta 5 argomenti totali, se non sono 5 argomenti, entra nell’if e restituisce un messaggio di errore.
        fprintf(stderr, "Uso: %s N M r1 r2\n", argv[0]);
        return 1;
    }*/

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