#include "physics.h"
#include "io.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

/*Funzione di lettura dei parametri*/
void read_params(GeneralData *gd) {
    /*Apertura del file*/
    FILE *fp = fopen("params.txt", "r");       
    /*Controllo sull'apertura del file*/
    if (!fp) {                                  
        printf("Errore: impossibile aprire params.txt\n");
        exit(1);
    }
    char name[64];                  //name conterrà: “UdL”, “UdM”, “UdT”
    double value;                   //value conterrà i valori numerici

    while (fscanf(fp, "%s %lf", name, &value) == 2) {           
        if (strcmp(name, "UnitLength") == 0) {                         
            gd->Lnorm = value;
        }
        else if (strcmp(name, "UnitMass") == 0) {
            gd->Mnorm = value;
        }
        else if (strcmp(name, "UnitVel") == 0) {
            gd->Vnorm = value;
        }
        else if (strcmp(name, "t_final") == 0) {
            gd->t_final = value;
        }
        else if (strcmp(name, "N_points") == 0) {
            gd->N_p = (int)value;
        }
        else if (strcmp(name, "N_grid") == 0) {
            gd->N_g = (int)value;
        }
        else if (strcmp(name, "N_tsteps_est") == 0) {
            gd->N_tsteps_est = (int)value;
        }
        else if (strcmp(name, "BoxLength") == 0) {
            gd->L_Box = value;
        }
        else if (strcmp(name, "A_deltaPar") == 0) {
            gd->delta_rho_max = value;
        }
        /*Controllo su eventuali parametri in più se ho fatto casino*/
        else {                                                      
            printf("Attenzione: parametro sconosciuto %s\n", name);
        }
    }

    /*Calcolo di UnitTime*/
    gd->Tnorm = gd->Lnorm / gd->Vnorm;

    /*Riscaliamo G*/
    const double G = 6.67430e-8;;
    gd->G_scaled = G * pow(gd->Tnorm,2) * gd->Mnorm / pow(gd->Lnorm,3);

    /*Chiusura file*/
    fclose(fp);                                 
}

/*Allocazione variabili*/
Particle *allocate_particles(int N) {
    Particle *p = malloc(N * sizeof(Particle));
    if (p == NULL) {
        printf("Errore allocazione particelle\n");
        exit(1);
    }
    return p;
}

/*Scrittura file testo condizioni iniziali*/
void write_ic_text(const char *filename, Particle *part, int N) {
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("Errore apertura file %s\n", filename);
        return;
    }

    /*Scrittura numero di particelle*/
    fprintf(fp, "%d\n", N);            
    /*Scrittura delle particelle*/
    for (int i = 0; i < N; i++) {           
        fprintf(fp, "%f %f %f\n", part[i].x, part[i].v, part[i].mass);                 
    }

    fclose(fp);
}

/*Scrittura file binario condizioni iniziali*/
void write_ic_binary(const char *filename, Particle *part, int N) {
    FILE *fp = fopen(filename, "wb");
    if (fp == NULL) {
        printf("Errore apertura file %s\n", filename);
        return;
    }

    /*Scrittura numero di particelle*/
    fwrite(&N, sizeof(int), 1, fp);             
    /*Scrittura particelle*/
    fwrite(part, sizeof(Particle), N, fp);     

    fclose(fp);
}

//Istogramma
void write_histogram(const char *filename, Particle *part, int N, double L, int nbins) {
    int *hist = calloc(nbins, sizeof(int));                         //allocazione dei contatori: creazione un array di nbins interi tutti inizializzati a zero
    FILE *fp;
    double dx;

    if (hist == NULL) {
        printf("Errore allocazione istogramma\n");
        return;
    }

    dx = L / nbins;                                                 //larghezza del bin

    /*Conteggio delle particelle che cadono in ciascun bin */
    for (int i = 0; i < N; i++) {
        int bin = (int)(part[i].x / dx);                            //trova il bin di ciascuna particella

        /*Controlli superflui in questo caso specifico ma utili in generale*/
        if (bin < 0) bin = 0;
        if (bin >= nbins) bin = nbins - 1;

        hist[bin]++;
    }

    fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("Errore apertura file %s\n", filename);
        free(hist);
        return;
    }

    /*Scrittura centro bin e conteggio*/
    for (int i = 0; i < nbins; i++) {
        double xcenter = (i + 0.5) * dx;
        fprintf(fp, "%f %d\n", xcenter, hist[i]);
    }

    fclose(fp);
    free(hist);
}

/*Lettura condizioni iniziali formato testo*/
Particle *read_ic_text(const char *filename, int *N) {
    FILE *fp = fopen(filename, "r");

    if (fp == NULL) {
        printf("Errore apertura file %s\n", filename);
        return NULL;
    }

    if (fscanf(fp, "%d", N) != 1) {
        printf("Errore lettura numero particelle da %s\n", filename);
        fclose(fp);
        return NULL;
    }

    Particle *part = malloc((*N) * sizeof(Particle));
    if (part == NULL) {
        printf("Errore allocazione particelle\n");
        fclose(fp);
        return NULL;
    }

    for (int i = 0; i < *N; i++) {
        if (fscanf(fp, "%lf %lf %lf",
                    &part[i].x,
                    &part[i].v,
                    &part[i].mass) != 3) {
            printf("Errore lettura particella %d da %s\n", i, filename);
            free(part);
            fclose(fp);
            return NULL;
        }
    }

    fclose(fp);
    return part;
}

/*Lettura condizioni iniziali formato binario*/
Particle *read_ic_bin(const char *filename, int *N) {
    FILE *fp  = fopen(filename, "rb");
    if (fp == NULL) {
        printf("Errore apertura file %s\n", filename);
        return NULL;
    }

    if (fread(N, sizeof(int), 1, fp) != 1) {
        printf("Errore lettura numero particelle da %s\n", filename);
        fclose(fp);
        return NULL;
    }

    Particle *part = malloc((*N) * sizeof(Particle));
    if (part == NULL) {
        printf("Errore allocazione particelle\n");
        fclose(fp);
        return NULL;
    }

    if (fread(part, sizeof(Particle), *N, fp) != (size_t)(*N)) {
        printf("Errore lettura particelle da %s\n", filename);
        free(part);
        fclose(fp);
        return NULL;
    }
    
    fclose(fp);
    return part;
}

/*Scrittura file andamento temporale*/
void write_snapshot(int step, Particle *part, int Np, double *rho, int Ng, double L) {
    char fname_part[256];
    char fname_rho[256];

    snprintf(fname_part, sizeof(fname_part), "particles_step_%05d.txt", step);
    snprintf(fname_rho, sizeof(fname_rho), "density_step_%05d.txt", step);

    /* Particelle */
    FILE *fp_part = fopen(fname_part, "w");
    if (fp_part == NULL) {
        printf("Errore apertura file %s\n", fname_part);
    } else {
        fprintf(fp_part, "%d\n", Np);
        for (int p = 0; p < Np; p++) {
            fprintf(fp_part, "%.15e %.15e %.15e\n",
                    part[p].x, part[p].v, part[p].mass);
        }
        fclose(fp_part);
    }

    /*Quantità su griglia*/
    double dx = L / Ng;

    FILE *fp_rho = fopen(fname_rho, "w");
    if (fp_rho == NULL) {
        printf("Errore apertura file %s\n", fname_rho);
    } else {
        for (int i = 0; i < Ng; i++) {
            double x = (i + 0.5) * dx;
            fprintf(fp_rho, "%.15e %.15e\n", x, rho[i]);
        }
        fclose(fp_rho);
    }

}

/*File output finale: particelle, densità, forza, potenziale*/
void write_final_data(Particle *part, int Np, double *rho, double *pot, double *force, int Ng, double L) {
    /*Particelle finali*/
    FILE *fp_part = fopen("final_particles.txt", "w");
    if (fp_part == NULL) {
        printf("Errore apertura file final_particles.txt\n");
    } else {
        fprintf(fp_part, "%d\n", Np);
        for (int p = 0; p < Np; p++) {
            fprintf(fp_part, "%.15e %.15e %.15e\n",
                    part[p].x, part[p].v, part[p].mass);
        }
        fclose(fp_part);
    }

    double dx = L / Ng;

    /*Densità finale*/
    FILE *fp_rho = fopen("density_final.txt", "w");
    if (fp_rho == NULL) {
        printf("Errore apertura file density_final.txt\n");
    } else {
        for (int i = 0; i < Ng; i++) {
            double x = (i + 0.5) * dx;
            fprintf(fp_rho, "%.15e %.15e\n", x, rho[i]);
        }
        fclose(fp_rho);
    }

    /*Potenziale finale*/
    FILE *fp_pot = fopen("potential_final.txt", "w");
    if (fp_pot == NULL) {
        printf("Errore apertura file potential_final.txt\n");
    } else {
        for (int i = 0; i < Ng; i++) {
            double x = (i + 0.5) * dx;
            fprintf(fp_pot, "%.15e %.15e\n", x, pot[i]);
        }
        fclose(fp_pot);
    }

    /*Forza finale*/
    FILE *fp_force = fopen("force_final.txt", "w");
    if (fp_force == NULL) {
        printf("Errore apertura file force_final.txt\n");
    } else {
        for (int i = 0; i < Ng; i++) {
            double x = (i + 0.5) * dx;
            fprintf(fp_force, "%.15e %.15e\n", x, force[i]);
        }
        fclose(fp_force);
    }
}
