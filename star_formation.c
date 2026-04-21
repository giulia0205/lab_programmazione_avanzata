/*
Write a midpoint integrator and apply it to solving the following ODE: d rho_star/ dt = c_star * rho_gas/t_gas ---- This is the simplest possible star formation law.
    rho_star: stellar density
    c_star: star formation efficiency ---- c_star=0.02
    rho_gas: gas density (constant) ------ rho_gas= 10 M_sun/kpc^3 ------ keeping rho_gas constant physically means that the gas is replenished while it form stars!
    t_gas: time consumption of gas ------- t_gas=1.5 Gyr
Integrate if for 2*t_gas and write output each 1/5t_gas. 
    How do you chose the integration interval h? 
    How do the accuracy of solution and the computing time vary with h?
This equation is analitically integrable and its solution is: rho_star= rho_star_0 + c_star*(rho_gas/t_gas)*(t-t_0) ---- coundary conditions: t_0=0, rho_star_0=0
You can use it to test your integration. 
    How the accuracy depend on h? Why?

What would change if the law were  d rho_star/ dt = c_star * (rho_gas/t_gas)*t ?
    What is the impact of the size of h, now?

If we do not suppose the gas to be replenished, a second, couple ODE for the gas is simply: d rho_gas / dt = - d rho_star / dt
Modify your code to integrate these two equations togheter (coupled system).
     What is changing now?
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*Definizione dei parametri*/
double c_star = 0.02;
double rho_gas_0 = 10.0;    //M_sun/kpc^3
double t_gas = 1.5;     //Gyr

/*Derivata costante*/
double ode_star_cost(double rho_gas) {
    return c_star * rho_gas/t_gas;
}

/*Derivata lineare*/
double ode_star_lin(double t) {
    return c_star * (rho_gas_0/t_gas) * t;
}

/*equazione del gas*/
double ode_gas(double rho_gas) {
    return -c_star * (rho_gas/t_gas);
}


int stamp(char *arg, char *mode) {
    return strcmp(arg, mode) == 0;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("Uso: %s COSTANTE | LINEARE | SISTEMA_ACCOPPIATO\n", argv[0]);
        return 1;
    }

    const double t_max = 2.0 * t_gas;         //tempo totale di integrazione
    const double h = 0.2 * t_gas;
    const int nsteps = (int)round(t_max / h);

    double rho_star[nsteps + 1];              //+1 perché vogliamo includere anche t=0
    double rho_gas[nsteps + 1];

    /*Metodo del punto medio cost*/
    if (stamp(argv[1], "COSTANTE")) {

        rho_star[0] = 0.0;

        for (int n = 0; n < nsteps; n++) {
            double k1    = ode_star_cost(rho_gas_0);
            double y_mid = rho_star[n] + 0.5 * h * k1;
            (void)y_mid;
            double k_mid = ode_star_cost(rho_gas_0);

            rho_star[n + 1] = rho_star[n] + h * k_mid;
        }

        FILE *f = fopen("starformation_cost.dat", "w");
        if (!f) {
            perror("Errore apertura file");
            return 1;
        }

        for (int n = 0; n <= nsteps; n++) {
            double t = n * h;
            fprintf(f, "%f %f\n", t, rho_star[n]);
        }
        fclose(f);
    }

    /*Metodo del punto medio lineare*/
    else if (stamp(argv[1], "LINEARE")) {

        rho_star[0] = 0.0;

        for (int n = 0; n < nsteps; n++) {
            double t = n * h;

            double k1    = ode_star_lin(t);
            double y_mid = rho_star[n] + 0.5 * h * k1;
            (void)y_mid; 
            double k_mid = ode_star_lin(t + 0.5 * h);

            rho_star[n + 1] = rho_star[n] + h * k_mid;
        }

        FILE *f = fopen("starformation_lin.dat", "w");
        if (!f) {
            perror("Errore apertura file");
            return 1;
        }

        for (int n = 0; n <= nsteps; n++) {
            double t = n * h;
            fprintf(f, "%f %f\n", t, rho_star[n]);
        }
        fclose(f);
    }

    /*Metodo del punto medio per il sistema accoppiato*/
    else if (stamp(argv[1], "SISTEMA_ACCOPPIATO")) {

        rho_star[0] = 0.0;
        rho_gas[0]  = rho_gas_0;

        for (int n = 0; n < nsteps; n++) {

            double k1_star = ode_star_cost(rho_gas[n]);
            double k1_gas  = ode_gas(rho_gas[n]);

            double y_mid_star = rho_star[n] + 0.5 * h * k1_star;
            double y_mid_gas  = rho_gas[n]  + 0.5 * h * k1_gas;
            (void)y_mid_star; 

            double k_mid_star = ode_star_cost(y_mid_gas);
            double k_mid_gas  = ode_gas(y_mid_gas);

            rho_star[n + 1] = rho_star[n] + h * k_mid_star;
            rho_gas[n + 1]  = rho_gas[n]  + h * k_mid_gas;
        }

        FILE *f = fopen("starformation_acc.dat", "w");
        if (!f) {
            perror("Errore apertura file");
            return 1;
        }

        for (int n = 0; n <= nsteps; n++) {
            double t = n * h;
            fprintf(f, "%f %f %f\n", t, rho_star[n], rho_gas[n]);
        }
        fclose(f);
    }

    return 0;
}




