#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define XTOL 3e-16
#define NMAX 200

void find_root_pair(double *poly, int deg, double *qpoly, double *r1, double *r2){
    int i;
    double c,d,g,h,u,v,det,err,vo,uo;
    double *rpoly = (double*)malloc(sizeof(double)*deg);
   
    u = poly[deg-1] / poly[deg];
    v = poly[deg-2] / poly[deg];

    int nsteps = 0;  err = 10*XTOL;
    while (err > XTOL && nsteps < NMAX){
        qpoly[deg] = qpoly[deg-1] = 0;
        for (i=deg-2; i>=0; i--)
            qpoly[i] = poly[i+2] - u*qpoly[i+1] - v*qpoly[i+2];
        c = poly[1] - u*qpoly[0] - v*qpoly[1];
        d = poly[0] - v*qpoly[0];

        rpoly[deg] = rpoly[deg-1] = 0;
        for (i=deg-2; i>=0; i--)
            rpoly[i] = qpoly[i+2] - u*rpoly[i+1] - v*rpoly[i+2];
        g = qpoly[1] - u*rpoly[0] - v*rpoly[1];
        h = qpoly[0] - v*rpoly[0];

        det = 1.0/(v*g*g + h*(h-u*g));
        uo = u; vo = v;
        u = u - det*(g*d - c*h);
        v = v - det*((g*u-h)*d - g*v*c);

        err = sqrt((u-uo)*(u-uo)/(uo*uo) + (v-vo)*(v-vo)/(vo*vo));
        nsteps++;
    }
    printf("n = %i\n", nsteps);
    // b^2 - 4*a*c
    double desc = u*u - 4*v; 

    // if we have imaginery roots, quit now
    if (desc < 0){
        *r1 = NAN; *r2 = NAN;
        return;
    }

    double sqt = 0.5*sqrt(desc);
    double boa = -0.5*u;
    
    *r1 = 1.0/(boa + sqt);
    *r2 = 1.0/(boa - sqt);
    free(rpoly);
}

double find_smallest_root(double *poly, int deg){
    int i;
    int nroots = 0;
    double r1, r2;
    double *realroots = (double*)malloc(sizeof(double)*deg);
    double *temppoly = (double*)malloc(sizeof(double)*deg);

    for (i=deg; i>=1; i-=2){
        find_root_pair(poly, i, temppoly, &r1, &r2);
        memcpy(poly, temppoly, sizeof(double)*deg);

        if (isnan(r1) || isnan(r2))
            continue;

        realroots[nroots] = r1; nroots++;
        realroots[nroots] = r2; nroots++;
        printf("%0.16f  %0.16f  \n", r1, r2);
    }

    // if we didn't find any roots, return a nan (only special number)
    if (nroots == 0){
        free(realroots);
        free(temppoly);
        return NAN;
    }

    // otherwise, find the root closest to zero
    double minroot = realroots[0];
    for (i=0; i<nroots; i++)
        if (minroot > realroots[i] && realroots[i] > 0) 
            minroot = realroots[i];

    free(realroots);
    free(temppoly);

    if (minroot < 0) 
        return NAN;
    return minroot;
}

int main(){
    //double poly[] = {1, -21, 175, -735, 1624, -1764, 720}; int deg=6;
    //double poly[] = {1, -59.0/6, 47./3, -49./6, 4./3}; int deg=4;
    double poly[] = {1, 559./6, -4183./6, -2689./2, 23161./6, -6938./3, 400}; int deg=6;
    find_smallest_root(poly, deg);
    return 0;
}
