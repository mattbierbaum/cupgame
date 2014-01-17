#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "cupgamelib.h"

#define XSCAN_UNITQUAD {2.0*i/N-1, 2.0*j/N-1, height}
#define XSCAN_UNITLINE {1.0*(i+j*N)/(N*N), 1.0*(i+j*N)/(N*N), height}
#define XSCAN_HEIGHT {1.0*i/N, 1.0*i/N, height+ 0.5*(double)j/N}
#define XSCAN_HEXAGON {6.0*i/N-3.0, 6.0*j/N-2.0, height} //0.75 + 5*(double)j/N}

int main(int argc, char **argv){
    if (argc != 5){
        printf("<N> <height> <restoring> <filename>\n");
        return 1;
    }

    int N;
    double h = 0.95;
    double r = 0.35;
    //double h = 0.95; double r = 0.3; restore = 0.75; parameters for a real game
    double height;
    double restore;
    char filename[1024];

    N = atoi(argv[1]);
    height = atof(argv[2]);
    restore = atof(argv[3]);
    strcpy(filename, argv[4]);

    int *bounces = malloc(sizeof(int)*N*N);

    double rate = 0.0;
    struct timespec start, end;
    clock_gettime(CLOCK_REALTIME, &start);

    int steps = 0;
    int userinformed = 0;
    #pragma omp parallel shared(userinformed)
    for (int i=0; i<N; i++){

        userinformed = 0;
        #pragma omp for nowait schedule(dynamic,N/16) reduction(+:steps)
        for (int j=0; j<N; j++){
            steps++;

            double vsq = 1.3*height;//30.0;
            double v = sqrt(vsq);
            double tmin = asin((height-4)/vsq)/2;
            double tmax = asin((height+4)/vsq)/2;
            double pmin = -asin(4./height);
            double pmax = asin(4./height);

            double theta = (tmax-tmin)*j/N + tmin;
            double phi = (pmax-pmin)*i/N + pmin;
            double xin[3] = {0, -height, 0.0}; //XSCAN_UNITQUAD;
            //double vin[3] = {0.0, 0.0, -1e-1};
            double vin[3] = {v*sin(theta)*sin(phi), v*sin(theta)*cos(phi), v*cos(theta)};

            //if (sqrt(xin[0]*xin[0] + xin[1]*xin[1]) - h > r)
            //    continue;
            //restore = 0.5 + 0.5*(double)j/N;

            int b = trackCollisions(3, xin, 3, vin, h, r, restore, 10000);
            bounces[i+j*N] = b;
        }

        if (!userinformed)
        {
            userinformed = 1;
            clock_gettime(CLOCK_REALTIME, &end);
            rate = steps/((end.tv_sec-start.tv_sec)+(end.tv_nsec-start.tv_nsec)/1e9);
            printf("done: %0.4f \t rate: %0.2f\r", (float)steps/(N*N), rate);
            fflush(stdout);
        }
    }
    printf("done: %0.4f \t rate: %0.2f\n", (float)steps/(N*N), rate);

    FILE *f = fopen(filename, "wb");
    fwrite(bounces, sizeof(int), N*N, f);
    fclose(f);

    free(bounces);
    return 0;
}
