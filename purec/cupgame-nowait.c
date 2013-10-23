#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "cupgamelib.h"

int main(int argc, char **argv){
    if (argc != 5){
        printf("<N> <height> <restoring> <filename>\n");
        return 1;
    }

    int N;
    double h = 0.7;
    double r = 0.3;
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
            double xin[3] = {2.0*i/N-1, 2.0*j/N-1, height};
            if (sqrt(xin[0]*xin[0] + xin[1]*xin[1]) - h > r)
                continue;
            //double xin[3] = {1.0*i/N, 1.0*i/N, 1}; //0.35 + 0.5*(double)j/N};
            //restore = 0.5 + 0.5*(double)j/N;
            //double xin[3] = {6.0*i/N-3.0, 6.0*i/N-2.0, 0.75 + 5*(double)j/N};
            double vin[3] = {0.0, 0.0, -1e-1};

            int b = trackCollisions(3, xin, 3, vin, h, r, restore, 1000);
            bounces[i+j*N] = b;
            steps++;
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
