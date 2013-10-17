#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include "cupgamelib.h"

int main(int argc, char **argv){
    if (argc != 5){
        printf("Incorrect arguments supplied, must be <N> <height> <restoring> <filename>\n");
        return 1;
    }

    int N;
    double h = 0.7;
    double r = 0.299;
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
    #pragma omp parallel shared(steps)
    for (int i=0; i<N; i++){

        #pragma omp for nowait
        for (int j=0; j<N; j++){
            double xin[3] = {1.0*i/N, 1.0*j/N, height};
            double vin[3] = {0.0, 0.0, -1e-1};

            if (sqrt(xin[0]*xin[0] + xin[1]*xin[1]) - h > r)
                continue;

            int b = trackCollisions(3, xin, 3, vin, h, r, restore, 1000);
            bounces[i+j*N] = b;
        }
        #pragma omp master
        {
            steps += N;
            printf("done: %0.4f \t rate: %0.2f\r", (float)steps/(N*N), rate);
            fflush(stdout);

            clock_gettime(CLOCK_REALTIME, &end);
            rate = steps/((end.tv_sec-start.tv_sec)+(end.tv_nsec-start.tv_nsec)/1e9);
        }
    }
    printf("done: %0.4f \t rate: %0.2f\n", (float)steps/(N*N), rate);

    FILE *f = fopen(filename, "wb");
    fwrite(bounces, sizeof(int), N*N, f);
    fclose(f);

    free(bounces);
    return 0;
}
