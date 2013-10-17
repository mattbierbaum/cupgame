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
    double r = 0.29;
    double height;
    double restore;
    char filename[1024];

    N = atoi(argv[1]);
    height = atof(argv[2]);
    restore = atof(argv[3]);
    strcpy(filename, argv[4]);

    int *bounces = malloc(sizeof(int)*N*N);

    FILE *f = fopen(filename, "wb");
    for (int i=0; i<N*N; i++){
        int zero = 0;
        fwrite(&zero, sizeof(int), 1, f);
    }
    fflush(f);
    fseek(f, 0, SEEK_SET);

    double rate = 0.0;
    struct timespec start, end;
    clock_gettime(CLOCK_REALTIME, &start);

    int subf = 8;
    for (int i=0; i<N/subf; i++){
        printf("slice %d \t rate: %0.2f\r", i*subf, rate);
        fflush(stdout);

        #pragma omp parallel for
        for (int k=0; k<subf; k++){
            for (int j=0; j<N; j++){
                int ind = i*subf + k;
                double xin[3] = {1.0*ind/N, 1.0*j/N, height};
                double vin[3] = {0.0, 0.0, -1e0};

                if (sqrt(xin[0]*xin[0] + xin[1]*xin[1]) - h > r)
                    continue;

                int b = trackCollisions(3, xin, 3, vin, h, r, restore, 1000);
                bounces[ind*N+j] = b;
            }
        }
        fwrite(bounces+i*subf*N, sizeof(int), N*subf, f);
        fflush(f);

        clock_gettime(CLOCK_REALTIME, &end);
        rate = (i+1)*N*subf/((end.tv_sec-start.tv_sec)+(end.tv_nsec-start.tv_nsec)/1e9);
    }
    fclose(f);
   
    free(bounces);
    return 0;
}
