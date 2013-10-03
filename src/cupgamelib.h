#ifndef __CUPGAME_H__
#define __CUPGAME_H__

int trackTrajectory(double *pos, double *vel, double h, double r, 
        double damp, int maxbounces, double *traj, int *clen, int mlen);
int trackCollisions(double *pos, double *vel, double h, double r, 
        double damp, int maxbounces);
double singleCollisions(int NP, double *pos, int NV, double *vel, double h, double r);


#endif
