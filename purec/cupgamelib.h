#ifndef __CUPGAME_H__
#define __CUPGAME_H__

int trackTrajectory(int NP, double *pos, int NV, double *vel, double h, double r, 
        double damp, int maxbounces, int NT, double *traj);
int trackCollisions(int NP, double *pos, int NV, double *vel, double h, double r, 
        double damp, int maxbounces);
double singleCollisions(int NP, double *pos, int NV, double *vel, double h, double r);


#endif
