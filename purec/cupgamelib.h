#ifndef __CUPGAME_H__
#define __CUPGAME_H__

//========================================================
// global constants for the calculation
//========================================================
#define XTOL 1e-14
#define NMAX 512

#define DEG     8
#define DEGSIZE 9
#define MAXCUPS 50
#define TSAMPLES 25

#define RESULT_NOTHING   0
#define RESULT_COLLISION 1
#define RESULT_INCUP     2

#define EPS 1e-10

#define M_PI 3.14159265358979323846

//========================================================
/* These are functions that should be called externally */
int trackTrajectory(int NP, double *pos, int NV, double *vel, double h, double r,
        double damp, int maxbounces, int NT, double *traj);
int trackCollisions(int NP, double *pos, int NV, double *vel, double h, double r,
        double damp, int maxbounces);
double singleCollisions(int NP, double *pos, int NV, double *vel, double h, double r);

//========================================================
/* internal use functions only */
double polyeval(double *poly, int deg, double x);
void find_root_pair(double *poly, int deg, double *qpoly,
        double *r1, double *r2);
void find_all_roots(double *poly, int deg, double *roots, int *nroots);
double find_smallest_root(double *poly, int deg);

void build_torus_poly(double *pos, double *vel, double h, double r,
        double *cup, double *poly);
void build_cylinder_poly(double *pos, double *vel, double r,
        double cylx, double cyly, double *poly);
void build_zero_poly(double *pos, double *vel, double *poly);
void build_plane_poly(double *pos, double *vel, double r, double *poly);

double mymod(double a, double b);
double dot(double *r1, double *r2);
void cross(double *a, double *b, double *out);
void position(double *x0, double *v0, double t, double *out);
void velocity(double *v0, double t, double *out);
double distance_to_torus(double *pos, double h, double r, double *cup);

void cups_near(double *pos, double *cups, int *ncups);
void cups_near_hex(double *pos, double *cups, int *ncups);
void cups_near_6cup(double *pos, double *cups, int *ncups);
void collision_normal(double *pos, double h, double *cup, double *out);
void rotate_about_axis(double *r, double *axis, double theta, double *out);
void reflect_velocity(double *pos, double *vel, double h,
        double f, double *cup, double *out);
int collides_with_cup(double *pos, double *vel, double h, double r,
        double *cup, double *tcoll);
int collision_near(double *pos, double *vel, double h, double r,
        double *tcoll, double *cup);
int collision_time(double *pos, double *vel, double h,
        double r, double *tcoll, double *cup);
int cup_hash(double *cup);

#endif
