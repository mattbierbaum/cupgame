#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <stdlib.h>
#include <math.h>

typedef double (scalarFunction)(double, void *);
double fminbound(scalarFunction *func, double a, double b, void *args, int max_iter, double xtol);

inline double mymod(double a, double b){
      return a - b*(int)(a/b) + b*(a<0);
}

inline int signum(double a){
    return (a>0) - (a<0);
}

double dot(double *r1, double *r2){
    return r1[0]*r2[0] + r1[1]*r2[1] + r1[2]*r2[2];
}

void cross(double *a, double *b, double *out){
    out[0] = a[1]*b[2]-a[2]*b[1];
    out[1] = a[2]*b[0]-a[0]*b[2];
    out[2] = a[0]*b[1]-a[1]*b[0];
}

void position(double *x0, double *v0, double t, double *out){
    double a[] = {0,0,-1};

    int i;
    for (i=0; i<3; i++)
        out[i] = x0[i] + v0[i]*t + 0.5*a[i]*t*t;
}

void velocity(double *v0, double t, double *out){
    double a[] = {0,0,-1};

    int i;
    for (i=0; i<3; i++)
        out[i] = v0[i] + a[i]*t;
}

void normpos2(double *r, double *out){
    int i;
    out[2] = r[2];
    for (i=0; i<2; i++){
        out[i] = r[i];
        out[i] = mymod(out[i]+1, 2) - 1; 
    }
}

void normpos(double *r, double *out){
    int i,j;
    double rt3 = sqrt(3.0);

    for (i=0; i<3; i++)
        out[i] = r[i];

    double mx[] = {mymod(out[0], 1), mymod(out[1],rt3)};
    double dx[] = {mx[0] - int(out[0]/1)*1, mx[1] - int(out[1]/rt3)*rt3};
    double pts[][2] = {{0,0},{1,0},{0,rt3},{1,rt3},{0.5,rt3/2}};

    double mindist = 1e10;
    int minindex = 0;

    for (i=0; i<5; i++){
        double dist = pow(pts[i][0]-mx[0],2) + pow(pts[i][1]-mx[1],2);
        if (dist < mindist){
            mindist = dist;
            minindex = i;
        }
    }
    dx[0] -= pts[minindex][0];
    dx[1] -= pts[minindex][1];
 
    out[0] = dx[0];
    out[1] = dx[1];
}

int getcup(double *r, double h){
    int out[2];
    int i;
    for (i=0; i<2; i++)
        out[i] = 2*(((int)(r[i])+(r[i]>0))/2);
    double r2[3];
    normpos(r, r2);
    double r2r = sqrt(r2[0]*r2[0] + r2[1]*r2[1]);
    return (out[0]*out[0] + out[1]*out[1] + 1)*(r2r < h); 
}

double distance_to_torus(double *pos, double hsq, double rsq){
    double r[3];
    normpos(pos, r);
    double term1 = (dot(r,r) + hsq - rsq);
    double term2 = 4*hsq*(r[0]*r[0] + r[1]*r[1]);
    double dist = term1*term1 - term2;
    return dist;
}

double distance_to_torus_t(double t, void *args){
    double *dargs = (double*)args;
    double x0[3]; double v0[3]; double hsq; double rsq;
    x0[0] = dargs[0]; x0[1] = dargs[1];  x0[2] = dargs[2];
    v0[0] = dargs[3]; v0[1] = dargs[4];  v0[2] = dargs[5];
    hsq = dargs[6]; rsq = dargs[7];

    double r[3];
    position(x0, v0, t, r);
    double dist = fabs(distance_to_torus(r, hsq, rsq));
    return dist;
}

// finds the collision time for an initial condition
// by finding the sign crossing and minimizing in those bounds
double collision_time(double *pos, double *vel, double hsq, double rsq, double *dist,
        double *traj, int *clen, int mlen, double xtol){
    double px0=pos[0]; double py0=pos[1]; double pz0=pos[2];
    double vx0=vel[0]; double vy0=vel[1]; double vz0=vel[2];
    double dt=5e-3;
    double t=0;
    
    double tdist = distance_to_torus(pos, hsq, rsq);
    double vlen = sqrt(vx0*vx0 + vy0*vy0 + (vz0-t)*(vz0-t));
    int sign = tdist > 0;
    int prevsign = sign;

    double x = px0 + vx0*t;
    double y = py0 + vy0*t;
    double z = pz0 + vz0*t - 0.5*t*t;

    double tmp[3]; 
    double vtmp[3];
    while (z > 0 && prevsign == sign){
        prevsign = sign;
        t = t + dt;//0.1*fabs(tdist/vlen) + 1e-3;
        x = px0 + vx0*t;
        y = py0 + vy0*t;
        z = pz0 + vz0*t - 0.5*t*t;

        /*if (mlen >= 0 && *clen < mlen){
            traj[*clen*3+0] = x;
            traj[*clen*3+1] = y;
            traj[*clen*3+2] = z;
            clen[0] += 1;
        }*/

        tmp[0] = x; tmp[1] = y; tmp[2] = z;
        tdist = distance_to_torus(tmp, hsq, rsq);
        vlen = sqrt(vx0*vx0 + vy0*vy0 + (vz0-t)*(vz0-t));

        sign = tdist > 0;
    }

    double topt, fopt;
    if (sign != prevsign){
        double args[8];
        args[0] = pos[0]; args[1] = pos[1]; args[2] = pos[2];
        args[3] = vel[0]; args[4] = vel[1]; args[5] = vel[2];
        args[6] = hsq;    args[7] = rsq;
        topt = fminbound(distance_to_torus_t, t-dt, t, (void*)args, 100, xtol);  //FIXME - magic number - related to xtol
    }
    else {
        topt = t;
    }

    position(pos, vel, topt, tmp);
    fopt = fabs(distance_to_torus(tmp, hsq, rsq));
    *dist = fopt;

    return topt;
}

void collision_normal(double *pos, double *vel, double h, double *out){
    double t[3];
    normpos(pos, t);
    double theta = atan2(t[1],t[0]);
    double yin = h*sin(theta);
    double xin = h*cos(theta);
  
    out[0] = pos[0]-xin;
    out[1] = pos[1]-yin;
    out[2] = pos[2];

    double len = sqrt(dot(out,out));
    out[0] /= len;
    out[1] /= len;
    out[2] /= len;
}

void rotate_about_axis(double *r, double *axis, double theta, double *out){
    int i;
    double crs[3]; cross(axis, r, crs);
    double dt = dot(axis, r);
    for (i=0; i<3; i++)
        out[i] = r[i]*cos(theta) + crs[i]*sin(theta) + axis[i]*dt*(1-cos(theta));
}

void reflect_velocity(double *pos, double *vel, double h, double f, double *out){
    double ax[3]; collision_normal(pos, vel, h, ax);
    rotate_about_axis(vel, ax, M_PI, out);
    out[0] *= -f; out[1] *= -f; out[2] *= -f;
}

double singleCollisions(double *x0, double *v0, double h, double r){
    int i;
    double hsq = h*h;
    double rsq = r*r;
    double xtol = 1e-4;

    double pos[3]; double vel[3]; 
    for (i=0; i<3; i++){ pos[i]=x0[i]; vel[i]=v0[i]; }

    double tcoll, closest_approach;
    tcoll = collision_time(pos, vel, hsq, rsq, &closest_approach, NULL, NULL, -1, xtol);

    position(pos, vel, tcoll, pos);
    if (closest_approach < 1e-2 && pos[2] > 0)
        return pos[2];
    else
        return 0;
}

int trackCollisions(double *x0, double *v0, double h, double r, double damp, int maxbounces,
        double *traj, int *clen, int mlen){
    int i;
    double hsq = h*h;
    double rsq = r*r;
    double factor = damp;
    double xtol = 1e-6;

    double pos[3]; double vel[3]; 
    for (i=0; i<3; i++){ pos[i]=x0[i]; vel[i]=v0[i]; }

    int tbounces = 0;
    while (tbounces < maxbounces){
        // get the next collision
        double tcoll, closest_approach;
        tcoll = collision_time(pos, vel, hsq, rsq, &closest_approach, traj, clen, mlen, xtol);

        // figure out where it hit and what speed
        position(pos, vel, tcoll, pos);
        velocity(vel, tcoll, vel);

        if (closest_approach < 1e-2 && pos[2] > 0){  //FIXME - was 2e-2 (this is good) should be smaller - error in collision
            tbounces += 1;
            reflect_velocity(pos, vel, h, factor, vel);
            position(pos, vel, 2e-2, pos); ///FIXME - another magic number 

            // FIXME - remove
            /*if (mlen >= 0 && *clen < mlen){
                traj[*clen*3+0] = pos[0];
                traj[*clen*3+1] = pos[1];
                traj[*clen*3+2] = pos[2];
                clen[0] += 1;
            }*/

            //tbounces = getcup(pos, h);
        }
        else {
            normpos(pos, pos);
            break;
        }
    }
    return tbounces;
}

void handler (const char * reason,
        const char * file,
        int line,
        int gsl_errno){
    //printf("%s:%i %s\n", file, line, reason);
}

double fminbound(scalarFunction *func, double a, double b, void *args, int max_iter, double xtol){
    gsl_set_error_handler(handler);

    int status;
    int iter = 0;
    double m = (b+a)/2;
    const gsl_min_fminimizer_type *T;
    gsl_min_fminimizer *s;
    gsl_function F;

    F.function = func;
    F.params = args;

    T = gsl_min_fminimizer_brent;
    T = gsl_min_fminimizer_quad_golden;
    s = gsl_min_fminimizer_alloc(T);
    gsl_min_fminimizer_set(s, &F, m, a, b);

    do {
        iter++;
        status = gsl_min_fminimizer_iterate(s);

        m = gsl_min_fminimizer_x_minimum(s);
        a = gsl_min_fminimizer_x_lower(s);
        b = gsl_min_fminimizer_x_upper(s);

        status = gsl_min_test_interval(a, b, xtol, 0.0);

        //if (status == GSL_SUCCESS)
        //    printf ("Converged:\n");

    }
    while (status == GSL_CONTINUE && iter < max_iter);

    gsl_min_fminimizer_free(s);
    return m;
}

