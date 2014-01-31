#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "cupgamelib.h"

// h is R, the radius of the hole
// r is the radius of the ring

//============================================================================
// Root finder.  These functions create and solve for the roots of an 8th
// order polynomial which describes the intersection of a torus and parabola
//============================================================================
double polyeval(double *poly, int deg, double x){
    int i;
    double out = 0.0;
    for (i=deg; i>=0; i--) out = out*x + poly[i];
    return out;
}

void find_root_pair(double *poly, int deg, double *qpoly, 
        double *r1, double *r2){
    int i;
    double c,d,g,h,u,v,det,err,vo,uo;
    double rpoly[DEGSIZE];
   
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

    // b^2 - 4*a*c
    double desc = u*u - 4*v; 

    // if we have imaginery roots, quit now
    if (desc < 0){
        *r1 = NAN; *r2 = NAN;
        return;
    }

    double sqt = 0.5*sqrt(desc);
    double boa = -0.5*u;
    
    *r1 = (boa + sqt);
    *r2 = (boa - sqt);
}

void find_all_roots(double *poly, int deg, double *roots, int *nroots){
    int i;
    double r1, r2;
    double temppoly[DEGSIZE];

    *nroots = 0;
    for (i=deg; i>=1; i-=2){
        find_root_pair(poly, i, temppoly, &r1, &r2);
        memcpy(poly, temppoly, sizeof(double)*deg);

        if (isnan(r1) || isnan(r2))
            continue;

        roots[*nroots] = r1; *nroots += 1;
        roots[*nroots] = r2; *nroots += 1;
    }
}

double find_smallest_root(double *poly, int deg){
    int i=0, nroots=0;
    double realroots[DEGSIZE];

    find_all_roots(poly, deg, realroots, &nroots);

    // if we didn't find any roots, return a nan (only special number)
    if (nroots == 0) return NAN;

    // otherwise, find the root closest to zero
    double minroot = NAN;
    for (i=0; i<nroots; i++)
        if ((isnan(minroot) || minroot > realroots[i]) && realroots[i] > 0) 
            minroot = realroots[i];

    return minroot;
}

void build_torus_poly(double *pos, double *vel, 
        double h, double r, double *cup, double *poly){
    double r2,h2,x2,y2,z2,xvx,yvy,zvz,xyzv,vx2,vy2,vz2,vx,vy,vz,x,y,z;
    vx = vel[0]; vy = vel[1]; vz = vel[2];
    x = pos[0]-cup[0];  
    y = pos[1]-cup[1];  
    z = pos[2];
    r2 = r*r; h2 = h*h; 
    x2 = x*x; y2 = y*y; z2 = z*z;
    vx2 = vx*vx; xvx = x*vx;
    vy2 = vy*vy; yvy = y*vy;
    vz2 = vz*vz; zvz = z*vz;
    xyzv = (xvx+yvy+zvz);
    
    poly[8] = 1.0/16;
    poly[7] = -0.5*vz;
    poly[6] = 0.5*(vx2 + vy2 + 3*vz2 - z);
    poly[5] = -2*(vx2 + vy2)*vz + (xvx + yvy) -2*vz*vz2 + 3*zvz;
    poly[4] = (vx2*vx2 + 2*(vy2 + vz2)*vx2 -2*z*vx2 - 4*vz*xvx
        + vy2*vy2 + vz2*vz2 + 0.5*(-r2 + h2 + x2 + y2 + 3*z2) 
        - 2*z*vy2 + 2*vy2*vz2 - 6*z*vz2 - 4*yvy*vz);
    poly[3] = (4*xyzv*vx2 + 4*xyzv*vy2 + 4*xyzv*vz2 - 4*z*xvx - 4*z*yvy + 2*r2*vz
        - 2*h2*vz - 2*x2*vz - 2*y2*vz - 6*z2*vz);
    poly[2] = (-2*z2*z + 2*vx2*z2 + 2*vy2*z2 + 6*vz2*z2 + 2*r2*z - 2*h2*z - 2*x2*z 
        - 2*y2*z + 8*xvx*zvz + 8*yvy*zvz - 2*r2*vx2 - 2*h2*vx2 + 6*x2*vx2 
        + 2*y2*vx2 - 2*r2*vy2 - 2*h2*vy2 + 2*x2*vy2 + 6*y2*vy2 - 2*r2*vz2 + 2*h2*vz2 
        + 2*x2*vz2 + 2*y2*vz2 + 8*xvx*yvy);
    poly[1] = (4*xyzv*x2 - 4*r2*xvx - 4*h2*xvx + 4*y2*xvx + 4*z2*xvx + 4*y2*yvy + 4*z2*yvy 
        - 4*r2*yvy - 4*h2*yvy + 4*z2*zvz - 4*r2*zvz + 4*h2*zvz + 4*y2*zvz);
    poly[0] = (r2*r2 + h2*h2 + x2*x2 + y2*y2 + z2*z2 - 2*r2*h2 - 2*r2*x2 - 2*h2*x2 - 2*r2*y2
        - 2*h2*y2 + 2*x2*y2 - 2*r2*z2 + 2*h2*z2 + 2*x2*z2 + 2*y2*z2);
}

void build_cylinder_poly(double *pos, double *vel, double r, 
        double cylx, double cyly, double *poly){
    double vx,vy,x,y;
    vx = vel[0];    vy = vel[1]; 
    x = pos[0]-cylx; y = pos[1]-cyly; 
    poly[2] = vx*vx + vy*vy;
    poly[1] = 2*vx*x + 2*vy*y;
    poly[0] = x*x + y*y - r*r;
}

void build_zero_poly(double *pos, double *vel, double *poly){
    poly[0] = pos[2]; poly[1] = vel[2]; poly[2] = -0.5;
}

void build_plane_poly(double *pos, double *vel, double r, double *poly){
    poly[0] = pos[2]-r-1e-6; poly[1] = vel[2]; poly[2] = -0.5;
}

//============================================================================
// These are helper functions that calculate positions of cups and 
// the positions of trajectories
//============================================================================
inline double mymod(double a, double b){
      return a - b*(int)(a/b) + b*(a<0);
}

inline double dot(double *r1, double *r2){
    return r1[0]*r2[0] + r1[1]*r2[1] + r1[2]*r2[2];
}

inline void cross(double *a, double *b, double *out){
    out[0] = a[1]*b[2]-a[2]*b[1];
    out[1] = a[2]*b[0]-a[0]*b[2];
    out[2] = a[0]*b[1]-a[1]*b[0];
}

inline void position(double *x0, double *v0, double t, double *out){
    int i;
    double a[] = {0,0,-1};
    for (i=0; i<3; i++) out[i] = x0[i] + v0[i]*t + 0.5*a[i]*t*t;
}

inline void velocity(double *v0, double t, double *out){
    int i;
    double a[] = {0,0,-1};
    for (i=0; i<3; i++) out[i] = v0[i] + a[i]*t;
}

double distance_to_torus(double *pos, double h, double r, double *cup){
    double dr[3];
    dr[0] = pos[0] - cup[0];
    dr[1] = pos[1] - cup[1];
    dr[2] = pos[2];
    double term1 = (h - sqrt(dr[0]*dr[0] + dr[1]*dr[1]));
    return (term1*term1 + dr[2]*dr[2] - r*r);
}

//============================================================================
// Cup neighborlist generator - finite and infinite sets
//============================================================================
// a square grid of cups at even grid points
void cups_near(double *pos, double *cups, int *ncups){
    int i, j;
    *ncups = 0;
    for (i=-1; i<=1; i++){
        for (j=-1; j<=1; j++){
            cups[2*(*ncups)+0] = pos[0]-(mymod(pos[0]+1, 2)-1) + 2*i;
            cups[2*(*ncups)+1] = pos[1]-(mymod(pos[1]+1, 2)-1) + 2*j;
            *ncups += 1;
        }
    }
}

void cups_near_hex(double *pos, double *cups, int *ncups){
    *ncups = 0;

    int i=0, j=0;
    double rt3 = sqrt(3.0);
    int mp = (int)round(pos[1]/rt3);
    int np = (int)round((pos[0]-mp)/2);

    for (i=-1; i<=1; i++){
        for (j=-1; j<=1; j++){
            cups[2*(*ncups)+0] = 2*(i+np) + (j+mp);
            cups[2*(*ncups)+1] = rt3*(j+mp);
            *ncups += 1;
        }
    }
}

void cups_near_6cup(double *pos, double *cups, int *ncups){
    if (fabs(pos[0]) > 6.0 || fabs(pos[1]) > 6.0){
        *ncups = 0;
        return;
    }

    *ncups = 6;
    double rt3 = sqrt(3.0);

    cups[2*0+0] = 0.0;  cups[2*1+0] = 0.0;
    cups[2*0+1] = 0.0;  cups[2*1+1] = 2.0*rt3;

    cups[2*2+0] = -2.0; cups[2*3+0] = 2.0;
    cups[2*2+1] = 0.0;  cups[2*3+1] = 0.0;

    cups[2*4+0] = -1.0; cups[2*5+0] = 1.0;
    cups[2*4+1] = rt3;  cups[2*5+1] = rt3;

    for (int i=0; i<*ncups; i++)
        cups[2*i+1] -= rt3/2;
}

void collision_normal(double *pos, double h, double *cup, double *out){
    double t[3];
    t[0] = pos[0] - cup[0];
    t[1] = pos[1] - cup[1];
    t[2] = pos[2];
    
    double theta = atan2(t[1],t[0]);
    double xin = h*cos(theta);
    double yin = h*sin(theta);
  
    out[0] = t[0]-xin;
    out[1] = t[1]-yin;
    out[2] = t[2];

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

void reflect_velocity(double *pos, double *vel, double h, 
        double f, double *cup, double *out){
    double ax[3]; 
    collision_normal(pos, h, cup, ax);
    rotate_about_axis(vel, ax, M_PI, out);
    out[0] *= -f; out[1] *= -f; out[2] *= -f;
}


//============================================================================
// finds the collision time for an initial condition
// by finding the roots of a poly and finding the nearest collision
//============================================================================
int collides_with_cup(double *pos, double *vel, double h, double r, 
        double *cup, double *tcoll){
    /* 
     * This functions determines whether a particular trajectory collides
     * with the cup specified by h, r, (cx, cy).  It returns:
     *  0 : There was no collision 
     *  1 : There was a collision and it occured at tcoll
     *  2 : The ball passed through the center of the cup
     * It does not modify the values of pos, vel; modified tcoll
    */
    double tpos[3], tvel[3], poly[DEGSIZE];
    memcpy(tpos, pos, sizeof(double)*3);
    memcpy(tvel, vel, sizeof(double)*3);
    
    // try to find a collision with the object
    build_torus_poly(tpos, tvel, h, r, cup, poly);
    *tcoll = find_smallest_root(poly, DEG);

    // we didn't collide with the torus.  did we fall through?
    if (isnan(*tcoll)){
        build_zero_poly(tpos, tvel, poly);
        *tcoll = find_smallest_root(poly, 2);

        position(tpos, tvel, *tcoll, tpos); 

        // there is no collision if we fall outside of the torus
        if ((tpos[0]-cup[0])*(tpos[0]-cup[0])+(tpos[1]-cup[1])*(tpos[1]-cup[1]) > (h+r)*(h+r))
            return RESULT_NOTHING;
        return RESULT_INCUP;
    } 
    return RESULT_COLLISION;
}

int collision_near(double *pos, double *vel, double h, double r,
        double *tcoll, double *cup){
    int i, result, ncups, event;
    double tevent, cups[MAXCUPS];

    event = RESULT_NOTHING;
    tevent = NAN;

    cups_near_hex(pos, cups, &ncups);
    for (i=0; i<ncups; i++){
        result = collides_with_cup(pos, vel, h, r, &cups[2*i], tcoll);
        if (result == RESULT_COLLISION || result == RESULT_INCUP){
            if ((isnan(tevent) || tevent > *tcoll) && *tcoll > 0){
                cup[0] = cups[2*i+0];
                cup[1] = cups[2*i+1];
                event = result;
                tevent = *tcoll;
            }
        }
    }
    *tcoll = tevent;
    return event;
}

int collision_time(double *pos, double *vel, double h, 
        double r, double *tcoll, double *cup){
    int result, napproaches, i;
    double tpos[3], tvel[3], poly[DEGSIZE], tcloseapproaches[2];

    int event = RESULT_NOTHING;
    double tevent = NAN;

    // find where the trajectory intersects the z=h plane
    build_plane_poly(pos, vel, r, poly);
    find_all_roots(poly, 2, tcloseapproaches, &napproaches);
   
    for (i=0; i<napproaches; i++){
        if (tcloseapproaches[i] < 0) continue;
        position(pos, vel, tcloseapproaches[i], tpos);
        velocity(vel, tcloseapproaches[i], tvel);
        result = collision_near(tpos, tvel, h, r, tcoll, cup);
        *tcoll = *tcoll + tcloseapproaches[i];

        if (result == RESULT_COLLISION || result == RESULT_INCUP){
            if ((isnan(tevent) || tevent > *tcoll) && *tcoll > 0){
                event = result;
                tevent = *tcoll;
            }
        }
    }

    result = collision_near(pos, vel, h, r, tcoll, cup);
    if (result == RESULT_COLLISION || result == RESULT_INCUP){
        if ((isnan(tevent) || tevent > *tcoll) && *tcoll > 0){
            event = result;
            tevent = *tcoll;
        }
    }

    if (event == RESULT_NOTHING){
        build_zero_poly(pos, vel, poly);
        find_all_roots(poly, 2, tcloseapproaches, &napproaches);
        tevent = tcloseapproaches[0] > 0 ? tcloseapproaches[0] : tcloseapproaches[1];
    }

    *tcoll = tevent;
    return event;
}


double singleCollisions(int NP, double *pos, int NV, double *vel, double h, double r){
    double tcoll, cup[2];
    if (NP != 3 || NV != 3)
        return -1;

    collision_time(pos, vel, h, r, &tcoll, cup);
    position(pos, vel, tcoll, pos);
    return tcoll;
}

int cup_hash(double *cup){
    return (int)(100*sqrt(cup[0]*cup[0]+cup[1]*cup[1])*atan2(cup[0],cup[1]));
}

int trackCollisions(int NP, double *pos, int NV, double *vel, 
        double h, double r, double damp, int maxbounces){
    if (NP != 3 || NV != 3)
        return -1;

    int result;
    double factor = damp;
    double tcoll, vlen;

    double tpos[3], tvel[3], cup[2]; 
    memcpy(tpos, pos, sizeof(double)*3);
    memcpy(tvel, vel, sizeof(double)*3);

    int tbounces = 0;
    while (tbounces < maxbounces){
        // get the next collision
        result = collision_time(tpos, tvel, h, r, &tcoll, cup);

        if (result == RESULT_NOTHING)
#ifdef RECORD_MISSES
            break;
#else
            return 0;
#endif
        if (result == RESULT_INCUP)
#ifdef RECORD_CUPINDEX
            return cup_hash(cup);
#else
            break;
#endif

        // figure out where it hit and what speed
        position(tpos, tvel, tcoll, tpos);
        velocity(tvel, tcoll, tvel);
        vlen = dot(tvel, tvel);

        if (tpos[2] < 0 || vlen < EPS) break;

        reflect_velocity(tpos, tvel, h, factor, cup, tvel);
        position(tpos, tvel, EPS, tpos); 
        velocity(tvel, EPS, tvel);
        tbounces++;
    }

#ifdef RECORD_CUPINDEX
    return 0;
#else
    return tbounces;
#endif
}

/*void trackSlice(int NP, double *pos, int NV, double *vel,
        double h, double r, double damp, int maxbounces, int NS, int *bounces){
    if (NP != NV || NP/3 != NS){
        printf("Incorrect dimensions for input arrays");
        return;
    }

    //#pragma omp parallel for
    //for (int i=0; i
}*/

int trackTrajectory(int NP, double *pos, int NV, double *vel, double h, double r, 
        double damp, int maxbounces, int NT, double *traj){
    if (NP != 3 || NV != 3)
        return -1;

    int result, clen=0;
    double factor = damp;
    double tcoll;

    double tpos[3], tvel[3], cup[2]; 
    memcpy(tpos, pos, sizeof(double)*3);
    memcpy(tvel, vel, sizeof(double)*3);

    double t = 0.0;
    double ttpos[3];
    memcpy(tpos, pos, sizeof(double)*3);
    memcpy(tvel, vel, sizeof(double)*3);

    int tbounces = 0;
    while (tbounces < maxbounces){
        // get the next collision
        result = collision_time(tpos, tvel, h, r, &tcoll, cup);

        // FIXME - only draw 2d lines between collision times
        for (t=0; t<tcoll; t+=(tcoll/TSAMPLES)){
            position(tpos, tvel, t, ttpos);
            if (NT >= 0 && clen < NT/3-3){
                memcpy(traj+3*clen, ttpos, sizeof(double)*3);
                clen += 1;
            }
        }

        if (result == RESULT_NOTHING || result == RESULT_INCUP)
            break;

        // figure out where it hit and what speed
        position(tpos, tvel, tcoll, tpos);
        velocity(tvel, tcoll, tvel);

        if (NT >= 0 && clen < NT/3-3){
            memcpy(traj+3*clen, tpos, sizeof(double)*3);
            clen += 1;
        }

        //double dist = distance_to_torus(tpos, h, r, cup);
        //double vlen = dot(tvel, tvel);
        //printf("%i %e %e %e\n", result, tcoll, dist, vlen);

        if (result == RESULT_COLLISION){
            reflect_velocity(tpos, tvel, h, factor, cup, tvel);

            // if display normals
            if (0){
                double ax[3];
                collision_normal(tpos, h, cup, ax);
                if (NT >= 0 && clen < NT/3-9){
                    memcpy(traj+3*clen, tpos, sizeof(double)*3);clen += 1;
                    ttpos[0] = tpos[0]+ax[0];
                    ttpos[1] = tpos[1]+ax[1];
                    ttpos[2] = tpos[2]+ax[2];
                    memcpy(traj+3*clen, ttpos, sizeof(double)*3);clen += 1;
                    memcpy(traj+3*clen, tpos, sizeof(double)*3);clen += 1;
                }
            }

            tbounces++;
        }

        position(tpos, tvel, EPS, tpos); 
        velocity(tvel, EPS, tvel);
    }
    return 3*clen;
}

