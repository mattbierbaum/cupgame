import numpy as np
import scipy as sp
import pylab as pl
from scipy.optimize import fminbound
import matplotlib as mpl
import mpl_toolkits.mplot3d as mp3d
import mpl_toolkits.mplot3d.axes3d as a3
from scipy import weave

# kinematic equations of motion
def position(x0, v0, t):
    a = np.array([0,0,-1])
    return x0+v0*t+0.5*a*t**2

def velocity(v0, t):
    a = np.array([0,0,-1])
    return v0 + a*t

# this is the transform which determines the lattice of toruses
def normpos_square(r):
    rt = r.copy()
    cc = rt[:2].astype('int')
    rt[:2] -= 2*((cc+(rt[:2]>0))/2)
    return rt

def normpos_triangle(r):
    rt3 = np.sqrt(3)

    rt = r.copy()
    mx,my = rt[0] % 1,           rt[1] % np.sqrt(3)      #shifted to corner
    dx,dy = mx - int(rt[0]/1)*1, my - int(rt[1]/rt3)*rt3 #shift in coordinates 

    pts = [[0,0],[1,0],[0,rt3],[1,rt3],[0.5,rt3/2]]

    dists = np.array([(l[0]-mx)**2 + (l[1]-my)**2 for l in pts])
    mi = np.argmin(dists)
    dx -= pts[mi][0]
    dy -= pts[mi][1]

    rt[:2] = np.array([dx, dy])
    return rt

normpos = normpos_square

# the distance functions for toruses
def distance_to_torus2(t, x0, v0, hsq, rsq):
    r = position(x0, v0, t)
    r = normpos(r)
    return abs((r.dot(r) + hsq - rsq)**2 - 4*hsq*(r[0]*r[0] + r[1]*r[1]))

def distance_to_torus(r, hsq, rsq):
    r = normpos(r)
    return (r.dot(r) + hsq - rsq)**2 - 4*hsq*(r[0]*r[0] + r[1]*r[1])

def inside_torus(r, hsq, rsq):
    return np.sign(distance_to_torus(r, hsq, rsq))

# finds the collision time for an initial condition
# by finding the sign crossing and minimizing in those bounds
def collision_time(pos, vel, hsq, rsq):
    px0,py0,pz0 = pos 
    vx0,vy0,vz0 = vel
    dt, t = 5e-2, 0
    
    sign = inside_torus(pos, hsq, rsq)
    prevsign = sign

    x = px0 + vx0*t
    y = py0 + vy0*t
    z = pz0 + vz0*t - 0.5*t**2

    while z > 0 and prevsign == sign:
        prevsign = sign
        t = t + dt
        x = px0 + vx0*t
        y = py0 + vy0*t
        z = pz0 + vz0*t - 0.5*t**2
        sign = inside_torus(np.array([x,y,z]), hsq, rsq)

    if sign != prevsign:
        topt, fopt, iter, funcalls = fminbound(distance_to_torus2, t-dt, t, 
            args=(pos, vel, hsq, rsq), disp=0, full_output=1, xtol=1e-5)
    else:
        topt = t
        fopt = distance_to_torus(np.array([x,y,z]), hsq, rsq) 
    return topt,fopt 

def collision_normal(pos, vel, h):
    x,y,z = normpos(pos)
    theta = np.arctan2(y,x)
    yin = h*np.sin(theta)
    xin = h*np.cos(theta)
   
    dr = pos-np.array([xin, yin, 0])
    dr = dr / np.sqrt(dr.dot(dr))
    return dr 

def rotate_about_axis(r, axis, theta):
    return r*np.cos(theta) + np.cross(axis, r)*np.sin(theta) + axis*(axis.dot(r))*(1-np.cos(theta))

def reflect_velocity(pos, vel, h):
    ax = collision_normal(pos, vel, h)
    return -rotate_about_axis(vel, ax, np.pi)

def plotPoints(x,y,z):
    fig = pl.figure()
    axs = a3.Axes3D(fig)
    axs.scatter(x,y,z)
    pl.show()
    
def plotPorcupine(lines, x, y, z):
    segs = mp3d.art3d.Line3DCollection(lines)
    fig = pl.figure()
    axs = a3.Axes3D(fig)
    axs.scatter(x,y,z)
    axs.add_collection3d(segs)
    pl.show()

def testscan():
    m = 1.0
    N = 40
    r = 0.15
    h = 0.5-r

    hsq, rsq = h*h, r*r

    x0 = np.array([0, 0, 1])
    v0 = np.array([0, 0, 0])
    out = np.zeros((N,N))
    angles = np.mgrid[-m:m:N*1j, -m:np.sqrt(3)*m:N*1j].T
   
    segs, xs, ys, zs = [], [], [], []
    for i in xrange(N):
        for j in xrange(N):
            xt = x0+np.array([angles[i,j,0], angles[i,j,1], 0])
            tcoll, closest_approach = collision_time(xt, v0, hsq=hsq, rsq=rsq)

            if abs(closest_approach) < 1e-4:
                pos = position(xt, v0, tcoll)
                vel = velocity(v0, tcoll)
                #segs.append(collision_resultant(pos, vel, h))
                xs.append(pos[0])
                ys.append(pos[1])
                zs.append(pos[2])
                out[i,j] = pos[2]
   
    #plotPorcupine(segs,xs,ys,zs)
    plotPoints(xs,ys,zs)
    #pl.figure()
    #pl.imshow(out)
    return out

def orbits():
    # tirangular
    m = 0.5
    N = 40
    r = 0.15
    h = 0.5-r
   
    #square
    m = 1.0
    N = 50
    r = 0.3
    h = 1-r

    hsq, rsq = h*h, r*r

    x0 = np.array([0, 0, 2])
    v0 = np.array([0, 0, 0])
    bounces = np.zeros((N,N))
    inout = np.zeros((N,N))
    angles = np.mgrid[0:m:N*1j, 0:m:N*1j].T
   
    xs, ys, zs = [], [], []
    for i in xrange(N):
        print i
        for j in xrange(i,N):
            offset = np.array([angles[i,j,0], angles[i,j,1], 0])
            if np.sqrt(offset[0]**2 + offset[1]**2)  - h > r:
                continue
            pos = x0+offset
            vel = v0
            tbounces = 0
            while 1:
                # get the next collision
                tcoll, closest_approach = collision_time(pos, vel, hsq=hsq, rsq=rsq)

                # figure out where it hit and what speed
                pos = position(pos, vel, tcoll)
                vel = velocity(vel, tcoll)

                if closest_approach < 2e-4:
                    tbounces += 1
                    xs.append(pos[0])
                    ys.append(pos[1])
                    zs.append(pos[2])
                    vel = 0.9*reflect_velocity(pos, vel, h)
                    pos = position(pos, vel, 2e-2)
                else:
                    pos = normpos(pos)
                    if pos.dot(pos) < hsq:
                        inout[i,j] = 1
                    else:
                        inout[i,j] = 0
                    break
                if tbounces > 5000:
                    print "too many bounces"
                    break
            bounces[i,j] = tbounces
   
    #plotPoints(xs,ys,zs)
    return bounces, inout


def corbits():
    #square
    m = 1.0
    N = 500
    r = 0.3
    h = 1-r

    x0 = np.array([0, 0, 2]).astype('double')
    v0 = np.array([0, 0, 0]).astype('double')
    bounces = np.zeros((N,N))
    angles = np.mgrid[0:m:N*1j, 0:m:N*1j].T
    angles = np.mgrid[-m:m:N*1j, -m:m:N*1j].T
  
    with open("cupgamefinal.c") as f:
        import hashlib
        clib = "\n".join(f.readlines())
        sha1 = hashlib.sha1(clib)
        uuid = sha1.hexdigest()

    for i in xrange(N):
        print i
        for j in xrange(N):
            offset = np.array([angles[i,j,0], angles[i,j,1], 0])
            if np.sqrt(offset[0]**2 + offset[1]**2)  - h > r:
                continue
            pos = x0+offset
            vel = v0
            tbounces = np.array([0])
            weave.inline("tbounces[0] = trackCollisions(pos, vel, h, r, 10000); char id[] = \"%s\";" % uuid, 
                    ["tbounces", "pos", "vel", "h", "r"], support_code=clib,
                    libraries=["m", "gsl", "blas"])
            bounces[i,j] = tbounces[0]
   
    return bounces

