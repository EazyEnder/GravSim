import numpy as np
from Particle import Particle
from Vector2 import Vector2
from math import *
from utils import printProgressBar
from QuadTree.QuadTree import QuadTree
from QuadTree.QRegion import QRegion
from plotsUtils import plotGraph

SIZE = 2500

P_NUMBER = 50000
P_MASS = 1

from random import random


dt = 1.

r_angle = np.random.rand(P_NUMBER)*2*np.pi
r_radius = np.sqrt(np.random.rand(P_NUMBER))*SIZE/4
particles_pos = np.array([[r_radius[i]*np.cos(r_angle[i])+SIZE/2,SIZE/2+r_radius[i]*np.sin(r_angle[i])] for i in range(P_NUMBER)])
BASE_VELOCITY = 0.010
p_vel_x = -np.sin(r_angle)*r_radius*BASE_VELOCITY
p_vel_y = np.cos(r_angle)*r_radius*BASE_VELOCITY

qt = QuadTree(QRegion(Vector2(SIZE/2,SIZE/2),SIZE))
for i in range(P_NUMBER):
    qt.insert(Particle(P_MASS,particles_pos[i][0],particles_pos[i][1],p_vel_x[i],p_vel_y[i],dt=dt))


#Integration de Verlet
a=0.3 #"amortissement"

def compute(qt):
    updated_particles = []
    for p1 in qt.getAllParticles():
        f_sum = p1.host.getForce(p1,precision=1.)

        #Verlet
        new_x = 2 * p1.r.x - p1.old_r.x + f_sum.x/p1.m*dt**2
        new_y = 2 * p1.r.y - p1.old_r.y + f_sum.y/p1.m*dt**2
        #p1.save()
        p = p1.clone()
        p.old_r = p1.r
        p.r = Vector2(new_x,new_y)
        updated_particles.append(p)
    return updated_particles

T_MAX = 1000
t = 0

import matplotlib.pyplot as plt
fig, ax = plt.subplots(1,1)
ax.set_aspect('equal', adjustable='box')

offset = SIZE/2
while(t < T_MAX):
    particles = compute(qt)
    qt = QuadTree(QRegion(Vector2(SIZE/2,SIZE/2),SIZE))
    for p in particles:
        qt.insert(p)
    printProgressBar(t,T_MAX)
    ax.cla()
    
    plotGraph([[p.r.x+offset, p.r.y+offset] for p in particles], SIZE+offset*2, frame_size=500 ,ax=ax, periodic=False)
    #plt.pause(0.1)
    plt.savefig(fname="exportOctTree/fig1_2D"+str(int(t))+".jpg")
    t += dt

plt.show()