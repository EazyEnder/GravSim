import numpy as np
from Particle import Particle
from Vector2 import Vector2
from math import *
from utils import printProgressBar
from QuadTree.QuadTree import QuadTree
from QuadTree.QRegion import QRegion
from plotsUtils import plotGraph

SIZE = 100

P_NUMBER = 50
P_MASS = 5

G = 0.5

from random import random


dt = 0.5

qt = QuadTree(QRegion(Vector2(SIZE/2,SIZE/2),SIZE))
for i in range(P_NUMBER):
    qt.insert(Particle(P_MASS,SIZE*random(),SIZE*random(),0,0,dt=dt))


#Integration de Verlet
a=0.3 #"amortissement"

def compute(qt):
    updated_particles = []
    for p1 in qt.getAllParticles():
        f_sum = p1.host.getForce(p1)

        #Verlet
        new_x = 2 * p1.r.x - p1.old_r.x + f_sum.x/p1.m*dt**2
        new_y = 2 * p1.r.y - p1.old_r.y + f_sum.y/p1.m*dt**2
        p1.save()
        p = p1.clone()
        p.old_r = p1.r
        p.r = Vector2(new_x,new_y)
        updated_particles.append(p)
    return updated_particles

T_MAX = 100
t = 0

import matplotlib.pyplot as plt
fig, ax = plt.subplots(1,1)
ax.set_aspect('equal', adjustable='box')

while(t < T_MAX):
    particles = compute(qt)
    qt = QuadTree(QRegion(Vector2(SIZE/2,SIZE/2),SIZE))
    for p in particles:
        qt.insert(p)
    printProgressBar(t,T_MAX)
    ax.cla()
    
    ax.scatter([p.r.x for p in particles], [p.r.y for p in particles],s=0.5)
    qt.drawLines(ax)
    plt.pause(0.1)
    t += dt

#qt.drawForces(ax)
#particles = qt.getAllParticles()
#for i in range(len(particles)):
#   ax.scatter([s[0].x for s in particles[i].states],[s[0].y for s in particles[i].states],s=0.3)

plt.show()