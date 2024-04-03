import numpy as np
from Particle import Particle
from Vector2 import Vector2
from math import *
from utils import printProgressBar
from QuadTree.QuadTree import QuadTree
from QuadTree.QRegion import QRegion

SIZE = 100

P_NUMBER = 1000
P_MASS = 1

G = 0.25

from random import random


dt = 0.5

qt = QuadTree(QRegion(Vector2(0,0),SIZE*5))
for i in range(P_NUMBER):
    qt.insert(Particle(P_MASS,SIZE*random(),SIZE*random(),0,0,dt=dt))


#Integration de Verlet
a=0.3 #"amortissement"

def compute(particles):
    updated_particles = []
    for p1 in particles:
        f_sum = Vector2(0,0)
        pot_e_sum = 0.
        for p2 in particles:
            if(p2.equals(p1)):
                continue
            f_sum = f_sum.add(p2.r.sub(p1.r).normalize().multiplyScalar(G * (p1.m*p2.m) * 1/(p1.distTo(p2)**2+a**2))) 
            pot_e_sum = pot_e_sum - G * (p1.m*p2.m) * 1/(np.sqrt(p1.distTo(p2)**2+a**2))
        #Verlet
        new_x = 2 * p1.r.x - p1.old_r.x + f_sum.x/p1.m*dt**2
        new_y = 2 * p1.r.y - p1.old_r.y + f_sum.y/p1.m*dt**2
        p1.save()
        p = p1.clone()
        p.old_r = p1.r
        p.potential_energy = pot_e_sum
        p.r = Vector2(new_x,new_y)
        updated_particles.append(p)
    return updated_particles

T_MAX = 100
t = 0
import copy
while(t < T_MAX):
    particles = compute(particles)
    printProgressBar(t,T_MAX)
    t += dt

import matplotlib.pyplot as plt
plt.figure(figsize=(10,10)) 

for i in range(len(particles)):
   plt.scatter([s[0].x for s in particles[i].states],[s[0].y for s in particles[i].states],s=0.3)

plt.show()