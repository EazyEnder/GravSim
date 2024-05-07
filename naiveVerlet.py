import numpy as np
from Particle import Particle
from Vector2 import Vector2
from math import *
from utils import printProgressBar

SIZE = 100

P_NUMBER = 50
P_MASS = 1

G = 0.25

from random import random


dt = 0.5

initial_particles = []
for i in range(P_NUMBER):
    initial_particles.append(Particle(P_MASS,SIZE*random(),SIZE*random(),0,0,dt=dt))

#Integration de Verlet
a=0. #"amortissement"

def compute(particles,a=0.):
    updated_particles = []
    for p1 in particles:
        f_sum = Vector2(0,0)
        pot_e_sum = 0.
        for p2 in particles:
            if(p2.equals(p1)):
                continue
            f_sum = f_sum.add(p2.r.sub(p1.r).normalize().multiplyScalar(G * (p1.m*p2.m) * 1/(p1.distTo(p2)**2+a**2))) 
            pot_e_sum = 0#pot_e_sum - G * (p1.m*p2.m) * 1/(np.sqrt(p1.distTo(p2)**2+a**2))
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

T_MAX = 200
t = 0

softening = np.linspace(0.3,1.2,5)
result = []
import copy
for a in range(len(softening)):
    particles = copy.deepcopy(initial_particles)
    while(t < T_MAX):
        particles = compute(particles,softening[a])
        printProgressBar(a*T_MAX+t,T_MAX*len(softening))
        t += dt
    t = 0
    result.append(particles)
    
    

import matplotlib.pyplot as plt

"""
fig, (ax, ax2) = plt.subplots(1,2)

ax2.hist([abs(p.old_r.sub(p.r).length())/dt for p in particles], color="black")
ax2.set_xlabel("Vitesse v$")
ax2.set_ylabel("Nombre de Particules")

for i in range(len(particles)):
   ax.scatter([s[0].x for s in particles[i].states],[s[0].y for s in particles[i].states],s=0.3)

ax.set_xlabel("pos $x$")
ax.set_ylabel("pos $y$")"""

from computeEnergy import plotEnergy
fig, (ax, ax2, ax3) = plt.subplots(3,1)
for i in range(len(softening)):
    plotEnergy(result[i], dt, G, softening[i], [fig,ax, ax2, ax3])
plt.show()