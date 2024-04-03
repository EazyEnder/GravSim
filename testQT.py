import matplotlib.pyplot as plt
import numpy as np
from Vector2 import Vector2
from Particle import Particle
from QuadTree.QuadTree import QuadTree
from QuadTree.QRegion import QRegion
from utils import printProgressBar
from timeit import default_timer as timer

import time

N = 500
positions = np.random.normal(0.5,0.3,(N,2))

qt = QuadTree(QRegion(Vector2(.5,.5),1.))
for p in positions:
    qt.insert(Particle(1.,p[0],p[1],0,0,1))

fig, (ax, ax2, ax3) = plt.subplots(1,3)

#print(qt.getForce(qt.getAllParticles()[0]).length())

ax.scatter(positions[:,0], positions[:,1],s=1.5, color="black")
ax.set_aspect('equal', adjustable='box')
ax.set_title("Forces experienced by each particle")
ax.set_xlabel("pos X")
ax.set_ylabel("pos Y")
ax2.scatter(positions[:,0], positions[:,1],s=1.5, color="black")
ax2.set_aspect('equal', adjustable='box')
ax2.set_title("Center of mass of each QuadTree")
ax2.set_xlabel("pos X")
ax2.set_ylabel("pos Y")

qt.computeMassCenter()
qt.drawLines(ax)
qt.drawForces(ax)
qt.drawLines(ax2)
qt.drawMassCenter(ax2)

x = np.linspace(0,100,100)
times = []
for i in x:
    start = timer()
    for p in qt.getAllParticles():
        qt.getForce(p, i)
    end = timer()
    times.append(end - start)
    printProgressBar(i,np.max(x))
ax3.scatter(x,times)

plt.show()