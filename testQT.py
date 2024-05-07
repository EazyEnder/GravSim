import matplotlib.pyplot as plt
import numpy as np
from Vector2 import Vector2
from Particle import Particle
from QuadTree.QuadTree import QuadTree
from QuadTree.QRegion import QRegion
from utils import printProgressBar
from timeit import default_timer as timer


N = 100
positions = np.random.normal(0.5,0.3,(N,2))

qt = QuadTree(QRegion(Vector2(.5,.5),1.))
for p in positions:
    qt.insert(Particle(1.,p[0],p[1],0,0,1))

fig, (ax, ax2) = plt.subplots(1,2)

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
start = timer()
qt.drawForces(ax)
end = timer()
qt.drawLines(ax2)
qt.drawMassCenter(ax2)
#print(end-start)

x = []
y = []
u = []
v = []
for p in qt.getAllParticles():
    x.append(p.r.x)
    y.append(p.r.y)
    f_sum = Vector2(0,0)
    for p1 in qt.getAllParticles():
        if(p1.equals(p)):
            continue
        f_sum = f_sum.add(p1.r.sub(p.r).normalize().multiplyScalar(0.005 * (p.m*p1.m) * 1/(p.distTo(p1)**2+0.3**2)))
    u.append(f_sum.x)
    v.append(f_sum.y)

ax.quiver(x,y, u, v, color="orange", width=0.0025)

plt.show()