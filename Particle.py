import copy
from Vector2 import Vector2
from math import sqrt

#Objet: Particle
class Particle():
    def __init__(self,m,x,y,vx,vy,dt=1.):
        self.m = m
        self.r = Vector2(x,y)
        self.old_r = Vector2(x-vx*dt,y-vy*dt)
        self.states = []
        self.potential_energy = 0
    def distTo(self,p2):
        diff = self.r.sub(p2.r)
        return sqrt(diff.dot(diff))
    def equals(self,p2):
        return p2.r.equals(self.r) and p2.old_r.equals(self.old_r) and p2.m == self.m
    def clone(self):
        return copy.deepcopy(self)
    def save(self):
        self.states.append((self.old_r,self.potential_energy))
        