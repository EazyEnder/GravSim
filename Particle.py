import copy
from Vector2 import Vector2
from math import sqrt

class Particle():
    """2D Particle object with a mass, position and velocity
        >Args: -m: mass
            -x/y: x/y coordinate\n
            -vx/vy: velocity on the x/y axis\n
            -dt: temp step of the simulation (used to calculate the velocity of the particle)"""
    def __init__(self,m,x,y,vx,vy,dt=1.):
        self.m = m
        self.r = Vector2(x,y)
        self.old_r = Vector2(x-vx*dt,y-vy*dt)
        self.states = []
        self.potential_energy = 0
        self.host = None
    def distTo(self,p2):
        """Calculate the distance between this particle and another named p2"""
        diff = self.r.sub(p2.r)
        return sqrt(diff.dot(diff))
    def equals(self,p2):
        """Check if this particle is equal to another"""
        return p2.r.equals(self.r) and p2.old_r.equals(self.old_r) and p2.m == self.m
    def clone(self):
        """Clone this particle"""
        clone = Particle(self.m,self.r.x,self.r.y,self.r.x,self.r.y)
        clone.old_r = self.old_r
        clone.states = self.states
        clone.potential_energy = self.potential_energy
        clone.host = self.host
        return clone
    def save(self):
        """Save this particle in the 'states' list"""
        self.states.append((self.old_r,self.potential_energy))
        