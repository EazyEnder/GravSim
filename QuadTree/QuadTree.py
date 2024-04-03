
QT_PARTICLE_CAPACITY = 10
"""How many particles a qt can contain"""

G = 0.25

from Vector2 import Vector2
from utils import isInBox2
from QuadTree.QRegion import QRegion

import matplotlib.pyplot as plt
import matplotlib.patches as patches

class QuadTree():
    def __init__(self, region):
        self.region = region
        #NE-NW-SE-SW
        self.subqt = [None,None,None,None]
        #tuple with qt parent & child id
        self.parent = None
        self.particles = []
        #masscenter : [mass, Vector2(x,y)]
        self.masscenter = None

    def insert(self,particle):
        """Insert a particle in the quadtree
            >Args: particle
            >Return: if it's a success
        """
        if(not(self.region.containsParticle(particle))):
            return False
        if(len(self.particles) < QT_PARTICLE_CAPACITY and self.subqt[0] == None):
            particle.host = self
            self.particles.append(particle)
            return True
        
        if(self.subqt[0] == None):
            self.subdivide()

        for qt in self.subqt:
            if(qt.insert(particle)):
                return True
        
        return False
    
    def subdivide(self):
        """Subdivide the quadtree into 4 smaller trees"""

        #NE
        self.subqt[1] = QuadTree(QRegion(self.region.center.add(Vector2(-self.region.half_length/2,+self.region.half_length/2)),self.region.half_length/2))
        #NW
        self.subqt[0] = QuadTree(QRegion(self.region.center.add(Vector2(-self.region.half_length/2,-self.region.half_length/2)),self.region.half_length/2))
        #SE
        self.subqt[3] = QuadTree(QRegion(self.region.center.add(Vector2(+self.region.half_length/2,+self.region.half_length/2)),self.region.half_length/2))
        #SW
        self.subqt[2] = QuadTree(QRegion(self.region.center.add(Vector2(+self.region.half_length/2,-self.region.half_length/2)),self.region.half_length/2))
        
        for i in range(len(self.subqt)):
            self.subqt[i].parent = (self,i)

        #Copy the list bcs always a bad idea to modify when we iterate it
        qt_particles = []
        qt_particles.extend(self.particles)
        for p in qt_particles:
            for qt in self.subqt:
                if(qt.region.containsParticle(p)):
                    qt.insert(p)
                    self.particles.remove(p)
                    break
        if(len(self.particles) > 0):
            print("QT node not in subtree")

    def query(self, min, max):
        """Get all the particles in a demarcated area
            >Args: -min: Vector2
                -max: Vector2
            >Return: list containing all the particles in the area
        """

        result = []

        if(self.region.center.x + self.region.half_length < min.x or self.region.center.x - self.region.half_length > max.x or self.region.center.y + self.region.half_length < min.y or self.region.center.y - self.region.half_length > max.y):
            return result
        
        for p in self.particles:
            if(isInBox2(p.r, min, max)):
                result.append(p)

        if(self.subqt[0] == None):
            return result
        
        for qt in self.subqt:
            result.extend(qt.query(min,max))

        return result
    
    def computeMassCenter(self):
        """Compute the inertial center of all the particles in the QuadTree and his children"""
        #mass, vector 2 position
        masscenter = [0,Vector2(0,0)]
        if(self.subqt[0] != None):
            for qt in self.subqt:
                m_c = qt.masscenter
                if(m_c == None):
                    m_c = qt.computeMassCenter()
                masscenter[0] += m_c[0]
                masscenter[1] = masscenter[1].add(m_c[1].multiplyScalar(m_c[0]))
            masscenter[1] = masscenter[1].multiplyScalar(1/masscenter[0])
        
        if(len(self.particles) > 0 and self.subqt[0] == None):
            p_masscenter = [0, Vector2(0,0)]
            for p in self.particles:
                p_masscenter[0] += p.m
                p_masscenter[1] = p_masscenter[1].add(p.r.multiplyScalar(p.m))
            p_masscenter[1] = p_masscenter[1].multiplyScalar(1/p_masscenter[0])

            masscenter[0] = masscenter[0] + p_masscenter[0]
            masscenter[1] = masscenter[1].add(p_masscenter[1].multiplyScalar(p_masscenter[0])).multiplyScalar(1/masscenter[0])

        self.masscenter = masscenter
        return masscenter
        
    def drawMassCenter(self,ax,layer=1):
        """Draw inertial center for each quadtree"""
        if(self.masscenter != None and self.masscenter[0] > 0):
            colorlist = ["green","lime","gold","darkorange","indianred","crimson","darkorchid","slateblue","blue"]
            ax.scatter(self.masscenter[1].x,self.masscenter[1].y,c=colorlist[layer-1 % len(colorlist)],s=50/layer)

        if(self.subqt[0] == None):
            return

        for qt in self.subqt:
            qt.drawMassCenter(ax,layer=layer+1)
    
    def drawLines(self,ax):
        """Draw the border of the Quadtree and his children"""
        square = patches.Rectangle((self.region.center.x-self.region.half_length, self.region.center.y-self.region.half_length), self.region.half_length*2, self.region.half_length*2, edgecolor='red', facecolor='none', linewidth=0.5)
        ax.add_patch(square)
        
        if(self.subqt[0] == None):
            return

        for qt in self.subqt:
            qt.drawLines(ax)

    def drawForces(self,ax):
        """"""

        x = []
        y = []
        u = []
        v = []
        for p in self.particles:
            x.append(p.r.x)
            y.append(p.r.y)
            force = self.getForce(p,precision=1)
            u.append(force.x)
            v.append(force.y)

        ax.quiver(x,y, u, v, color="green", width=0.0025)
        
        if(self.subqt[0] == None):
            return

        for qt in self.subqt:
            qt.drawForces(ax)

    def getAllParticles(self):
        """Return a flatten list of all the particles in the QuadTree and his children"""
        part = []
        part.extend(self.particles)

        if(self.subqt[0] == None):
            return part

        for qt in self.subqt:
            part.extend(qt.getAllParticles())

        return part
    
    def getForce(self, particle, precision = 1):
        
        if(self.masscenter == None):
            self.computeMassCenter()

        force = Vector2(0,0)

        layer = self
        parent_index = -1
        while(layer.parent != None):
            for i in range(len(layer.parent[0].subqt)):
                if(i == parent_index):
                    continue
                force = force.add(layer.parent[0].subqt[i].getSubForce(particle, precision))

            layer,parent_index = layer.parent

        return force
    
    def getMinSize(self):
        min_size = self.region.half_length*2.
        if(self.subqt[0] != None):
            for qt in self.subqt:
                min_size = min(min_size,qt.getMinSize())
        return min_size
    
    def getClosestDistance(self, pos):
        c_distance = pos.sub(self.region.center).length()
        if(self.subqt[0] != None):
            for qt in self.subqt:
                c_distance = min(c_distance,qt.getClosestDistance(pos))
        return c_distance    
    
    def getSubForce(self, particle, precision=1):

        distToQt = self.getClosestDistance(particle.r)
        sizeQt = self.getMinSize()

        if(distToQt/sizeQt < precision):
            f_sum = Vector2(0,0)

            if(self.subqt[0] != None):
                f_sum = Vector2(0,0)
                for qt in self.subqt:
                    f_sum = f_sum.add(qt.getSubForce(particle, precision))
                return f_sum
            
            for p in self.particles:
                f_sum = f_sum.add(p.r.sub(particle.r).normalize().multiplyScalar(G * (p.m*particle.m) * 1/(p.distTo(particle)**2+0.3**2)))
            return f_sum
        else:
            if(self.masscenter == None):
                self.computeMassCenter()
            return Vector2(0,0).add(self.masscenter[1].sub(particle.r).normalize().multiplyScalar(G * (self.masscenter[0]*particle.m) * 1/(self.masscenter[1].sub(particle.r).length()**2+0.3**2)))
