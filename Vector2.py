from math import sqrt

class Vector2():
    """2D Vector
        >Args: x and y coordinates"""
    def __init__(self,x,y):
        self.x = x
        self.y = y
    def dot(self,v2):
        """dot product with an other vector v2"""
        return self.x*v2.x + self.y*v2.y
    def sub(self,v2):
        """sub this vector with an other vector v2"""
        return Vector2(self.x-v2.x,self.y-v2.y)
    def multiplyScalar(self,scalar):
        """multiply this vector to a scalar"""
        return Vector2(self.x*scalar,self.y*scalar)
    def add(self,v2):
        """add this vector to an other vector v2"""
        return Vector2(self.x+v2.x,self.y+v2.y)
    def length(self):
        """Return the length of this vector"""
        return sqrt(self.dot(self))
    def normalize(self):
        """Normalize this vector"""
        if(self.length() != 0.):
            return self.multiplyScalar(1/self.length())
        else:
            return Vector2(0.,0.)
    def equals(self,v2):
        """Check if this vector is equal to another"""
        return self.x == v2.x and self.y == v2.y