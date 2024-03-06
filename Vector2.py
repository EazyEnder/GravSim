from math import sqrt

class Vector2():
    def __init__(self,x,y):
        self.x = x
        self.y = y
    def dot(self,v2):
        return self.x*v2.x + self.y*v2.y
    def sub(self,v2):
        return Vector2(self.x-v2.x,self.y-v2.y)
    def multiplyScalar(self,scalar):
        return Vector2(self.x*scalar,self.y*scalar)
    def add(self,v2):
        return Vector2(self.x+v2.x,self.y+v2.y)
    def length(self):
        return sqrt(self.dot(self))
    def normalize(self):
        return self.multiplyScalar(1/self.length())
    def equals(self,v2):
        return self.x == v2.x and self.y == v2.y