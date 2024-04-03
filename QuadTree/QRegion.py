from utils import isInBox

class QRegion():

    def __init__(self, center, half_length):
        self.center = center
        self.half_length = half_length

    def containsParticle(self,particle):
        return isInBox(particle.r, self.center, self.half_length)