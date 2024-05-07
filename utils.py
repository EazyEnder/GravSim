import numpy as np
from math import floor

def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """Print a progress bar"""
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    if iteration == total: 
        print()

def isInBox(pos,center,h_length):
    return (pos.x <= center.x + h_length and pos.x >= center.x - h_length) and  (pos.y <= center.y + h_length and pos.y >= center.y - h_length)

def isInBox2(pos,min,max):
    return (pos.x <= max.x and pos.x >= min.x) and  (pos.y <= max.y and pos.y >= min.y)

def countParticles(p_pos,total_size=15,frame_size=480,periodic=True):
    """Count how many particles has a cell in a grid without interpolations, i.e a particle is only in one cell
        >Args: -p_pos: particles positions
            -size: size of the simulation box\n
            -frame_size: size of the grid
        >Return: 2D tensor of particle number in each cell of the grid"""
    rho = np.zeros((frame_size,frame_size))
    for pos in p_pos:
        i0 = floor(frame_size*pos[0]/total_size)
        i1 = floor(frame_size*pos[1]/total_size)
        if(periodic):
            if i1 >= frame_size:
                i1 = frame_size -1
            if i0 >= frame_size:
                i0 = frame_size -1
        try:
            rho[i0,i1] += 1
        except:
            "particle outbound"
    return rho

def countParticlesUsingBilinearInterp(p_pos,total_size=15,frame_size=480):
    """Not completed"""
    weights = []
    for pos in p_pos:
        p0 = frame_size*pos[0]/total_size
        p1 = frame_size*pos[1]/total_size

        i0 = floor(p0)
        i1 = floor(p1)
        
        center = np.array([i0+.5,i1+.5])
        dist_to_center = pos-center
        points = []
        morph = np.zeros((2,2))
        if(dist_to_center[0] > 0 and dist_to_center[1] > 0):
            points = [
                center,
                [center[0],center[1]+1],
                [center[0]+1,center[1]],
                [center[0]+1,center[1]+1]
            ]
            morph[1][1] = 1
            #"bas droit"
        elif(dist_to_center[0] > 0 and dist_to_center[1] < 0):
            points = [
                [center[0]-1,center[1]],
                [center[0]-1,center[1]+1],
                center,
                [center[0],center[1]+1]
            ]
            morph[0][1] = 1
            #"haut droit"
        elif(dist_to_center[0] < 0 and dist_to_center[1] > 0):
            points = [
                [center[0],center[1]-1],
                center,
                [center[0]+1,center[1]-1],
                [center[0]+1,center[1]]
            ]
            #"bas gauche"
            morph[1][0] = 1
        elif(dist_to_center[0] < 0 and dist_to_center[1] < 0):
            points = [
                [center[0]-1,center[1]-1],
                [center[0]-1,center[1]],
                [center[0],center[1]-1],
                center
            ]
            #"haut gauche"
            morph[0][0] = 1
        x = p0 - points[0][0]
        y = p1 - points[0][1]
        w0 = np.array([
            [(1-x)*(1-y),((1-x)*y)],
            [((1-y)*x),x*y]
        ])
        w = np.kron(morph,w0)



def findMassCenter(rho,center=(0.,0.),distance=10.):
    """Find the inertia center of a circle localized area
        >Args: -rho: 2D tensor of density
            -center: center of the localized area\n
            -distance: distance from 0. to 1. of how far we take particles/ radius of the circle area
        >Return: coordinates of the inertia center """
    norm_matrix = []
    for i in range(len(rho)):
        norm_matrix.append([])
        for j in range(len(rho[i])):
            r = rho[i,j]
            x = i / len(rho)
            y = j / len(rho[i])
            if(np.sqrt((x-center[0])**2 + (y-center[1])**2) <= distance):
                norm_matrix[i].append(r)
            else:
                norm_matrix[i].append(0.)
    norm_matrix = np.array(norm_matrix)
    xg = 0.
    yg = 0.
    for i in range(len(norm_matrix)):
        for j in range(len(norm_matrix[i])):
            x = i / len(norm_matrix)
            y = j / len(norm_matrix[i])
            xg += norm_matrix[i,j]*x
            yg += norm_matrix[i,j]*y
    m_sum = np.sum(norm_matrix)
    return (xg/m_sum,yg/m_sum)