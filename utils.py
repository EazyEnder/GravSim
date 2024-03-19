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

def countParticles(p_pos,total_size=15,frame_size=480):
    """Count how many particles has a cell in a grid without interpolations, i.e a particle is only in one cell
        >Args: -p_pos: particles positions
            -size: size of the simulation box\n
            -frame_size: size of the grid
        >Return: 2D tensor of particle number in each cell of the grid"""
    rho = np.zeros((frame_size,frame_size))
    for pos in p_pos:
        i0 = floor(frame_size*pos[0]/total_size)
        i1 = floor(frame_size*pos[1]/total_size)
        if i1 >= frame_size:
            i1 = frame_size -1
        if i0 >= frame_size:
            i0 = frame_size -1
        rho[i0,i1] += 1
    return rho

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