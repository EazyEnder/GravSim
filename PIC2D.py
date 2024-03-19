import scipy.sparse as sp
import scipy.sparse.linalg as linalg
import numpy as np
from utils import printProgressBar, countParticles
from serializePositions import serializePositions, deserializeArray

#Executable script

G_CST = 0.00025

def buildLaplacian2DKernel(N=10,dx=1):
        i = np.ones(N)
        diags = np.array([-1,0,1])
        vals  = np.vstack((i,-2*i,i))
        M = sp.spdiags(vals, diags, N, N)
        M = sp.lil_matrix(M)
        M[0,N-1] = 1
        M[N-1,0] = 1
        M /= dx**2
        M = sp.csr_matrix(M)
        return M

def buildLaplacianMatrix2D(N=10,dx=1):
    id = np.identity(N)
    kernel = buildLaplacian2DKernel(N,dx)
    M = sp.kron(id,kernel) + sp.kron(kernel,id)
    M = sp.csr_matrix(M)
    return M

def buildGradient2DKernel(N=10,dx=1):
    i = np.ones(N)
    diags = np.array([-1,1])
    vals  = np.vstack((-i,i))
    M = sp.spdiags(vals, diags, N, N)
    M = sp.lil_matrix(M)
    M[0,N-1] = -1
    M[N-1,0] = 1
    M /= 2.*dx
    M = sp.csr_matrix(M)
    return M

def buildGradientMatrix2D(N=10,dx=1):
    id = np.identity(N)
    kernel = buildGradient2DKernel(N,dx)
    M = sp.kron(id,kernel) + sp.kron(kernel,id)
    M = sp.csr_matrix(M)
    return M

from math import floor
def getAcc(particles_pos, Nx, size, rho0, cst_G=G_CST):
    N = particles_pos.shape[0]
    dx = size / Nx

    rho = np.zeros((Nx,Nx))
    for pos in particles_pos:
        rho[floor(pos[0]/dx),floor(pos[1]/dx)] += 1

    rho = rho.flatten()
        
    rho *= rho0*size**2 / N / dx**2 
    phi_grid = linalg.spsolve(buildLaplacianMatrix2D(Nx,dx), 4*np.pi*(rho)*cst_G)

    phi_grid = np.reshape(phi_grid,(Nx,Nx))

    (Gx_grid,Gy_grid) = np.gradient(phi_grid)
    Gx_grid = -Gx_grid
    Gy_grid = -Gy_grid

    x = np.floor(particles_pos[:,0]/dx).astype(int)
    y = np.floor(particles_pos[:,1]/dx).astype(int)
    
    Gx = Gx_grid[x,y]
    Gy = Gy_grid[x,y]

    a = (Gx,Gy)
    
    return a

import random


N = 100000
Nt = 500
Nx = 500
size = 15
r_angle = np.random.rand(N)*2*np.pi
#sqrt to have an uniform density in the disk
r_radius = np.sqrt(np.random.rand(N))*size/6
particles_pos = np.array([[r_radius[i]*np.cos(r_angle[i])+size/2,size/2+r_radius[i]*np.sin(r_angle[i])] for i in range(N)])


BASE_VELOCITY = 0.020

p_vel_x = -np.sin(r_angle)*r_radius*BASE_VELOCITY
p_vel_y = np.cos(r_angle)*r_radius*BASE_VELOCITY
dt = 1

def loadFile():
    f = open("export/sim1.txt","r")
    lines = f.readlines()
    l_last_line = lines[2*Nt-2]
    last_line = lines[2*Nt-1]
    array_pos = deserializeArray(last_line)
    array_velocity = (array_pos-deserializeArray(l_last_line))/dt
    f.close()

    global particles_pos
    global p_vel_x
    global p_vel_y
    particles_pos = array_pos
    p_vel_x = array_velocity[:,0]
    p_vel_y = array_velocity[:,1]

#loadFile()


import matplotlib.pyplot as plt

f = open("export/sim1.txt", "a")

def plotGraph(p_pos,size):
     plt.clf()
     rho = np.zeros((480,480))
     for pos in p_pos:
          rho[floor(480.*pos[0]/size),floor(480.*pos[1]/size)] += 1
     plt.imshow(np.power(rho,0.33), cmap="afmhot")
     plt.colorbar()

def launchSim(dt,Nt,Nx,size,p_pos,p_vel_x,p_vel_y):
    (p_ac_x,p_ac_y) = getAcc( p_pos, Nx, size, 1.)
    t = 0
    for i in range(Nt):
         printProgressBar(i,Nt)
         p_vel_x += p_ac_x * dt/2.
         p_vel_y += p_ac_y * dt/2.
         p_pos[:,0] += p_vel_x * dt
         p_pos[:,1] += p_vel_y * dt
         p_pos = np.mod(p_pos, size)
         (p_ac_x,p_ac_y) = getAcc( p_pos, Nx, size, 1, G_CST)
         p_vel_x += p_ac_x * dt/2.0
         p_vel_y += p_ac_y * dt/2.0
         plotGraph(p_pos,size)
         f.write(serializePositions(p_pos)+"\n")
         plt.savefig(fname="export/fig1_2D"+str(i)+"")
         #plt.pause(0.01)
         t += dt

launchSim(dt,Nt,Nx,size,particles_pos,p_vel_x,p_vel_y)

f.close()
print("\nDone")