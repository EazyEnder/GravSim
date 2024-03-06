import scipy.sparse as sp
import scipy.sparse.linalg as linalg
import numpy as np

G_CST = 0.001

def buildLaplacianMatrix1D(N=10,dx=1):
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

def buildGradientMatrix1D(N=10,dx=1):
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


def leapfrog(particles_pos, Nx, size, rho0):
	N = particles_pos.shape[0]
	dx = size / Nx
	j = np.floor(particles_pos/dx).astype(int)
	w_j   = ( (j+1)*dx - particles_pos  )/dx
	w_j1 = ( particles_pos - j*dx )/dx
	j1 = np.mod(j+1, Nx)
	rho  = np.bincount(j[:,0],   weights=w_j[:,0],   minlength=Nx);
	rho += np.bincount(j1[:,0], weights=w_j1[:,0], minlength=Nx);
	rho *= rho0*size / N / dx 
	
	phi_grid = linalg.spsolve(buildLaplacianMatrix1D(Nx,dx), 4*np.pi*(rho)*G_CST)
	
	G_grid = - buildGradientMatrix1D(Nx,dx) @ phi_grid
	
	G = w_j * G_grid[j] + w_j1 * G_grid[j1]
	
	a = G

	return a

from math import floor

def plotGraph(p_pos,p_vel,Nx=1000):
    plt.clf()
    plt.subplot(2,1,1)
    plt.scatter(p_pos, p_vel, s=.25,color='black')
    plt.axis([0,size,-6,6])
    plt.title("Portrait de phase")
    plt.xlabel("Position")
    plt.ylabel("Vitesse")

    plt.subplot(2,1,2)
    rho  = np.bincount(np.floor(p_pos*Nx/size).astype(int)[:,0],   minlength=Nx);
    plt.imshow([rho for i in range(100)])
    plt.title("Densité")
    plt.xlabel("Position")
    plt.colorbar()

     


N = 10000
Nt = 1000
Nx = 1000
size = 100
particles_pos = np.random.rand(N,1) * size
particles_vel = np.random.randn(N,1) * 0.
particles_vel[floor(N/2):] *= -1
dt = 1

import matplotlib.pyplot as plt

fig = plt.figure(figsize=(5,5), dpi=100)

def launchSim(dt,Nt,Nx,size,p_pos,p_vel):
    p_ac = leapfrog( p_pos, Nx, size, 1.)
    t = 0
    for i in range(Nt):
         p_vel += p_ac * dt/2.0
         p_pos += particles_vel * dt
         p_pos = np.mod(p_pos, size)
         p_ac = leapfrog( p_pos, Nx, size, 1)
         p_vel += p_ac * dt/2.0
         plt.cla()
         plotGraph(p_pos,p_vel)
         plt.pause(0.01)
         t += dt

launchSim(dt,Nt,Nx,size,particles_pos,particles_vel)
