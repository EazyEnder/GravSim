import matplotlib.pyplot as plt
import numpy as np
from serializePositions import deserializeArray
from utils import countParticles, findMassCenter, printProgressBar

def plotGraph(p_pos,size):
     plt.clf()
     rho = countParticles(p_pos,total_size=size)
     plt.imshow(np.power(rho,0.5), cmap="afmhot")
     plt.colorbar()

def centerFilm(sim_file_name="sim.txt",export_file_prefix="modified_export",size=15,exceptions=[]):
     """Plot density graph and keep to the center the 'main body' using an approximation of the inertia center
            >Args: -sim_file_name: simulation  file name   
                -size: size of the simulation box\n
                -exceptions: skip lines in the sim file; [ interval1, interval2, ...] with an interval : [x0,x1]
            >Return: Save plots as a list of png images"""
     f = open("export/"+sim_file_name,"r")
     lines = f.readlines()
     previous_center = (np.mod(0.5,1),np.mod(0.5,1))
     k = 0
     for i in range(len(lines)):
          l = lines[i]
          pos = deserializeArray(l)
          previous_center = findMassCenter(countParticles(pos),previous_center,0.2)
          pos[:,0] -= (previous_center[0]-.5)*size
          pos[:,1] -= (previous_center[1]-.5)*size
          pos = np.mod(pos, size)
          plotGraph(pos,size)
          printProgressBar(i,len(lines))

          flag = False
          for e in exceptions:
            if(i >= e[0] and i <= e[1] ):
                continue
          if flag:
              continue
          else:
               k += 1
          plt.savefig(fname="export/"+export_file_prefix+str(k)+"")
     f.close()

centerFilm("sim1.txt",exceptions=[[1000,1017]])