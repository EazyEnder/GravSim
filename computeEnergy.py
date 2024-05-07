import numpy as np

def computeKineticEnergy(particles,dt):
    #matrix : time - particle
    k_E = []
    for p in particles:
        p_velocity = np.sqrt(np.power(np.gradient([s[0].x for s in p.states]),2.)+np.power(np.gradient([s[0].y for s in p.states]),2.))/dt
        p_cinetic_E = .5 * p.m * np.power(p_velocity,2.)
        k_E.append(p_cinetic_E)
    k_E = np.array(k_E)
    return k_E

def getTotalKineticEnergy(particles,dt):
    k_E = computeKineticEnergy(particles,dt)
    return np.sum(k_E,axis=0)

def computePotentialEnergy(particles,dt):
    p_E = []
    for p in particles:
        p_potential_E = np.array([s[1] for s in p.states])
        p_E.append(p_potential_E)
    return p_E

def getTotalPotentialEnergy(particles,dt,G,a=0):
    #p_E = computePotentialEnergy(particles,dt)
    #return np.sum(p_E,axis=0)

    E = []
    for s in range(len(particles[0].states)):
        E_sum = 0
        for i in range(len(particles)):
            for j in range(i+1,len(particles)):
                E_sum -= G * particles[i].m * particles[j].m / np.sqrt(particles[i].states[s][0].sub(particles[j].states[s][0]).length()**2 + a**2)
        E.append(E_sum)

    return np.array(E)

import matplotlib.pyplot as plt
def plotEnergy(particles,dt,G,a=0,axs=None):
    Pot_E = getTotalPotentialEnergy(particles,dt,G,a)
    Kin_E = getTotalKineticEnergy(particles,dt)
    if(axs == None):
        plt.plot(Kin_E, label="$E_c$")
        plt.plot(Pot_E,label="$E_p$")
        plt.plot(Kin_E+Pot_E, label="$E_{tot}$")
        plt.legend()
        plt.xlabel("Time")
        plt.title("Energy vs time")
    else:
        fig = axs[0]
        fig.suptitle(f"Energy vs time for different softening factors $a$, $N_{'{par}'}={len(particles)}$")
        ax1 = axs[1]
        ax1.set_title("Kinetic Energy")
        ax2 = axs[2]
        ax2.set_title("Total Energy")
        ax3 = axs[3]
        ax3.set_title("Potential Energy")
        ax1.plot(Kin_E, label=f"$a={np.round(a,2)}$")
        ax3.plot(Pot_E,label=f"$a={np.round(a,2)}$")
        ax2.plot(Kin_E+Pot_E, label=f"$a={np.round(a,2)}$")
        ax1.legend()
        ax2.legend()
        ax3.legend()
        plt.xlabel("Time")
