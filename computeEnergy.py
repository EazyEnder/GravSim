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

def getTotalPotentialEnergy(particles,dt):
    p_E = computePotentialEnergy(particles,dt)
    return np.sum(p_E,axis=0)

import matplotlib.pyplot as plt
def plotEnergy(particles,dt):
    Pot_E = -getTotalPotentialEnergy(particles,dt)
    Kin_E = getTotalKineticEnergy(particles,dt)
    plt.plot(Kin_E, label="$E_c$")
    plt.plot(Pot_E,label="$E_p$")
    plt.plot(Kin_E+Pot_E, label="$E_{tot}$")
    plt.legend()
    plt.xlabel("Time")
    plt.title("Energy vs time")