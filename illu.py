import matplotlib.pyplot as plt
import numpy as np

X = np.linspace(-1,2,100)
Ysimu = [1 if x > 0 and x < 1 else 0.015 for x in X]
a = 0.1
Yreel =  0.014 / ((X+0.03)**2+a**2)

plt.plot(X,Ysimu,color="black",label="Simulation PP")
plt.plot(X,Yreel,color="orange",label="Réalité")
plt.legend()
plt.xlabel("Temps t")
plt.ylabel("Force exercé sur $p_1$ de $p_2$")
plt.title("Plus proche à t=0, pas de temps = 1")
plt.gca().set_aspect('equal', adjustable='box')

plt.show()