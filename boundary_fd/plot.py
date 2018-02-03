import numpy as np
import matplotlib.pyplot as plt


x = np.loadtxt("x.dat")
t = np.loadtxt("T.dat")

plt.plot(x,t)
plt.show()


