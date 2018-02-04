import numpy as np
import matplotlib.pyplot as plt


y = np.loadtxt("y.dat")
p = np.loadtxt("p.dat")

plt.plot(y,p)
plt.show()


