import numpy as np
import matplotlib.pyplot as plt


y = np.loadtxt("y.dat")
px = np.loadtxt("px.dat")
py = np.loadtxt("py.dat")

plt.plot(y,px,label="px")
plt.plot(y,py,label="py")

plt.legend()
plt.show()


