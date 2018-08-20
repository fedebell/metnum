import math
import numpy
import numpy
import pylab
from scipy.optimize import curve_fit
import math
import scipy.stats
import matplotlib.pyplot as plt


mag, abs_mag, bc = pylab.loadtxt("C:\cygwin64\metnum\pureising.txt", unpack =True)
n = range (len(mag))
mag_per = []

for i in range(len(mag)):
	if bc[i]==1:
		mag_per.append(mag[i])
		
n_per = range(len(mag_per))

pylab.title('magnetization in markov chain')
pylab.xlabel('N')
pylab.ylabel('m')
pylab.grid(color = "gray")
plt.plot(n, mag)
plt.show()

pylab.title("boundary condition")
pylab.xlabel("N")
pylab.ylabel("bc")
pylab.grid(color = "gray")
plt.plot(n, bc)
plt.show()

pylab.title("magnetization in periodic bc only")
pylab.xlabel("N_per")
pylab.ylabel("m")
pylab.grid(color = "gray")
plt.plot(n_per, mag_per)
plt.show()
