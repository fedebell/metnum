import math
import numpy
import numpy
import pylab
import math
import scipy.stats
import matplotlib.pyplot as plt

all_mag, all_abs_mag, all_bc = pylab.loadtxt("miscData//0.223_5_5_15.txt", unpack =True)

if len(all_mag) % 2 == 1:
	print ("error!")


mag = numpy.empty(len(all_mag)//2)
abs_mag = numpy.empty(len(all_abs_mag)//2)
bc = numpy.empty(len(all_bc)//2)

#creo array solo con cose pari (supersloppy)
for i in range(len(all_mag)):
	if i % 2 == 0:
		mag[(i)//2] = all_mag[i]
		abs_mag[i//2] = all_abs_mag[i]
		bc[i//2] = all_bc[i]


totall_mag = all_mag.sum()
totmag = mag.sum()

totall_abs_mag = all_abs_mag.sum()
totabs_mag = abs_mag.sum()

totall_bc = all_bc.sum()
totbc = bc.sum()

print("AVG MAG = ", totall_mag/len(all_mag))
print("AVG MAG REDUCED = ", totmag/len(mag))

print("AVG ABSMAG = ", totall_abs_mag/len(all_abs_mag))
print("AVG ABSMAG REDUCED = ", totabs_mag/len(abs_mag))

print("AVG BC = ", totall_bc/len(all_bc))
print("AVG BC REDUCED = ", totbc/len(bc))
