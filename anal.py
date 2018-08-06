import numpy
import pylab
#from scipy.optimize import curve_fit
import math
#from uncertainties import ufloat, unumpy

beta, per, aper = pylab.loadtxt("data.txt", unpack = True)

ratio = aper/per

free = - numpy.log(ratio)

print(beta)
print(free)

pylab.errorbar(beta, free, marker = 'o', color = 'blue', linestyle = '')
pylab.title('Energia libera della superficie in funzione della temperatura.')
pylab.ylabel('$F$', fontsize = 16)
pylab.xlabel('$\\beta$', fontsize = 16)
pylab.grid()
pylab.show()


