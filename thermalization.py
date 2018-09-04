import math
import numpy
import numpy
import pylab
from scipy.optimize import curve_fit
import math
from scipy import stats
#from scipy.stats import scipy.stats.ks_2samp

#apro file
mag, abs_mag, bc = pylab.loadtxt("thermalizationtrybig.txt", unpack =True)

n = len(bc)
#definisco le delte aperiodiche
bc = -(bc - 1)/2

#definisco percentuale che si butta via per la termalizzazione
t = 5 * (n//100)
#definisco la lunghezza >> tau
l = 30
#definisco numero eventi che voglio testare usando KS, ossia la lunghezza degli array di test
evt = 1000

testmagA = numpy.empty(evt)
testmagB = numpy.empty(evt)

testbcA = numpy.empty(evt)
testbcB = numpy.empty(evt)
#creo gli array su cui faccio ks per la mag e per le bc
for i in range(evt):
	testmagA[i] = mag[t + l*i]
	testmagB[i] = mag[t + l*(evt+i)]
	
	testbcA[i] = bc[t + l*i]
	testbcB[i] = bc[t + l*(evt+i)]

#faccio KS sui due tipi di osservabili
print("MAGNETIZZAZIONE")
print("(statistica, p value)", stats.ks_2samp(testmagA, testmagB))

print("BOUNDARY CONDITION")
print("(statistica, p value)", stats.ks_2samp(testbcA, testbcB))
