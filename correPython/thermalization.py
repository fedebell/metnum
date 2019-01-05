import math
import numpy
import numpy
import pylab
from scipy.optimize import curve_fit
import math
from scipy import stats
#from scipy.stats import scipy.stats.ks_2samp

#apro file
mag, abs_mag, bc = pylab.loadtxt("miscData//thermalizationtrybig.txt", unpack =True)

n = len(bc)
#definisco le delte aperiodiche
bc = -(bc - 1)/2

#definisco percentuale che si butta via per la termalizzazione
t = int(round(0 * (n//100)))
#definisco la lunghezza >> tau
l = 10
#definisco numero eventi che voglio testare usando KS, ossia la lunghezza degli array di test
evt = 5000

testmagA = numpy.empty(evt)
testmagB = numpy.empty(evt)

testbcA = numpy.empty(evt)
testbcB = numpy.empty(evt)

testmixA = numpy.empty(evt)
testmixB = numpy.empty(evt)

#creo gli array su cui faccio ks per la mag e per le bc
for i in range(evt):
	testmagA[i] = abs_mag[t + l*i]
	testmagB[i] = abs_mag[t + l*(evt+i)]
	
	testbcA[i] = bc[t + l*i]
	testbcB[i] = bc[t + l*(evt+i)]
	
	testmixA[i] = abs_mag[t + l*i] * (-2*bc[t + l*i] + 1)
	testmixB[i] = abs_mag[t + l*(evt+i)] * (-2*bc[t + l*(evt+i)] + 1)
	
#faccio KS sui due tipi di osservabili
print("MAGNETIZZAZIONE")
print("(statistica, p value)", stats.ks_2samp(testmagA, testmagB))

print("BOUNDARY CONDITION")
print("(statistica, p value)", stats.ks_2samp(testbcA, testbcB))

print("MAGNETIZZAZIONE * BC")
print("(statistica, p value)", stats.ks_2samp(testmixA, testmixB))
