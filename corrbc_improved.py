#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import math
import numpy
import numpy
import pylab
import math
import scipy.stats
import matplotlib.pyplot as plt

mag, abs_mag, bc = pylab.loadtxt("0.2391_4_4_12.txt", unpack =True)

i = 1
varb = []
n = []

#riscrivo le bc in termini della delta aperiodica
bc = (-1*bc + 1)/2

totsum = bc.sum()

#blocking 
while i < numpy.sqrt(len(bc)): #scorro sui posssibli divisori del numero di eventi
	#creo array i cui elementi sono le medie dei blocchetti desiderati
	n.append(i) #salvo la dimensione dei blocchi
	block = []
	avg = totsum
	for cnt in range(len(bc) - (len(bc)//i) * i):
		avg -= bc[(len(bc)//i) * i + cnt] #tolgo gli elementi in eccesso, che comunque non vengono passati dal for successivo...
	avg = avg / (i* (len(bc)//i) ) #true avg	
	for j in range(len(bc)//i): #scelgo il blocco
		temp_sum = 0
		for k in range(i): #mi muovo nel blocco e metto le somme ina una variabile temporanea
			temp_sum +=  bc[j*i+k] 
		temp_sum = temp_sum/i #medio la somma
		block.append(temp_sum)			
	#calcolo varianza nelle nuove variabili 
	varb.append(((numpy.array(block)-avg)**2).sum()/(len(block)-1))
	i += 1

#print utili per vedere se funziona

print ("varb = ", varb)
print ("n = ", n)

#creo variabili utili per computo tau

var = varb[0] 

print("var_single = ", var)

#creo la 2*tau con la formula 2tau=sigma_b^2/sigma_O_i^2 * len(block)

varb = numpy.array(varb)
#n = numpy.array(n)

tau = 0.5 * varb * n / var


print ("tau = ", tau)

pylab.title('tau')
pylab.xlabel('block_length')
pylab.ylabel('tau')
pylab.grid(color = "gray")
plt.plot(n, tau)
plt.show()
