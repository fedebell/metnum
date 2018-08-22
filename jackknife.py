#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import math
import numpy
import numpy
import pylab
import math
import scipy.stats
import matplotlib.pyplot as plt

mag, abs_mag, bc = pylab.loadtxt("0.5_10_10_10.txt", unpack = True)

varb = []
N = 500

#definisco la lunghezza dei blocchi su cui fare jackknife
#lng = numpy.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 400, 500])
lng = numpy.array([i for i in range(N)])
print(lng)


#riscrivo le bc in termini della delta aperiodica
bc = (-1*bc + 1)/2
totsum = bc.sum()

#ciclo sugli elementi dell'array delle lunghezze
for i in range(len(lng)):
	block = [] #resetto la variabile block
	#calcolo media (escludo gli ultimi elementi se lng non divide len(bc)
	avg = totsum
	for cnt in range(len(bc) - (len(bc)//lng[i]) * lng[i]):
		avg -= bc[(len(bc)//lng[i]) * lng[i] + cnt] #tolgo gli elementi in eccesso ,che comunque non vengono passati dal for successivo...
	avg = avg / (lng[i] * (len(bc)//lng[i]) ) #true avg
	#print("avg = ", avg)

	#creo array il cui elemento j-esimo è la media sui blocchi lunghi lng tranne il j-esimo blocco	
	Nb = len(bc)//lng[i] #Numero di blocchi
	for j in range(Nb): #scelgo il blocco DA ESCLUDERE 
		temp_sum = 0
		for k in range(lng[i]): #mi muovo nel blocco e metto le somme ina una variabile temporanea
			temp_sum +=  bc[j*lng[i]+k] 
		avg_exc = (totsum - temp_sum)/((len(bc) // lng[i]) * lng[i] - lng[i]) #media su tutti i  gli elementi tranne quelli appartenenti al j-esimo blocco
		block.append(avg_exc)

	block = numpy.array(block)
	
	print(i)
	#print(block)
	
				
	#calcolo varianza nelle nuove variabili 
	varb.append(((block-block.sum()/len(block))**2).sum()*((Nb-1.0)/Nb)) #len block è pari al numero di blocchi: escluderne uno conta come contarne uno
	
	#print(Nb) 
	
	#print(((block-block.sum()/len(block))**2).sum() * ((Nb-1.0)/Nb) )

varb = numpy.array(varb)
#print utili per vedere se funziona
print ("varb = ", varb)
print("error = ", numpy.sqrt(varb))

somma = 0.0
for i in range(N-100):
	somma += numpy.sqrt(varb)[i+100]
	
somma = somma/(N-100)

print(somma)

#plotto l'errore sulla <delta_aper> in funzione della blok length
pylab.title('errore su <delta_aper>')
pylab.xlabel('block_length')
pylab.ylabel('error')
pylab.grid(color = "gray")
plt.plot(lng, numpy.sqrt(varb))
plt.show()

