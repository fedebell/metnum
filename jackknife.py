import math
import numpy
import numpy
import pylab
import math
import scipy.stats
import matplotlib.pyplot as plt

mag, abs_mag, bc = pylab.loadtxt("C://cygwin64//metnum//ciao.txt", unpack =True)

varb = []

#definisco la lunghezza dei blocchi su cui fare jackknife
lng = numpy.array([400, 500, 600, 700, 800, 900, 1000])

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
	print("avg = ", avg)

	#creo array il cui elemento j-esimo è la media sui blocchi lunghi lng tranne il j-esimo blocco	
	for j in range(len(bc)//lng[i]): #scelgo il blocco DA ESCLUDERE 
		temp_sum = 0
		for k in range(lng[i]): #mi muovo nel blocco e metto le somme ina una variabile temporanea
			temp_sum +=  bc[j*lng[i]+k] 
		avg_exc = (totsum - temp_sum)/((len(bc) // lng[i]) * lng[i] - lng[i]) #media su tutti i  gli elementi tranne quelli appartenenti al j-esimo blocco
		block.append(avg_exc)	

	block = numpy.array(block)
				
	#calcolo varianza nelle nuove variabili 
	varb.append(((numpy.array(block)-avg)**2).sum()/len(block) * (len(block)-1)) #len block è pari al numero di blocchi: escluderne uno conta come contarne uno

varb = numpy.array(varb)
#print utili per vedere se funziona
print ("varb = ", varb)
print("error = ", numpy.sqrt(varb))

#plotto l'errore sulla <delta_aper> in funzione della blok length
pylab.title('errore su <delta_aper>')
pylab.xlabel('block_length')
pylab.ylabel('error')
pylab.grid(color = "gray")
plt.plot(lng, numpy.sqrt(varb))
plt.show()

