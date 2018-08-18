import math
import numpy
import numpy
import pylab
import math
import scipy.stats
import matplotlib.pyplot as plt

mag, abs_mag, bc = pylab.loadtxt("C:\cygwin64\metnum\ciao.txt", unpack =True)
n = range (len(mag))
mag_per = []

for i in range(len(mag)):
	if bc[i]==1:
		mag_per.append(mag[i])
		
n_per = range(len(mag_per))

avg = bc.sum()/len(bc) #questa media è la stessa per tutti i dati blockati

i = 1
c = 0
var = []
n = []
#blocking sulla mag
while i < numpy.sqrt(len(bc)): #scorro sui posssibli divisori del numero di eventi
	if len(bc)%i == 0: #seleziono un divisore buono
		#creo array i cui elementi sono le medie dei blocchetti desiderati
		n.append([i]) #salvo la dimensione dei blocchi
		block = []
		for j in range(len(mag)//i): #scelgo il blocco
			temp_sum = 0
			for k in range(i): #mi muovo nel blocco e metto le somme ina una variabile temporanea
				temp_sum +=  bc[j*i+k] 
			temp_sum = temp_sum/i #medio la somma
			block.append(temp_sum)			
		#calcolo varianza nelle nuove variabili 
		var.append(((numpy.array(block)-avg)**2).sum()/(len(block)*(len(block)-1)))
	i += 1
#plotto la varianza in funzione del numero di elementi nel blocchetto
pylab.title('magnetization in markov chain')
pylab.xlabel('N')
pylab.ylabel('var')
pylab.grid(color = "gray")
plt.plot(n, numpy.sqrt(var))
plt.show()
