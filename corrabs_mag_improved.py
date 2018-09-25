import math
import numpy
import numpy
import pylab
import math
import scipy.stats
import matplotlib.pyplot as plt

all_mag, all_abs_mag, all_bc = pylab.loadtxt("miscData//0.2275_18_18_54.txt", unpack =True)

varb = []
n = []

if len(all_abs_mag) % 2 == 1:
	print ("error!")
	 
abs_mag = numpy.empty(len(all_abs_mag)//2)


#creo array solo con cose pari (supersloppy)
for i in range(len(all_mag)):
	if i % 2 == 0:
		abs_mag[i//2] = all_abs_mag[i]


totsum = abs_mag.sum()

#refresh: se non parte da 1 fa casino a dividere in blocchi da 0
i = 1
#blocking 
while i < numpy.sqrt(len(abs_mag)): #scorro sui possibli divisori del numero di eventi
	#creo array i cui elementi sono le medie dei blocchetti desiderati
	n.append(i) #salvo la dimensione dei blocchi
	block = []
	avg = totsum
	for cnt in range(len(abs_mag) - (len(abs_mag)//i) * i):
		avg -= abs_mag[(len(abs_mag)//i) * i + cnt] #tolgo gli elementi in eccesso ,che comunque non vengono passati dal for successivo...
	avg = avg / (i* (len(abs_mag)//i) ) #true avg	
	for j in range(len(abs_mag)//i): #scelgo il blocco
		temp_sum = 0
		for k in range(i): #mi muovo nel blocco e metto le somme ina una variabile temporanea
			temp_sum +=  abs_mag[j*i+k] 
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
