import math
import numpy
import numpy
import pylab
import math
import scipy.stats
import matplotlib.pyplot as plt

mag, abs_mag, bc = pylab.loadtxt("C://cygwin64//metnum//ciao.txt", unpack =True)
n = range (len(mag))
mag_per = []

for i in range(len(mag)):
	if bc[i]==1:
		mag_per.append(mag[i])
		
n_per = range(len(mag_per))


i = 1
c = 0
varb = []
n = []

#riscrivo le bc in termini della delta aperiodica
bc = (-1*bc + 1)/2

avg = bc.sum()/len(bc) #questa media è la stessa per tutti i dati blockati
print ("avg = ", avg)

#blocking 
while i < len(bc): #scorro sui posssibli divisori del numero di eventi
	if len(bc)%i == 0: #seleziono un divisore buono
		#creo array i cui elementi sono le medie dei blocchetti desiderati
		n.append(i) #salvo la dimensione dei blocchi
	i += 1

i = 0
varb = numpy.zeros(len(n)) #creo numpy.array di lunghezza pari al numero di tipi di blocco

for i in range(len(n)):
	for j in range(len(mag)//int(n[i])):
		temp_sum = 0
		for k in range(n[i]):
			temp_sum += bc[j*n[i]+k]
		temp_sum = temp_sum/n[i]
		varb[i] = temp_sum

#print utili per vedere se funziona


print ("varb = ", varb)
print ("n = ", n)

#creo variabili utili per computo tau

var = varb[0] / len(mag) 

print("var_single = ", var)

#creo la 2*tau con la formula 2tau=sigma_b^2/sigma_O_i^2 * len(block)

varb = numpy.array(varb)
n = numpy.array(n)

tau = 0.5 * varb * n / var


print ("tau = ", tau)

pylab.title('2*tau')
pylab.xlabel('block_length')
pylab.ylabel('tau')
pylab.grid(color = "gray")
plt.plot(n, tau)
plt.show()
