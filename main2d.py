#Versione 2D dell'algortimo a cluster di wolff implementato per provare la tecnica. Funziona e in effetti faceno un po' di prove e misurando la magnetizzazione si osserva che vi è una tenperatura critica al valore predetto da onsager.

import numpy as np
import random
import sys

sys.setrecursionlimit(20000)

#Successione delle dimensioni spaziali
#Beta = np.array([0.2, 0.2, 0.2, 0.2, 0.2])
#L = np.array([10, 20, 30, 40, 50])
Beta = np.array([0.10])
L = np.array([100])

def singleCluster(latt, l, beta):
	#scrivere piu elegantemente
	#Scelgo una posizione casuale del reticolo
	i = np.random.randint(0, l)
	j = np.random.randint(0, l)
	spin = latt[i][j]
	#print(i, j, k)
	latt[i][j] = latt[i][j]*(-1) #Cambio subito valore al primo sito
	clusterize(latt, l, beta, i, j, spin)


#in c tutti questi passaggi di vettore dovrebbero essere implementati come puntatori, spero il python lo faccia
#e possibile pttimizzare il codice non passando i, j, k ?? 
#potrei ottimizzare non dovendo passare la variabile spin ma ricordando che la chiamata precedente cambia il segno a link, ma cosi facendo il codice sarebbe poco leggibile
#perche' dovrei scrivere latt[i+1][j][k] != latt[i][j][k]
def clusterize(latt, l , beta, i, j, spin):
	#print(i, j, k)
	for i1 in [-1, +1]:
		i2 = i
		if(i + i1 == l): i2 = -1
		if(i + i1 == -1): i2 = l
		if (latt[i2+i1][j] == spin and np.random.random_sample() < 1.0 - np.exp(-2*beta)):
			#print(np.random.random_sample(), 1.0 - np.exp(-2*beta))
			latt[i2+i1][j] = latt[i2+i1][j]*(-1) #cambia segno
			clusterize(latt, l, beta, i2+i1, j, spin)

	for j1 in [-1, +1]:
		j2 = j
		if(j + j1 == l): j2 = -1
		if(j + j1 == -1): j2 = l
		if (latt[i][j2+j1] == spin and np.random.random_sample() < 1.0 - np.exp(-2*beta)):
			latt[i][j2+j1] = latt[i][j2+j1]*(-1) #cambia segno
			clusterize(latt, l, beta, i, j2+j1, spin)


for l in L:
	for beta in Beta:

		#TODO inizializza lattice di numeri casuali tra 0 e 1
		latt = np.random.rand(l, l)
		for i in range(0, l):
			for j in range(0, l):
					latt[i][j] = 1
		#arrotonda all'intero piu vicino
		#np.round(latt, 0)
		#creazione dei link

		N = 20000
		magn = np.zeros(N)
		link = np.zeros((l, l , 3))
		for rep in range(0, N):
			somma = 0
			print(rep)
			singleCluster(latt, l, beta)
			for i in range(0, l):
				for j in range(0, l):
						somma += latt[i][j]
			magn[rep] = somma/(l*l)
		#print(latt)

		M = 0.0
		for i in range(10000, N):
			M += np.absolute(magn[i])		
		M = M/(N-10000)

		print(M)
