import numpy as np
import random
import sys

#Vorrei rimuovere il limite di ricorsione di Python
sys.setrecursionlimit(20000)

#Successione delle dimensioni spaziali
#Beta = np.array([0.2, 0.2, 0.2, 0.2, 0.2])
#L = np.array([10, 20, 30, 40, 50])
Beta = np.array([0.10])
L = np.array([5])
T = 5

def singleCluster(latt, boundary, l, T, beta):
	#Scelgo una posizione casuale del reticolo

	#Reticolo dei cluster
	cluster = np.zeros((l, l, T))
	for i in range(0, l):
		for j in range(0, l):
			for k in range(0, T):
				cluster[i][j][k] = 1

	i = np.random.randint(0, l)
	j = np.random.randint(0, l)
	k = np.random.randint(0, T)
	spin = latt[i][j][k]

	cluster[i][j][k] = -1 #Cambio subito valore al primo sito
	#Il parametro 1 indica che la funzione clusterize è stata chiamata in modalità cluster singolo, 0 significa che è stata chiamata in modalità surfaceCluster, la differenza è che ci sono due diverse probabilità di produrre un cluster (approfondire il perché)
	clusterize(latt, cluster, boundary, l, T, beta, i, j, k, spin, 1)
	
	#Effettivo inserimento del cluster nel reticolo
	for i in range(0, l):
		for j in range(0, l):
			for k in range(0, T):
				latt[i][j][k] = latt[i][j][k] * cluster[i][j][k]


def surfaceCluster(latt, boundary, l, T, beta):
	#Inizializzo le varibili che tengono trccia dei cluster

	cluster = np.random.rand(l, l, T)
	for i in range(0, l):
		for j in range(0, l):
			for k in range(0, T):
				cluster[i][j][k] = 1

	#E molto difficile far partire i cluster in maniera randomica, perche bisogna prendere una posizione causale su un una superficie dalla forma irregolare. Non lo implementiamo
	#chissa se funziona comunque.
	for i in range(0, l):
		for j in range(0, l):
			#Cioe se non appartiene a nessun cluster
			if(cluster[i][j][T-1] == 1):
				cluster[i][j][T-1] = -1 #Inizializza il primo elemento del cluster
				clusterize(latt, cluster, boundary, l, T, beta, i, j, T-1, latt[i][j][T-1], 0)


	#print("cluster: ", cluster)


	#TODO Bisogna implementare questo controllo in maniera piu efficiente

	for k in range(0, T):
		sheet = 1
		for j in range(0, l):
			for i in range(0, l):
				if(cluster[i][j][k] == -1): sheet = 0 
		if(sheet ==  1): break

	if(sheet == 1): 
		for i in range(0, l):
			for j in range(0, l):
				for k in range(0, T):
					latt[i][j][k] = latt[i][j][k] * cluster[i][j][k]
					
		boundary = boundary * (-1)
	
	return boundary


	

#TODO: Si puo ottimizzare e eliminare la variabile sheet.

			
#TODO: Implementare le routine in c all'interno di python. Costruire uno stack interno al programma per trasformare la ricorsione (lenta) in iterazione (piu veloce).

#e possibile pttimizzare il codice non passando i, j, k ?? 
#potrei ottimizzare non dovendo passare la variabile spin ma ricordando che la chiamata precedente cambia il segno a link, ma cosi facendo il codice sarebbe poco leggibile
#perche' dovrei scrivere latt[i+1][j][k] != latt[i][j][k]
def clusterize(latt, cluster, boundary, l , T, beta, i, j, k, spin, setting):

	#TODO Questi tre blocchi di codice uguali si possono scrivere una volta sola introducendo un ciclo e un vettore coord[0], coord[1], coord[3] inizalizzato
	#ai valori i, j, k e modificato in un ciclo sul numero di dimensioni dello spazio.

	#print(i, j, k)
	for i1 in [-1, +1]:

		#Queste condizioni servondo perche lo spazio e periodico e perche alcuni accoppiamenti J possono essere antiferromagnetici.
		bound = 1
		i2 = i
		if(i + i1 == l): 
			i2 = -1
			bound = boundary
		if(i + i1 == -1): 
			i2 = l
			bound = boundary

		if(setting == 1): p = 1.0 - np.exp(-bound*2*beta)
		else: p = 1-np.exp(-beta * (1 + bound *  latt[i][j][k]*latt[i2+i1][j][k]))


		if (latt[i2+i1][j][k] == spin and cluster[i2+i1][j][k] == 1 and np.random.random_sample() < p):
			#print(np.random.random_sample(), 1.0 - np.exp(-2*beta))
			#latt[i2+i1][j][k] = latt[i2+i1][j][k]*(-1) #cambia segno
			cluster[i2+i1][j][k] = -1
			clusterize(latt, cluster, boundary, l , T, beta, i2+i1, j, k, spin, setting)

	for j1 in [-1, +1]:
		bound = 1
		j2 = j
		if(j + j1 == l): 
			j2 = -1
			bound = boundary
		if(j + j1 == -1): 
			j2 = l
			bound = boundary

		if(setting == 1): p = 1.0 - np.exp(-bound*2*beta)
		else: p = 1-np.exp(-beta * (1 + bound *  latt[i][j][k]*latt[i][j2+j1][k]))

		if (latt[i][j2+j1][k] == spin and cluster[i][j2+j1][k] == 1 and np.random.random_sample() < p):
			cluster[i][j2+j1][k] = -1 #cambia segno
			clusterize(latt, cluster, boundary, l , T, beta, i, j2+j1, k, spin, setting)

	for k1 in [-1, +1]:
		bound = 1
		k2 = k
		if(k + k1 == T): 
			k2 = -1
			bound = boundary
		if(k + k1 == -1): 
			k2 = T
			bound = boundary

		if(setting == 1): p = 1.0 - np.exp(-bound*2*beta)
		else: p = 1-np.exp(-beta * (1 + bound *  latt[i][j][k]*latt[i][j][k2+k1]))


		if (latt[i][j][k2+k1] == spin and cluster[i][j][k2+k1] == 1 and np.random.random_sample() < p):
			cluster[i][j][k2+k1] =-1 #cambia segno
			clusterize(latt, cluster, boundary, l , T, beta, i, j, k2+k1, spin, setting)

#Funzioni autoesplicative

def coldStart(latt, l, T, spin):
	for i in range(0, l):
		for j in range(0, l):
			for k in range(0, T):
				latt[i][j][k] = spin
	return latt

def hotStart(latt, l, T):
	for i in range(0, l):
		for j in range(0, l):
			for k in range(0, T):
				if(np.random.randint(2) == 0): latt[i][j][k] = -1
				else: latt[i][j][k] = +1
	return latt



for l in L:
	for beta in Beta:

		#Inizializzo le boundary condition periodiche
		boundary = + 1
		latt = np.zeros((l, l, T))
		#latt = coldStart(latt, l, T, 1)
		latt = hotStart(latt, l, T)
		
		
		N = 20
		magn = np.zeros(N)

		per = 0
		aper = 0

		for rep in range(0, N):
			somma = 0

			#print(rep) 

			singleCluster(latt, boundary, l, T, beta)
			
			if(rep%10 == 0): boundary = surfaceCluster(latt, boundary, l, T, beta)

			if(boundary == 1): per += 1
			else: aper += 1

			for i in range(0, l):
				for j in range(0, l):
					for k in range(0, T):
						somma += latt[i][j][k] 
			magn[rep] = somma/(l*l*T)

		#print("prima:", latt)
		#print("dopo:", latt)

		M = 0.0
		for i in range(N/2, N):
			M += magn[i]		
		M = M/(N/2)

		print(M, per, aper)

