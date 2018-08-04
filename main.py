# coding=utf-8
#Questo codice riproduce i risultati dell'articolo tranne quando T è piccolo. Per qualche ragione non riesce mai a flippare gli spin con il surfaceCluster (ma il risultato della magnetizzazione media è corretto). Per esempio se nell'articolo sono riportate 70 configurazioni con condizioni antiperiodiche io ne trovo 0. Questo problema si ha solo per piccoli T (circa 5) perché non appena T diventa più grande questo algoritmo prouce risultati in buon accordo con quelli dell'articolo. Ovviamente deve essere sottoposto a ottimizzazione e riscritto in c++.

import numpy as np
import random
import sys

#Vorrei rimuovere il limite di ricorsione di Python
sys.setrecursionlimit(200000)

#Successione delle dimensioni spaziali
#Beta = np.array([0.2, 0.2, 0.2, 0.2, 0.2])
#L = np.array([10, 20, 30, 40, 50])
Beta = np.array([0.2391])
L = np.array([4])
T = 12

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

	cluster[i][j][k] = -1 #Cambio subito valore al primo sito
	clusterize(latt, cluster, boundary, l, T, beta, i, j, k)
	
	#Effettivo inserimento del cluster nel reticolo
	for i in range(0, l):
		for j in range(0, l):
			for k in range(0, T):
				latt[i][j][k] = latt[i][j][k] * cluster[i][j][k]


def surfaceCluster(latt, boundary, l, T, beta):
	#Inizializzo le varibili che tengono traccia dei cluster

	cluster = np.random.rand(l, l, T)
	for i in range(0, l):
		for j in range(0, l):
			for k in range(0, T):
				cluster[i][j][k] = 1

	#E' molto difficile far partire i cluster in maniera randomica, perchè bisogna prendere una posizione causale su un una superficie dalla forma irregolare. Non lo implementiamo
	#chissà se funziona comunque. FORSE NON FUNZIONA PER QUESTO
	#for i in range(0, l):
	#	for j in range(0, l):
	#		#Cioe se non appartiene a nessun cluster
	#		if(cluster[i][j][T-1] == 1):
	#			cluster[i][j][T-1] = -1 #Inizializza il primo elemento del cluster
	#			clusterize(latt, cluster, boundary, l, T, beta, i, j, T-1)

	#RANDOMIZED CLUSTER EXTRACTION


	#Basta settare questa a 1 per farlo ripartire
	nonclusterized = 1
	while(nonclusterized == 1):

		#Si tratta di un do while che cerca un sito non clusterizzato
		i = np.random.randint(0, l)
		j = np.random.randint(0, l)
		while(cluster[i][j][T-1] == -1):
			i = np.random.randint(0, l)
			j = np.random.randint(0, l)

		cluster[i][j][T-1] = -1 #Inizializza il primo elemento del cluster
		clusterize(latt, cluster, boundary, l, T, beta, i, j, T-1)

		nonclusterized = 0

		#Cerca se ce ne sono ancora
		for i in range(0, l):
			for j in range(0, l):
				if(cluster[i][j][T-1] == 1): 
					nonclusterized = 1
					break

	print("cluster: ", cluster)


	#TODO Bisogna implementare questo controllo in maniera più efficiente

	position = 0
	sheet = 0
	for k in range(0, T):
		sheet = 1
		for j in range(0, l):
			for i in range(0, l):
				if(cluster[i][j][k] == -1): sheet = 0 
		if(sheet ==  1): 
			position = k
			break
	
	
	#Trovato l'errore devo invertire solo la parte di sotto del cluster

	if(sheet == 1): 
		for i in range(0, l):
			for j in range(0, l):
				for k in range(position, T):
					latt[i][j][k] = latt[i][j][k] * cluster[i][j][k]
					
		boundary = boundary * (-1)
	
	return boundary


#TODO: Si può ottimizzare e eliminare la variabile sheet.

			
#TODO: Implementare le routine in c all'interno di python. Costruire uno stack interno al programma per trasformare la ricorsione (lenta) in iterazione (più veloce).

#e possibile pttimizzare il codice non passando i, j, k ?? 
#potrei ottimizzare non dovendo passare la variabile spin ma ricordando che la chiamata precedente cambia il segno a link, ma così facendo il codice sarebbe poco leggibile
#perchè dovrei scrivere latt[i+1][j][k] != latt[i][j][k]
def clusterize(latt, cluster, boundary, l , T, beta, i, j, k):

	#TODO Questi tre blocchi di codice uguali si possono scrivere una volta sola introducendo un ciclo e un vettore coord[0], coord[1], coord[3] inizalizzato
	#ai valori i, j, k e modificato in un ciclo sul numero di dimensioni dello spazio.

	#print(i, j, k)
	for i1 in [-1, +1]:

		#Queste condizioni servondo perchè lo spazio e periodico e perche alcuni accoppiamenti J possono essere antiferromagnetici.
		i2 = i
		if(i + i1 == l): 
			i2 = -1
			
		if(i + i1 == -1): 
			i2 = l

		p = 1.0-np.exp(-beta * (1.0 + latt[i][j][k]*latt[i2+i1][j][k]))

		#Se la posizione considerata già non appartiene al cluster allora la inserisco
		if (cluster[i2+i1][j][k] == 1 and np.random.random_sample() < p):
			#print(np.random.random_sample(), 1.0 - np.exp(-2*beta))
			#latt[i2+i1][j][k] = latt[i2+i1][j][k]*(-1) #cambia segno
			cluster[i2+i1][j][k] = -1
			clusterize(latt, cluster, boundary, l , T, beta, i2+i1, j, k)

	for j1 in [-1, +1]:

		j2 = j
		if(j + j1 == l): 
			j2 = -1
		if(j + j1 == -1): 
			j2 = l

		p = 1.0-np.exp(-beta * (1.0 + latt[i][j][k]*latt[i][j2+j1][k]))

		if (cluster[i][j2+j1][k] == 1 and np.random.random_sample() < p):
			cluster[i][j2+j1][k] = -1 #cambia segno
			clusterize(latt, cluster, boundary, l , T, beta, i, j2+j1, k)

	for k1 in [-1, +1]:
		bound = 1
		k2 = k
		if(k + k1 == T): 
			k2 = -1
			bound = boundary
		if(k + k1 == -1): 
			k2 = T
			bound = boundary

		p = 1.0-np.exp(-beta * (1.0 + bound *  latt[i][j][k]*latt[i][j][k2+k1]))

		if (cluster[i][j][k2+k1] == 1 and np.random.random_sample() < p):
			cluster[i][j][k2+k1] = -1 #cambia segno
			clusterize(latt, cluster, boundary, l , T, beta, i, j, k2+k1)

#Funzioni autoesplicative

def coldStart(latt, l, T, spin):
	for i in range(0, l):
		for j in range(0, l):
			for k in range(0, T):
				latt[i][j][k] = spin

def hotStart(latt, l, T):
	for i in range(0, l):
		for j in range(0, l):
			for k in range(0, T):
				if(np.random.randint(2) == 0): latt[i][j][k] = -1
				else: latt[i][j][k] = +1


def energy(latt, l, T, boundary):
	return

for l in L:
	for beta in Beta:

		#Inizializzo le boundary condition periodiche
		boundary = + 1
		latt = np.zeros((l, l, T))
		coldStart(latt, l, T, 1)
		#hotStart(latt, l, T)
		
		N = 1000

		magn = np.zeros(N)

		per = 0
		aper = 0

		for rep in range(0, N):
			somma = 0.0

			print(rep) 
			
			#print("prima:", latt)
			boundary = surfaceCluster(latt, boundary, l, T, beta)
			#print("dopo:", latt)

			singleCluster(latt, boundary, l, T, beta)

			if(boundary == 1): 
				per += 1
				for i in range(0, l):
					for j in range(0, l):
						for k in range(0, T):
							somma += latt[i][j][k] 
				magn[rep] = somma/(l*l*T)
			else: aper += 1

			

		M2 = 0.0
		M4 = 0.0
		for i in range(0, N):
			M2 += magn[i]**2		
		M2 = M2/per

		for i in range(0, N):
			M4 += magn[i]**4		
		M4 = M4/per

		
		print(M2, per, aper)

