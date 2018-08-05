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
Beta = np.array([0.2275])
L = np.array([10])
T = 10

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

	cluster = np.zeros((l, l, T))
	for i in range(0, l):
		for j in range(0, l):
			for k in range(0, T):
				cluster[i][j][k] = 1

	#Questo pezzo di codice non può essere più ottimizzato di così.
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
	#Questo ciclo moralmente è un do while.
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

	#print("cluster: ", cluster)


	#TODO Bisogna implementare questo controllo in maniera più efficiente

	sheet = 0
	for k in range(0, T):
		sheet = 1
		for j in range(0, l):
			for i in range(0, l):
				if(cluster[i][j][k] == -1): 
					sheet = 0
					break #Risarmia di finire il ciclo su l ma non riesco a risparmiare quello su j
			if(sheet == 0): break 
		if(sheet ==  1): 
			for i in range(0, l):
				for j in range(0, l):
					for s in range(k, T):
						latt[i][j][s] = latt[i][j][s] * cluster[i][j][s]
			boundary = boundary * (-1)
			break

	return boundary


#TODO: Si può ottimizzare e eliminare la variabile sheet.

			
#TODO: Implementare le routine in c all'interno di python. Costruire uno stack interno al programma per trasformare la ricorsione (lenta) in iterazione (più veloce).

#e possibile pttimizzare il codice non passando i, j, k ?? 
#potrei ottimizzare non dovendo passare la variabile spin ma ricordando che la chiamata precedente cambia il segno a link, ma così facendo il codice sarebbe poco leggibile
#perchè dovrei scrivere latt[i+1][j][k] != latt[i][j][k]

def clusterize(latt, cluster, boundary, l , T, beta, i, j, k):
	#TODO Questi tre blocchi di codice uguali si possono scrivere una volta sola introducendo un ciclo e un vettore coord[0], coord[1], coord[3] inizalizzato
	#ai valori i, j, k e modificato in un ciclo sul numero di dimensioni dello spazio.
	size = np.array([l, l, T])
	#print(size)
	#Ci sono tre dimensioni
	for d in range(0, 3):
		#print(d)
		#print(coord)
		for a in [-1, +1]:
			bound = 1
			coord = np.array([i, j, k])
			#Queste condizioni servondo perchè lo spazio e periodico e perche alcuni accoppiamenti J possono essere antiferromagnetici.
			b = coord[d]
			if((coord[d] + a == size[d]) or (coord[d] + a == -1)): 
				b = size[d] - 1 - (coord[d] + a)
				if(d == 2): bound = boundary
			coord[d] = a+b
			p = 1.0-np.exp(-beta * (1.0 + bound*latt[i][j][k]*latt[coord[0]][coord[1]][coord[2]]))
			#Se la posizione considerata già non appartiene al cluster allora la inserisco
			if (cluster[coord[0]][coord[1]][coord[2]] == 1 and np.random.random_sample() < p):
				cluster[coord[0]][coord[1]][coord[2]] = -1
				clusterize(latt, cluster, boundary, l , T, beta, coord[0], coord[1], coord[2])

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
		boundary = - 1
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
			#boundary = surfaceCluster(latt, boundary, l, T, beta)
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

			print(latt)

		

		M2 = 0.0
		M4 = 0.0
		for i in range(0, N):
			M2 += magn[i]**2		
		M2 = M2/aper

		for i in range(0, N):
			M4 += magn[i]**4		
		M4 = M4/aper

		
		print(M2, per, aper)

