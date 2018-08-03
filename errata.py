#TODO: attenzione ignora l'ultimo layer dove gli accoppiamenti potrebbero essere antiferromagnetici. Questo non dovrebbe fare la differenza, ma 
	#eventualmente provare a correggere. fatto
	flip = 0
	for k in range(0, T):
		sheet = 1
		for j in range(0, l):
			for i in range(0, l):
				#Se non viene eliminato il bond allora setta sheet a false
				#Non sono arrivato sul bordo
				if(k < T-1):
					if(flip == 0): 
						if(np.random.random_sample() > np.exp(-beta * (1 + latt[i][j][k]*latt[i][j][k+1]))):
							print(np.exp(-beta * (1 + latt[i][j][k]*latt[i][j][k+1])))
							sheet = 0
							continue
						#Questo if entra a partire dal layer successivo a quello in cui si 
					else: latt[i][j][k] = latt[i][j][k] *(-1)
				#Sono arrivato sul bordo
				else: #if(k = T -1):
					if(flip == 0): 
						if(np.random.random_sample() > np.exp(-beta * (1 + boundary * latt[i][j][k]*latt[i][j][1]))):
							sheet = 0
							continue
					else:
						latt[i][j][k] = latt[i][j][k] *(-1)
						boundary = boundary* (-1)

		if(sheet == 1): 
			flip = 1
			print("flip")
			print(latt)
