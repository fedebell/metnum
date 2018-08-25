#Eseguo il fit

#error = unumpy.std_devs(I_D)
init = numpy.array([1.2, 0.221, 1.2])
#Errori tutti statistici
par, cov = curve_fit(surfaceTension, surf, temp, init, err, absolute_sigma = "true")
#trattazione statistica degli errori
print(par, cov)
	
#Di nuovo co capisco il chi quadro, non cambia nulla se cambio da true a false
sigma_0 = par[0]
t_c = par[1]
exp = par[2]

chisq = ((surf_surfaceTension(temp, a, b))/err)**2
somma = sum(chisq)
	
ndof = len(temp) - 2 #Tolgo due parametri estratti dal fit
	
p=1.0-scipy.stats.chi2.cdf(somma, ndof)
	
print("Chisquare/ndof = %f/%d" % (somma, ndof))
print("p = ", p)
	
div = 1000
bucket = numpy.array([0.0 for i in range(div)])
retta = numpy.array([0.0 for i in range(div)])
inc = (l_temp.max()-l_temp.min())/div 
for k in range(len(bucket)):
        bucket[k]=float(k)*inc + temp.min()
	retta[k] = surfaceTension(bucket[k], par[0], par[1])
	
pylab.plot(bucket, retta, color = "red")
		
#pylab.savefig("plot" + str(temp[i]) + ".png", dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format=None, transparent=False, bbox_inches=None, pad_inches=0.1, frameon=None)

#Plotto un grafico:

t_c = 0.221

div = 1000
bucket = numpy.array([0.0 for i in range(div)])
retta = numpy.array([0.0 for i in range(div)])
inc = (temp.max()-temp.min())/div 
for k in range(len(bucket)):
        bucket[k]=float(k)*inc + t_c
	retta[k] = 1.3*(bucket[k]-t_c)**(2*0.629971)
	
pylab.plot(bucket, retta, color = "red")
		
#pylab.savefig("plot" + str(temp[i]) + ".png", dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format=None, transparent=False, bbox_inches=None, pad_inches=0.1, frameon=None)
