import math
import numpy
import numpy
import pylab
from scipy.optimize import curve_fit
import math
import scipy.stats
#carico array di beta, dimensioni, misura free energy e relativo errore
beta, l, F, error = pylab.loadtxt("newData.txt~", unpack = "True")

temp = []
surf = []
err = []

def freeEnergy(l, a, b):
	#return b + a*l**2
	#return c + d*(1/l)**2 + b*l + a*l**2
	#return c + d*(1/l)  + a*l**2
	#return c + b*(1/l)**2  + a*l**2
	#return c + b * numpy.log(l) + a*l**2
	#return c + d*(1/l)**2 + b*(1/l)**4 + e*(1/l)**6 + a*l**2
	#return c + d*(1/l)**2 + b*(1/l)**4 + e*(1/l)**6 + f*l + a*l**2
	#return a*l**2 + c - numpy.log( 1.0 + b*(1/l)**2) 
	#return a*(l**2) + c + b/(a*(l**2))
	#return a*(l**2) + c + b/(a*(l**2)) + d/(a*(l**2))
	return b + a * l**2 - numpy.log(1.0 + 0.25/( a * l**2))
	
def surfaceTension(b, sigma_0, a_theta, a):
	return sigma_0*numpy.float_power(abs(1 - 0.2218/b), 1.260)*( 1.0 + a_theta * numpy.float_power(abs(1 - 0.2216544/b), 0.51) + a*abs(1 - 0.2216544/b)) 
	
def search(element, array):
	flag = 0
	for i in range(len(array)):
		if(array[i] == element): flag = 1
	return flag

#scorro il file .txt e appendo a temp solo le beta di quelli con beta giusto
for i in range(len(beta)):
	if(search(beta[i], temp) == 0):
		temp.append(beta[i])
		
temp = numpy.array(temp)

temp.sort()

for i in range(len(temp)):
	l_temp = []
	F_temp = []
	error_temp = []
	#riscorro tutto il .txt e metto negli array temporanei i valori di l, F, err con beta giusta 
	for j in range(len(l)):
		if(beta[j] == temp[i]):
			if(search(l[j], l_temp) == 0):
				l_temp.append(l[j])
				F_temp.append(F[j])
				error_temp.append(error[j])
	
	l_temp = numpy.array(l_temp)
	F_temp = numpy.array(F_temp)
	error_temp = numpy.array(error_temp)
	
	#Eseguo il grafico per ogni temperatura, indicizzata proprio da temp
	
	pylab.figure(num=None, figsize=(12, 6), dpi=80, facecolor='w', edgecolor='k')
	pylab.rc('font',size=13)
	pylab.title('Energia libera alla temperatura ' + str(temp[i]), fontsize = "16")
	pylab.xlabel('L', size = "14")
	pylab.ylabel('freeEnergy', size = "14")
	pylab.grid(color = "gray")
	pylab.errorbar(l_temp, F_temp, error_temp, 0, "o", color="black")
	
	 
	#Eseguo il fit per ogni temperatura

	#error = unumpy.std_devs(I_D)
	init = numpy.array([0.2, 0.3])
	#Errori tutti statistici
	par, cov = curve_fit(freeEnergy, l_temp, F_temp, init, error_temp, absolute_sigma = "false")
	#trattazione statistica degli errori
	print(par, cov)
	
	a = par[0]
	b = par[1]
	#c = par[2]
	#d = par[3]
	#e = par[4]
	#f = par[5]

	#creo array surf[i] = surface tension ricavata dal fit a beta[i]
	surf.append(a)
	err.append(math.sqrt(cov[0][0]))
	
	#printo i chi^2
	chisq = ((F_temp-freeEnergy(l_temp, a, b))/error_temp)**2
	somma = sum(chisq)
	
	ndof = len(F_temp) - 2 #Tolgo due parametri estratti dal fit
	
	p=1.0-scipy.stats.chi2.cdf(somma, ndof)
	
	print("beta = %f, Chisquare/ndof = %f" % (temp[i], somma/ndof))
	#print("p = ", p)
	
	
"""	#plotto per ogni temperatura le funzioni di fit coi parametri ricavati
	div = 1000
	bucket = numpy.array([0.0 for i in range(div)])
	retta = numpy.array([0.0 for i in range(div)])
	inc = (l_temp.max()-l_temp.min())/div 
	for k in range(len(bucket)):
	        bucket[k]=float(k)*inc + l_temp.min()
	        retta[k] = freeEnergy(bucket[k], a, b)	
	
	pylab.plot(bucket, retta, color = "red")
		
	#pylab.savefig("plot" + str(temp[i]) + ".png", dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format=None, transparent=False, bbox_inches=None, pad_inches=0.1, frameon=None)

	pylab.show()"""

#print("beta = ", temp)
#print("sigma = ", surf)
	
surf = numpy.array(surf)
err = numpy.array(err)

#taglio i primi valori

#temp = numpy.delete(temp, numpy.array([6, 7]))
#surf = numpy.delete(surf, numpy.array([6, 7]))
#err = numpy.delete(err, numpy.array([6, 7]))

#plotto l'andamento in funzione di beta delle surftension fittate
pylab.figure(num=None, figsize=(12, 6), dpi=80, facecolor='w', edgecolor='k')
pylab.rc('font',size=13)
pylab.title('Surface tension', fontsize = "16")
pylab.xlabel('beta', size = "14")
pylab.ylabel('sigma', size = "14")
pylab.grid(color = "gray")
#err = err*10
pylab.errorbar(temp, surf, err, 0, "o", color="black")
#err = err/10

#Eseguo il fit

#error = unumpy.std_devs(I_D)
init = numpy.array([1.5157, 0.0, 0.0])
#Errori tutti statistici
par, cov = curve_fit(surfaceTension, temp, surf, init, err, absolute_sigma = "false")
#trattazione statistica degli errori
#print(par, cov)

#par = init
	
#Di nuovo co capisco il chi quadro, non cambia nulla se cambio da true a false
sigma_0 = par[0]
a_theta = par[1]
a = par[2]

print(surfaceTension(temp, par[0], par[1], par[2]))

chisq = ((surf - surfaceTension(temp, par[0], par[1], par[2]))/err)**2
somma = sum(chisq)
	
ndof = len(temp) - 3 #Tolgo tre parametri estratti dal fit
	
p=1.0-scipy.stats.chi2.cdf(somma, ndof)
	
print("Chisquare/ndof = %f" % (somma/ndof))
print("p = ", p)

print("Surface tension = ", sigma_0, "pm", numpy.sqrt(cov[0][0]))
print("a_theta = ", a_theta, "pm", numpy.sqrt(cov[1][1]))
print("a = ", a, "pm", numpy.sqrt(cov[2][2]))

div = 1000
bucket = numpy.array([0.0 for i in range(div)])
retta = numpy.array([0.0 for i in range(div)])
inc = (temp.max()-temp.min())/div 
for k in range(len(bucket)):
	bucket[k]=float(k)*inc + temp.min()
	retta[k] = surfaceTension(bucket[k], par[0], par[1], par[2])
	
pylab.plot(bucket, retta, color = "red")
		
#pylab.savefig("plot" + str(temp[i]) + ".png", dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format=None, transparent=False, bbox_inches=None, pad_inches=0.1, frameon=None)
	 

pylab.show()
