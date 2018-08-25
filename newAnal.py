import math
import numpy
import numpy
import pylab
from scipy.optimize import curve_fit
import math
import scipy.stats
beta, l, F, error = pylab.loadtxt("data.txt", unpack = "True")

temp = []
surf = []
err = []

def freeEnergy(l, a, c):
	return c + a*l**2
	#return c + d*(1/l) + b*l + a*l**2
	#return c + d*(1/l)  + a*l**2
	
def surfaceTension(t, sigma_0, t_c, exp):
	return sigma_0*numpy.float_power(abs(t - t_c)/t_c, exp)
	
	
def search(element, array):
	flag = 0
	for i in range(len(array)):
		if(array[i] == element): flag = 1
	return flag

for i in range(len(beta)):
	if(search(beta[i], temp) == 0):
		temp.append(beta[i])
		
temp = numpy.array(temp)


for i in range(len(temp)):
	l_temp = []
	F_temp = []
	error_temp = []
	for j in range(len(l)):
		if(beta[j] == temp[i]):
			l_temp.append(l[j])
			F_temp.append(F[j])
			error_temp.append(error[j])
	
	l_temp = numpy.array(l_temp)
	F_temp = numpy.array(F_temp)
	error_temp = numpy.array(error_temp)
	
	#Eseguo il grafico
	
	#pylab.figure(num=None, figsize=(12, 6), dpi=80, facecolor='w', edgecolor='k')
	#pylab.rc('font',size=13)
	#pylab.title('Energia libera alla temperatura ' + str(temp[i]), fontsize = "16")
	#pylab.xlabel('L', size = "14")
	#pylab.ylabel('freeEnergy', size = "14")
	#pylab.grid(color = "gray")
	#pylab.errorbar(l_temp, F_temp, error_temp, 0, "o", color="black")
	 
	#Eseguo il fit

	#error = unumpy.std_devs(I_D)
	init = numpy.array([0.2, 0.3])
	#Errori tutti statistici
	par, cov = curve_fit(freeEnergy, l_temp, F_temp, init, error_temp, absolute_sigma = "true")
	#trattazione statistica degli errori
	print(par, cov)
	
	#Di nuovo co capisco il chi quadro, non cambia nulla se cambio da true a false
	a = par[0]
	c = par[1]
	#d = par[2]
	#b = par[3]

	surf.append(a)
	err.append(math.sqrt(cov[0][0]))
	
	chisq = ((F_temp-freeEnergy(l_temp, a, c))/error_temp)**2
	somma = sum(chisq)
	
	ndof = len(F_temp) - 2 #Tolgo due parametri estratti dal fit
	
	p=1.0-scipy.stats.chi2.cdf(somma, ndof)
	
	print("beta = %f, Chisquare/ndof = %f/%d" % (temp[i], somma, ndof))
	print("p = ", p)
	
	div = 1000
	bucket = numpy.array([0.0 for i in range(div)])
	retta = numpy.array([0.0 for i in range(div)])
	inc = (l_temp.max()-l_temp.min())/div 
	for k in range(len(bucket)):
	        bucket[k]=float(k)*inc + l_temp.min()
	        retta[k] = freeEnergy(bucket[k], par[0], par[1])
	
	
	#pylab.plot(bucket, retta, color = "red")
		
	#pylab.savefig("plot" + str(temp[i]) + ".png", dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format=None, transparent=False, bbox_inches=None, pad_inches=0.1, frameon=None)

	#pylab.show()
	
surf = numpy.array(surf)
err = numpy.array(err)

#taglio i primi valori

#temp = numpy.delete(temp, numpy.array([6, 7]))
#surf = numpy.delete(surf, numpy.array([6, 7]))
#err = numpy.delete(err, numpy.array([6, 7]))

pylab.figure(num=None, figsize=(12, 6), dpi=80, facecolor='w', edgecolor='k')
pylab.rc('font',size=13)
pylab.title('Surface tension', fontsize = "16")
pylab.xlabel('beta', size = "14")
pylab.ylabel('sigma', size = "14")
pylab.grid(color = "gray")
err = err*10
pylab.errorbar(temp, surf, err, 0, "o", color="black")
err = err/10

#Eseguo il fit

#error = unumpy.std_devs(I_D)
init = numpy.array([1.42, 0.221, 2*0.629971])
#Errori tutti statistici
par, cov = curve_fit(surfaceTension, surf, temp, init, err, absolute_sigma = "true")
#trattazione statistica degli errori
print(par, cov)

par = init
	
#Di nuovo co capisco il chi quadro, non cambia nulla se cambio da true a false
sigma_0 = par[0]
t_c = par[1]
exp = par[2]

print(surfaceTension(temp, par[0], par[1], par[2]))

chisq = ((surf - surfaceTension(temp, par[0], par[1], par[2]))/err)**2
somma = sum(chisq)
	
ndof = len(temp) - 2 #Tolgo due parametri estratti dal fit
	
p=1.0-scipy.stats.chi2.cdf(somma, ndof)
	
print("Chisquare/ndof = %f/%d" % (somma, ndof))
print("p = ", p)
	
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



