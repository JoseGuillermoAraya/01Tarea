import numpy as np
from pylab import *
from scipy import *
import matplotlib.pyplot as plt
#importar datos
data=np.loadtxt("sun_AM0.dat")
#transformar a unidades cgs (erg*s-1*cm-2*cm-1) para flujo y angstrom para longitud de onda
#nm->armstrong. *10
#(W*m-2*nm-1)->(erg*s-1*cm-2*cm-1). 10^7*(10^2)^-2*(10^2)^-1
w_length=data[:,0]
w_length_armstrong=data[:,0]*10.0
flux=data[:,1]
flux_cgs=data[:,1]*10.0

#graficar en escala logaritmica
plt.yscale('log')
plt.xscale('log')

plt.plot(w_length_armstrong,flux_cgs,label='Espectro Solar')
plt.xlabel('Longitud de onda (Armstrong)')
plt.ylabel('Flujo ($erg*s^{-1} cm^{-3}$)')
plt.title('Espectro del Sol')
plt.legend(loc='lower left')
grid(True)
savefig("espectro.png")
show()

#Calculo luminosidad total
##calculo de la luminosidad por unidad de area integrando el espectro
#----------------------------------
#el problema estaba en no escalar bien el arreglo a integrar.
#al implementar el flujo en un arreglo, el salto es de 1, f[0], f[1], etc. Entonces la diferencia en el eje x entre x[0] y x[1] debe ser 1 también
#al intentar integrar el flujo en cgs, el salto en x entre x[0] y x[1] era 10^-7
#para solucionarlo, ocupare flujo en las unidades iniciales (salto en x calza con 1) y luego transformaré las unidades para obtener luminosidad en cgs
#----------------------------------
#error 2, asumir que el salto en los datos era constante
#-----------------------------------------------------

#metodo de los trapecios para integrar
#-------------------------------------
#hacer la integración paso por paso, y en cada momento usar el salto en x calculado directamente del arreglo como w_length[i+1]-w_length[i]
#------------------------------------
ini=0
fin=len(flux)
suma=0
for i in xrange(ini,fin-1):
    suma+=(flux[i]+flux[i+1])*(w_length[i+1]-w_length[i])/2

luminosidad_area=suma*10**3 #transformar unidades a cgs
print('luminosidad por unidad de area del Sol (erg*s^-1*cm^-2)= '+str(luminosidad_area))
##calculo de la luminosidad total, multiplicando por la superficie de la esfera con radio 1UA(1.496*10^13 cm)
luminosidad_total=4.0*np.pi*(1.496*10.0**13.0)**2*luminosidad_area
print('luminosidad total del Sol (erg*s^-1)= '+str(luminosidad_total))

#-----------------------------------------------------------------------
