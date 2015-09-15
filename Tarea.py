import numpy as np
from pylab import *
import scipy.integrate as spi
from astropy import constants as const
import matplotlib.pyplot as plt
#---------------------------------------------------------------------------------------------------------------------------------
##metodo para calcular la integral con la tolerancia deseada
##dado un paso dx, tomamos un rectangulo cuya altura será f(a+dx/4), y de ancho dx/2 (desde a hasta a+dx/2)
def refinar_integral(f,a,b,tolerancia):#funcion,inicio,final,tolerancia
    #calculo integral inicial
    #evitar singularidades
    e=000.1
    a=a+e
    b=b-e

    dif=1000000.0 #iniciar diferencia entre initegrales
    dx=b-a#salto inicial
    #fin=((b-a)/dx) +1
    suma_1=f(a)+f(b)
    integral_1=dx/2*suma_1
    while dif>=tolerancia:
        dx=dx/2#reducir salto
        fin=(b-a)/dx+1#cantidad de iteraciones
        suma=suma_1
        #agregar solamente los términos faltantes (impares) (2i+1)
        for i in xrange(0,int((fin-1)/2)):
            suma+=2*(f(a+dx*(2*i+1)))
        integral=integral_1
        integral=dx/2*suma
        dif=np.fabs(integral-integral_1)
        integral_1=integral
        suma_1=suma
    '''
    for i in xrange(0,int(fin)):#forzar que fin sea int
        suma_1+=(f(a+i*dx+dx/4)+f(a+i*dx+3*dx/4))#para cada punto a, calculamos los dos rectangulos entre a y a+dx
    integral_1=dx/2*suma_1
    #iterar hasta que se alcance la tolerancia esperada
    while dif>=tolerancia:
        dx=dx/2#reducir salto
        fin=(b-a)/dx+1#cantidad de iteraciones
        suma=0

        #agregar solamente los términos faltantes (impares)
        for i in xrange(0,int(fin/2)):
            suma+=(f(a+(2*i+1)*dx+dx/4)+f(a+(2*i+1)*dx+dx/4+dx/2))#para cada punto a, calculamos los dos rectangulos entre a y a+dx
        integral=dx/2*suma

        for i in xrange(0,int(fin)):#forzar que fin sea int
            suma+=(f(a+i*dx+dx/4)+f(a+i*dx+3*dx/4))#para cada punto a, calculamos los dos rectangulos entre a y a+dx
        integral=dx/2*suma
        dif=np.fabs(integral-integral_1)
        integral_1=integral
        suma_1=suma
        '''
    return integral
#---------------------------------------------------------------------------------------------------------------------------------
#importar datos
data=np.loadtxt("sun_AM0.dat")
#transformar a unidades cgs (erg*s-1*cm-2*cm-1) para flujo y angstrom para longitud de onda
#nm->angstrom. *10
#(W*m-2*nm-1)->(erg*s-1*cm-2*cm-1). 10^7*(10^2)^-2*(10^2)^-1
w_length=data[:,0]
w_length_angstrom=data[:,0]*10.0
flux=data[:,1]
flux_cgs=data[:,1]*10.0

#graficar en escala logaritmica
plt.yscale('log')
plt.xscale('log')

plt.plot(w_length_angstrom,flux_cgs,label='Espectro Solar')
plt.xlabel('Longitud de onda (Angstrom)')
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
#error 2, asumir que el salto en los datos era constantsante
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
#Integral numerica de la funcion de planck
T=5778#temperatura solar kelvins
planck_sinctes=lambda x:tan(x)**3/(cos(x)**2*(exp(x)-1))
constantes=2*pi*const.h.cgs/(const.c.cgs)**2 * (const.k_B.cgs*T/const.h.cgs)**4
tolerancia=1000.0
planck_numerico=refinar_integral(planck_sinctes,0.0,pi/2.0,tolerancia)*constantes
print('Integral funcion de planck= '+str(planck_numerico))



#-----------------------------------------------------------------------
#Integracion mediante librerias de scipy
lum_area=trapz(flux,w_length,axis=0)*10**3#factor 10^3 para convertir de W*m^-2 a erg*s^-1*cm^-2
lum_sol=4.0*np.pi*(1.496*10.0**13.0)**2*lum_area
print('luminosidad solar (scipy)= '+str(lum_sol))


integral_planck=spi.quad(planck_sinctes,0,np.pi/2)[0]*constantes
print('integral funcion de planck(scipy)= '+str(integral_planck))
#-----------------------------------------------------------------------
