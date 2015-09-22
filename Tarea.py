import numpy as np
from pylab import *
import scipy.integrate as spi
from astropy import constants as const
import matplotlib.pyplot as plt
import astropy.units as u
#------------------------------------------------------
##método de los trapecios para dos arreglos
def calcular_trapecio(x,y,ini,fin):
    suma=0
    for i in xrange(ini,fin-1):
        suma+=(y[i]+y[i+1])*(x[i+1]-x[i])/2
    return suma

#---------------------------------------------------------------------------------------------------------------------------------
##metodo para calcular la integral con la tolerancia deseada
##
def refinar_integral(f,a,b,tolerancia):#funcion,inicio,final,tolerancia
    #calculo integral inicial
    dif=tolerancia+1 #iniciar diferencia entre initegrales
    dx=b-a#salto inicial
    suma_1=0#debería ser f(a)+f(b), pero en este caso la función a integrar tiende a 0 en a y en b.
    integral_1=dx/2*suma_1
    while dif>=tolerancia:
        dx=dx/2#reducir salto
        fin=(b-a)/dx+1#cantidad de iteraciones
        suma=suma_1
        #agregar solamente los términos faltantes (impares) (2i+1)
        for i in xrange(0,int((fin-1)/2)):
            suma+=2*(f(a+dx*(2*i+1)))
        integral=dx/2*suma
        dif=np.fabs(integral-integral_1)
        integral_1=integral
        suma_1=suma
    return integral
#---------------------------------------------------------------------------------------------------------------------------------
#importar datos
data=np.loadtxt("sun_AM0.dat")
#transformar a unidades cgs (erg*s-1*cm-2*cm-1) para flujo y angstrom para longitud de onda
#nm->angstrom. *10
#(W*m-2*nm-1)->(erg*s-1*cm-2*Angstrom-1). 10^7*(10^2)^-2*(10)^-1
w_length=data[:,0]
w_length_angstrom=data[:,0]*10.0
flux=data[:,1]
flux_cgs=data[:,1]*10.0**2

#graficar en escala logaritmica
plt.yscale('log')
plt.xscale('log')

plt.plot(w_length_angstrom,flux_cgs,label='Espectro Solar')
plt.xlabel('Longitud de onda (Angstrom)')
plt.ylabel('Espectro ($erg*s^{-1} cm^{-2}*angstrom^{-1}$)')
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
l=calcular_trapecio(w_length,flux,ini,fin)
luminosidad_area=l*10**3 #transformar unidades a cgs
print('luminosidad por unidad de area del Sol (erg*s^-1*cm^-2)= '+str(luminosidad_area))
##calculo de la luminosidad total, multiplicando por la superficie de la esfera con radio 1UA(1.496*10^13 cm)
luminosidad_total=4.0*np.pi*(1.496*10.0**13.0)**2*luminosidad_area
print('luminosidad total del Sol (erg*s^-1)= '+str(luminosidad_total))
print('-------------------------------------------------')

#-----------------------------------------------------------------------
#Integral numerica de la funcion de planck
T=5778#temperatura solar kelvins
planck_sinctes=lambda x:tan(x)**3/(cos(x)**2*(exp(tan(x))-1))
constantes=(2*pi*const.h.cgs/(const.c.cgs)**2) * (const.k_B.cgs*(T*u.K)/const.h.cgs)**4
tolerancia=0.01#notar que el rango de tolerancia debe ser algo entre 0.1 y 0.01. para tolerancias mayores el error es demasiado
planck_numerico=refinar_integral(planck_sinctes,0.0,pi/2.0,tolerancia)*constantes

print('Integral funcion de planck= '+str(planck_numerico))
print('valor analitico integral= '+str((pi**4/15)*constantes))

error=np.fabs(planck_numerico-(pi**4/15)*constantes)
print('error= '+str(error))
print('-------------------------------------------------')

##para un R(radio efectivo del sol), se tiene que cumplir que la integral de la función de planck P, por el area de la esfera debe ser la luminosidad total
##4*pi*R**2*P=l0, donde L0 es la luminosidad del sol calculada anteriormente
##R=(L0/(4*pi*P))**(1/2)
radio_efectivo=(luminosidad_total*(u.erg/u.s)/(4.0*pi*planck_numerico))**(1.0/2.0)
print('radio efectivo del sol= '+str(radio_efectivo))
print('-------------------------------------------------')

#-----------------------------------------------------------------------
#Integracion mediante librerias de scipy
lum_area=trapz(flux,w_length,axis=0)*10**3#factor 10^3 para convertir de W*m^-2 a erg*s^-1*cm^-2
lum_sol=4.0*np.pi*(1.496*10.0**13.0)**2*lum_area
error_2= np.fabs(lum_area-luminosidad_area)
print('luminosidad solar por unidad de area (scipy)= '+str(lum_area))
print('error= '+str(error_2))
print('-------------------------------------------------')


integral_planck=spi.quad(planck_sinctes,0.0,np.pi/2.0)[0]*constantes
error_3=np.fabs(integral_planck-planck_numerico)
print('integral funcion de planck(scipy)= '+str(integral_planck))
print('error= '+str(error_3))
#-----------------------------------------------------------------------
##tiempo de integracion con libreria python:
###L_solar=10000 loops, best of 3: 47.4 microseg por loop
###Planck=100 loops, 2.03 ms por loop

##tiempo de integracion numerica:
###L_solar=
###Planck=1000, 357 microseg
#------------------------
