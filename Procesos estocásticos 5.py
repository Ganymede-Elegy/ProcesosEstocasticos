import numpy as np;
import math;
import random as rnd; 
import matplotlib.pyplot as plt;
from scipy.integrate import odeint; 
import matplotlib; 
matplotlib.rc('xtick', labelsize=18);
matplotlib.rc('ytick', labelsize=18);
matplotlib.rc('axes', labelsize=18);

#Condiciones iniciales
t0 = 0.0; 
r10=1.5;
r20=0.5;

#Parámetros generales
alfa = 0.3;
beta = 10.1;
Omega = 1000;
pasos=100000;

#Parámetros de la forma de Hill
f = 6;
kb = 6.0;
n = 3;

#Función para la forma de Hill
def gr(x):
    gr = (1.0+f*((x/kb)**n))/(1.0+(x/kb)**n);
    return gr;

#Funciones para armar los números aleatorios.
def dist_exp(a):
    r = rnd.random();
    return -(1.0/a)*math.log(r);

def dist_reaccion(ni, A):
    r = rnd.random();
    if(r < ni[0]/A):
        return 0;
    elif(r < (ni[0]+ni[1])/A):
        return 1;
    elif(r < (ni[0]+ni[1] + ni[2])/A):
        return 2;
    elif(r < (ni[0]+ni[1] + ni[2] + ni[3])/A):
        return 3;

#Matriz estequiométrica
S = [[1.0, 0.0, -1.0, 0.0],[0.0, 1.0, 0.0, -1.0]];

#Vectores de cantidad de sustancia y del tiempo
Y = np.zeros([2,pasos+1]);
t = np.zeros(pasos+1)

for i in range(pasos):
    ni = np.array([alfa*gr(Y[1][i]/Omega), alfa*gr(Y[0][i]/Omega), beta*Y[0][i]/Omega, beta*Y[1][i]/Omega]);
    a = sum(ni);
    tau = dist_exp(a);
    mu = dist_reaccion(ni, a);
    Y[0][i+1] = Y[0][i] + S[0][mu];
    Y[1][i+1] = Y[1][i] + S[1][mu];
    t[i+1] = t[i] + tau;
    
plt.plot(Y[0],Y[1]);
plt.show();

#Grafica el sistema determinista (Copia adaptada del código del profe)
def sistema(y,t):
    x1=y[0];
    x2=y[1];
    
    dx1 = alfa*(1.0 + f*np.power(x2,n))/(1.0+np.power(x2,n)) - x1;
    dx2 = alfa*(1.0 + f*np.power(x1,n))/(1.0+np.power(x1,n)) - x2;
    return np.array([dx1, dx2]);

x1i = 1.5;
x2i = 0.5;

tspan = np.linspace(0,600,10000);
ys = odeint(sistema, [x1i,x2i],tspan);
plt.plot(ys[:,0],ys[:,1],'b-');
plt.show();
