import numpy as np;
import math;
import random as rnd; 
import matplotlib.pyplot as plt;

#Condiciones iniciales
t0 = 0.0;
x0 = 0.0;
pasos=150;
Omega = 1000.0 #Tamaño del sistema

#Parámetros para las tasas
alfa = 4.0; 
beta = 0.1; 

#Matriz estequiométrica
S=[1, -1];

Y=np.zeros(pasos+1);
D=np.zeros(pasos+1);

t = np.zeros(pasos+1);

#Le ponemos la condición inicial al vector de partículas.
Y[0] = x0;
D[0] = 0;
t[0] = t0;

def dist_exp(a):
    r = rnd.random();
    return -(1.0/a)*math.log(r);

def dist_reaccion(ni, A):
    r = rnd.random();
    if(r < ni[0]/A):
        return 0;
    if(r < (ni[0]+ni[1])/A):
        return 1;
    
for i in range(pasos):
    ni = np.array([alfa, beta*Y[i]]);
    a = sum(ni);
    tau = dist_exp(a);
    mu = dist_reaccion(ni, a);
    D[i] = (alfa/beta)*(1-math.exp(-beta*t[i]));
    Y[i+1] = Y[i] + S[mu];
    t[i+1] = t[i] + tau;

plt.xlabel("t");
plt.ylabel("Y");
plt.plot(t,Y);
plt.plot(t,D);
plt.show();