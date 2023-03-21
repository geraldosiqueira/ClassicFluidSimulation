# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 15:55:18 2022

@author: Geraldo Siqueira
"""

import numpy as np
import matplotlib.pyplot as plt
import random
#from numba import jit


k = 1
Kb = 1
a = 1
e = 1

N = 63
Lx = 8.81
Lx_real = Lx + a
Ly = np.sqrt(3)*Lx/2
Ly_real = Ly + a
Nt = 60000
delta = 0.1
rho = N/(Lx*Ly)

#jit
def Inicio(N, Nt):
    j = 0
    i = 0
    x = np.zeros( [N, Nt] )
    y = np.zeros( [N, Nt] )
    x0 = np.arange(0, Lx, Lx/(np.sqrt(N) + 1))
    y0 = np.arange(0, Ly+a, Ly/(np.sqrt(N) - 1))
    for n in range(0,N):
        if i%2 == 0: 
            x[n, 0] = x0[i]
            y[n, 0] = Ly/(2*np.sqrt(N)) + y0[j]
            i += 1
        elif i%2 != 0:
            x[n, 0] = x0[i]
            y[n, 0] = y0[j]
            i += 1
        if i == len(x0):
            i = 0
            j += 1
    for n in range(0,N):
        if x[n, 0] > Lx:
            x[n, 0] -= Lx
        if x[n, 0] < 0:
            x[n, 0] += Lx
        if y[n, 0] > Ly:
            y[n, 0] -= Ly
        if y[n, 0] < 0:
            y[n, 0] += Ly   
    return x, y, x0, y0

#@jit
def MCStep(x, y, n, Lx, Ly):
    i = 0 
    dif = 0
    x[:, n+1] = x[:, n]
    y[:, n+1] = y[:, n]
    Part = random.randint(0,N-1)
    varx = random.uniform(-delta, delta)
    vary = random.uniform(-delta, delta)
    x[Part, n+1] += varx
    y[Part, n+1] += vary
    for k in range(N):
        if k != Part:
            dif = (x[Part,n+1]-x[k,n+1])*(x[Part,n+1]-x[k,n+1]) + (y[Part,n+1]-y[k,n+1])*(y[Part,n+1]-y[k,n+1])
            if dif < a*a:
                x[Part, n+1] -= varx
                y[Part, n+1] -= vary
                break

            if x[Part, n+1] + a/2 > Lx:
                dif = ((x[Part,n+1]- Lx) -x[k,n+1])*((x[Part,n+1]- Lx) -x[k,n+1]) + (y[Part,n+1]-y[k,n+1])*(y[Part,n+1]-y[k,n+1])
                if dif < a*a  :
                    x[Part, n+1] -= varx
                    y[Part, n+1] -= vary
                    break
      
            if x[Part, n+1] - a/2 < 0:
                dif = (x[Part,n+1] + Lx -x[k,n+1])*(x[Part,n+1]+ Lx -x[k,n+1]) + (y[Part,n+1]-y[k,n+1])*(y[Part,n+1]-y[k,n+1])
                if dif < a*a  :
                    x[Part, n+1] -= varx
                    y[Part, n+1] -= vary
                    break

            if y[Part, n+1] + a/2 > Ly:
                dif = (x[Part,n+1]-x[k,n+1])*(x[Part,n+1]-x[k,n+1]) + (y[Part,n+1]- Ly -y[k,n+1])*(y[Part,n+1]- Ly -y[k,n+1])
                if dif < a*a  :
                    x[Part, n+1] -= varx
                    y[Part, n+1] -= vary
                    break  

            if y[Part, n+1] - a/2 < 0:
                dif = (x[Part,n+1]-x[k,n+1])*(x[Part,n+1]-x[k,n+1]) + (y[Part,n+1]+ Ly -y[k,n+1])*(y[Part,n+1]+ Ly -y[k,n+1])
                if dif < a*a  :
                    x[Part, n+1] -= varx
                    y[Part, n+1] -= vary
                    break               

    if x[Part, n+1] > Lx:
        x[Part, n+1] -= Lx
    if x[Part, n+1] < 0:
        x[Part, n+1] += Lx
    if y[Part, n+1] > Ly:
        y[Part, n+1] -= Ly
    if y[Part, n+1] < 0:
        y[Part, n+1] += Ly

    if x[Part, n+1] != x[Part, n] or y[Part, n+1] != y[Part, n]:
        i = 1
    return x, y, i

#@jit
def CalculoG(x, y, Nhist, t, Lt):
    
    dx = (2*Lt)/Nhist
    Xh = np.arange(0, 2*Lt, dx)
    Hist = np.zeros(int(Nhist)+1)
    
    for i in range(len(x)):
        for j in range(len(x)):
            if j != i: 
                r = np.sqrt((x[i,t]-x[j,t])*(x[i,t]-x[j,t]) + (y[i,t]-y[j,t])*(y[i,t]-y[j,t]))
                if r < 2*Lt:
                    l = int(r) - 1 
                    Hist[l] += 1
    Hist = Hist/(2)
    return Xh, Hist

aux = 0
x_aux, y_aux, x0, y0= Inicio(N, Nt)
Parametro = np.sqrt(2)*Lx
Nhist = Parametro

for n in range(0,Nt-1):
    x_aux, y_aux, i = MCStep(x_aux, y_aux, n, Lx, Ly)
    aux += i
print("Probabilidade de aceitação:", (aux)/Nt)
print("rho:", (N)/(Lx*Ly))

#Calculando para vários valores de rho

Lx_valores = [8.99, 9.25, 9.53, 10.19, 11.01, 15.57]
x = np.zeros( [N, Nt] )
y = np.zeros( [N, Nt] )
x_auxv = np.zeros(N)
y_auxv = np.zeros(N)
g_aux = np.zeros(len(Lx_valores))
rho_aux = np.zeros(len(Lx_valores))
P = np.zeros(len(Lx_valores))
V = np.zeros(len(Lx_valores))
T = np.zeros(len(Lx_valores))
ratio = np.zeros(len(Lx_valores))
x[:, 0] = x_aux[:, Nt-1]
y[:, 0] = y_aux[:, Nt-1]
k = 0

for Lx in Lx_valores:
    Ly = np.sqrt(3)*Lx/2
    rho1 = N/(Lx*Ly)

    if x_auxv[3] != 0 and y_auxv[3] != 0:
        x[:, 0] = x_auxv[:]
        y[:, 0] = y_auxv[:]

    aux = 0
    Parametro = np.sqrt(2)*Lx
    Nhist = Parametro

    for n in range(0,Nt-1):
        x, y, i = MCStep(x, y, n, Lx, Ly)
        aux += i
    print("Probabilidade de aceitação:", (aux)/Nt)
    print("rho:", (N)/(Lx*Ly))
    rho_aux[k] = rho1 
    V[k] = Lx*Ly
    
    Xh2, Hist2 = CalculoG(x, y, Nhist, 1, Parametro)
    plt.figure(1)
    plt.plot(Xh2, Hist2, label=f"ρ = {rho1}")

    Xh3, Hist3 = CalculoG(x, y, Nhist, Nt-1, Parametro)
    g_aux[k] = Hist3[0]
    plt.figure(2)
    plt.plot(Xh3, Hist3, label="ρ = %.2f" % rho1)
    
    x[:, Nt-1] = x_auxv[:]
    y[:, Nt-1] = y_auxv[:]

    k += 1
    
for j in range (len(g_aux)):
    P[j] = rho_aux[j] + (np.pi/2)*rho_aux[j]*rho_aux[j]*a*a*g_aux[j]
    T[j] = P[j]*V[j]/N    
    ratio[j] = P[j]*V[j]/(N*T[j]) 
    
plt.figure(1)
plt.legend()

plt.figure(2)
plt.legend()

plt.figure(3)
plt.plot(rho_aux, P, 'r-')
plt.title("Pressão em função da densidade")
plt.xlabel("ρ")
plt.ylabel("Pressão Média")

plt.figure(4)
plt.plot(rho_aux, T, 'r-')
#plt.title("Temperatura em função da densidade")
plt.xlabel("ρ")
plt.ylabel("1/T")