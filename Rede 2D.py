# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 09:25:30 2022

@author: geral
"""

import numpy as np
import matplotlib.pyplot as plt
import random
from matplotlib.animation import FuncAnimation
from numba import jit

k = 1
Kb = 1
a = 1
e = 1

Naux = 70000
Lx = 6
Ly = np.sqrt(3)*Lx/2
Nt = 60000
delta = 0.1
#rho = N/(Lx*Ly)

def Inicio(Lx, Ly):
    N = 0
    x_aux = np.zeros( [Naux, 1] )
    y_aux = np.zeros( [Naux, 1] )
    for i in range(0, Naux-1):
        x_aux[i, 0] = random.uniform(0, Lx)
        y_aux[i, 0] = random.uniform(0, Ly)
        for j in range(Naux):
            if j != i and (x_aux[j, 0] != 0) and (y_aux[j, 0] != 0):
                aux = ((x_aux[i,0]-x_aux[j,0])*(x_aux[i,0]-x_aux[j,0]) + (y_aux[i,0]-y_aux[j,0])*(y_aux[i,0]-y_aux[j,0]))**2 
                if aux < a*a:
                    x_aux[i, 0] = 0
                    y_aux[i, 0] = 0
                    #i -= 1
                    break        
        if (x_aux[i, 0] != 0) or (y_aux[i, 0] != 0):    
            N += 1
    
    x_vetor = []
    y_vetor = []
    for i in range(Naux):
        if (x_aux[i, 0] != 0) or (y_aux[i, 0] != 0):
            x_vetor.append(x_aux[i, 0])
            y_vetor.append(y_aux[i, 0])
    x = np.zeros( [N, Nt] )
    y = np.zeros( [N, Nt] )
    x[:,0] = x_vetor[:]
    y[:,0] = y_vetor[:]
    return x, y, N

def MCStep(x, y, n):
    i = 0 
    dif2 = 0
    x[:, n+1] = x[:, n]
    y[:, n+1] = y[:, n]
    Part = random.randint(0,N-1)
    varx = random.uniform(-delta, delta)
    vary = random.uniform(-delta, delta)
    x[Part, n+1] += varx
    y[Part, n+1] += vary
    for k in range(N):
        if k != Part:
            dif2 = ((x[Part,n+1]-x[k,n+1])*(x[Part,n+1]-x[k,n+1]) + (y[Part,n+1]-y[k,n+1])*(y[Part,n+1]-y[k,n+1]))**2
            if dif2 < a*a:
                x[Part, n+1] -= varx
                y[Part, n+1] -= vary
                break
            '''
            if x[Part, n+1] + a/2 > L:
                d = x[Part,n+1] + a/2 - L
                difup = (-(L+a/2) + d) - x[(Part + 1)%N ,n+1] #(-(L+a/2) + d) simula o centro de outra partícula fora de L
                if difup*difup < a*a  :
                    x[Part, n+1] -= varx
                    break
                elif x[Part, n+1] > L:
                    x[Part, n+1] -= 2*L
            if x[Part,n+1] - a/2 < -L:
                d = x[Part,n+1] - a/2 + L
                difdown = -((L+(a/2)) + d) + x[(Part - 1)%N ,n+1] 
                if difdown*difdown < a*a  :
                    x[Part, n+1] -= varx
                    break
                elif x[Part, n+1] < -L:
                    x[Part, n+1] += 2*L
            '''
    if x[Part, n+1] != x[Part, n] or y[Part, n+1] != y[Part, n]:
        i = 1
    return x, y, i


aux = 0
x, y, N = Inicio(Lx, Ly)
'''
for j in range(N):
    if (x[j, 0] != 0) and (y[j, 0] != 0):
        #plt.plot(x[j, 0], y[j, 0], '.', ms = 90)            
        plt.plot((x[j, 0] - a/2, x[j, 0] + a/2), (y[j,0], y[j,0]), '-')
        plt.plot((x[j, 0] , x[j, 0] ), (y[j,0]- a/2, y[j,0]+ a/2), '-')
'''
for n in range(0,Nt-1):
    x, y, i = MCStep(x, y, n)
    aux += i
print("Probabilidade de aceitação:", (aux)/Nt)

plt.plot(x[:, Nt-1], y[:, Nt-1], '.')
