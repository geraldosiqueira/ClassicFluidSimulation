# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 15:41:35 2022

@author: Geraldo Siqueira
"""

import numpy as np
import matplotlib.pyplot as plt
import random
from numba import jit

k = 1
Kb = 1
a = 1

L = 18
Nt = 5000000
N = 30
delta = 0.1
rho = N/(2*L)

@jit
def Inicio(N, Nt):
    x = np.zeros( [N, Nt] )
    x0 = np.arange(-L, L, 2*L/N)
    for j in range(0,N):
        x[j,0] = a/2 + x0[j] + random.uniform(0, (2*L/N) - (a))
    return x

@jit
def MCStep(x, n):
    i = 0 
    dif = 0
    x[:, n+1] = x[:, n]
    Part = random.randint(0,N-1)
    varx = random.uniform(-delta, delta)
    x[Part, n+1] += varx
    for k in [-1, 1]:
        dif = x[Part, n+1] - x[(Part + k)%N, n+1]
        if dif*dif < a*a:
            x[Part, n+1] -= varx
            break
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
    if x[Part, n+1] != x[Part, n]:
        i = 1
    return x, i
        

def CalculoG(x, Nhist, t, Lt):
    
    dx = (2*Lt)/Nhist
    Xh = np.arange(0, 2*Lt, dx)
    Hist = np.zeros(Nhist)
    
    for i in range(len(x)):
        for j in range(len(x)):
            if j != i: 
                r = np.sqrt((x[i,t]-x[j,t])*(x[i,t]-x[j,t]))
                if r < 2*Lt:
                    l = int(r) - 1 
                    Hist[l] += 1
    Hist = Hist/(2)
    return Xh, Hist

x = Inicio(N, Nt)
aux = 0
Nhist = 2*L
for n in range(Nt-1):
    x, i = MCStep(x, n)
    aux += i
print(rho)
print(rho + (rho*rho*29))
    
Xh1, Hist1 = CalculoG(x, Nhist, 1, 18)
plt.figure(1)
plt.plot(Xh1, Hist1, '-')

Xh2, Hist2 = CalculoG(x, Nhist, Nt//256, 18)
plt.figure(2)
plt.plot(Xh2, Hist2, '-')

Xh3, Hist3 = CalculoG(x, Nhist, Nt//128, 18)
plt.figure(3)
plt.plot(Xh3, Hist3, '-')

Xh4, Hist4 = CalculoG(x, Nhist, Nt-1, 18)
plt.figure(4)
plt.plot(Xh4, Hist4, '-')
