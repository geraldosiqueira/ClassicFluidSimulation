"""
Created on Sat Oct  1 10:55:49 2022

@author: Geraldo SIqueira
"""

import numpy as np
import matplotlib.pyplot as plt
import random
from numba import jit

k = 1
Kb = 1
a = 1

L = 18
Nt = 500000
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


x = Inicio(N, Nt)
aux = 0
Nhist = 2*L
for n in range(Nt-1):
    x, i = MCStep(x, n)
    aux += i
print("Probabilidade de aceitação:", (aux)/Nt)

y = np.zeros(N)
for l in range(0, N):
    if x[l, Nt-1] - a/2 > -L and x[l, Nt-1] + a/2 < L:
        plt.plot((x[l, Nt-1] - a/2, x[l, Nt-1] + a/2), (y[l], y[l]), '-')
    elif x[l, Nt-1] - a/2 < -L:
        plt.plot((-L, x[l, Nt-1] + a/2), (y[l], y[l]), '-')
        d = x[l,Nt-1] - a/2 + L
        plt.plot((L + d, L), (y[l], y[l]), '-')
    elif x[l, Nt-1] + a/2 > L:
        plt.plot((x[l, Nt-1] - a/2, L), (y[l], y[l]), '-')
        d = x[l,Nt-1] + a/2 - L
        plt.plot((-L, -L + d), (y[l], y[l]), '-')
plt.xlim(-L, L)
plt.title("Sistema de varas rígidas - ρ = %.2f" %rho)


print(x[:,Nt-1])
for k in range(0,N-1):
    print(x[k,Nt-1] - x[k+1,Nt-1])
