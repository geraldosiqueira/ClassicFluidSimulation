# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 15:15:09 2022

@author: Pichau
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
#rho = N/(Lx*Ly)


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

def MCStep(x, y, n):
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
                '''
                if y[Part, n+1] + a/2 > Ly:
                    dif = (x[Part,n+1]- Lx -x[k,n+1])*(x[Part,n+1]- Lx -x[k,n+1]) + (y[Part,n+1]- Ly -y[k,n+1])*(y[Part,n+1]- Ly -y[k,n+1])
                    dif1 = (x[Part,n+1]- Lx -x[k,n+1])*(x[Part,n+1]- Lx -x[k,n+1]) + (y[Part,n+1]-y[k,n+1])*(y[Part,n+1]-y[k,n+1])
                    if dif < a*a or dif1 < a*a :
                        x[Part, n+1] -= varx
                        y[Part, n+1] -= vary
                        break
                elif y[Part, n+1] - a/2 < 0:
                    dif = (x[Part,n+1]- Lx -x[k,n+1])*(x[Part,n+1]- Lx -x[k,n+1]) + (y[Part,n+1]+ Ly -y[k,n+1])*(y[Part,n+1]+ Ly -y[k,n+1])
                    dif1 = (x[Part,n+1]- Lx -x[k,n+1])*(x[Part,n+1]- Lx -x[k,n+1]) + (y[Part,n+1]-y[k,n+1])*(y[Part,n+1]-y[k,n+1])
                    if dif < a*a or dif1 < a*a :
                        x[Part, n+1] -= varx
                        y[Part, n+1] -= vary
                        break
                else:
                    '''
                dif = ((x[Part,n+1]- Lx) -x[k,n+1])*((x[Part,n+1]- Lx) -x[k,n+1]) + (y[Part,n+1]-y[k,n+1])*(y[Part,n+1]-y[k,n+1])
                if dif < a*a  :
                    x[Part, n+1] -= varx
                    y[Part, n+1] -= vary
                    break
      
            if x[Part, n+1] - a/2 < 0:
                '''
                if y[Part, n+1] + a/2 > Ly:
                    dif = (x[Part,n+1]+ Lx -x[k,n+1])*(x[Part,n+1]+ Lx -x[k,n+1]) + (y[Part,n+1]- Ly -y[k,n+1])*(y[Part,n+1]- Ly -y[k,n+1])
                    dif1 = (x[Part,n+1]+ Lx -x[k,n+1])*(x[Part,n+1]+ Lx -x[k,n+1]) + (y[Part,n+1]-y[k,n+1])*(y[Part,n+1]-y[k,n+1])
                    if dif < a*a or dif1 < a*a :
                        x[Part, n+1] -= varx
                        y[Part, n+1] -= vary
                        break
                elif y[Part, n+1] - a/2 < 0:
                    dif = (x[Part,n+1]+ Lx -x[k,n+1])*(x[Part,n+1]+ Lx -x[k,n+1]) + (y[Part,n+1]+ Ly -y[k,n+1])*(y[Part,n+1]+ Ly -y[k,n+1])
                    dif1 = (x[Part,n+1]+ Lx -x[k,n+1])*(x[Part,n+1]+ Lx -x[k,n+1]) + (y[Part,n+1]-y[k,n+1])*(y[Part,n+1]-y[k,n+1])
                    if dif < a*a or dif1 < a*a :
                        x[Part, n+1] -= varx
                        y[Part, n+1] -= vary
                        break
                else:
                    '''
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

aux = 0
x, y, x0, y0= Inicio(N, Nt)

for n in range(0,Nt-1):
    x, y, i = MCStep(x, y, n)
    aux += i
print("Probabilidade de aceitação:", (aux)/Nt)
print((N+1)/(Lx*Ly))
'''
for j in range(N):
    if (x[j, 0] != 0) and (y[j, 0] != 0):
        #plt.plot(x[j, 0], y[j, 0], '.', ms = 90)            
        plt.plot((x[j, Nt-1] - a/2, x[j, Nt-1] + a/2), (y[j,Nt-1], y[j,Nt-1]), '-')
        plt.plot((x[j, Nt-1] , x[j, Nt-1] ), (y[j,Nt-1]- a/2, y[j,Nt-1]+ a/2), '-')


'''
plt.figure(1)
plt.plot(x[:, 0], y[:, 0], '.')
plt.plot((0 , Lx), (0, 0), '-')
plt.plot((0 , Lx), (Ly, Ly), '-')
plt.plot((0 , 0), (0, Ly), '-')
plt.plot((Lx , Lx), (0, Ly), '-')

plt.figure(2)
plt.plot(x[:, Nt-1], y[:, Nt-1], '.')
plt.plot((0 , Lx), (0, 0), '-')
plt.plot((0 , Lx), (Ly, Ly), '-')
plt.plot((0 , 0), (0, Ly), '-')
plt.plot((Lx , Lx), (0, Ly), '-')
