"""
Created on Mon Oct 10 16:16:38 2022

@author: Geraldo Siqueira
"""

import numpy as np
import matplotlib.pyplot as plt
import random
from numba import jit
from matplotlib.animation import FuncAnimation

k = 1
Kb = 1
a = 1
e = 1

N = 63
Lx = 10.19
Lx_real = Lx + a
Ly = np.sqrt(3)*Lx/2
Ly_real = Ly + a
Nt = 90000
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
    return x, y

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


aux = 0
x, y = Inicio(N, Nt)

for n in range(0,Nt-1):
    x, y, i = MCStep(x, y, n)
    aux += i
print("Probabilidade de aceitação:", (aux)/Nt)

dt = 0.01

#Configurando a figura onde será feita a animação
fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(autoscale_on=False, xlim=(0, Lx), ylim=(0, Ly))
ax.set_aspect('equal')
    
# Objetos que receberão os dados (configurações) e texto a serem mostrados em cada quadro
confs, = ax.plot([], [], 'ko', ms=27, alpha=0.4)
texto = ax.text(0.05, 0.9, '', transform=ax.transAxes)

# Converte os dados recebidos em cada iteração em objetos a serem mostrados na figura
def animate(n):
    confs.set_data(x[:,200*n], y[:,200*n])
    texto.set_text('Tempo = %.2f' % (n*dt))
    return confs, texto

# Constrói a animação
ani = FuncAnimation(fig, animate, frames=Nt//200,interval=50)


plt.show()
#ani.save("animacao fluidos 2d rho070.gif")