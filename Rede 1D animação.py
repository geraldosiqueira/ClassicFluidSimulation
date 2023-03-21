"""
Created on Sun Oct  2 13:12:52 2022

@author: Geraldo Siqueira
"""

import numpy as np
import matplotlib.pyplot as plt
import random
from matplotlib.animation import FuncAnimation
from numba import jit

k = 1
Kb = 1
a = 1

L = 18
Nt = 100000
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
    return x
        

x = Inicio(N, Nt)
for n in range(Nt-1):
    x = MCStep(x, n)
y = np.zeros(N)

dt = 0.01

#Configurando a figura onde será feita a animação
fig = plt.figure(figsize=(L, L))
ax = fig.add_subplot(autoscale_on=False, xlim=(-L, L), ylim=(-L, L))
ax.set_aspect('equal')
    
# Objetos que receberão os dados (configurações) e texto a serem mostrados em cada quadro
confs, = ax.plot([], [], 'ko', ms=30, alpha=0.8)
texto = ax.text(0.05, 0.9, '', transform=ax.transAxes)

# Converte os dados recebidos em cada iteração em objetos a serem mostrados na figura
def animate(n):
    confs.set_data(x[:,200*n], y[:])
    texto.set_text('Tempo = %.2f' % (n*dt))
    return confs, texto

# Constrói a animação
ani = FuncAnimation(fig, animate, frames=Nt//200,interval=50)

plt.title("Sistema de varas rígidas - ρ = %.2f" %rho)
plt.show()
ani.save("teste.gif")