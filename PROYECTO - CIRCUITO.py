# -*- coding: utf-8 -*-
"""
Universidad del Valle de Guatemala
Métodos Numéricos
Sección 40
Sergio Alejandro Vasquez Marroquin - 161259
31/10/2023

METODO DE SOLUCIONES NUMERICAS DE EDO - PROYECTO 4

"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def RK41(f1, f2, f3, h, t0, i1_0, i2_0, i3_0, tf):
    t_lista = [t0]
    i1_lista= [i1_0]
    i2_lista = [i2_0]
    i3_lista = [i3_0]

    t = t0
    i1 = i1_0
    i2 = i2_0
    i3 = i3_0

    while t < tf:
        i1_actual = i1
        i2_actual = i2
        i3_actual = i3
        
        k11 = eval(f1)
        k12 = eval(f2)
        k13 = eval(f3)
        
        t = t + h/2
        i1 = i1_actual + k11 * h/2
        i2 = i2_actual + k12 * h/2
        i3 = i3_actual + k13 * h/2

        k21 = eval(f1)
        k22 = eval(f2)
        k23 = eval(f3)
        
        i1 = i1_actual + k21 * h/2
        i2 = i2_actual + k22 * h/2
        i3 = i3_actual + k23 * h/2

        k31 = eval(f1)
        k32 = eval(f2)
        k33 = eval(f3)
        
        t = t + h/2
        
        i1 = i1_actual + k31 * h
        i2 = i2_actual + k32 * h
        i3 = i3_actual + k33 * h

        k41 = eval(f1)
        k42 = eval(f2)
        k43 = eval(f3)

        i1 = i1_actual + (k11 + 2*k21 + 2*k31 + k41) * h/6
        i2 = i2_actual + (k12 + 2*k22 + 2*k32 + k42) * h/6
        i3 = i3_actual + (k13 + 2*k23 + 2*k33 + k43) * h/6

        t_lista.append(t)
        i1_lista.append(i1)
        i2_lista.append(i2)
        i3_lista.append(i3)

    return t_lista, i1_lista, i2_lista, i3_lista

def puntoMedio1(f1, f2, f3, h, t0, i1_0, i2_0, i3_0, tf):
    t_lista = [t0]
    i1_lista= [i1_0]
    i2_lista = [i2_0]
    i3_lista = [i3_0]

    t = t0
    i1 = i1_0
    i2 = i2_0
    i3 = i3_0

    while t < tf:
        i1_actual = i1
        i2_actual = i2
        i3_actual = i3
        
        k11 = eval(f1)
        k12 = eval(f2)
        k13 = eval(f3)
        
        t = t + h/2
        i1 = i1_actual + h / 2 * k11
        i2 = i2_actual + h / 2 * k12
        i3 = i3_actual + h / 2 * k13

        t = t + h

        t_lista.append(t)
        i1_lista.append(i1)
        i2_lista.append(i2)
        i3_lista.append(i3)

    return t_lista, i1_lista, i2_lista, i3_lista

# Nuevos valores de R y L (1 ohmios y 1 H)
R1 = 1  # Ohmios
R2 = 2  # Ohmios
R3 = 3  # Ohmios
L1 = 1  # Henrys
L2 = 2  # Henrys
L3 = 3  # Henrys

f1 = "-(R2/L1)*i2 - (R3/L1)*i3"
f2 = "-(R2/L2 + R2/L1)*i2 + (R1/L2 - R3/L1)*i3"
f3 = "(R1/L3 - R2/L1)*i2 - (R1/L3 + R3/L1 + R3/L3)*i3"

h = 0.1
t0 = 0
i1_0 = 0  # Amperios
i2_0 = 1  # Amperios
i3_0 = 0  # Amperios
tf = 10

# Condiciones iniciales enlistadas
condiciones_iniciales = [0, 1, 0]

# Sistema de ecuaciones
def REAL(X, t):
    ode1 = -(R2/L1) * X[1] - (R3/L1) * X[2]
    ode2 = -(R2/L2 + R2/L1) * X[1] + (R1/L2 - R3/L1) * X[2]
    ode3 = (R1/L3 - R2/L1) * X[1] - (R1/L3 + R3/L1 + R3/L3) * X[2]
    return [ode1, ode2, ode3]

# Tiempo de integración
t = np.linspace(0, 10, 1000)

# Resolver el sistema de ecuaciones diferenciales
solucion = odeint(REAL, condiciones_iniciales, t)


# Solución con Runge-Kutta de orden 4
RK4 = RK41(f1, f2, f3, h, t0, i1_0, i2_0, i3_0, tf)

# Solución con método de punto medio
puntoMedio = puntoMedio1(f1, f2, f3, h, t0, i1_0, i2_0, i3_0, tf)


# Gráfica
plt.figure(figsize=(10, 8))

# Graficas las corrientes reales
plt.plot(t, solucion[:, 0], label='i1 Real')
plt.plot(t, solucion[:, 1], label='i2 Real')
plt.plot(t, solucion[:, 2], label='i3 Real')

# Grafica de las corrientes con el Met. de Runge-Kutta [RK4]
plt.plot(RK4[0], RK4[1], label="i1 RK4")
plt.plot(RK4[0], RK4[2], label="i2 RK4")
plt.plot(RK4[0], RK4[3], label="i3 RK4")

# Grafica de las corrientes con el Met. de RK2 o Punto Medio
plt.plot(puntoMedio[0], puntoMedio[1], label="i1 Pt. Med.")
plt.plot(puntoMedio[0], puntoMedio[2], label="i2 Pt. Med.")
plt.plot(puntoMedio[0], puntoMedio[3], label="i3 Pt. Med.")
plt.xlabel('Tiempo')
plt.ylabel('Corriente')
plt.legend()