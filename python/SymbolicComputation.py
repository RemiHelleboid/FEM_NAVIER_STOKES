#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 15:33:25 2020

@author: remi
"""

from sympy import *
import numpy as np


x, y, z = symbols("x, y ,z")

l1 = 1 - x - y
l2 = x
l3 = y 

phi1 = l1 * (2*l1 - 1)
phi2 = l2 * (2*l2 - 1)
phi3 = l3 * (2*l3 - 1)
phi4 = 4 * l1 * l2
phi5 = 4 * l2 * l3
phi6 = 4 * l3 * l1

I1 = integrate(phi1*phi1, (y,0, 1-x), (x, 0, 1))
I2 = integrate(l2*l3, (y,0, 1-x), (x, 0, 1))
I3 = integrate(l3*l1, (y,0, 1-x), (x, 0, 1))


Basis_func = [phi1, phi2, phi3, phi4, phi5, phi6]

mass_matrix = np.zeros((6,6))
V=[]
M = zeros(6,6)

for k in range(6):
    for j in range(6):
        expression = Basis_func[k] * Basis_func[j]
        Integral = integrate(expression, (y,0, 1-x), (x, 0, 1))
        print(Integral)
        mass_matrix[k][j] = Integral
        V.append(Integral)

print(mass_matrix)
print(V)

variables = [x, y]


GradPhi1 = [phi1.diff(var) for var in variables]
GradPhi2 = [phi2.diff(var) for var in variables]
GradPhi3 = [phi3.diff(var) for var in variables]
GradPhi4 = [phi4.diff(var) for var in variables]
GradPhi5 = [phi5.diff(var) for var in variables]
GradPhi6 = [phi6.diff(var) for var in variables]

Grad = [GradPhi1, GradPhi2, GradPhi3, GradPhi4, GradPhi5, GradPhi6]
Div = [np.sum(g) for g in Grad]

print("Div  ", Div)

print("Grad", Grad)

stiff_matrix = np.zeros((6,6))

for k in range(6):
    for j in range(6):
        expression = Grad[k][0]*Grad[j][0] + Grad[k][1]*Grad[j][1]
        # print(expression)
        Integral = integrate(expression, (y,0, 1-x), (x, 0, 1))
        print(Integral)
        stiff_matrix[k, j] = Integral
print(stiff_matrix)



L = [l1, l2, l3]
for k in range(3):
    for j in range(3):
        expression = L[k]*L[j]
        # print(expression)
        Integral = integrate(expression, (y,0, 1-x), (x, 0, 1))
        print(Integral)
        stiff_matrix[k, j] = Integral
print(stiff_matrix)



div_matrix = np.zeros((6,3))


for k in range(6):
    for j in range(3):
        expression = L[j] * Div[k]
        Integral = integrate(expression, (y,0, 1-x), (x, 0, 1))
        print(Integral)
        div_matrix[k,j] = Integral
print(div_matrix)











