import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from numpy import sqrt
from numpy import exp
from scipy.integrate import quad, dblquad, nquad

from sympy.vector import Vector
from sympy import *
x, y, z, x1, y1, x2, y2, x3, y3 = symbols('x y z x1 y1 x2 y2 x3 y3')

J = np.array([[x2-x1, x3-x1],
              [y2-y1, y3-y1]]).reshape(2,2)

detJ = (x2-x1) * (y3-y1) - (y2-y1) * (x3-x1)
print(detJ)

B =  np.array([[y3-y1, y1-y2],
              [x1-x3, x2-x1]]).reshape(2,2)
BT = B.transpose()

BBT = BT @ B

def P1_basis_func(index):
    if index == 1:
        return(x)
    elif index == 2:
        return(y)
    elif index == 3:
        return(1 - x - y)
    else:
        return "Erreur d'indice"

def P1_basis_func_grad(index):
    if index == 1:
        return(np.array([1, 0]))
    elif index == 2:
        return(np.array([0, 1]))
    elif index == 3:
        return(np.array([-1, -1]))
    else:
        return "Erreur d'indice"

def P2_basis_func(index):
    if index <= 3:
        return P1_basis_func(index) * (2 * P1_basis_func(index) - 1)
    if index == 4:
        return 4 * P1_basis_func(2) * P1_basis_func(3)
    if index == 5:
        return 4 * P1_basis_func(3) * P1_basis_func(1)
    if index == 6:
        return 4 * P1_basis_func(1) * P1_basis_func(2)

def P2_basis_func_grad(index):
    if index <= 3:
        return P1_basis_func_grad(index) * ( 4 * P1_basis_func(index) - 1)
    elif index == 4:
        return 4 * (P1_basis_func(2) * P1_basis_func_grad(3) + P1_basis_func(3) * P1_basis_func_grad(2))
    elif index == 5:
        return 4 * (P1_basis_func(3) * P1_basis_func_grad(1) + P1_basis_func(1) * P1_basis_func_grad(3))
    elif index == 6:
        return 4 * (P1_basis_func(1) * P1_basis_func_grad(2) + P1_basis_func(2) * P1_basis_func_grad(1))


def integrate_ref_triangle_P1_stiff():
    M = []
    for index_1 in range(1, 4):
        L = []
        for index_2 in range(1, 4):
            # print(index_1, index_2)
            L.append(integrate(np.dot(P1_basis_func_grad(index_1), P1_basis_func_grad(index_2)),
                            (y,0,1-x), (x,0,1)))
        M.append(L)
    M = np.array(M)
    M.reshape(3, 3)
    return M

def integrate_triangle_P1_stiff():
    M = []
    for index_1 in range(1, 4):
        L = []
        for index_2 in range(1, 4):
            # print(index_1, index_2, BBT @ P1_basis_func_grad(index_2))
            m = (integrate(np.dot(P1_basis_func_grad(index_1), BBT @ P1_basis_func_grad(index_2)),
                            (y,0,1-x), (x,0,1)))
            L.append(m)
        M.append(L)
    M = np.array(M)

    M.reshape(3, 3)
    return M

def integrate_triangle_P1_stiff_2():
    M = []
    for index_1 in range(1, 4):
        L=[]
        for index_2 in range(1, 4):
            print(index_1, index_2  )
            m = np.dot(P1_basis_func_grad(index_1), BBT @ P1_basis_func_grad(index_2))
            L.append(m)
        M.append(L)
    M = np.array(M)
    M.reshape(3, 3)
    return M

def integrate_ref_triangle_P2_stiff():
    M = []
    for index_1 in range(1, 7):
        L = []
        for index_2 in range(1, 7):
            # print(np.dot(P2_basis_func_grad(index_1), P2_basis_func_grad(index_2)))
            L.append(integrate(np.dot(P2_basis_func_grad(index_1), P2_basis_func_grad(index_2)),
                            (y,0,1-x), (x,0,1)))
        M.append(L)
    M = np.array(M)
    M.reshape(6, 6)
    return M


def integrate_triangle_P2_stiff():
    M = []
    for index_1 in range(1, 7):
        L = []
        for index_2 in range(1, 7):
            # print(np.dot(P2_basis_func_grad(index_1), P2_basis_func_grad(index_2)))
            L.append(integrate(np.dot(P2_basis_func_grad(index_1), BBT @ P2_basis_func_grad(index_2)),
                            (y,0,1-x), (x,0,1)))
        M.append(L)
    M = np.array(M)
    M.reshape(6, 6)
    for k in range(6):
        for j in range(6):
            # print("m", M[k, j])
            M[k, j] = simplify(M[k, j])
    return M

def integrate_triangle_P2_stiff_2():
    M = []
    for index_1 in range(1, 7):
        L = []
        for index_2 in range(1, 7):
            # print(np.dot(P2_basis_func_grad(index_1), P2_basis_func_grad(index_2)))
            L.append(integrate(np.dot( BBT @ P2_basis_func_grad(index_1), P2_basis_func_grad(index_2)),
                            (y,0,1-x), (x,0,1)))
        M.append(L)
    M = np.array(M)
    M.reshape(6, 6)
    for k in range(6):
        for j in range(6):
            # print("m", M[k, j])
            M[k, j] = simplify(M[k, j])
    return M


def integrate_ref_triangle_P1_mass():
    M = []
    for index_1 in range(1, 4):
        L = []
        for index_2 in range(1, 4):
            # print(index_1, index_2)
            L.append(integrate(P1_basis_func(index_1) * P1_basis_func(index_2),
                            (y,0,1-x), (x,0,1)))
        M.append(L)
    M = np.array(M)
    M.reshape(3, 3)
    return M

def integrate_ref_triangle_P2_mass():
    M = []
    for index_1 in range(1, 7):
        L = []
        for index_2 in range(1, 7):
            # print(index_1, index_2)
            L.append(integrate(P2_basis_func(index_1) * P2_basis_func(index_2),
                            (y,0,1-x), (x,0,1)))
        M.append(L)
    M = np.array(M)
    M.reshape(6, 6)
    return(M)


MP1 = integrate_ref_triangle_P1_mass()
print("Mass Matrix P1 \n", MP1, "\n")

MP2 = integrate_ref_triangle_P2_mass()
print("Mass Matrix P2 \n" ,MP2, "\n")

# SP1 = integrate_ref_triangle_P1_stiff()
# print("Stiff Matrix  \n", SP1, "\n")
#
# SPT1 = integrate_triangle_P1_stiff()
# SPT11 = integrate_triangle_P1_stiff_2()
# print("Stiff Matrix Triangle quelconque P1 \n", SPT1-SPT11, "\n")

SP2 = integrate_ref_triangle_P2_stiff()
print("Stiff Matrix P2 \n", SP2, "\n")

SPT2 = integrate_triangle_P2_stiff()
SPT22 = integrate_triangle_P2_stiff_2()
print("Stiff Matrix Triangle quelconque P2 \n", SPT2, "\n")
