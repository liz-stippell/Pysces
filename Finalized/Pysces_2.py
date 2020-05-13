#!/usr/bin/env python
# coding: utf-8

# In[3]:


from sympy import *
from sympy.physics.quantum import *
from sympy.abc import *
from sympy.plotting import *
from sympy import init_printing
init_printing() 

a, b, c, d, e, f, g, h, h_b, j, k, l, m, n, p, p_x, p_y, p_z, x, y, z, alpha, omega = symbols("a, b, c, d, e, f, g, h, h_b, j, k, l, m, n, p, p_x, p_y, p_z, x, y, z, alpha, omega")


def COMM(A, B, f):
    return expand((Commutator(A, B)*f).doit())

def Hamiltonian(x):
    return Operator(((-(h_b)**2)/(2*m))*((Derivative("1", x, x))))

def OP(x):
    return Operator(x)

def p_OP(x):
    return Operator(-I*h_b*(Derivative("1", x)))

def f(x):
    return Operator(Function('f')(x))

def Prod_Rule_Sp(A, B, C):
    return A*diff_sp(f, C) + B*diff(A, C)

def diff_sp(F, x):
    return Operator(Derivative(Function('f')(x)))

def Lx_Ly():
    return -OP(y)*OP(p_x)*Commutator(OP(z), OP(p_z)) + OP(x)*OP(p_y)*Commutator(OP(p_z), OP(z))

def expand_COMM_1(x, y, y_2):
    return Commutator(x, y) - Commutator(x, y_2)
   
def expand_COMM_2(a, b, c, d):
    return a*d*Commutator(b, c)

def expand_COMM_3(a, b, c, d):
    return b*c*Commutator(a, d)

def HO():
    return (p**2)/(2*m) + (1/2)*k*(x**2)

def Planewave(x):
    return exp(I*(k*x-omega*t))

def PiB(x, L, n):
    return sin((n*pi*x)/L)

def PiB_Norm(x, L, n):
    return sqrt(2/L)*sin((n*pi*x)/L)

def Conj(x):
    return x.replace(I, -I)

def Norm(A, x, y, z):
    return 1/sqrt(Integral(A*Conj(A), (x, y, z)))

def Expect(A, B, x, y, z):
    return Integral(A*Operator(B)*Operator(Conj(A)), (x, y, z))

def Laguerre(r, n):
    return (exp(r)*Derivative(exp(-r)*(r**n), (r, n)))

def Laguerre_Assoc(n, l):
    return ((-1)**l)*Derivative(Laguerre(r, n+l), (r, l))

def Radial(n, l):
    return (r**l)*Laguerre_Assoc(n, l)*exp((-Z/n)*r)


# In[ ]:




