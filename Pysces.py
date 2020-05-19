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

""" The above are a list of commonly used variables for the Pysces library. """


def COMM_1(A, B, f):
    return expand((Commutator(A, B)*f).doit())
"""

This function is not used directly, it is only used in the COMM() function below. 

Parameters:

A: The first operator in the commutator
B: The second operator in the commutator
f: The auxiliary function. This is defined below, as F(x).

Returns:

The commutator of A and B with respect to F(x).

"""



def COMM(x, y, a, b = None):
            if b == None:
                return COMM_1(OP(x), OP(y), a)
            elif b == 0:
                return Commutator(OP(x), OP(a)) + Commutator(OP(y), OP(a))
            else:       
                return Commutator(OP(x), OP(a)) - Commutator(OP(x), OP(b)) + Commutator(OP(y), OP(a)) + Commutator(OP(y), OP(b))
"""

This function has a few different outputs, depending on the parameters.

Parameters:
x: the first term of the first operator
y: the second term of the firts operator
a: the first term of the second operator
b: the second term of the second operator

The commutator looks similar to this: [x + y, a + b]
An example would be the angular momentum operators, [L_x, L_y], or any similar combination.

Returns:

If b == None (this means there is no b parameter listed), the output will be COMM_1(), as defined above.
If b == 0 (0 would be written in the argument), the output will be [x, a] + [y, a], such as used in the case of [L_x, x]
Else, [x, a] - [x, b] + [y, a] + [y, b] would be returned, such as the case of [L_x, L_y]

"""



def FACTOR_COMM(A, B, C, D):
    if A == p_y and C == y:
        return (OP(B)*OP(D)*Commutator(OP(A), OP(C)) )
    if A == y and C == p_y:
        return (OP(B)*OP(D)*Commutator(OP(A), OP(C)) )
    if A == p_x and C == x:
        return (OP(B)*OP(D)*Commutator(OP(A), OP(C)) )
    if A == x and C == p_x:
        return( OP(B)*OP(D)*Commutator(OP(A), OP(C)) )
    if A == p_z and C == z:
        return( OP(B)*OP(D)*Commutator(OP(A), OP(C)) )
    if A == z and C == p_z:
        return( OP(B)*OP(D)*Commutator(OP(A), OP(C)) )
        
    if A == p_y and D == y:
        return( OP(B)*OP(C)*Commutator(OP(A), OP(D)) )
    if A == y and D == p_y:
        return( OP(B)*OP(C)*Commutator(OP(A), OP(D)) )
    if A == p_x and D == x:
        return( OP(B)*OP(C)*Commutator(OP(A), OP(D)) )
    if A == x and D == p_x:
        return( OP(B)*OP(C)*Commutator(OP(A), OP(D)) )
    if A == p_z and D == z:
        return( OP(B)*OP(C)*Commutator(OP(A), OP(D)) )
    if A == z and D == p_z:
        return( OP(B)*OP(C)*Commutator(OP(A), OP(D)) )
   
    if B == p_y and C == y:
         return( OP(A)*OP(D)*Commutator(OP(B), OP(C)) )
    if B == y and C == p_y:
         return( OP(A)*OP(D)*Commutator(OP(B), OP(C)) )
    if B == p_x and C == x:
         return( OP(A)*OP(D)*Commutator(OP(B), OP(C)) )
    if B == x and C == p_x:
         return( OP(A)*OP(D)*Commutator(OP(B), OP(C)) )
    if B == p_z and C == z:
         return( OP(A)*OP(D)*Commutator(OP(B), OP(C)) )
    if B == z and C == p_z:
          return( OP(A)*OP(D)*Commutator(OP(B), OP(C)) )
    
    if B == p_y and D == y:
        return(OP(A)*OP(C)*Commutator(OP(B), OP(D)))
    if B == y and D == p_y:
        return(OP(A)*OP(C)*Commutator(OP(B), OP(D)))
    if B == p_x and D == x:
        return(OP(A)*OP(C)*Commutator(OP(B), OP(D)))
    if B == x and D == p_x:
        return(OP(A)*OP(C)*Commutator(OP(B), OP(D)))
    if B == p_z and D == z:
         return(OP(A)*OP(C)*Commutator(OP(B), OP(D)))
    if B == z and D == p_z:
          return(OP(A)*OP(C)*Commutator(OP(B), OP(D)))

    if A == y and C == y:
        return 0
    if A == p_y and C == p_y:
        return 0
    if A == x and C == x:
        return 0
    if A == p_x and C == p_x:
        return 0
    if A == z and C == z:
        return 0
    if A == p_z and C == p_z: 
        return 0
    
    if A == y and D == y:
        return 0
    if A == p_y and D == p_y:
         return 0
    if A == x and D == x:
         return 0
    if A == p_x and D == p_x:
         return 0
    if A == z and D == z:
        return 0
    if A == p_z and D == p_z:
         return 0
        
    if B == y and C == y:
         return 0
    if B == p_y and C == p_y:
         return 0
    if B == x and C == x:
         return 0
    if B == p_x and C == p_x:
         return 0
    if B == z and C == z:
         return 0
    if B == p_z and C == p_z:
         return 0
    
    if B == y and D == y:
        return 0
    if B == p_y and D == p_y:
        return 0
    if B == x and D == x:
        return 0
    if B == p_x and D == p_x:
        return 0
    if B == z and D == z:
        return 0
    if B == p_z and D == p_z:
        return 0
"""

Parameters:

A: The first term of the first operator
B: The second term of the first operator
C: The first term of the second operator
D: The second term of the second operator

Returns:

This function factors out terms based on the criteria set by the "if" statements. If a term in the first operator matches a term in the second operator, the commutator will be reported as zero ("0").
The criterion set examines the variables of each operator, and factors out operators that have "mismatched" variables (such as p with respect to z and y), while keeping the terms that are with respect to the same variable within a commutator, such as [p_x, x].

Note that if there are multiplte commutators, the FACTOR_COMM() functions can be added or subtracted ("+" or "-") in a single line for ease.


"""



def OP(x):
    return Operator(x)

"""

Parameters:

x: what the position operator is with respect to.

Returns:

The chosen parameter as an operator. Example: OP(y) will return "y"

"""



def P_OP(x):
    return Operator(-I*h_b*(Derivative("1", x)))
"""

Parameters:

x: what the linear momentum operator is with respect to.

Returns:

The linear momentum operator, with respect to the parameter.

Note that this equation is included in the Hamiltonian() function. 
Note that the "1" in the derivative is a placeholder, which will be manually replaced.

"""



def HAMILTONIAN(x):
    return P_OP(x)*(1/(2*m))
"""

This is subject to change.

Parameters:

x: The variable in which the derivative is with respect to.

Returns:

The Hamiltonian for a Particle in a Box, with respect to the chosen parameter. 

Note that (-(h_b)**2)*Derivative("1", x, x) == P_OP(x)
Note that the "1" in the derivative is a placeholder, which will be manually replaced.

"""



def EXPRESSION_REPLACE(R, K):
    if K == x:
        return sympify(str(R).replace(str(Derivative(1, K)*F(K)), str(Derivative(F(K), K).doit())).replace(str('Derivative(1, x)*x*F(x)'), str(Derivative(K*F(K), K).doit())))
    if K == y:
        return sympify(str(R).replace(str(Derivative(1, K)*F(K)), str(Derivative(F(K), K).doit())).replace(str('Derivative(1, y)*y*F(y)'), str(Derivative(K*F(K), K).doit())))
    if K == z:
        return sympify(str(R).replace(str(Derivative(1, K)*F(K)), str(Derivative(F(K), K).doit())).replace(str('Derivative(1, z)*z*F(z)'), str(Derivative(K*F(K), K).doit())))
"""

Parameters:

R: The expanded commutator to be replaced
K: The variable/parameter with respect to the chosen commutator.

Returns:

This replaces the Derivative(1, x) present in the expanded commutator with either Derivative(F(x), x) or Derivative(x*F(x), x).
Note that the above states "x", but can be done for x, y, or z variables.
Note that the R parameter is a str()


"""



def F(x):
    return Operator(Function('F')(x))
"""

Parameters:

x: This is what the auxiliary function is with respect to. It should also match the parameters of the other arguments in the COMM() function. (example: if one operator is P_OP(y), the auxiliary function should be F(y).

Returns:

This is the auxiliary function, commonly used in the COMM() function.
This simply returns "F(x)", x being with chosen parameter.

"""



def HO():
    return (p**2)/(2*m) + (1/2)*k*(x**2)
"""

Parameters:

Currently, there are no needed parameters for this equation.

Returns:

The Harmonic Oscillator Hamiltonian. Note that "p" is the linear momentum operator, and so the first term is the Hamiltonian() mentioned above.

"""



def PLANEWAVE(x):
    return exp(I*(k*x-omega*t))
"""

Parameters:

x: What the function is with respect to.

Returns:

The PlaneWave function with respect to the parameter.

"""



def PIB(x, L, n):
    return sin((n*pi*x)/L)
"""

Parameters:

x: a variable.
L: Length of the box.
n: an integer.

Returns:

The WaveFunction for Particle in a Box with respect to the chosen parameters.

"""



def PIB_NORMALIZE(x, L, n):
    return sqrt(2/L)*sin((n*pi*x)/L)
"""

Parameters:

x: a variable.
L: Length of the box.
n: an integer.

Returns:

The normalized WaveFunction for Particle in a Box, with respect to the chosen variables. This answer can also be calculated via the NORMALIZE() function.

"""



def CONJUGATE(x):
    return x.replace(I, -I)
"""

Parameters:

x: The term of interest. This is commonly a WaveFunction.

Returns:

The complex conjugate of a WaveFunction. If there are imaginary terms, they get negated, and if there are no imaginary terms, the WaveFunction is not affected.

"""



def NORMALIZE(A, x, y, z):
    return 1/sqrt(Integral(A*CONJUGATE(A), (x, y, z)))
"""

Parameters:

A: The WaveFunction/expression of interest
x: What the integral is taken with respect to
y: The lower bound of the integral. If bounds are not listed, this is -oo
z: The upper bound of the integral. If bounds are not listed, this is oo

Returns:

The normalization constant, with respect to the given parameters. To find the normalized WaveFunction, the output for NORMALIZE() must be multiplied by the original WaveFunction. A normalized WaveFunction indicates that the probability of finding a particle within certain bounds must be equal to one.

"""



def EXPECTATION(A, B, x, y, z):
    return Integral(A*B*CONJUGATE(A), (x, y, z))
"""

Parameters:

A: The normalized WaveFunction/expression of interest
B: The operator of interest.
x: What the integral is taken with respect to
y: The lower bound of the integral. If bounds are not listed, this is -oo
z: The upper bound of the integral. If bounds are not listed, this is oo

Returns:

The expectation value for the given operator and normalized WaveFunction. An expectation value is the average value of an operator for a given WaveFunction. 

"""



def PLOT(A, B, C, D):
    return plot(sympify(str(A)), (B, C, D))
"""

Parameters:

A: The function/Normalized WaveFunction of interest
B: This is "x" usually (x-axis)
C: The lower bound of the x-axis domain (for Particle in a Box, 0)
D: The upper bound of the x-axis domain (for Particle in a Box, 1)

Returns:

A plotted function of the function/Normalized WaveFunction of interest.

Note that the cell often has to be run TWICE in order to print the plot/graph.


"""



def LAGUERRE(r, n):
    return simplify((exp(r)*Derivative(exp(-r)*(r**n), (r, n))).doit())
"""

Parameters:

r: What the equation is with respect to. This is commonly "r" or "x"
n: The principle Quantum Number

Returns:

The Laguerre polynomial. This is commonly used to solve the LAGUERRE_ASSOC() function.

"""

def LAGUERRE_2(r, n):
    return (exp(r)*Derivative(exp(-r)*(r**n), (r, n)))

def LAGUERRE_ASSOC(n, l):
    return (((-1)**l)*Derivative(LAGUERRE_2(r, n+l), (r, l)).doit())
"""

Parameters:

n: The principle Quantum Number
l: The angular momentum Quantum Number

Returns:

The Laguerre Associated polynomial, commonly used to solve the RADIAL() function.

"""




# In[ ]:




