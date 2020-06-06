#!/usr/bin/env python
# coding: utf-8

# In[3]:


from sympy import *
from sympy.physics.quantum import *
from sympy.abc import *
from sympy.plotting import *
from sympy import init_printing
init_printing() 


I = sqrt(-1)


def comm_1(A, B, f):
    """

    This function is not used directly, it is only used in the COMM() function below. 

    Parameters:

    A: The first operator in the commutator
    B: The second operator in the commutator
    f: The auxiliary function. This is defined below, as F(x).

    Returns:

    The commutator of A and B with respect to F(x).

    """
    
    return expand((Commutator(Operator(A), Operator(B))*f).doit())




def comm(x, y, a, b = None):
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
    
    if b == None:
        return comm_1(Operator(x), Operator(y), a)
    elif b == 0:
        return Commutator(Operator(x), Operator(a)) + Commutator(Operator(y), Operator(a))
    else:       
        return Commutator(Operator(x), Operator(a)) - Commutator(Operator(x), Operator(b)) + Commutator(Operator(y), Operator(a)) + Commutator(Operator(y), Operator(b))




def factor_comm(A, B, C, D):
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
    
    x, y, z, p_x, p_y, p_z = symbols("x y z p_x p_y p_z")
    if A == p_y and C == y:
        return (Operator(B)*Operator(D)*Commutator(Operator(A), Operator(C)) )
    if A == y and C == p_y:
        return (Operator(B)*Operator(D)*Commutator(Operator(A), Operator(C)) )
    if A == p_x and C == x:
        return (Operator(B)*Operator(D)*Commutator(Operator(A), Operator(C)) )
    if A == x and C == p_x:
        return (Operator(B)*Operator(D)*Commutator(Operator(A), Operator(C)) )
    if A == p_z and C == z:
        return (Operator(B)*Operator(D)*Commutator(Operator(A), Operator(C)) )
    if A == z and C == p_z:
        return (Operator(B)*Operator(D)*Commutator(Operator(A), Operator(C)) )
        
    if A == p_y and D == y:
        return( Operator(B)*Operator(C)*Commutator(Operator(A), Operator(D)) )
    if A == y and D == p_y:
        return( Operator(B)*Operator(C)*Commutator(Operator(A), Operator(D)) )
    if A == p_x and D == x:
        return( Operator(B)*Operator(C)*Commutator(Operator(A), Operator(D)) )
    if A == x and D == p_x:
        return( Operator(B)*Operator(C)*Commutator(Operator(A), Operator(D)) )
    if A == p_z and D == z:
        return( Operator(B)*Operator(C)*Commutator(Operator(A), Operator(D)) )
    if A == z and D == p_z:
        return( Operator(B)*Operator(C)*Commutator(Operator(A), Operator(D)) )
   
    if B == p_y and C == y:
         return( Operator(A)*Operator(D)*Commutator(Operator(B), Operator(C)) )
    if B == y and C == p_y:
         return( Operator(A)*Operator(D)*Commutator(Operator(B), Operator(C)) )
    if B == p_x and C == x:
         return( Operator(A)*Operator(D)*Commutator(Operator(B), Operator(C)) )
    if B == x and C == p_x:
         return( Operator(A)*Operator(D)*Commutator(Operator(B), Operator(C)) )
    if B == p_z and C == z:
         return( Operator(A)*Operator(D)*Commutator(Operator(B), Operator(C)) )
    if B == z and C == p_z:
         return( Operator(A)*Operator(D)*Commutator(Operator(B), Operator(C)) )
    
    if B == p_y and D == y:
        return(Operator(A)*Operator(C)*Commutator(Operator(B), Operator(D)))
    if B == y and D == p_y:
        return(Operator(A)*Operator(C)*Commutator(Operator(B), Operator(D)))
    if B == p_x and D == x:
        return(Operator(A)*Operator(C)*Commutator(Operator(B), Operator(D)))
    if B == x and D == p_x:
        return(Operator(A)*Operator(C)*Commutator(Operator(B), Operator(D)))
    if B == p_z and D == z:
        return(Operator(A)*Operator(C)*Commutator(Operator(B), Operator(D)))
    if B == z and D == p_z:
         return(Operator(A)*Operator(C)*Commutator(Operator(B), Operator(D)))

    if A == C:
        return 0
    if A == D:
        return 0
    if B == C:
         return 0
    if B == D:
        return 0




def p_operator(A = None):
    """

    Parameters:

    A: what the linear momentum operator is with respect to.

    Returns:

    The linear momentum operator, with respect to the parameter.
 
    Note that the "1" in the derivative is a placeholder, which will be replaced.
    If A == None, the general linear operator is printed, with all three positional arguments "x", "y", and "z"

    """
    
    h_b, m = symbols("h_b m")
    if A == None:
        return Operator(-I*h_b*(Derivative("1", x))) + Operator(-I*h_b*(Derivative("1", y))) + Operator(-I*h_b*(Derivative("1", z)))
    else:
        return Operator(-I*h_b*(Derivative("1", A)))




def kinetic_energy(A = None):
    """

    Parameters:

    x: The variable in which the derivative is with respect to.

    Returns:

    The kinetic energy operator, with respect to the chosen parameter. 

    Note that the "1" in the derivative is a placeholder, which will be replaced.
    If A == None, the general linear operator is printed, with all three positional arguments "x", "y", and "z"

    """
    
    h_b, m = symbols("h_b m")
    if A == None:
        return (Operator((-I*h_b*(Derivative("1", x)))**2)*(1/(2*m))) + (Operator((-I*h_b*(Derivative("1", y)))**2)*(1/(2*m))) + (Operator((-I*h_b*(Derivative("1", z)))**2)*(1/(2*m)))
    else:
        return (Operator((-I*h_b*(Derivative("1", A)))**2)*(1/(2*m)))



def v(x):
    """
    
    Parameters:
    
    x: The variable for the given function. This is usually "x", "y" or "z".
    
    Returns:
    
    The general potential energy operator, V(x)
    
    """
    
    return Operator(Function("v")(x))



def hamiltonian(A):
    """
    
    Parameters:
    
    A: The variable for the given function. This is usually "x", "y" or "z".
    
    Returns:
    
    The Hamiltonian operator, made up of the kinetic energy operator and the general potential energy operator.
    
    """
    
    return kinetic_energy(A) + v(A)



def expression_replace(expr, K):
    """

    Parameters:

    expr: The expanded commutator to be replaced
    K: The variable/parameter with respect to the chosen commutator.

    Returns:

    This replaces the Derivative(1, x) present in the expanded commutator with either Derivative(F(x), x) or Derivative(x*F(x), x).
    Note that the above states "x", but can be done for x, y, or z variables.
    This will also replace [p_x, x] and similar commutators, which are usually computed when using angular momentum operators or similar operators.
    Note that the "expr" parameter is a str()


    """
    
    L_x, L_z, L_y, p_x, p_y, p_z, x, y, z = symbols("L_x, L_z, L_y, p_x, p_y, p_z, x, y, z")
    
    if str('[p_x,x]') in str(expr):
        return sympify(str(sympify(str(expr).replace(str('[p_x,x]'), str(expression_replace(comm(p_operator(K), Operator(K), f(K)), K))))).replace(str(p_y*z - p_z*y), str(-L_x)))
    if K == x:
        return sympify(str(expr).replace(str(Derivative(1, K)*f(K)), str(Derivative(f(K), K).doit())).replace(str('Derivative(1, x)*x*f(x)'), str(Derivative(K*f(K), K).doit())))
    
    if str('[p_y, y]') in str(expr):
        return sympify(str(sympify(str(R).replace(str('[p_y,y]'), str(expression_replace(comm(p_operator(K), Operator(K), f(K)), K))))).replace(str(p_x*z - p_z*x), str(-L_y)))
    if K == y:
        return sympify(str(expr).replace(str(Derivative(1, K)*f(K)), str(Derivative(f(K), K).doit())).replace(str('Derivative(1, y)*y*f(y)'), str(Derivative(K*f(K), K).doit())))
    
    if str('[p_z,z]') in str(expr):
        return sympify(str(sympify(str(expr).replace(str('[p_z,z]'), str(expression_replace(comm(p_operator(K), Operator(K), f(K)), K))))).replace(str(p_x*y - p_y*x), str(-L_z)))
    elif K == z:
        return sympify(str(expr).replace(str(Derivative(1, K)*f(K)), str(Derivative(f(K), K).doit())).replace(str('Derivative(1, z)*z*f(z)'), str(Derivative(K*f(K), K).doit())))





def f(x):
    """

    Parameters:

    x: This is what the auxiliary function is with respect to. It should also match the parameters of the other arguments in the COMM() function. (example: if one operator is P_OP(y), the auxiliary function should be F(y).

    Returns:

    This is the auxiliary function, commonly used in the COMM() function.
    This simply returns "F(x)", x being with chosen parameter.

    """
    
    return Operator(Function('f')(x))




def HO():
    """

    Parameters:

    Currently, there are no needed parameters for this equation.

    Returns:

    The Harmonic Oscillator Hamiltonian. Note that "p" is the linear momentum operator, and so the first term is the KINETIC_ENERGY() function mentioned above.

    """
    
    p, k, m, x = symbols("p k m x")
    return (p**2)/(2*m) + (1/2)*k*(x**2)




def planewave(x):
    """

    Parameters:

    x: What the function is with respect to.

    Returns:

    The PlaneWave function with respect to the parameter.

    """

    k, omega, t = symbols("k omega t")
    return exp(I*(k*x-omega*t))



def PIB(x, L, n):
    """

    Parameters:

    x: a variable.
    L: Length of the box.
    n: an integer.

    Returns:

    The WaveFunction for Particle in a Box with respect to the chosen parameters.

    """
    
    return sin((n*pi*x)/L)




def PIB_normalize(x, L, n):
    """

    Parameters:

    x: a variable.
    L: Length of the box.
    n: an integer.

    Returns:

    The normalized WaveFunction for Particle in a Box, with respect to the chosen variables. This answer can also be calculated via the NORMALIZE() function.

    """
    
    return sqrt(2/L)*sin((n*pi*x)/L)




def gaussian(alpha, x_0, p):
    """
    
    Parameters:
    
    alpha:alpha parameter of the gaussian
    x_0: x_0 parameter of the gaussian
    p: p parameter of the gaussian
    
    Returns:
    
    The moving gaussian wave function.
    
    """
    
    alpha, x, x_0, p, h_b = symbols("alpha x x_0 p h_b")
    return exp(-((alpha/2)*(x-x_0)**2 + (I*p)/h_b*(x-x_0)))




def gaussian_normalize(alpha, gamma, x_0, p_0):
    """
    
    Parameters:
    
    alpha: alpha parameter of the normalized gaussian
    gamma: gamma parameter of the normalized gaussian
    x_0: x_0 parameter of the normalized gaussian
    p_0: p_0 parameter of the normalized gaussian
    
    Returns:
    
    The normalized gaussian WaveFunction, with respect to the chosen variables.
    
    """
    
    alpha, gamma, x, x_0, p_0, h_b = symbols("alpha gamma x x_0 p_0 h_b")
    return ((2*alpha)/pi)**(1/4) * exp(-alpha*(x - x_0)**2 + ((I*p_0)/h_b)*(x - x_0) + (I*gamma)/h_b)




def conjugate(expr):
    """

    Parameters:

    x: The term of interest. This is commonly a WaveFunction.

    Returns:

    The complex conjugate of a WaveFunction. If there are imaginary terms, they get negated, and if there are no imaginary terms, the WaveFunction is not affected.

    """
    
    return expr.replace(I, -I)




def normalize_constant(WaveFunc, var, lower, upper):
    """

    Parameters:

    WaveFunc: The WaveFunction/expression of interest
    var: What the integral is taken with respect to
    lower: The lower bound of the integral. If bounds are not listed, this is -oo
    upper: The upper bound of the integral. If bounds are not listed, this is oo

    Returns:

    The normalization constant, with respect to the given parameters. To find the normalized WaveFunction, the output for NORMALIZE() must be multiplied by the original WaveFunction. A normalized WaveFunction indicates that the probability of finding a particle within certain bounds must be equal to one.

    """
    
    return 1/sqrt(Integral(WaveFunc*conjugate(WaveFunc), (var, lower, upper)))




def expectation_value(WaveFunc_1, Operator, WaveFunc_2, var, lower, upper):
    """

    Parameters:
    
    WaveFunc_1: The "bra" normalized WaveFunction
    Operator: The operator of interest
    WaveFunc_2: The "ket" normalized WaveFunction
    var: What the integral is taken with respect to
    lower: The lower bound of the integral. If bounds are not listed, this is -oo
    upper: The upper bound of the integral. If bounds are not listed, this is oo
    
    For reference:
    <bra|
    |ket>
    so:
    <WaveFunc_1|Operator|WaveFunc_2>

    Returns:

    The expectation value for the given operator and normalized WaveFunction. An expectation value is the average value of an operator for a given WaveFunction. 

    """
    
    
    if B == kinetic_energy(x):
        return sympify(str(sympify(str(Integral(conjugate(A)*B, (x, y, z))).replace(str(Derivative("1", x)**2), str(Derivative(A, x, x)))).doit()).replace(str('sin(pi*n)'), str(0)).replace(str('cos(pi*n)'), str(0)))
    if B == p_operator(x):
        return Integral(conjugate(A)*B, (x, y, z)).replace(Derivative("1", x), Derivative(A, x).doit()) 
    else:
        return Integral(conjugate(WaveFunc_1)*Operator*WaveFunc_2, (var, lower, upper))



def overlap(WaveFunc_1, WaveFunc_2, var, lower, upper):
    """
    
    Parameters:
    
    WaveFunc_1: The "bra" normalized WaveFunction
    WaveFunc_2: The "ket" normalized WaveFunction
    var: What the integral is taken with respect to
    lower: The lower bound of the integral. If bounds are not listed, this is -oo
    upper: The upper bound of the integral. If bounds are not listed, this is oo 
    
    For reference:
    <bra|
    |ket>
    so:
    <WaveFunc_1|WaveFunc_2>
    
    Returns:
    
    The overlap of the two WaveFunctions of interest over given bounds. Note that if these are the same WaveFunctions, the overlap is 1.
    
    """
    
    return Integral(WaveFunc_1*WaveFunc_2, (var, lower, upper))




def plot_function(func, B, lower, upper):
    """

    Parameters:

    func: The function/Normalized WaveFunction of interest
    B: This is "x" usually (x-axis)
    lower: The lower bound of the x-axis domain (for Particle in a Box, 0)
    upper: The upper bound of the x-axis domain (for Particle in a Box, 1)

    Returns:

    A plotted function of the function/Normalized WaveFunction of interest.

    Note that the cell often has to be run TWICE in order to print the plot/graph.


    """
        
    return plot(sympify(str(func)), (B, lower, upper))




def laguerre(r, n):
    """

    Parameters:

    r: What the equation is with respect to. This is commonly "r" or "x"
    n: The principle Quantum Number

    Returns:

    The Laguerre polynomial. This is commonly used to solve the LAGUERRE_ASSOC() function.

    """
    
    return simplify((exp(r)*Derivative(exp(-r)*(r**n), (r, n))).doit())




def laguerre_2(r, n):
    """

    Parameters:

    r: What the equation is with respect to. This is commonly "r" or "x"
    n: The principle Quantum Number

    Returns:

    The Laguerre polynomial, without simplification and without the derivative computed.

    """
    
    return (exp(r)*Derivative(exp(-r)*(r**n), (r, n)))




def laguerre_assoc(n, l):
    """

    Parameters:

    n: The principle Quantum Number
    l: The angular momentum Quantum Number

    Returns:

    The Laguerre Associated polynomial, commonly used to solve the RADIAL() function.

    """
    
    return simplify(((-1)**l)*Derivative(laguerre_2(r, n+l), (r, l)).doit())



def kronecker(i, j):
    """
    
    Parameters:
    
    i: this is usually "i" and the first variable of the kronecker delta.
    j: this is usually given a numerical value, and is the second variable of the kronecker delta.
    
    Returns:
    
    This function does not print anything, but instead is used in conjunction with the code found in the kronecker_delta file on:
    https://github.com/liz-stippell/Pysces
    
    """
    if i == j:
        return 1
    else:
        return 0
    
    

""" The following are for Spherical Harmonics """


def spherical(expr):
    """
    
    Parameters:
    
    expr: The expression of interest to be changed into spherical coordinates
    
    Returns:
    
    The expression of interest, A, in terms of spherical coordinates
    
    """
    
    x, r, theta, phi, y, z = symbols("x r theta phi y z")
    return expr.replace(x, r*sin(theta)*cos(phi)).replace(y, r*sin(theta)*sin(phi)).replace(z, r*cos(theta))



def L_2(j, m):
    """
    
    Parameters:
    
    j: 
    m:
    
    Returns:
    
    The L^2 vector magnitude eigenvalue for spherical harmonics.
    
    """
    
    h_b = Symbol("h_b")
    return Bra(str(j), str(","), str(m))*j*(j+1)*h_b**2*Ket(str(j), str(","), str(m))



def L_z(j, m):
    """
    
    Parameters:
    
    j:
    m:
    
    Returns:
    
    The L_z projection (in the z direction) eigenvalue for spherical harmonics.
    
    """
    
    h_b = Symbol("h_b")
    return Bra(str(j), str(","), str(m))*m*h_b*Ket(str(j), str(","), str(m))



def L_raising_operator(j = None, m = None):
    """
    
    Parameters:
    
    j:
    m:
    
    Returns:
    
    If j == None and m == None, the general formula for the raising operator for spherical harmonics is returned.
    
    Else, the formula for the raising operator is computed using Dirac notation
    
    """
    
    L_x, L_y, h_b = symbols("L_x L_y h_b")
    if j == None and m == None:
        return Operator(L_x) + I*Operator(L_y)
    else:
        return Bra(str(j), str(","), str(m))*h_b*sqrt(j*(j+1)-m*(m+1))*Ket(str(j), str(','), str(m+1))

    
    
def L_lowering_operator(j = None, m = None):
     """
    
    Parameters:
    
    j:
    m:
    
    Returns:
    
    If j == None and m == None, the general formula for the lowering operator for spherical harmonics is returned.
    
    Else, the formula for the lowering operator is computed using Dirac notation
    
    """
        
     L_x, L_y, h_b = symbols("L_x L_y h_b")
     if j == None and m == None:
        return Operator(L_x) - I*Operator(L_y)
     else:
        return Bra(str(j), str(","), str(m))*h_b*sqrt(j*(j+1)-m*(m-1))*Ket(str(j), str(','), str(m-1))


    
def L_x(j = None, m = None):
     """
    
    Parameters:
    
    j:
    m:
    
    Returns:
    
    If j == None and m == None, the general formula for the L_x operator for spherical harmonics is returned.
    
    Else, the formula for the L_x operator is computed using Dirac notation
    
    """
        
     if j == None and m == None:
        L_R, L_L = symbols("L_+ L_-")
        return (1/2)*(L_R + L_L)
     else:
        return (1/2)*(L_raising_operator(j, m) + L_lowering_operator(j, m))
    
# In[ ]:




