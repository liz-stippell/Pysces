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


def comm_1(commutator_1, commutator_2, aux):
    """
    
    comm_1(commutator_1, commutator_2, aux)

    This function is not used directly, it is only used in the comm() function below. 

    Parameters:

    commutator_1: The first operator in the commutator
    commutator_2: The second operator in the commutator
    aux: The auxiliary function. This is defined below, as F(x).

    Returns:

    The commutator of commutator_1 and commutator_2 with respect to the auxiliary function.

    """
    
    return expand((Commutator(Operator(commutator_1), Operator(commutator_2))*aux).doit())




x, p_y, y, p_x, z, p_z, L_z, L_y, L_x = symbols("x, p_y, y, p_x, z, p_z, L_z, L_y, L_x")
L_z = x*p_y - y*p_x
L_y = z*p_x - x*p_z
L_x = y*p_z - z*p_y

def comm(commutator_1, commutator_2, aux):
    """
    
    comm(commutator_1, commutator_2, aux)
    
    This function has a few different outputs, depending on the parameters.

    Parameters:

    commutator_1: the first operator
    commutator_2: the second operator
    aux: the auxiliary function

    Returns:

    This function automatically returns a solved commutator. For more complicated commutators involving the angular momentum operators (L_z, L_y, L_x), please make sure that the correct notation is being used (ex: if you want the angular momentum operator in the "x" direction, please use "L_x"
    
    
    Please note the following may be helpful:
    
    x, p_y, y, p_x, z, p_z, L_z, L_y, L_x = symbols("x, p_y, y, p_x, z, p_z, L_z, L_y, L_x")
    L_z = x*p_y - y*p_x
    L_y = z*p_x - x*p_z
    L_x = y*p_z - z*p_y

    """
    
    if A == L_z and B == L_y:
        return (Commutator(Operator(x*p_y), Operator(z*p_x))*aux + Commutator(Operator(y*p_x), Operator(x*p_z))*aux - Commutator(Operator(x*p_y), Operator(x*p_z))*aux - Commutator(Operator(y*p_x), Operator(z*p_x))*aux)
    if A == L_z and B == L_x:
        return (Commutator(Operator(x*p_y), Operator(y*p_z))*aux + Commutator(Operator(y*p_x), Operator(z*p_y))*aux - Commutator(Operator(x*p_y), Operator(z*p_y))*aux - Commutator(Operator(y*p_x), Operator(y*p_z))*aux)
    if A == L_y and B == L_z:
        return (Commutator(Operator(z*p_x), Operator(x*p_y))*aux + Commutator(Operator(x*p_z), Operator(y*p_x))*aux - Commutator(Operator(z*p_x), Operator(y*p_x))*aux - Commutator(Operator(x*p_z), Operator(x*p_y))*aux)
    if A == L_y and B == L_x:
        return (Commutator(Operator(z*p_x), Operator(y*p_z))*aux + Commutator(Operator(x*p_z), Operator(z*p_y))*aux - Commutator(Operator(z*p_x), Operator(z*p_y))*aux - Commutator(Operator(x*p_z), Operator(y*p_z))*aux)
    if A == L_x and B == L_z:
        return (Commutator(Operator(y*p_z), Operator(x*p_y))*aux + Commutator(Operator(z*p_y), Operator(y*p_x))*aux - Commutator(Operator(y*p_z), Operator(y*p_x))*aux - Commutator(Operator(z*p_y), Operator(x*p_y))*aux)
    if A == L_x and B == L_y:
        return (Commutator(Operator(y*p_z), Operator(z*p_x))*aux + Commutator(Operator(z*p_y), Operator(x*p_z))*aux - Commutator(Operator(y*p_z), Operator(x*p_z))*aux - Commutator(Operator(z*p_y), Operator(z*p_x))*aux)
    
    if A == L_z:
        return (Commutator(Operator(x*p_y), Operator(commutator_2))*aux - Commutator(Operator(y*p_x), Operator(commutator_2))*aux)
    if A == L_y:
        return (Commutator(Operator(z*p_x), Operator(commutator_2))*aux - Commutator(Operator(x*p_z), Operator(commutator_2))*aux)
    if A == L_x:
        return (Commutator(Operator(y*p_z), Operator(commutator_2))*aux - Commutator(Operator(z*p_y), Operator(commutator_2))*aux)
    if B == L_z:
        return (Commutator(Operator(commutator_1), Operator(x*p_y))*aux - Commutator(Operator(commutator_1), Operator(y*p_x))*aux)
    if B == L_y:
        return (Commutator(Operator(commutator_1), Operator(z*p_x))*aux - Commutator(Operator(commutator_1), Operator(x*p_z))*aux)
    if B == L_x:
        return (Commutator(Operator(commutator_1), Operator(y*p_z))*aux - Commutator(Operator(commutator_1), Operator(z*p_y))*aux)
    else:
        return (expression_replace(comm_1(Operator(commutator_1), Operator(commutator_2), aux), sympify(str('x'))) or expression_replace(comm_1(Operator(commutator_1), Operator(commutator_2), aux), sympify(str('y'))) or expression_replace(comm_1(Operator(commutator_1), Operator(commutator_2), aux), sympify(str('z')))) 
    
    
    
    
def factorization(expr, var):
    """
    
    factorization(expr, var)
    
    Parameters:
    
    expr: the expression of interest
    var: the variable of interest. This is most likely going to be the same variable used in the auxiliary function.
    
    
    Returns:
    
    The simplified commutator. This is only used when both of the operators in the commutator are angular momentum operators.
    
    """
    
    L_z, L_y, L_x = symbols("L_z, L_y, L_x")
    if var == z:
        return comm(Operator(z), p_operator(z), f(z))*L_z
    if var == y:
        return comm(Operator(y), p_operator(y), f(y))*L_y
    if var == x:
        return comm(Operator(x), p_operator(x), f(x))*L_x




def p_operator(var = None):
    """
    
    p_operator(var = None)

    Parameters:

    var: The variable that the linear momentum operator is with respect to.

    Returns:

    The linear momentum operator, with respect to the parameter.
 
    Note that the "1" in the derivative is a placeholder, which will be replaced.
    If var == None, the general linear operator is printed, with all three positional arguments "x", "y", and "z"

    """
    
    h_b, m = symbols("h_b m")
    if var == None:
        return Operator(-I*h_b*(Derivative("1", x))) + Operator(-I*h_b*(Derivative("1", y))) + Operator(-I*h_b*(Derivative("1", z)))
    else:
        return Operator(-I*h_b*(Derivative("1", var)))




def kinetic_energy(var = None):
    """
    
    kinetic_energy(var = None)

    Parameters:

    var: The variable in which the derivative is with respect to.

    Returns:

    The kinetic energy operator, with respect to the chosen parameter. 

    Note that the "1" in the derivative is a placeholder, which will be replaced.
    If var == None, the general linear operator is printed, with all three positional arguments "x", "y", and "z"

    """
    
    h_b, m = symbols("h_b m")
    if var == None:
        return (Operator((-I*h_b*(Derivative("1", x)))**2)*(1/(2*m))) + (Operator((-I*h_b*(Derivative("1", y)))**2)*(1/(2*m))) + (Operator((-I*h_b*(Derivative("1", z)))**2)*(1/(2*m)))
    else:
        return (Operator((-I*h_b*(Derivative("1", var)))**2)*(1/(2*m)))



def v(var):
    """
    
    v(var)
    
    Parameters:
    
    var: The variable for the given function. This is usually "x", "y" or "z".
    
    Returns:
    
    The general potential energy operator, V(x)
    
    """
    
    return Operator(Function("v")(var))



def hamiltonian(var):
    """
    
    hamiltonian(var)
    
    Parameters:
    
    var: The variable for the given function. This is usually "x", "y" or "z".
    
    Returns:
    
    The Hamiltonian operator, made up of the kinetic energy operator and the general potential energy operator.
    
    """
    
    return kinetic_energy(var) + v(var)



def expression_replace(expr, var):
    """
    
    This is only used within the "comm()" function. 
    
    expression_replace(expr, var)

    Parameters:

    expr: The expanded commutator to be replaced
    var: The variable/parameter with respect to the chosen commutator.

    Returns:

    This replaces the Derivative(1, x) present in the expanded commutator with either Derivative(F(x), x) or Derivative(x*F(x), x).
    Note that the above states "x", but can be done for x, y, or z variables.
    This will also replace [p_x, x] and similar commutators, which are usually computed when using angular momentum operators or similar operators.
    Note that the "expr" parameter is a str()


    """
    
    L_x, L_z, L_y, p_x, p_y, p_z, x, y, z = symbols("L_x, L_z, L_y, p_x, p_y, p_z, x, y, z")
    
    if str('[p_x,x]') in str(expr):
        return sympify(str(sympify(str(expr).replace(str('[p_x,x]'), str(expression_replace(comm(p_operator(var), Operator(var), f(var)), K))))).replace(str(p_y*z - p_z*y), str(-L_x)))
    if var == x:
        return sympify(str(expr).replace(str(Derivative(1, var)*f(var)), str(Derivative(f(var), var).doit())).replace(str('Derivative(1, x)*x*f(x)'), str(Derivative(K*f(var), var).doit())))
    
    if str('[p_y, y]') in str(expr):
        return sympify(str(sympify(str(R).replace(str('[p_y,y]'), str(expression_replace(comm(p_operator(var), Operator(var), f(var)), K))))).replace(str(p_x*z - p_z*x), str(-L_y)))
    if var == y:
        return sympify(str(expr).replace(str(Derivative(1, var)*f(var)), str(Derivative(f(var), var).doit())).replace(str('Derivative(1, y)*y*f(y)'), str(Derivative(var*f(var), var).doit())))
    
    if str('[p_z,z]') in str(expr):
        return sympify(str(sympify(str(expr).replace(str('[p_z,z]'), str(expression_replace(comm(p_operator(var), Operator(var), f(var)), K))))).replace(str(p_x*y - p_y*x), str(-L_z)))
    elif var == z:
        return sympify(str(expr).replace(str(Derivative(1, var)*f(var)), str(Derivative(f(var), var).doit())).replace(str('Derivative(1, z)*z*f(z)'), str(Derivative(var*f(var), var).doit())))




def f(var):
    """
    
    f(var)

    Parameters:

    var: This is what the auxiliary function is with respect to. It should also match the parameters of the other arguments in the comm() function. (example: if one operator is p_operator(y), the auxiliary function should be f(y).

    Returns:

    This is the auxiliary function, commonly used in the comm() function.
    This simply returns "f(x)", x being with chosen parameter.

    """
    
    return Operator(Function('f')(var))




def HO():
    """
    
    HO()

    Parameters:

    Currently, there are no needed parameters for this equation.

    Returns:

    The Harmonic Oscillator Hamiltonian. Note that "p" is the linear momentum operator, and so the first term is the KINETIC_ENERGY() function mentioned above.

    """
    
    p, k, m, x = symbols("p k m x")
    return (p**2)/(2*m) + (1/2)*k*(x**2)




def planewave(x):
    """
    
    planewave(x)

    Parameters:

    x: What the function is with respect to.

    Returns:

    The PlaneWave function with respect to the parameter.

    """

    k, omega, t = symbols("k omega t")
    return exp(I*(k*x-omega*t))



def PIB(x, L, n):
    """
    
    PIB(x, L, n)

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
    
    PIB_normalize(x, L, n)

    Parameters:

    x: a variable.
    L: Length of the box.
    n: an integer.

    Returns:

    The normalized WaveFunction for Particle in a Box, with respect to the chosen variables. This answer can also be calculated using the normalize_constant() function.

    """
    
    return sqrt(2/L)*sin((n*pi*x)/L)




def gaussian(alpha, x_0, p):
    """
    
    gaussian(alpha, x_0, p)
    
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
    
    gaussian_normalize(alpha, gamma, x_0, p_0)
    
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
    
    conjugate(expr)

    Parameters:

    x: The term of interest. This is commonly a WaveFunction.

    Returns:

    The complex conjugate of a WaveFunction. If there are imaginary terms, they get negated, and if there are no imaginary terms, the WaveFunction is not affected.

    """
    
    return expr.replace(I, -I)




def normalize_constant(WaveFunc, var, lower, upper):
    """
    
    normalize_constant(WaveFunc, var, lower, upper)

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
    
    expectation_value(WaveFunc_1, Operator, WaveFunc_2, var, lower, upper)

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
    
    overlap(WaveFunc_1, WaveFunc_2, var, lower, upper)
    
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
    
    plot_function(func, B, lower, upper)

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
    
    laguerre(r, n)

    Parameters:

    r: What the equation is with respect to. This is commonly "r" or "x"
    n: The principle Quantum Number

    Returns:

    The Laguerre polynomial. This is commonly used to solve the LAGUERRE_ASSOC() function.

    """
    
    return simplify((exp(r)*Derivative(exp(-r)*(r**n), (r, n))).doit())




def laguerre_2(r, n):
    """
    
    laguerre_2(r, n)

    Parameters:

    r: What the equation is with respect to. This is commonly "r" or "x"
    n: The principle Quantum Number

    Returns:

    The Laguerre polynomial, without simplification and without the derivative computed.

    """
    
    return (exp(r)*Derivative(exp(-r)*(r**n), (r, n)))




def laguerre_assoc(n, l):
    """
    
    laguerre_assoc(n, l)

    Parameters:

    n: The principle Quantum Number
    l: The angular momentum Quantum Number

    Returns:

    The Laguerre Associated polynomial, commonly used to solve the RADIAL() function.

    """
    
    return simplify(((-1)**l)*Derivative(laguerre_2(r, n+l), (r, l)).doit())



def kronecker(i, j):
    """
    
    kronecker(i, j)
    
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
    
    spherical(expr)
    
    Parameters:
    
    expr: The expression of interest to be changed into spherical coordinates
    
    Returns:
    
    The expression of interest, A, in terms of spherical coordinates
    
    """
    
    x, r, theta, phi, y, z = symbols("x r theta phi y z")
    return expr.replace(x, r*sin(theta)*cos(phi)).replace(y, r*sin(theta)*sin(phi)).replace(z, r*cos(theta))



def L_2(j, m):
    """
    
    L_2(j, m)
    
    Parameters:
    
    j: The total angular momentum quantum number
    m: The magnetic quantum number
    
    Returns:
    
    The L^2 vector magnitude eigenvalue for spherical harmonics.
    
    """
    
    h_b = Symbol("h_b")
    return Bra(str(j), str(","), str(m))*j*(j+1)*h_b**2*Ket(str(j), str(","), str(m))



def L_z(j, m):
    """
    
    L_z(j, m)
    
    Parameters:
    
    j: The total angular momentum quantum number
    m: The magnetic quantum number
    
    Returns:
    
    The L_z projection (in the z direction) eigenvalue for spherical harmonics.
    
    """
    
    h_b = Symbol("h_b")
    return Bra(str(j), str(","), str(m))*m*h_b*Ket(str(j), str(","), str(m))



def L_raising_operator(j = None, m = None):
    """
    
    L_raising_operator(j = None, m = None)
    
    Parameters:
    
    j: The total angular momentum quantum number
    m: The magnetic quantum number
    
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
     
     L_lowering_operator(j = None, m = None)
    
     Parameters:
    
    j: The total angular momentum quantum number
    m: The magnetic quantum number
    
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
    
    L_x(j = None, m = None)
    
    Parameters:
    
    j: The total angular momentum quantum number
    m: The magnetic quantum number
    
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




