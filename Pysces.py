#!/usr/bin/env python
# coding: utf-8

# In[3]:


from sympy import *
from sympy.physics.quantum import *
from sympy.abc import *
from sympy.plotting import *
from sympy import init_printing
init_printing() 


""" The above are a list of commonly used variables for the Pysces library. """


def COMM_1(A, B, f):
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




def COMM(x, y, a, b = None):
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
        return COMM_1(Operator(x), Operator(y), a)
    elif b == 0:
        return Commutator(Operator(x), Operator(a)) + Commutator(Operator(y), Operator(a))
    else:       
        return Commutator(Operator(x), Operator(a)) - Commutator(Operator(x), Operator(b)) + Commutator(Operator(y), Operator(a)) + Commutator(Operator(y), Operator(b))




def FACTOR_COMM(A, B, C, D):
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




def P_OPERATOR(A = None):
    """

    Parameters:

    x: what the linear momentum operator is with respect to.

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




def KINETIC_ENERGY(A = None):
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



def V(x):
    """
    
    Parameters:
    
    x: The variable for the given function. This is usually "x", "y" or "z".
    
    Returns:
    
    The general potential energy operator, V(x)
    
    """
    
    return Operator(Function("V")(x))



def HAMILTONIAN(A):
    """
    
    Parameters:
    
    A: The variable for the given function. This is usually "x", "y" or "z".
    
    Returns:
    
    The Hamiltonian operator, made up of the kinetic energy operator and the general potential energy operator.
    
    """
    
    return KINETIC_ENERGY(A) + V(A)



def EXPRESSION_REPLACE(R, K):
    """

    Parameters:

    R: The expanded commutator to be replaced
    K: The variable/parameter with respect to the chosen commutator.

    Returns:

    This replaces the Derivative(1, x) present in the expanded commutator with either Derivative(F(x), x) or Derivative(x*F(x), x).
    Note that the above states "x", but can be done for x, y, or z variables.
    Note that the R parameter is a str()


    """
    
    if K == x:
        return sympify(str(R).replace(str(Derivative(1, K)*F(K)), str(Derivative(F(K), K).doit())).replace(str('Derivative(1, x)*x*F(x)'), str(Derivative(K*F(K), K).doit())))
    if K == y:
        return sympify(str(R).replace(str(Derivative(1, K)*F(K)), str(Derivative(F(K), K).doit())).replace(str('Derivative(1, y)*y*F(y)'), str(Derivative(K*F(K), K).doit())))
    if K == z:
        return sympify(str(R).replace(str(Derivative(1, K)*F(K)), str(Derivative(F(K), K).doit())).replace(str('Derivative(1, z)*z*F(z)'), str(Derivative(K*F(K), K).doit())))




def F(x):
    """

    Parameters:

    x: This is what the auxiliary function is with respect to. It should also match the parameters of the other arguments in the COMM() function. (example: if one operator is P_OP(y), the auxiliary function should be F(y).

    Returns:

    This is the auxiliary function, commonly used in the COMM() function.
    This simply returns "F(x)", x being with chosen parameter.

    """
    
    return Operator(Function('F')(x))




def HO():
    """

    Parameters:

    Currently, there are no needed parameters for this equation.

    Returns:

    The Harmonic Oscillator Hamiltonian. Note that "p" is the linear momentum operator, and so the first term is the KINETIC_ENERGY() function mentioned above.

    """
    
    p, k, m, x = symbols("p k m x")
    return (p**2)/(2*m) + (1/2)*k*(x**2)




def PLANEWAVE(x):
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




def PIB_NORMALIZE(x, L, n):
    """

    Parameters:

    x: a variable.
    L: Length of the box.
    n: an integer.

    Returns:

    The normalized WaveFunction for Particle in a Box, with respect to the chosen variables. This answer can also be calculated via the NORMALIZE() function.

    """
    
    return sqrt(2/L)*sin((n*pi*x)/L)




def CONJUGATE(x):
    """

    Parameters:

    x: The term of interest. This is commonly a WaveFunction.

    Returns:

    The complex conjugate of a WaveFunction. If there are imaginary terms, they get negated, and if there are no imaginary terms, the WaveFunction is not affected.

    """
    
    return x.replace(I, -I)




def NORMALIZE(A, x, y, z):
    """

    Parameters:

    A: The WaveFunction/expression of interest
    x: What the integral is taken with respect to
    y: The lower bound of the integral. If bounds are not listed, this is -oo
    z: The upper bound of the integral. If bounds are not listed, this is oo

    Returns:

    The normalization constant, with respect to the given parameters. To find the normalized WaveFunction, the output for NORMALIZE() must be multiplied by the original WaveFunction. A normalized WaveFunction indicates that the probability of finding a particle within certain bounds must be equal to one.

    """
    
    return 1/sqrt(Integral(A*CONJUGATE(A), (x, y, z)))




def EXPECTATION(A, B, x, y, z):
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
    
    if B == KINETIC_ENERGY(x):
        return sympify(str(sympify(str(Integral(CONJUGATE(A)*B, (x, y, z))).replace(str(Derivative("1", x)**2), str(Derivative(A, x, x)))).doit()).replace(str('sin(pi*n)'), str(0)).replace(str('cos(pi*n)'), str(0)))
    if B == P_OPERATOR(x):
        return Integral(CONJUGATE(A)*B, (x, y, z)).replace(Derivative("1", x), Derivative(A, x).doit())
    
    if B == KINETIC_ENERGY(y):
        return sympify(str(sympify(str(Integral(CONJUGATE(A)*B, (x, y, z))).replace(str(Derivative("1", y)**2), str(Derivative(A, y, y)))).doit()).replace(str('sin(pi*n)'), str(0)).replace(str('cos(pi*n)'), str(0)))
    if B == P_OPERATOR(y):
        return Integral(CONJUGATE(A)*B, (x, y, z)).replace(Derivative("1", y), Derivative(A, y).doit())
    
    if B == KINETIC_ENERGY(z):
        return sympify(str(sympify(str(Integral(CONJUGATE(A)*B, (x, y, z))).replace(str(Derivative("1", z)**2), str(Derivative(A, z, z)))).doit()).replace(str('sin(pi*n)'), str(0)).replace(str('cos(pi*n)'), str(0)))
    if B == P_OPERATOR(z):
        return Integral(CONJUGATE(A)*B, (x, y, z)).replace(Derivative("1", z), Derivative(A, z).doit())
    
    else:
        return Integral(CONJUGATE(A)*B*A, (x, y, z))




def PLOT(A, B, C, D):
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
        
    return plot(sympify(str(A)), (B, C, D))




def LAGUERRE(r, n):
    """

    Parameters:

    r: What the equation is with respect to. This is commonly "r" or "x"
    n: The principle Quantum Number

    Returns:

    The Laguerre polynomial. This is commonly used to solve the LAGUERRE_ASSOC() function.

    """
    
    return simplify((exp(r)*Derivative(exp(-r)*(r**n), (r, n))).doit())




def LAGUERRE_2(r, n):
    """

    Parameters:

    r: What the equation is with respect to. This is commonly "r" or "x"
    n: The principle Quantum Number

    Returns:

    The Laguerre polynomial, without simplification and without the derivative computed.

    """
    
    return (exp(r)*Derivative(exp(-r)*(r**n), (r, n)))




def LAGUERRE_ASSOC(n, l):
    """

    Parameters:

    n: The principle Quantum Number
    l: The angular momentum Quantum Number

    Returns:

    The Laguerre Associated polynomial, commonly used to solve the RADIAL() function.

    """
    
    return simplify(((-1)**l)*Derivative(LAGUERRE_2(r, n+l), (r, l)).doit())





# In[ ]:




