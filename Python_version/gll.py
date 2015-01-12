# -*- coding: utf-8 -*-
"""
Functions for working with GLL points.
"""

import numpy as np
from scipy.interpolate import interp1d


# Gauss Lobatto Legendre points and integration weights
GLL_POINTS = {
    4: np.array([-1, -0.6546536707, 0, 0.6546536707, 1]),
    5: np.array([-1, -0.7650553239, -0.2852315164, 0.2852315164, 0.7650553239,
                 1]),
    6: np.array([-1, -0.8302238962, -0.4688487934, 0, 0.4688487934,
                 0.8302238962, 1]),
    7: np.array([-1, -0.8717401485, -0.5917001814, -0.2092992179, 0.2092992179,
                 0.5917001814, 0.8717401485, 1]),
}

GLL_WEIGHTS = {
    4: np.array([0.1, 0.5444444444, 0.7111111111, 0.5444444444, 0.1]),
    5: np.array([0.0666666666, 0.3784749562, 0.5548583770, 0.5548583770,
                 0.3784749562, 0.0666666666]),
    6: np.array([0.0476190476, 0.2768260473, 0.4317453812, 0.4876190476,
                 0.4317453812, 0.2768260473, 0.0476190476]),
    7: np.array([0.0357142857, 0.2107042271, 0.3411226924, 0.4124587946,
                 0.4124587946, 0.3411226924, 0.2107042271, 0.0357142857]),
}

# Gauss Lobatto Jacobi points and integration weights

GLJ_POINTS = {
    4: np.array([-1.0, -0.5077876295, 0.1323008207, 0.7088201421, 1.0]),
    5: np.array([-1.0, -0.6507788566, -0.1563704318, 0.3734893787,
                 0.7972962733, 1.0]),
    6: np.array([-1.0, -0.7401236486, -0.3538526341, 0.09890279315,
                 0.5288423045, 0.8508465697, 1.0]),
    7: np.array([-1.0, -0.7993818545, -0.4919057913, -0.1117339354,
                 0.2835397724, 0.6337933270, 0.8856884817, 1.0]),
}

GLJ_WEIGHTS = {
    4: np.array([0.01333333333, 0.2896566946, 0.7360043695, 0.794338936,
                 0.1666666667]),
    5: np.array([0.006349206349, 0.1503293754, 0.452292685, 0.6858215721,
                 0.5909214468, 0.1142857143]),
    6: np.array([0.003401360544, 0.08473655296, 0.2803032119, 0.5016469619,
                 0.5945754451, 0.4520031342, 0.08333333333]),
    7: np.array([0.001984126984, 0.0510689152, 0.1792187805, 0.3533996177,
                 0.4909749105, 0.5047839706, 0.3550776151, 0.06349206349]),
}


def lagrange_derivative(ksiGLL):
    """Calculates the values of the derivative of the Lagrange polynomials
    at the GLL points"""
    N = len(ksiGLL) - 1
    deriv=np.zeros((N+1,N+1),dtype='d')
    for i in range(N+1):
        for j in range(N+1):
            if i == 0 and j == 0:
                deriv[i,j] = -N*(N+1.)/4.
            elif i == N and j == N:
                deriv[i,j] = N*(N+1.)/4.
            elif i == j:
                deriv[i,j] = 0.
            else:
                prod1=1.
                for m in range(N+1):
                    if m!=i:
                        prod1 *= 1/(ksiGLL[i]-ksiGLL[m])
                sum1=0.
                for m in range(N+1):
                    prod2=1.
                    for k in range(N+1):
                        if k!=m and k!=i:
                            prod2 *= (ksiGLL[j]-ksiGLL[k])
                    sum1 += prod2
                deriv[i,j]=prod1*sum1
    return deriv


def glj_derivative(ksiGLJ):
    """Calculates the values of the derivative of the polynomials associated
    to the Gauss-Lobatto-Jacobi quadrature at the GLJ points"""
    N=len(ksiGLJ)-1
    deriv=np.zeros((N+1,N+1),dtype='d')
    lenX=100000
    xi = np.linspace(-1,1,lenX)
    xi2=(xi[1:len(xi)]+xi[0:len(xi)-1])/2
    dxi=xi[1]-xi[0]
    P = np.polynomial.legendre.legval(xi, np.identity(8))
    Pb=np.zeros_like(P[N])
    Pb[1:]=(P[N][1:]+P[N+1][1:])/(1.+xi[1:])
    Pb[0]=round(Pb[1]) # =+-(N+1)
    PbInterp = interp1d(xi, Pb)
    for i in np.arange(N+1):
        for j in np.arange(N+1):
            if i == 0 and j == 0:
                deriv[i,j] = -N*(N+2.)/6.
            elif i == 0 and 0 < j < N:
                deriv[i,j] = 2*(-1)**N*PbInterp(ksiGLJ[j])/((1+ksiGLJ[j])*(N+1))
            elif i == 0 and j == N:
                deriv[i,j] = (-1)**N/(N+1)
            elif 0 < i < N and j == 0:
                deriv[i,j] = (-1)**(N+1)*(N+1)/(2*PbInterp(ksiGLJ[i])*(1+ksiGLJ[i]))
            elif 0 < i < N and 0 < j < N and i != j:
                deriv[i,j] = 1/(ksiGLJ[j]-ksiGLJ[i])*PbInterp(ksiGLJ[j])/PbInterp(ksiGLJ[i])
            elif 0 < i < N and i == j:
                deriv[i,j] = -1/(2*(1+ksiGLJ[i]))
            elif 0 < i < N and j == N:
                deriv[i,j] = 1/(1-ksiGLJ[i])*1/PbInterp(ksiGLJ[i])
            elif i == N and j == 0:
                deriv[i,j] = (-1)**(N+1)*(N+1)/4.
            elif i == N and 0 < j < N:
                deriv[i,j] = -1/(1-ksiGLJ[j])*PbInterp(ksiGLJ[j])
            elif i == N and j == N:
                deriv[i,j] = (N*(N+2)-1.)/4.
    return deriv


def jacobian(ticks, param):
    """Calculates the jacobian dx/dksi of the substitution on the GLL points
       (and GLJ points for the first element in axisymmetric).
       Returns a matrix nSpec*(N+1) containing its value for each element and
       each points"""
    dE=ticks[1:]-ticks[:len(ticks)-1] # dE[i] : length of element number i
    dXdKsi=np.zeros((len(dE),len(param.ksiGLL)),dtype='d')
    for e in np.arange(len(dE)):
        for k in np.arange(len(param.ksiGLL)):
            dXdKsi[e,k] = dE[e]/2 # Here it does not depend on ksi
    return dXdKsi


def jacobian_inverse(ticks, param):
    """Calculate the jacobian dksi/dx of the substitution on the GLL points
       (and GLJ points for the first element in axisymmetric).
       Returns a matrix nSpec*(N+1) containing its value for each element and
       each points"""
    dE=ticks[1:]-ticks[:len(ticks)-1] # dE[i] : length of element number i
    dKsiDx=np.zeros((len(dE),len(param.ksiGLL)),dtype='d')
    for e in np.arange(len(dE)):
        for k in np.arange(len(param.ksiGLL)):
            dKsiDx[e,k] = 2/dE[e]   # Here it does not depend on ksi
    return dKsiDx
