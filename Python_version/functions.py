# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 10:36:57 2013

This file gathers all the fundamental functions used for 1D spectral
elements simulations.

@author: Alexis Bottero (alexis.bottero@gmail.com)
"""

# --- MODULES AND PACKAGES ---
import numpy as np  # NumPy (multidimensional arrays, linear algebra, ...)
from numpy import diff
import scipy as sp  # SciPy (signal and image processing library)
from scipy import misc
from scipy.interpolate import interp1d
import matplotlib as mpl         # Matplotlib (2D/3D plotting library)
import matplotlib.pyplot as plt  # Matplotlib's pyplot: MATLAB-like syntax
from pylab import *              # Matplotlib's pylab interface

# --- FUNCTIONS ---
import defines       # Contains all the constants and parameters

def globalArray(nSpec,nGLL):
    """Returns a matrix A. A[element_number,GLL_considered] -> index in the
    global array, if we work in axisym and that the element number is 0 the points
    are GLJ points"""
    ibool = np.zeros((nSpec,nGLL),dtype='d')
    for e in np.arange(nSpec):
        for i in np.arange(nGLL):
            ibool[e,i] = (nGLL-1)*e+i
    return ibool

def estimateDt(grid,param):
    dzMin=(grid.z[1:]-grid.z[:len(grid.z)-1]).min()
    dh=param.length/(len(grid.z)-1)
    vMax = np.sqrt(param.meanMu/param.meanRho).max()
    dt = param.cfl*dh/vMax # defines.CFL*dzMin/vMax # TODO : See why...?!
    return dt
#    z0=((ksi+1)/(1-ksi)*ticks[elt_number+1]+ticks[elt_number])*1/(1+(ksi+1)/(1-ksi))
def invProjection(ksi, elt_number, ticks):
    """From the abscissa of a point (ksi) onto the reference interval,
       the number of the element to which it belongs (elt_number) and
       the abscissa of the borders of the elements (ticks) returns
       its abscissa onto the interval [zmin, zmax]"""
    z0=1./2.*((ksi+1.)*ticks[elt_number+1]+(1.-ksi)*ticks[elt_number])
    return z0

def lagrangeDeriv(ksiGLL):
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

def GLJderiv(ksiGLJ):
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
                deriv[i,j] = (-1)**(N+1)*(N+1)/4
            elif i == N and 0 < j < N:
                deriv[i,j] = -1/(1-ksiGLJ[j])*PbInterp(ksiGLJ[j])
            elif i == N and j == N:
                deriv[i,j] = (N*(N+2)-1.)/4.
    return deriv

def jacob(ticks,param):
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

def jacobInv(ticks,param):
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

def makeStiffness(grid,param):
    """Computation of stiffness matrices"""
    Ke=np.zeros((param.nSpec,param.nGLL,param.nGLL),dtype='d')
    for e in np.arange(param.nSpec):
        if param.axisym and e == 0 : # The first element
            Ke[e,:,:]=makeStiffness0(grid,param) # Build the stiffness matrix for first element
        else:
            for i in np.arange(param.nGLL):
                for j in np.arange(param.nGLL):
                    sumk=0
                    for k in np.arange(param.nGLL):
                        if param.axisym is not True:
                            sumk+=param.wGLL[k]*grid.mu[e,k]*param.deriv[i,k]* \
                            param.deriv[j,k]*(grid.dKsiDx[e,k]**2)* \
                            grid.dXdKsi[e,k]
                        else:
                            sumk+=param.wGLL[k]*grid.mu[e,k]*param.deriv[i,k]* \
                            param.deriv[j,k]*(grid.dKsiDx[e,k]**2)* \
                            grid.dXdKsi[e,k]*invProjection(param.ksiGLL[k],e,grid.ticks)/invProjection(param.ksiGLL[i],e,grid.ticks)
                    Ke[e,i,j]= sumk
    return Ke

def makeStiffness0(grid,param):
    """Computation of stiffness matrix for the first element"""
    K0=np.zeros((param.nGLJ,param.nGLJ),dtype='d')
    for i in np.arange(param.nGLJ):
        for j in np.arange(param.nGLJ):
            sumk=0
            for k in np.arange(param.nGLJ):
                sumk+=param.wGLJ[k]*grid.mu[0,k]*param.derivGLJ[i,k]* \
                param.derivGLJ[j,k]*(grid.dKsiDx[0,k]**2)*grid.dXdKsi[0,k]
            K0[i,j]= sumk
    return K0

def makeMass(grid,param):
    """Computation of global mass matrix"""
    M=np.zeros(param.nGlob,dtype='d')
    for e in np.arange(param.nSpec):
        for i in np.arange(param.nGLL):
            if param.axisym and e != 0:
                Me=param.wGLL[i]*grid.rho[e,i]*grid.dXdKsi[e,i]#*invProjection(param.ksiGLL[i],e,grid.ticks) # Axisym
            if param.axisym and e == 0:
                Me=param.wGLJ[i]*grid.rho[e,i]*grid.dXdKsi[e,i]
            else:
                Me=param.wGLL[i]*grid.rho[e,i]*grid.dXdKsi[e,i]#*invProjection(param.ksiGLL[i],e,grid.ticks)
            M[param.ibool[e,i]]+=Me
    return M
