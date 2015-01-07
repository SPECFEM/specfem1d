# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 10:36:57 2013

This file gathers all the fundamental functions used for 1D spectral
elements simulations.

@author: Alexis Bottero (alexis.bottero@gmail.com)
"""

import numpy as np


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
