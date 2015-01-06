# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 10:51:04 2013

This file gathers all the constants, classes and parameters used for 1D
spectral elements simulations.
This file is the only one that should be edited. Feel free to change the
constants defined at the beginning.
For the moment just one type of source (ricker) and one type of medium
is implemented.

@author: Alexis Bottero (alexis.bottero@gmail.com)
"""

try:
    # Python 3
    from configparser import SafeConfigParser
except ImportError:
    # Python 2
    from ConfigParser import SafeConfigParser

import numpy as np  # NumPy (multidimensional arrays, linear algebra, ...)
import scipy as sp  # SciPy (signal and image processing library)
import matplotlib as mpl         # Matplotlib (2D/3D plotting library)
import matplotlib.pyplot as plt  # Matplotlib's pyplot: MATLAB-like syntax
from pylab import *              # Matplotlib's pylab interface

# --- FUNCTIONS --- #
import functions        # Contains fundamental functions


# Gauss Lobatto Legendre points and integration weights
ksiGLL = {
    4: np.array([-1, -0.6546536707, 0, 0.6546536707, 1]),
    5: np.array([-1, -0.7650553239, -0.2852315164, 0.2852315164, 0.7650553239,
                 1]),
    6: np.array([-1, -0.8302238962, -0.4688487934, 0, 0.4688487934,
                 0.8302238962, 1]),
    7: np.array([-1, -0.8717401485, -0.5917001814, -0.2092992179, 0.2092992179,
                 0.5917001814, 0.8717401485, 1]),
}

wGLL = {
    4: np.array([0.1, 0.5444444444, 0.7111111111, 0.5444444444, 0.1]),
    5: np.array([0.0666666666, 0.3784749562, 0.5548583770, 0.5548583770,
                 0.3784749562, 0.0666666666]),
    6: np.array([0.0476190476, 0.2768260473, 0.4317453812, 0.4876190476,
                 0.4317453812, 0.2768260473, 0.0476190476]),
    7: np.array([0.0357142857, 0.2107042271, 0.3411226924, 0.4124587946,
                 0.4124587946, 0.3411226924, 0.2107042271, 0.0357142857]),
}

# Gauss Lobatto Jacobi points and integration weights

ksiGLJ = {
    4: np.array([-1.0, -0.5077876295, 0.1323008207, 0.7088201421, 1.0]),
    5: np.array([-1.0, -0.6507788566, -0.1563704318, 0.3734893787,
                 0.7972962733, 1.0]),
    6: np.array([-1.0, -0.7401236486, -0.3538526341, 0.09890279315,
                 0.5288423045, 0.8508465697, 1.0]),
    7: np.array([-1.0, -0.7993818545, -0.4919057913, -0.1117339354,
                 0.2835397724, 0.6337933270, 0.8856884817, 1.0]),
}

wGLJ = {
    4: np.array([0.01333333333, 0.2896566946, 0.7360043695, 0.794338936,
                 0.1666666667]),
    5: np.array([0.006349206349, 0.1503293754, 0.452292685, 0.6858215721,
                 0.5909214468, 0.1142857143]),
    6: np.array([0.003401360544, 0.08473655296, 0.2803032119, 0.5016469619,
                 0.5945754451, 0.4520031342, 0.08333333333]),
    7: np.array([0.001984126984, 0.0510689152, 0.1792187805, 0.3533996177,
                 0.4909749105, 0.5047839706, 0.3550776151, 0.06349206349]),
}


class FakeGlobalSectionHead(object):
    def __init__(self, fp):
        self.fp = fp
        self.sechead = '[global]\n'
    def readline(self):
        if self.sechead:
            try:
                return self.sechead
            finally:
                self.sechead = None
        else:
            return self.fp.readline()


class Parameter(object):
    """Contains all the constants and parameters necessary for 1D spectral
    element simulation"""

    def __init__(self):
        """Init"""
        cp = SafeConfigParser(defaults={
            # True if axial symmetry
            'axisym': True,
            # "Physical" length of the domain (in meters)
            'LENGTH': 3000,
            # Number of elements
            'NSPEC': 250,
            # Degree of the basis functions
            'N': 4,
            # Degree of basis functions in the first element
            'NGLJ': 4,
            # Number of time steps
            'NTS': 2,
            # Courant CFL number
            'CFL': 0.45,
            # Grid description
            'GRID_TYPE': 'homogeneous',
            'GRID_FILE': 'grid_homogeneous.txt',
            'TICKS_FILE': 'ticks_homogeneous.txt',
            # kg/m^3
            'DENSITY': 2500,
            # Pa
            'RIGIDITY': 30000000000,
            # Duration of the source in dt
            'TSOURCE': 100,
            # GLL point number on which the source is situated
            'ISOURCE': 0,
            # Maximum amplitude
            'MAX_AMPL': 1e7,
            # Source's type
            'SOURCE_TYPE': 'ricker',
            # Decay rate for the ricker
            'DECAY_RATE': 2.628,
            # One image is displayed each DPLOT time step
            'DPLOT': 10,
        })
        with open('Par_file') as f:
            cp.readfp(FakeGlobalSectionHead(f))

        self.axisym = cp.getboolean('global', 'AXISYM')
        self.length = cp.getfloat('global', 'LENGTH')
        self.nSpec = cp.getint('global', 'NSPEC')
        self.N = cp.getint('global', 'N')
        self.NGLJ = cp.getint('global', 'NGLJ')
        self.nts = cp.getint('global', 'NTS')
        self.cfl = cp.getfloat('global', 'CFL')
        self.gridType = cp.get('global', 'GRID_TYPE').strip("'\"")
        self.gridFile = cp.get('global', 'GRID_FILE').strip("'\"")
        self.ticksFile = cp.get('global', 'TICKS_FILE').strip("'\"")
        self.meanRho = cp.getfloat('global', 'DENSITY')
        self.meanMu = cp.getfloat('global', 'RIGIDITY')
        self.tSource = cp.getfloat('global', 'TSOURCE')
        self.iSource = cp.getint('global', 'ISOURCE')
        self.maxAmpl = cp.getfloat('global', 'MAX_AMPL')
        self.sourceType = cp.get('global', 'SOURCE_TYPE').strip("'\"")
        self.decayRate = cp.getfloat('global', 'DECAY_RATE')
        self.dplot = cp.getfloat('global', 'DPLOT')

        self.nGLL = self.N + 1              # Number of GLL points per elements
        self.nGLJ = self.NGLJ + 1           # Number of GLJ in the first element
        self.nGlob = (self.nSpec - 1) * self.N + self.NGLJ + 1  # Number of points in the array
        self.ibool = functions.globalArray(self.nSpec, self.nGLL)  # Global array TODO add GLJ
        self.dt = 0                       # Time step (will be updated)

        # Gauss Lobatto Legendre points and integration weights :
        try:
            # Position of the GLL points in [-1,1]
            self.ksiGLL = ksiGLL[self.N]
            # Integration weights
            self.wGLL = wGLL[self.N]
        except KeyError:
            raise ValueError('N = %d is invalid!' % (self.N, ))
        try:
            # Position of the GLJ points in [-1,1]
            self.ksiGLJ = ksiGLJ[self.NGLJ]
            # Integration weights
            self.wGLJ = wGLJ[self.NGLJ]
        except KeyError:
            raise ValueError('NGLJ = %d is invalid!' % (self.NGLJ, ))

        # Derivatives of the Lagrange polynomials at the GLL points
        self.deriv = functions.lagrangeDeriv(self.ksiGLL)
        self.derivGLJ = functions.GLJderiv(self.ksiGLJ)


class OneDgrid:
    """Contains the grid properties"""

    def __init__(self,param):
        """Init"""
        self.param=param
        self.z=np.zeros(param.nGlob)
        self.rho=np.zeros((param.nSpec,param.nGLL))
        self.mu=np.zeros((param.nSpec,param.nGLL))
        self.ticks=np.zeros(param.nSpec+1)
        if param.gridType == 'homogeneous':
            self.ticks=linspace(0,param.length,param.nSpec+1)
            for e in arange(param.nSpec):
                for i in np.arange(param.nGLL):
                    self.rho[e,i] = param.meanRho
                    self.mu[e,i] = param.meanMu
            for i in np.arange(param.nGlob-1)+1:
                if i < param.nGLJ:
                    self.z[i]=functions.invProjection(param.ksiGLJ[i],0,self.ticks)
                if i >= param.nGLL and i < param.nGlob-1:
                    self.z[i]=functions.invProjection(param.ksiGLL[i%param.N],i//param.N,self.ticks)
                else:
                    self.z[param.nGlob-1]=self.ticks[len(self.ticks)-1]
        elif param.gridType == 'gradient':
            print "typeOfGrid == 'gradient' Has not been implemented yet"
            raise
        elif param.gridType == 'miscellaneous':
            print "typeOfGrid == 'miscellaneous' Has not been implemented yet"
            raise
        elif param.gridType == 'file':
            self.z, self.rho, self.mu = np.loadtxt(param.gridFile, unpack=True)
            self.ticks = np.loadtxt(param.ticksFile)
        else :
            print "Unknown grid's type"
            raise
        # Jacobians at the GLL (and GLJ for the first element in axisym)
        # points (arrays nSpec*(N+1) elements)
        self.dXdKsi=functions.jacob(self.ticks,param)
        self.dKsiDx=functions.jacobInv(self.ticks,param)

    def plotGrid(self,fig=0):
        """Plot the grid
        my_ticks gives the abscissa of the borders
        TODO I should test : _the types of the parameters
                             _their sizes"""
        plt.figure(fig)
        sub1=plt.subplot(211)
        plt.hold(True)
        for e in arange(self.param.nSpec):
            for i in np.arange(self.param.nGLL):
                plt.plot(self.z[self.param.ibool[e,i]],self.rho[e,i],'b+')
        sub1.set_title(r'$\rho(z)$')
        plt.xticks(self.ticks)
        plt.grid(True)
        sub2=plt.subplot(212)
        for e in arange(self.param.nSpec):
            for i in np.arange(self.param.nGLL):
                plt.plot(self.z[self.param.ibool[e,i]],self.mu[e,i],'r+')
        sub2.set_title(r'$\mu(z)$')
        plt.suptitle(" Grid ")
        plt.xticks(self.ticks)
        plt.grid(True)
        plt.show()
        plt.hold(False)

class Source:
    """Contains the source properties"""

    def __init__(self,param):
        """Init"""
        self.typeOfSource=param.sourceType
        self.ampl=param.maxAmpl
        if self.typeOfSource == 'ricker':
            self.hdur = param.tSource*param.dt # Duration of the source (s)
            self.decayRate = param.decayRate #2.628
            self.alpha = self.decayRate/self.hdur
        else:
            print "Unknown source's type"
            raise

    def __getitem__(self,t):
        """What happens when we do source[t]"""
        t-=self.hdur
        return self.ampl*-2.*(self.alpha**3)*t*np.exp(-self.alpha*self.alpha*t*t)/np.sqrt(np.pi)

    def plotSource(self,fig=1):
        """Plot the source"""
        t=linspace(0,self.hdur,1000)
        plt.figure(fig)
        plt.plot(t,self[t],'b')
        plt.title('Source(t)')
        plt.grid(True)
        plt.show()
