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

### THIS PART CAN BE MODIFIED -> ###
AXISYM=True
# Grid description
GRID_TYPE='homogeneous' # Grid's type
LENGTH=3000 # "Physical" length of the domain (in meters)
DENSITY=2500 #kg/m^3
RIGIDITY=30000000000 #Pa
GRID_FILE='grid_homogeneous.txt'
TICKS_FILE='ticks_homogeneous.txt'

# Source description
TSOURCE=100 #500.         # Duration of the source in dt
ISOURCE=0 #501 #501 #501 #501 #501     # GLL point number on which the source is situated
MAX_AMPL=1e7         # Maximum amplitude
SOURCE_TYPE='ricker' # Source's type
DECAY_RATE=2.628 #10 #2.628

# Other constants
NSPEC=250        # Number of elements
N=4              # Degree of the basis functions
NGLJ=N           # Degree of basis functions in the first element
# For the moment NGLJ need to be = N
NTS=20000        # Number of time steps
CFL=0.45         # Courant CFL number

# Plots
DPLOT=10        # One image is displayed each DPLOT time step
### <- THIS PART CAN BE MODIFIED ###

### THIS PART IS NOT SUPPOSED TO BE MODIFIED -> ###
# --- MODULES AND PACKAGES ---
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


try:
    ksiGLL = ksiGLL[N]
    wGLL = wGLL[N]
except KeyError:
    print >> sys.stderr, 'ERROR: N = %d is invalid!' % (N, )
    sys.exit()

try:
    ksiGLJ = ksiGLJ[NGLJ]
    wGLJ = wGLJ[NGLJ]
except KeyError:
    print >> sys.stderr, 'ERROR: NGLJ = %d is invalid!' % (NGLJ, )
    sys.exit()


class Parameter:
    """Contains all the constants and parameters necessary for 1D spectral
    element simulation"""

    def __init__(self):
        """Init"""
        self.axisym=AXISYM              # True if axisymmetry
        self.length=LENGTH              # "Physical" length of the domain (in meters)
        self.nSpec=NSPEC                # Number of elements
        self.N=N                        # Degree of the basis functions
        self.NGLJ=NGLJ                  # Degree of the basis functions in the first element
        self.nGLL=self.N+1              # Number of GLL points per elements
        self.nGLJ=self.NGLJ+1           # Number of GLJ in the first element
        self.nGlob=(self.nSpec-1)*self.N+self.NGLJ+1 # Number of points in the array
        self.ibool=functions.globalArray(self.nSpec,self.nGLL)  # Global array TODO add GLJ
        self.nts=NTS                    # Number of time steps
        self.cfl=CFL                    # Courant CFL number
        self.dt=0                       # Time step (will be updated)
        #Grid description
        self.gridType=GRID_TYPE         # Grid's type
        self.meanRho=DENSITY
        self.meanMu=RIGIDITY
        # Source description :
        self.tSource=TSOURCE            # Duration of the source in dt
        self.iSource=ISOURCE            # GLL point number on which the source is situated
        self.maxAmpl=MAX_AMPL           # Maximum amplitude
        self.sourceType=SOURCE_TYPE     # Source's type
        self.decayRate=DECAY_RATE       # Decay rate for the ricker
        # Gauss Lobatto Legendre points and integration weights :
        self.ksiGLL=ksiGLL              # Position of the GLL points in [-1,1]
        self.wGLL=wGLL                  # Integration weights
        self.ksiGLJ=ksiGLJ              # Position of the GLJ points in [-1,1]
        self.wGLJ=wGLJ                  # Integration weights
        # Derivatives of the Lagrange polynomials at the GLL points
        self.deriv=functions.lagrangeDeriv(self.ksiGLL)
        self.derivGLJ=functions.GLJderiv(self.ksiGLJ)

class OneDgrid:
    """Contains the grid properties"""

    def __init__(self,param):
        """Init"""
        self.param=param
        self.gridFile=GRID_FILE
        self.ticksFile=TICKS_FILE
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
            [self.z,self.rho,self.mu]=np.loadtxt(defines.GRID_FILE).transpose()
            self.ticks=np.loadtxt(defines.TICKS_FILE)
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

### <- THIS PART IS NOT SUPPOSED TO BE MODIFIED ###

