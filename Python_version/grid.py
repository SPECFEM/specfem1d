# -*- coding: utf-8 -*-
'''
Definitions of the grid.
'''

import numpy as np
import matplotlib.pyplot as plt

import functions
import gll


class OneDimensionalGrid(object):
    """Contains the grid properties"""

    def __init__(self,param):
        """Init"""
        self.param=param
        self.z=np.zeros(param.nGlob)
        self.rho=np.zeros((param.nSpec,param.nGLL))
        self.mu=np.zeros((param.nSpec,param.nGLL))
        self.ticks=np.zeros(param.nSpec+1)
        if param.gridType == 'homogeneous':
            self.ticks = np.linspace(0, param.length, param.nSpec + 1)
            for e in np.arange(param.nSpec):
                for i in np.arange(param.nGLL):
                    self.rho[e,i] = param.meanRho
                    self.mu[e,i] = param.meanMu
            for i in np.arange(param.nGlob-1)+1:
                if i < param.nGLJ:
                    self.z[i] = functions.project_inverse(param.ksiGLJ[i], 0,
                                                          self.ticks)
                if i >= param.nGLL and i < param.nGlob-1:
                    self.z[i] = functions.project_inverse(
                        param.ksiGLL[i % param.N],
                        i // param.N,
                        self.ticks)
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
        self.dXdKsi = gll.jacobian(self.ticks, param)
        self.dKsiDx = gll.jacobian_inverse(self.ticks, param)

    def plot(self, fig=0):
        """Plot the grid
        my_ticks gives the abscissa of the borders
        TODO I should test : _the types of the parameters
                             _their sizes"""
        plt.figure(fig)
        sub1=plt.subplot(211)
        plt.hold(True)
        for e in np.arange(self.param.nSpec):
            for i in np.arange(self.param.nGLL):
                plt.plot(self.z[self.param.ibool[e,i]],self.rho[e,i],'b+')
        sub1.set_title(r'$\rho(z)$')
        plt.xticks(self.ticks)
        plt.grid(True)
        sub2=plt.subplot(212)
        for e in np.arange(self.param.nSpec):
            for i in np.arange(self.param.nGLL):
                plt.plot(self.z[self.param.ibool[e,i]],self.mu[e,i],'r+')
        sub2.set_title(r'$\mu(z)$')
        plt.suptitle(" Grid ")
        plt.xticks(self.ticks)
        plt.grid(True)
        plt.show()
        plt.hold(False)
