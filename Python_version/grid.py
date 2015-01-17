# -*- coding: utf-8 -*-
'''
Definitions of the grid.
'''

from __future__ import (absolute_import, division, print_function)

import numpy as np

import functions
import gll


class OneDimensionalGrid(object):
    """Contains the grid properties"""

    def __init__(self, param):
        """Init"""
        self.param = param
        self.z = np.zeros(param.nGlob)
        self.rho = np.zeros((param.nSpec, param.nGLL))
        self.mu = np.zeros((param.nSpec, param.nGLL))
        self.ticks = np.zeros(param.nSpec + 1)

        if param.gridType == 'homogeneous':
            self.ticks = np.linspace(0, param.length, param.nSpec + 1)
            self.rho.fill(param.meanRho)
            self.mu.fill(param.meanMu)

            self.z[1:param.nGLJ] = functions.project_inverse(
                param.ksiGLJ[1:param.nGLJ],
                0,
                self.ticks)

            ksiGLL = param.ksiGLL[1:]
            for i in range(param.nGLL, param.nGlob, param.N):
                self.z[i:i + param.N] = functions.project_inverse(ksiGLL,
                                                                  i // param.N,
                                                                  self.ticks)

            self.z[-1] = self.ticks[-1]

        elif param.gridType == 'gradient':
            msg = "typeOfGrid == 'gradient' has not been implemented yet"
            raise NotImplementedError(msg)
        elif param.gridType == 'miscellaneous':
            msg = "typeOfGrid == 'miscellaneous' has not been implemented yet"
            raise NotImplementedError(msg)
        elif param.gridType == 'file':
            self.z, self.rho, self.mu = np.loadtxt(param.gridFile, unpack=True)
            self.ticks = np.loadtxt(param.ticksFile)
        else:
            raise ValueError('Unknown grid type: %s' % (param.gridType, ))
        # Jacobians at the GLL (and GLJ for the first element in axisym)
        # points (arrays nSpec*(N+1) elements)
        self.dXdKsi = gll.jacobian(self.ticks, param)
        self.dKsiDx = gll.jacobian_inverse(self.ticks, param)

    def plot(self):
        """Plot the grid
        my_ticks gives the abscissa of the borders
        TODO I should test : _the types of the parameters
                             _their sizes"""
        import matplotlib.pyplot as plt
        from matplotlib.ticker import FixedLocator

        fig, ax = plt.subplots(2, 1, sharex=True)

        ax[0].plot(self.z[self.param.ibool].flat, self.rho.flat, 'b+')
        ax[0].set_title(r'$\rho(z)$')
        ax[0].xaxis.set_minor_locator(FixedLocator(self.ticks))
        ax[0].xaxis.grid(True, which='minor', alpha=0.5)
        ax[0].yaxis.grid(True)

        ax[1].plot(self.z[self.param.ibool].flat, self.mu.flat, 'r+')
        ax[1].set_title(r'$\mu(z)$')
        ax[1].xaxis.set_minor_locator(FixedLocator(self.ticks))
        ax[1].xaxis.grid(True, which='minor', alpha=0.5)
        ax[1].yaxis.grid(True)

        plt.suptitle('Grid')

        plt.show()
