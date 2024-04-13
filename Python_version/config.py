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

from __future__ import (absolute_import, division, print_function)

import argparse
try:
    # Python 3
    from configparser import ConfigParser
except ImportError:
    # Python 2
    from ConfigParser import SafeConfigParser

import numpy as np

import gll


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
        cp = ConfigParser(defaults={
            # True if axial symmetry
            'axisym': True,
            # Boundary's type
            'BOUND_TYPE':'NONE',
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
            # Time steps between snapshots (0 == disabled)
            'SNAPSHOT': 0,
            # Plot grid, source, and periodic results
            'PLOT': False,
            # One image is displayed each DPLOT time step
            'DPLOT': 10,
        })
        with open('Par_file') as f:
            try:
                # Python 3
                cp.read_string('[global]\n' + f.read(), source='Par_file')
            except AttributeError:
                # Python 2
                cp.readfp(FakeGlobalSectionHead(f))

        self.axisym = cp.getboolean('global', 'AXISYM')
        self.boundType = cp.get('global','BOUND_TYPE').strip("'\"")
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
        self.snapshot = cp.getint('global', 'SNAPSHOT')
        self.plot = cp.getboolean('global', 'PLOT')
        self.dplot = cp.getfloat('global', 'DPLOT')

        parser = argparse.ArgumentParser(
            description='Spectral element method in a 1D medium')
        parser.add_argument('--no-plot', action='store_true',
                            help='Force disable plotting')
        args = parser.parse_args()
        self.plot = self.plot and not args.no_plot

        # Number of GLL points per elements
        self.nGLL = self.N + 1
        # Number of GLJ in the first element
        self.nGLJ = self.NGLJ + 1
        # Number of points in the array
        self.nGlob = (self.nSpec - 1) * self.N + self.NGLJ + 1
        self.ibool = self.make_global_index()
        # Time step (will be updated)
        self.dt = 0
        # dh (will be updated)
        self.dh = 0


        # Gauss Lobatto Legendre points and integration weights:
        try:
            # Position of the GLL points in [-1,1]
            self.ksiGLL = gll.GLL_POINTS[self.N]
            # Integration weights
            self.wGLL = gll.GLL_WEIGHTS[self.N]
        except KeyError:
            raise ValueError('N = %d is invalid!' % (self.N, ))
        try:
            # Position of the GLJ points in [-1,1]
            self.ksiGLJ = gll.GLJ_POINTS[self.NGLJ]
            # Integration weights
            self.wGLJ = gll.GLJ_WEIGHTS[self.NGLJ]
        except KeyError:
            raise ValueError('NGLJ = %d is invalid!' % (self.NGLJ, ))

        # Derivatives of the Lagrange polynomials at the GLL points
        self.deriv = gll.lagrange_derivative(self.ksiGLL)
        self.derivGLJ = gll.glj_derivative(self.ksiGLJ)

    def make_global_index(self):
        """Returns a matrix A. A[element_number,GLL_considered] -> index in the
        global array, if we work in axisym and that the element number is 0 the
        points are GLJ points"""
        # TODO: add GLJ
        ibool = np.zeros((self.nSpec, self.nGLL), dtype=np.intp)
        for e in np.arange(self.nSpec):
            for i in np.arange(self.nGLL):
                ibool[e, i] = (self.nGLL - 1) * e + i
        return ibool


class Source(object):
    """Contains the source properties"""

    def __init__(self, param):
        """Init"""
        self.typeOfSource = param.sourceType
        self.ampl = param.maxAmpl
        if self.typeOfSource == 'ricker':
            self.hdur = param.tSource * param.dt  # Duration of the source (s)
            self.decayRate = param.decayRate
            self.alpha = self.decayRate / self.hdur
        else:
            raise ValueError('Unknown source type: %s' % (self.typeOfSource, ))

    def __getitem__(self, t):
        """What happens when we do source[t]"""
        t -= self.hdur
        return -2 * self.ampl * self.alpha**3 / np.sqrt(np.pi) * (
            t * np.exp(-(self.alpha * t)**2))

    def plotSource(self):
        """Plot the source"""
        import matplotlib.pyplot as plt
        t = np.linspace(0, self.hdur, 1000)
        plt.figure()
        plt.plot(t, self[t], 'b')
        plt.title('Source(t)')
        plt.grid(True)
        plt.show()
