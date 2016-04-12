# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 10:36:57 2013

This file gathers all the fundamental functions used for 1D spectral
elements simulations.

@author: Alexis Bottero (alexis.bottero@gmail.com)
"""

from __future__ import (absolute_import, division, print_function)

import numpy as np


def estimate_timestep(grid, param):
    dh = param.length / (len(grid.z) - 1)
    vMax = np.sqrt(param.meanMu / param.meanRho)
    dt = param.cfl * dh / vMax
    # TODO : See why...?!
    # dzMin = np.diff(grid.z).min()
    # dt = param.CFL * dzMin / vMax
    return dt, dh


def project_inverse(ksi, elt_number, ticks):
    """From the abscissa of a point (ksi) onto the reference interval,
       the number of the element to which it belongs (elt_number) and
       the abscissa of the borders of the elements (ticks) returns
       its abscissa onto the interval [zmin, zmax]"""
    z0 = 1. / 2. * ((ksi + 1) * ticks[elt_number + 1] +
                    (1 - ksi) * ticks[elt_number])
    return z0


def make_stiffness_matrix(grid, param):
    """Computation of stiffness matrices"""
    Ke = np.zeros((param.nSpec, param.nGLL, param.nGLL))
    for e in range(param.nSpec):
        if param.axisym and e == 0:
            gllj = param.wGLJ
            ngllj = param.nGLJ
            deriv = param.derivGLJ
        else:
            gllj = param.wGLL
            ngllj = param.nGLL
            deriv = param.deriv

        tmpe = gllj * grid.mu[e, :] * grid.dKsiDx[e, :]**2 * grid.dXdKsi[e, :]
        if param.axisym and e != 0:
            tmpe *= project_inverse(param.ksiGLL, e, grid.ticks)

        for i in range(ngllj):
            tmpi = deriv[i, :] * tmpe
            if param.axisym and e != 0:
                tmpi /= project_inverse(param.ksiGLL[i], e, grid.ticks)

            Ke[e, i, :] = np.inner(tmpi, deriv)

    return Ke


def make_mass_matrix(grid, param):
    """Computation of global mass matrix"""
    M = np.zeros(param.nGlob)
    # NOTE: We cannot simply drop this outer loop because param.ibool contains
    # repeated indices into M and the addition would not work correctly if we
    # tried to slice it in.
    for e in range(param.nSpec):
        if param.axisym and e == 0:
            Me = param.wGLJ * grid.rho[e, :] * grid.dXdKsi[e, :]
        else:
            Me = param.wGLL * grid.rho[e, :] * grid.dXdKsi[e, :]
        # We can still slice along the second index because the individual GLL
        # points are assumed to not repeat a point within an spectral element.
        M[param.ibool[e, :]] += Me
    return M
