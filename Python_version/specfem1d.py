#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 10:30:04 2013

Main script for 1D spectral elements.
    u(z,t) = sum(ui(t)li(z)) for i=1...nSpec*(N+1)
    with  nSpec  : number of elements:
         N+1  : number of basis functions per element
         li(z): Lagrange polynomial of degree N
    ->   ui(t): Time dependent expansion coeffs (what we are looking for)

   The parameters of the run can be edited in the file Par_file

@author: Alexis Bottero, CNRS Marseille, France (alexis.bottero@gmail.com)
Modified on Wed Apr 6 14:52 2016
 Add First-order Mur boundary
@author: Zeming Su, Xi'an Jiaotong University, China (suzeming1992@gmail.com)
"""

from __future__ import (absolute_import, division, print_function)

import numpy as np

from config import Parameter
from config import Source
import functions
from grid import OneDimensionalGrid

# --- MAIN BODY --- #
# Initialization
param = Parameter()
grid = OneDimensionalGrid(param)
if param.plot:
    grid.plot()

param.dt, param.dh = functions.estimate_timestep(grid, param)
source = Source(param)
if param.plot:
    source.plotSource()

# Computation of stiffness matrices
Ke = functions.make_stiffness_matrix(grid, param)

# Computation of global mass matrices
M = functions.make_mass_matrix(grid, param)

# Time integration
u = np.zeros(param.nGlob)
vel = np.zeros_like(u)
acc = np.zeros_like(u)

if param.plot:
    import matplotlib.pyplot as plt
    plt.ion()
    fig = plt.figure()
    plt.hold(False)

if param.axisym and (param.plot or param.snapshot):
    cz = np.concatenate((-grid.z[:0:-1], grid.z))

# Compute the start index and end index of ABC
staInx = 0
endInx = param.nGlob
if param.boundType == 'ABC':
    if param.iSource > param.N:
        staInx = 1
    if param.iSource < param.nGlob-param.N:
        endInx = -1
# Main time loop :
for it in range(param.nts):
    print('it = %d (t = %f s)' % (it, it * param.dt))
    if it > 0:
        if param.boundType == 'ABC':
            u[staInx:endInx] += param.dt * vel[staInx:endInx] + acc[staInx:endInx] * param.dt**2 / 2
            vel += param.dt / 2 * acc
            acc.fill(0)
        elif param.boundType == 'NONE':
            u += param.dt * vel + acc * param.dt**2 / 2
            vel += param.dt / 2 * acc
            acc.fill(0)
        else:
             raise ValueError('Unknown type of boundary conditon.')

    for e in range(param.nSpec):
        acc[param.ibool[e, :]] -= np.dot(Ke[e, :, :], u[param.ibool[e, :]])

    acc[param.iSource] += source[it*param.dt]
    acc /= M
    vel += param.dt / 2 * acc
    # process ABC boundary condition
    if param.boundType == 'ABC' :
        if staInx == 1:
            u[0] = u[1] + (param.cfl*param.dh-param.dh)/(param.cfl*param.dh+param.dh) * (u[1]-u[0])
        if endInx == -1:
            u[-1] = u[-2] + (param.cfl*param.dh-param.dh)/(param.cfl*param.dh+param.dh) * (u[-2]-u[-1])

    if param.snapshot and (it % param.snapshot == 0 or it == param.nts - 1):
        name = 'snapshot_forward_normal%05d' % (it, )
        if param.axisym:
            c = np.concatenate((u[:0:-1], u))
            np.savetxt(name, np.column_stack((cz, c)))
        else:
            np.savetxt(name, np.column_stack((grid.z, u)))

    if param.plot and it % param.dplot == 0:
        if param.axisym:
            c = np.concatenate((u[:0:-1], u))
            plt.plot(cz, c)
            plt.xlim([-max(grid.z), max(grid.z)])
            plt.ylim([-10, 10])
            plt.grid(True)
        else:
            plt.plot(grid.z, u)
            # plt.ylim([-0.10, 0.10])
        plt.title("it : {}".format(it))
        plt.draw()
