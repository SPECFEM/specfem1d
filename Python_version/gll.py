# -*- coding: utf-8 -*-
"""
Functions for working with GLL points.
"""

from __future__ import (absolute_import, division, print_function)

import numpy as np


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
    deriv = np.zeros((N + 1, N + 1))

    # AKA: diffs[i,j] = ksiGLL[i] - ksiGLL[j]
    diffs = ksiGLL[:, np.newaxis] - ksiGLL[np.newaxis, :]

    for i in range(N + 1):
        # Exclude i since ksiGLL[i] - ksiGLL[i] is obviously 0.
        prod1 = np.product(diffs[i, :i]) * np.product(diffs[i, i + 1:])
        prod1 = 1 / prod1

        for j in range(N + 1):
            if i == 0 and j == 0:
                deriv[i, j] = -N * (N + 1) / 4.0

            elif i == N and j == N:
                deriv[i, j] = N * (N + 1) / 4.0

            elif i < j:
                # Just like prod1, but excluding i and j.
                prod2 = (np.product(diffs[j, :i]) *
                         np.product(diffs[j, i + 1:j]) *
                         np.product(diffs[j, j + 1:]))

                deriv[i, j] = prod1 * prod2

            elif i > j:
                # Just like prod1, but excluding j and i.
                prod2 = (np.product(diffs[j, :j]) *
                         np.product(diffs[j, j + 1:i]) *
                         np.product(diffs[j, i + 1:]))

                deriv[i, j] = prod1 * prod2

    return deriv


def glj_derivative(ksiGLJ):
    """Calculates the values of the derivative of the polynomials associated
    to the Gauss-Lobatto-Jacobi quadrature at the GLJ points"""
    N = len(ksiGLJ) - 1
    # Only need interior points. Exterior points need special treatment.
    ksiGLJ = ksiGLJ[1:N]

    coef = np.zeros(N + 2)
    coef[N:] = 1.0
    L = np.polynomial.legendre.Legendre(coef)
    Ld = L.deriv()
    glj_poly = Ld(ksiGLJ)

    deriv = np.zeros((N + 1, N + 1))

    # Top row
    deriv[0, 0] = -N * (N + 2) / 6.0
    deriv[0, 1:N] = 2 * (-1)**N * glj_poly / ((1 + ksiGLJ) * (N + 1))
    deriv[0, N] = (-1)**N / (N + 1)

    # Left column
    deriv[1:N, 0] = (-1)**(N + 1) * (N + 1) / (2 * glj_poly * (1 + ksiGLJ))

    # Middle block
    for i in range(1, N):
        # The index for ksiGLJ and glj_poly is shifted because we dropped the
        # first element earlier.
        deriv[i, i] = -1 / (2 * (1 + ksiGLJ[i - 1]))
        deriv[i, i + 1:N] = (glj_poly[i:] / glj_poly[i - 1] /
                             (ksiGLJ[i:] - ksiGLJ[i - 1]))
        deriv[i + 1:N, i] = (glj_poly[i - 1] / glj_poly[i:] /
                             (ksiGLJ[i - 1] - ksiGLJ[i:]))

    # Right column
    deriv[1:N, N] = 1 / (1 - ksiGLJ) / glj_poly

    # Bottom row
    deriv[N, 0] = (-1)**(N + 1) * (N + 1) / 4.0
    deriv[N, 1:N] = -1 / (1 - ksiGLJ) * glj_poly
    deriv[N, N] = (N * (N + 2) - 1) / 4.0

    return deriv


def jacobian(ticks, param):
    """Calculates the jacobian dx/dksi of the substitution on the GLL points
       (and GLJ points for the first element in axisymmetric).
       Returns a matrix nSpec*(N+1) containing its value for each element and
       each points"""
    Np1 = len(param.ksiGLL)
    dE = np.diff(ticks) / 2
    dXdKsi = np.repeat(dE, Np1).reshape((-1, Np1))
    return dXdKsi


def jacobian_inverse(ticks, param):
    """Calculate the jacobian dksi/dx of the substitution on the GLL points
       (and GLJ points for the first element in axisymmetric).
       Returns a matrix nSpec*(N+1) containing its value for each element and
       each points"""
    Np1 = len(param.ksiGLL)
    dE = 2 / np.diff(ticks)
    dKsiDx = np.repeat(dE, Np1).reshape((-1, Np1))
    return dKsiDx
