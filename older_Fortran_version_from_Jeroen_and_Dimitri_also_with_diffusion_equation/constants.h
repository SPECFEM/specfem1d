
!=====================================================================
!
!               S p e c f e m 1 D  V e r s i o n  1 . 0
!               ---------------------------------------
!
!                 Jeroen Tromp and Dimitri Komatitsch
!    Seismological Laboratory - California Institute of Technology
!         (c) California Institute of Technology November 2007
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

! number of spectral elements
  integer, parameter :: NSPEC = 250

! number of GLL points (polynomial degree plus one)
  integer, parameter :: NGLL = 5

! Courant-Friedrichs-Lewy (CFL) stability value
! see e.g. http://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition
  double precision, parameter :: courant_CFL = 0.45d0

! number of time steps to compute
  integer, parameter :: NSTEP = 20000

! model parameters  (SI)
  double precision, parameter :: LENGTH = 3.0d+03 ! m
  double precision, parameter :: DENSITY = 2.5d+03 ! kg/m^3
  double precision, parameter :: RIGIDITY = 3.0d+10 ! Pa

! number of global points
  integer, parameter :: NGLOB = (NGLL-1)*NSPEC + 1

! for the Gauss-Lobatto-Legendre points and weights
  double precision, parameter :: GAUSSALPHA = 0.d0,GAUSSBETA = 0.d0

! pi
  double precision, parameter :: PI = 3.141592653589793d0

