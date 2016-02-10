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

!========= Logicals =========

! Axisymmetric run?
  logical,parameter :: AXISYM = .false.
! Calculate synthetics? Just work when AXISYM = .true. so far
  logical,parameter :: SYNTHETICS = .false.
! Run forward or backward?
  logical, parameter :: RUN_BACKWARDS = .false.
! Use a plane wave source or a force point source?
  logical, parameter :: USE_PLANE_WAVE_SOURCE = .false.
! Fixed boundary conditions?
  logical, parameter :: FIXED_BC = .false.
! Assemble global siffness matrix?
  logical :: ASSEMBLE_GLOBAL_STIFFNESS_MATRIX = .false.

!========= Spectral elements =========

! Number of spectral elements
  integer, parameter :: NSPEC = 500
! Number of GLL points (polynomial degree plus one)
  integer, parameter :: NGLL = 4
! Courant-Friedrichs-Lewy (CFL) stability value
! see e.g. http://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition
  double precision, parameter :: courant_CFL = 0.40d0 !0.45d0
! Number of time steps to compute
  integer, parameter :: NSTEP = 1641

!========= Outputs =========

! Every how many time steps do we compute synthetics? (SYNTHETICS must be .true.)
  integer,parameter :: NSYNTH = 100 !NSTEP
! Snapshots every how many steps ?
  integer,parameter :: NSNAP = 200 !NSTEP

!========= Model parameters  (SI) =========

  double precision, parameter :: LENGTH = 50.0d+03 ! m
  double precision, parameter :: DENSITY = 1.5d+03 ! kg/m^3
  double precision, parameter :: RIGIDITY = 3.0d+10 ! Pa
! -> c = sqrt(RIGIDITY/DENSITY)  m/s
! Add C*v damping, off by default
  double precision, parameter :: C = 0 !15000

!========= Source parameters =========

  integer, parameter :: ISPEC_OF_THE_SOURCE = 1 !NSPEC/2
  integer, parameter :: I_OF_THE_SOURCE = 1! 2
  double precision, parameter :: HDUR_DIVIDED_BY_DT_OF_THE_SOURCE = 82.035395797065604
  double precision, parameter :: AMPLITUDE_OF_THE_SOURCE = 1.0d7

!==================DO=NOT=CHANGE=WHAT=IS=BELOW========================

! Number of global points
  integer, parameter :: NGLOB = (NGLL-1)*NSPEC + 1
! Number of GLJ points (For axisymmetric runs)
  integer, parameter :: NGLJ = NGLL
! For the Gauss-Lobatto-Legendre points and weights
  double precision, parameter :: GAUSSALPHA = 0.d0,GAUSSBETA = 0.d0
! Numbers
  double precision, parameter :: ZERO=0.d0,ONE=1.d0,TWO=2.d0
  double precision, parameter :: TINYVAL = 1.d-15
! Pi
  double precision, parameter :: PI = 3.141592653589793d0

