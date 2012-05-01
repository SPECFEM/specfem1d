
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

! If you use this code for your own research, please cite:
!
! @ARTICLE{KoTr99,
! author={D. Komatitsch and J. Tromp},
! title={Introduction to the spectral-element method for 3-{D} seismic wave propagation},
! journal={Geophys. J. Int.},
! year=1999,
! volume=139,
! number=3,
! pages={806-822},
! doi={10.1046/j.1365-246x.1999.00967.x}}
!
! @ARTICLE{KoVi98,
! author={D. Komatitsch and J. P. Vilotte},
! title={The spectral-element method: an efficient tool to simulate the seismic response of 2{D} and 3{D} geological structures},
! journal={Bull. Seismol. Soc. Am.},
! year=1998,
! volume=88,
! number=2,
! pages={368-392}}

  program wave

  implicit none

  include "constants.h"

! use a plane wave source or a point source (a force)
  logical, parameter :: USE_PLANE_WAVE_SOURCE = .false.

! fixed boundary conditions
  logical, parameter :: FIXED_BC = .true.

  integer ispec,i,j,iglob,it

! Gauss-Lobatto-Legendre points of integration
  double precision, dimension(NGLL) :: xigll

! weights
  double precision, dimension(NGLL) :: wgll

! array with derivatives of Lagrange polynomials
  double precision, dimension(NGLL,NGLL) :: hprime

! anchors
  double precision, dimension(NSPEC) :: x1,x2

! global grid points
  double precision, dimension(NGLOB) :: x

! material properties
  double precision, dimension(NGLL,NSPEC) :: rho,mu

! Jacobian `matrix' and Jacobian
  double precision, dimension(NGLL,NSPEC) :: dxi_dx,jacobian

! local mass matrix
  double precision mass_local

! global mass matrix
  double precision, dimension(NGLOB) :: mass_global

! displacement, velocity and acceleration
  double precision, dimension(NGLOB) :: displ,veloc,accel

! local to global numbering
  integer, dimension(NGLL,NSPEC) :: ibool

! time marching
  double precision dh,v
  double precision deltat,deltatover2,deltatsqover2

! source
  integer ispec_source,i_source,iglob_source
  double precision hdur,source_amp,stf
  double precision, external :: source_time_function

! receiver
  integer ireceiver
  double precision seismogram(NSTEP)

! derivatives
  double precision du_dxi,epsilon,sigma,templ,temp(NGLL)

! movie
  character(len=50) moviefile

!++++++++++++++++++++++++++++++++++++++++++++++++++

  call define_derivative_matrix(xigll,wgll,hprime)

! evenly spaced achors between 0 and 1
  do ispec = 1,NSPEC
    x1(ispec) = LENGTH*dble(ispec-1)/dble(NSPEC)
    x2(ispec) = LENGTH*dble(ispec)/dble(NSPEC)
  enddo

! set up the mesh properties
  do ispec = 1,NSPEC
    do i = 1,NGLL
      rho(i,ispec) = DENSITY
      mu(i,ispec) = RIGIDITY
      dxi_dx(i,ispec) = 2. / (x2(ispec)-x1(ispec)) ! this is d(xi) / dx
      jacobian(i,ispec) = (x2(ispec)-x1(ispec)) / 2.
    enddo
  enddo

! set up local to global numbering
  iglob = 1
  do ispec = 1,NSPEC
    do i = 1,NGLL
      if(i > 1) iglob = iglob+1
      ibool(i,ispec) = iglob
    enddo
  enddo

! compute the position of the global grid points
  do ispec = 1,NSPEC
    do i = 1,NGLL
      iglob = ibool(i,ispec)
      x(iglob) = 0.5*(1.-xigll(i))*x1(ispec)+0.5*(1.+xigll(i))*x2(ispec)
    enddo
  enddo

! calculate the assembled global mass matrix
  mass_global(:) = 0.
  do ispec = 1,NSPEC
    do i = 1,NGLL
      mass_local = wgll(i)*rho(i,ispec)*jacobian(i,ispec)
      iglob = ibool(i,ispec)
      mass_global(iglob) = mass_global(iglob) + mass_local
    enddo
  enddo

! estimate the time step in seconds
  dh = LENGTH/dble(NGLOB-1)
  v = dsqrt(RIGIDITY/DENSITY)
  deltat = courant_CFL*dh/v
  print *,'time step estimate: ',deltat,' seconds'

! set the source
  ispec_source = NSPEC/2
  i_source = 2
  hdur = 500.*deltat
  source_amp = 10000000.d0

! set the receiver
  ireceiver = 2*NGLL-2

! time marching parameters
  deltatover2 = deltat/2.
  deltatsqover2 = deltat*deltat/2.

! initialize the fields to zero
  displ(:) = 0.
  veloc(:) = 0.
  accel(:) = 0.

  if(USE_PLANE_WAVE_SOURCE) then
    do iglob = 1,NGLOB
      displ(iglob) = dsin(PI*x(iglob)/LENGTH)
    enddo
  endif

! main time loop
  do it = 1,NSTEP

! `predictor' update displacement using explicit finite-difference time scheme (Newmark)
    displ(:) = displ(:) + deltat*veloc(:) + deltatsqover2*accel(:)
    veloc(:) = veloc(:) + deltatover2*accel(:)
    accel(:) = 0.

    do ispec = 1,NSPEC

      do i = 1,NGLL
! compute d(u) / d(xi)
        du_dxi = 0.
        do j = 1,NGLL
          iglob = ibool(j,ispec)
          du_dxi = du_dxi + displ(iglob)*hprime(i,j)
        enddo

! strain (i.e., d(u) / dx in 1D)
        epsilon = du_dxi*dxi_dx(i,ispec)

! stress
        sigma = mu(i,ispec)*epsilon

        temp(i) = jacobian(i,ispec)*sigma*dxi_dx(i,ispec)
      enddo ! first loop over the GLL points

      do i = 1,NGLL
        templ = 0.
        do j = 1,NGLL
          templ = templ + temp(j)*hprime(j,i)*wgll(j)
        enddo

! `corrector' update of acceleration in the Newmark scheme
! the minus sign comes from the integration by part done in the weak formulation of the equations
        iglob = ibool(i,ispec)
        accel(iglob) = accel(iglob) - templ

      enddo ! second loop over the GLL points

    enddo ! end loop over all spectral elements

! add source at global level
  if(.not. USE_PLANE_WAVE_SOURCE) then
    iglob_source = ibool(i_source,ispec_source)
    stf = source_time_function(dble(it-1)*deltat-hdur,hdur)
    accel(iglob_source) = accel(iglob_source) + stf*source_amp
  endif

! fixed boundary conditions at global level
    if(FIXED_BC) then
      accel(1) = 0.
      accel(NGLOB) = 0.
    endif

! divide by the mass matrix, which is strictly (i.e. perfectly) diagonal
    accel(:) = accel(:)/mass_global(:)

! `corrector' update velocity
    veloc(:) = veloc(:) + deltatover2*accel(:)

! write out snapshots
    if(mod(it,100) == 0) then
      print *,'time step ',it,' out of ',NSTEP
      write(moviefile,"('snapshot',i5.5)") it
      open(unit=10,file=moviefile,status='unknown')
      do iglob = 1,NGLOB
        write(10,*) sngl(x(iglob)),sngl(displ(iglob))
      enddo
      close(10)
    endif

    seismogram(it) = displ(ireceiver)

  enddo ! end time loop

  open(unit=12,file='seismogram',status='unknown')
  do it = 1,NSTEP
    write(12,*) sngl(dble(it-1)*deltat-hdur),sngl(seismogram(it))
  enddo
  close(12)

  end program wave

