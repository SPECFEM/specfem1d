
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

  program diffusion

  implicit none

  include "constants.h"

! fixed boundary conditions
  logical, parameter :: FIXED_BC = .true.

! model parameters  (SI)
  double precision, parameter :: THERMALCONDUCTIVITY = 10.0d-01 ! cal/m/s/K
  double precision, parameter :: HEATCAPACITY = 0.3d+03 ! cal/kg/K

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
  double precision, dimension(NGLL,NSPEC) :: rho,heat_capacity,thermal_conductivity

! Jacobian `matrix' and Jacobian
  double precision, dimension(NGLL,NSPEC) :: dxi_dx,jacobian

! local mass matrix
  double precision mass_local

! global mass matrix
  double precision, dimension(NGLOB) :: mass_global

! temperature and temperature gradient
  double precision, dimension(NGLOB) :: temperature,grad_temperature

! local to global numbering
  integer, dimension(NGLL,NSPEC) :: ibool

! time marching
  double precision deltat,deltatover2,deltatsqover2
  double precision dh,diffusivity

! end fluxes
  double precision flux_1,flux_NGLOB

! end temperatures
  double precision temperature_1,temperature_NGLOB

! derivatives
  double precision dt_dx,flux,templ,temp(NGLL)

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
      thermal_conductivity(i,ispec) = THERMALCONDUCTIVITY
      heat_capacity(i,ispec) = HEATCAPACITY
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
      mass_local = wgll(i)*rho(i,ispec)*heat_capacity(i,ispec)*jacobian(i,ispec)
      iglob = ibool(i,ispec)
      mass_global(iglob) = mass_global(iglob) + mass_local
    enddo
  enddo

! estimate the time step
  dh = LENGTH/dble(NGLOB-1)
  diffusivity = THERMALCONDUCTIVITY/(HEATCAPACITY*DENSITY)
!!!!!!  deltat = 0.5*dh*dh/diffusivity
  deltat = 800000.
  print *,'time step estimate: ',deltat,' seconds'

  if(FIXED_BC) then
! set up the temperatures at the ends
    temperature_1 = 10.
    temperature_NGLOB = 0.
  else
! set up the fluxes at the ends
    flux_1 = 0.
    flux_NGLOB = 0.
  endif

! time marching parameters
  deltatover2 = deltat/2.
  deltatsqover2 = deltat*deltat/2.

! initialize
  temperature(:) = 0.
  grad_temperature(:) = 0.

  do it = 1,NSTEP

! update temperature
    temperature(:) = temperature(:) + deltatover2*grad_temperature(:)
    grad_temperature(:) = 0.

    do ispec = 1,NSPEC

      do i = 1,NGLL
! compute dt_dx
        templ = 0.
        do j = 1,NGLL
          iglob = ibool(j,ispec)
          templ = templ + temperature(iglob)*hprime(i,j)
        enddo
        dt_dx = templ*dxi_dx(i,ispec)

! heat flux
        flux = -thermal_conductivity(i,ispec)*dt_dx

        temp(i) = jacobian(i,ispec)*flux*dxi_dx(i,ispec)
      enddo ! first loop over the GLL points

      do i = 1,NGLL
        templ = 0.
        do j = 1,NGLL
          templ = templ + temp(j)*hprime(j,i)*wgll(j)
        enddo

! update the gradient of temperature
        iglob = ibool(i,ispec)
        grad_temperature(iglob) = grad_temperature(iglob) + templ
      enddo ! second loop over the GLL points

    enddo ! end loop over all spectral elements

! fixed boundary conditions at global level
    if(FIXED_BC) then
! left side
      temperature(1) = temperature_1
      templ = 0.
      do i = 1,NGLL
        iglob = ibool(i,1)
        templ = templ + temperature(iglob)*hprime(1,i)
      enddo
      dt_dx = templ*dxi_dx(1,1)
      flux_1 = -thermal_conductivity(1,1)*dt_dx
! right side
      temperature(NGLOB) = temperature_NGLOB
      templ = 0.
      do i = 1,NGLL
        iglob = ibool(i,NSPEC)
        templ = templ + temperature(iglob)*hprime(NGLL,i)
      enddo
      dt_dx = templ*dxi_dx(NGLL,NSPEC)
      flux_NGLOB = -thermal_conductivity(NGLL,NSPEC)*dt_dx
    endif

! add in the end fluxes
    grad_temperature(1) = grad_temperature(1) + flux_1
    grad_temperature(NGLOB) = grad_temperature(NGLOB) - flux_NGLOB

! divide by the mass matrix
    grad_temperature(:) = grad_temperature(:)/mass_global(:)

! update temperature
    temperature(:) = temperature(:) + deltatover2*grad_temperature(:)

! write out snapshots
    if(mod(it,100) == 0) then
      print *,'time step ',it,' out of ',NSTEP
      write(moviefile,"('snapshot',i5.5)") it
      open(unit=10,file=moviefile,status='unknown')
      do iglob = 1,NGLOB
        write(10,*) sngl(x(iglob)),sngl(temperature(iglob))
      enddo
      close(10)
    endif

  enddo ! end time loop

  end program diffusion

