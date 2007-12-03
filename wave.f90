
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

! number of timesteps
  integer, parameter :: NSTEP = 10000

! time step in seconds
  double precision, parameter :: DT = 0.0002 ! s

! fixed boundary conditions
  logical, parameter :: FIXED_BC = .true.

! model parameters  (SI)
  double precision, parameter :: LENGTH = 3.0d+03 ! m
  double precision, parameter :: DENSITY = 2.5d+03 ! kg/m^3
  double precision, parameter :: RIGIDITY = 3.0d+10 ! Pa

! pi
  double precision, parameter :: PI = 3.141592653589793d0

  integer ispec,i,j,iglob,itime

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
  double precision, dimension(NGLL,NSPEC) :: dxidx,jacobian

! local mass matrix
  double precision mass_local

! global mass matrix
  double precision, dimension(NGLOB) :: mass_global

! displacement, velocity and acceleration
  double precision, dimension(NGLOB) :: displ,veloc,accel

! local to global numbering
  integer, dimension(NGLL,NSPEC) :: ibool

! time marching
  double precision dh,v,courant,time_step
  double precision deltat,deltatover2,deltatsqover2

! source
  integer ispec_source,i_source
  double precision hdur,source_amp
  double precision, external :: source_time_function

! receiver
  integer ireceiver
  double precision seismogram(NSTEP)

! derivatives
  double precision dudx,sigma,templ,temp(NGLL)

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
      dxidx(i,ispec) = 2. / (x2(ispec)-x1(ispec))
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

! get the global grid points
  do ispec = 1,NSPEC
    do i = 1,NGLL
      iglob = ibool(i,ispec)
      x(iglob) = 0.5*(1.-xigll(i))*x1(ispec)+0.5*(1.+xigll(i))*x2(ispec)
    enddo
  enddo

! calculate the global mass matrix
  mass_global(:) = 0.
  do ispec = 1,NSPEC
    do i = 1,NGLL
      mass_local = wgll(i)*rho(i,ispec)*jacobian(i,ispec)
      iglob = ibool(i,ispec)
      mass_global(iglob) = mass_global(iglob) + mass_local
    enddo
  enddo

! estimate the time step
  dh = LENGTH/dble(NGLOB-1)
  v = dsqrt(RIGIDITY/DENSITY)
  courant = 0.2
  time_step = courant*dh/v
!  print *,'time step estimate: ',time_step,' seconds'

! set the source
  ispec_source = (NSPEC+1)/2
  i_source = (NGLL+1)/2
  hdur = 50.*DT
  source_amp = 0.001

! set the receiver
  ireceiver = 2*NGLL-2

! time marching parameters
  deltat = DT
  deltatover2 = deltat/2.
  deltatsqover2 = deltat*deltat/2.

! initialize
  displ(:) = 0.
  veloc(:) = 0.
  accel(:) = 0.

  do iglob = 1,NGLOB
    displ(iglob) = dsin(PI*x(iglob)/LENGTH)
  enddo

  do itime = 1,NSTEP

! `predictor' update displacement using finite-difference time scheme (Newmark)
    displ(:) = displ(:) + deltat*veloc(:) + deltatsqover2*accel(:)
    veloc(:) = veloc(:) + deltatover2*accel(:)
    accel(:) = 0.

    do ispec = 1,NSPEC

      do i = 1,NGLL
! get dudx
        templ = 0.
        do j = 1,NGLL
          iglob = ibool(j,ispec)
          templ = templ + displ(iglob)*hprime(i,j)
        enddo
        dudx = templ*dxidx(i,ispec)

! stress
        sigma = mu(i,ispec)*dudx

        temp(i) = jacobian(i,ispec)*sigma*dxidx(i,ispec)
      enddo ! first loop over the GLL points

      do i = 1,NGLL
        templ = 0.
        do j = 1,NGLL
          templ = templ + temp(j)*hprime(j,i)*wgll(j)
        enddo

! `corrector' update acceleration
        iglob = ibool(i,ispec)
        accel(iglob) = accel(iglob) - templ
      enddo ! second loop over the GLL points

    enddo ! end loop over all spectral elements

! add source at global level
!    iglob_source = ibool(ispec_source,i_source)
!    stf = source_time_function(dble(itime-1)*DT-hdur,hdur)
!    accel(iglob_source) = accel(iglob_source) + stf*source_amp

! fixed boundary conditions at global level
    if(FIXED_BC) then
      accel(1) = 0.
      accel(NGLOB) = 0.
    endif

! divide by the mass matrix
    accel(:) = accel(:)/mass_global(:)

! `corrector' update velocity
    veloc(:) = veloc(:) + deltatover2*accel(:)

! write out snapshots
    if(mod(itime-1,1000) == 0) then
      write(moviefile,"('snapshot',i5.5)") itime
      open(unit=10,file=moviefile,status='unknown')
      do iglob = 1,NGLOB
        write(10,*) sngl(x(iglob)),sngl(displ(iglob))
      enddo
      close(10)
    endif

    seismogram(itime) = displ(ireceiver)

  enddo ! end time loop

  open(unit=12,file='seismogram',status='unknown')
  do itime = 1,NSTEP
    write(12,*) sngl(dble(itime-1)*DT-hdur),sngl(seismogram(itime))
  enddo
  close(12)

  end program wave

