  program diffusion

  implicit none

  include "constants.h"

! number of timesteps
  integer, parameter :: NSTEP = 30000

! time step in seconds
  double precision, parameter :: DT = 100000000. ! s

! fixed boundary conditions
  logical, parameter :: FIXED_BC = .true.

! model parameters  (SI)
  double precision, parameter :: LENGTH = 3.0d+03 ! m
  double precision, parameter :: DENSITY = 2.5d+03 ! kg/m^3
  double precision, parameter :: THERMALCONDUCTIVITY = 10.0d-01 ! cal/m/s/K
  double precision, parameter :: HEATCAPACITY = 0.3d+03 ! cal/kg/K

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
  double precision, dimension(NGLL,NSPEC) :: rho,heat_capacity,thermal_conductivity

! Jacobian `matrix' and Jacobian
  double precision, dimension(NGLL,NSPEC) :: dxidx,jacobian

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
  double precision dh,diffusivity,time_step

! end fluxes
  double precision flux_1,flux_NGLOB

! end temperatures
  double precision temperature_1,temperature_NGLOB

! derivatives
  double precision dtdx,flux,templ,temp(NGLL)

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
      mass_local = wgll(i)*rho(i,ispec)*heat_capacity(i,ispec)*jacobian(i,ispec)
      iglob = ibool(i,ispec)
      mass_global(iglob) = mass_global(iglob) + mass_local
    enddo
  enddo

! estimate the time step
  dh = LENGTH/dble(NGLOB-1)
  diffusivity = THERMALCONDUCTIVITY/(HEATCAPACITY*DENSITY)
  time_step = 0.5*dh*dh/diffusivity
!  print *,'time step estimate: ',time_step,' seconds'

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
  deltat = DT
  deltatover2 = deltat/2.
  deltatsqover2 = deltat*deltat/2.

! initialize
  temperature(:) = 0.
  grad_temperature(:) = 0.

  do itime = 1,NSTEP

! update temperature
    temperature(:) = temperature(:) + deltatover2*grad_temperature(:)
    grad_temperature(:) = 0.

    do ispec = 1,NSPEC

      do i = 1,NGLL
! get dtdx
        templ = 0.
        do j = 1,NGLL
          iglob = ibool(j,ispec)
          templ = templ + temperature(iglob)*hprime(j,i)
        enddo
        dtdx = templ*dxidx(i,ispec)

! heat flux
        flux = -thermal_conductivity(i,ispec)*dtdx

        temp(i) = jacobian(i,ispec)*flux*dxidx(i,ispec)
      enddo ! first loop over the GLL points

      do i = 1,NGLL
        templ = 0.
        do j = 1,NGLL
          templ = templ + temp(j)*hprime(i,j)*wgll(j)
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
        templ = templ + temperature(iglob)*hprime(i,1)
      enddo
      dtdx = templ*dxidx(1,1)
      flux_1 = -thermal_conductivity(1,1)*dtdx
! right side
      temperature(NGLOB) = temperature_NGLOB
      templ = 0.
      do i = 1,NGLL
        iglob = ibool(i,NSPEC)
        templ = templ + temperature(iglob)*hprime(i,NGLL)
      enddo
      dtdx = templ*dxidx(NGLL,NSPEC)
      flux_NGLOB = -thermal_conductivity(NGLL,NSPEC)*dtdx
    endif

! add in the end fluxes
    grad_temperature(1) = grad_temperature(1) + flux_1
    grad_temperature(NGLOB) = grad_temperature(NGLOB) - flux_NGLOB

! divide by the mass matrix
    grad_temperature(:) = grad_temperature(:)/mass_global(:)

! update temperature
    temperature(:) = temperature(:) + deltatover2*grad_temperature(:)

! write out snapshots
    if(mod(itime-1,1000) == 0) then
      write(moviefile,"('snapshot',i5.5)") itime
      open(unit=10,file=moviefile,status='unknown')
      do iglob = 1,NGLOB
        write(10,*) sngl(x(iglob)),sngl(temperature(iglob))
      enddo
      close(10)
    endif

  enddo ! end time loop

  end program diffusion
