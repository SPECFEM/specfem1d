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
! This file contains the modules : grid, otherVariables
! and globalStiffMat
!=====================================================================

module grid

  implicit none

  include "constants.h"

  ! Gauss-Lobatto-Legendre points of integration
    double precision, dimension(NGLL) :: xigll
  ! Gauss-Lobatto-Jacobi points of integration
    double precision, dimension(NGLJ) :: xiglj
  ! Weights
    double precision, dimension(NGLL) :: wgll
  ! GLJ weights
    double precision, dimension(NGLJ) :: wglj
  ! Array with derivatives of Lagrange polynomials
    double precision, dimension(NGLL,NGLL) :: hprime
  ! Array with derivatives of GLJ polynomials
    double precision, dimension(NGLJ,NGLJ) :: hprimeBar,hprimeBarwglj
  ! Anchors
    double precision, dimension(NSPEC) :: x1,x2
  ! Global grid points
    double precision, dimension(NGLOB) :: x
  ! Material properties
    double precision, dimension(NGLL,NSPEC) :: rho,mu
  ! Jacobian `matrix' and Jacobian
    double precision, dimension(NGLL,NSPEC) :: dxi_dx,jacobian
  ! Local to global numbering
    integer, dimension(NGLL,NSPEC) :: ibool
  ! Movie
    character(len=50) moviefile
  ! OUTPUT_FILES directory
    character*15 :: path='./OUTPUT_FILES/'

contains

  ! ====================== Implementation part ===============

  subroutine defineDerivationMatrices()
    ! Define the matrices l'(xi) (and l_bar'(xi) for axisym)
    call define_derivative_matrix(xigll,wgll,hprime)
    if (AXISYM) call define_GLJ_derivation_matrix(xiglj,wglj,hprimeBar,hprimeBarwglj)
  end subroutine defineDerivationMatrices

  subroutine makeGrid()
    ! Define the grid (evenly spaced anchors for now -> TODO implement irregular grid)
    integer ispec,i,j,iglob,it

    ! Evenly spaced achors between 0 and 1
    do ispec = 1,NSPEC
      x1(ispec) = LENGTH*dble(ispec-1)/dble(NSPEC)
      x2(ispec) = LENGTH*dble(ispec)/dble(NSPEC)
    enddo
   ! Set up the mesh properties
    do ispec = 1,NSPEC
      do i = 1,NGLL
        rho(i,ispec) = DENSITY
        mu(i,ispec) = RIGIDITY
        dxi_dx(i,ispec) = 2. / (x2(ispec)-x1(ispec)) ! This is d(xi) / dx
        jacobian(i,ispec) = (x2(ispec)-x1(ispec)) / 2.
      enddo
    enddo
   ! Set up local to global numbering
    iglob = 1
    do ispec = 1,NSPEC
      do i = 1,NGLL
        if(i > 1) iglob = iglob+1
        ibool(i,ispec) = iglob
      enddo
    enddo
   ! Compute the position of the global grid points
    do ispec = 1,NSPEC
      do i = 1,NGLL
          iglob = ibool(i,ispec)
        if ((ispec == 1) .and. AXISYM) then
          x(iglob) = 0.5*(1.-xiglj(i))*x1(ispec)+0.5*(1.+xiglj(i))*x2(ispec)
        else
          x(iglob) = 0.5*(1.-xigll(i))*x1(ispec)+0.5*(1.+xigll(i))*x2(ispec)
        endif
      enddo
    enddo
  end subroutine makeGrid

  subroutine writeOutSnapshots(it,displ)
    ! Write out snapshots of the wavefield
    integer it,iglob
    double precision, dimension(NGLOB) :: displ

    if(mod(it,NSNAP) == 0) then
      print *,'time step ',it,' out of ',NSTEP
      if(RUN_BACKWARDS) then
        write(moviefile,"('snapshot_backward',i5.5)") it
      else
        write(moviefile,"('snapshot_forward_normal',i5.5)") it
      endif
      open(unit=10,file=path//moviefile,status='unknown')
      do iglob = 1,NGLOB
        write(10,*) sngl(x(iglob)),sngl(displ(iglob))
      enddo
      close(10)
    endif
  end subroutine writeOutSnapshots

  subroutine checkReconstruction(it,displ)
    ! To check that wavefield reconstruction is accurate when we run backwards
    integer it,iglob
    double precision, dimension(NGLOB) :: displ

    if(mod(it,100) == 1) then ! Value of 1 is set purposely, because we use NSTEP-it+1
      if(RUN_BACKWARDS) then
        write(moviefile,"('snapshot_backward_inverted_time',i5.5)") NSTEP-it+1
        open(unit=10,file=path//moviefile,status='unknown')
        do iglob = 1,NGLOB
          write(10,*) sngl(x(iglob)),sngl(displ(iglob))
        enddo
        close(10)
      endif
    endif
  end subroutine checkReconstruction

end module grid

!==============================================================================

module otherVariables

  use grid
  use gaussm3

  implicit none

  ! Local mass matrix
    double precision mass_local
  ! Global mass matrix
    double precision, dimension(NGLOB) :: mass_global
  ! Time marching
    double precision dh,v
    double precision deltat,deltatover2,deltatsqover2
  ! Source
    integer ispec_source,i_source,iglob_source
    double precision hdur,source_amp,stf
    double precision, external :: source_time_function
  ! Receiver
    integer ireceiver
    double precision seismogram(NSTEP)
  ! Derivatives
    double precision du_dxi,epsilon,sigma,templ,temp(NGLL)
  ! Displacement, velocity and acceleration
    double precision, dimension(NGLOB) :: displ,veloc,accel
  ! Synthetic displacement
    double precision, dimension(NGLOB) :: syntheticDispl
  ! Movie synthetics
    character(len=50) moviefileSynth,moviefileRatio

contains

  ! ====================== Implementation part ===============

  subroutine estimateTimeStep()

    ! estimate the time step in seconds
    if(AXISYM) then
      dh = x(2) ! In axisymetric formulation the smallest interval is the first one
    else
      dh = LENGTH/dble(NGLOB-1)
    endif

    v = dsqrt(RIGIDITY/DENSITY)
    deltat = courant_CFL*dh/v

    print *,'time step estimate: ',deltat,' seconds'
    ! simply inverting delta_t is fine even when there is C*v damping in the case of an explicit Newmark scheme

    if(RUN_BACKWARDS) deltat = - deltat

  end subroutine estimateTimeStep

  subroutine assembleGlobalMassMatrix()
  ! Calculate the assembled global mass matrix
    integer ispec,i,iglob
    mass_global(:) = 0.
    do ispec = 1,NSPEC
      do i = 1,NGLL
        if ((ispec == 1) .and. AXISYM) then
          mass_local = wglj(i)*(rho(i,ispec) + C*deltat/2)*jacobian(i,ispec) !TODO for AXISYM add the r
        else
          mass_local = wgll(i)*(rho(i,ispec) + C*deltat/2)*jacobian(i,ispec)
        endif
        iglob = ibool(i,ispec)
        mass_global(iglob) = mass_global(iglob) + mass_local
      enddo
    enddo
  end subroutine assembleGlobalMassMatrix

  subroutine setSourceReceivers()
    ! Set the source and receivers
    ! Set the source
    ispec_source = ISPEC_OF_THE_SOURCE !NSPEC/2
    i_source = I_OF_THE_SOURCE !2
    hdur = HDUR_DIVIDED_BY_DT_OF_THE_SOURCE*deltat
    source_amp = AMPLITUDE_OF_THE_SOURCE !1.0d7
    ! Set the receiver
    ireceiver = 2*NGLL-2
   end subroutine setSourceReceivers

  subroutine initializeTimeMarching()
    integer iglob
      ! Time marching parameters
    deltatover2 = deltat/2.
    deltatsqover2 = deltat*deltat/2.
    ! Tnitialize the fields to zero
    displ(:) = 0.
    veloc(:) = 0.
    accel(:) = 0.
    if(USE_PLANE_WAVE_SOURCE) then
      do iglob = 1,NGLOB
        displ(iglob) = dsin(PI*x(iglob)/LENGTH)
      enddo
    endif
    call system('mkdir -p OUTPUT_FILES')
  end subroutine initializeTimeMarching

  subroutine readBack()
    ! Read back last time frame of forward run

    open(unit=12,file=path//'dump_last_snapshot',status='old',form='unformatted')
    read(12) displ
    read(12) veloc
    read(12) accel
    close(12)
  end subroutine readBack

  subroutine mainTimeLoop()
  ! Main time loop

    integer ispec,i,j,k,iglob,it

    do it = 1,NSTEP ! Main time loop

      ! `Predictor' update displacement using explicit finite-difference time scheme (Newmark)
      displ(:) = displ(:) + deltat*veloc(:) + deltatsqover2*accel(:)
      veloc(:) = veloc(:) + deltatover2*accel(:)
      accel(:) = 0.

      do ispec = 1,NSPEC
        do k = 1,NGLL
          ! Compute d(u) / d(xi)
          du_dxi = 0.
          do i = 1,NGLL
            iglob = ibool(i,ispec)
            if(AXISYM .and. (ispec == 1)) then ! If axisym and first element
              du_dxi = du_dxi + displ(iglob)*hprimeBar(k,i)
            else
              du_dxi = du_dxi + displ(iglob)*hprime(k,i)
            endif
          enddo
          ! Strain (i.e., d(u) / dx in 1D)
          epsilon = du_dxi*dxi_dx(k,ispec)
          ! stress
          sigma = mu(k,ispec)*epsilon
          temp(k) = jacobian(k,ispec)*sigma*dxi_dx(k,ispec)
        enddo ! first loop over the GLL points
        do j = 1,NGLL
          templ = 0.
          do k = 1,NGLL
            if(AXISYM) then
              if (ispec > 1) then ! Not first element
                templ = templ + temp(k)*hprime(k,j)*wgll(k)*x(ibool(k,ispec))/x(ibool(j,ispec))
              else ! First element
                templ = templ + temp(k)*hprimeBarwglj(k,j) ! hprimeBarwglj(k,j) = hprime_bar(k,j)*wglj(k)
              endif
            else ! Not axisym
              templ = templ + temp(k)*hprime(k,j)*wgll(k)
            endif
          enddo
          ! `Corrector' update of acceleration in the Newmark scheme
          ! The minus sign comes from the integration by part done in the weak formulation of the equations
          iglob = ibool(j,ispec)
          accel(iglob) = accel(iglob) - templ
        enddo ! Second loop over the GLL points
      enddo ! End loop over all spectral elements

      if(.not. USE_PLANE_WAVE_SOURCE) call addSourceAtGlobalLevel(it)

      ! Fixed boundary conditions at global level
      if(FIXED_BC) then
        if(.not. AXISYM) accel(1) = 0.
        accel(NGLOB) = 0.
      endif

      ! Add C damping
      accel(:) = accel(:) - C*veloc(:)
      ! Divide by the mass matrix, which is strictly (i.e. perfectly) diagonal
      accel(:) = accel(:)/mass_global(:)
      ! `Corrector' update velocity
      veloc(:) = veloc(:) + deltatover2*accel(:)

      if(AXISYM .and. SYNTHETICS .and. (mod(it,NSYNTH) == 0)) then
        call calculateSynthetics(it,syntheticDispl) ! To Compare with the real solution
        call writeOutSyntheticSnapshots(it,displ,syntheticDispl)
      endif

      call writeOutSnapshots(it,displ)
      if(RUN_BACKWARDS) call checkReconstruction(it,displ)

      seismogram(it) = displ(ireceiver)

    enddo ! End time loop

  end subroutine mainTimeLoop

  subroutine addSourceAtGlobalLevel(it)
    ! Add source at global level
    integer it
    iglob_source = ibool(i_source,ispec_source)
    ! The source is not at the same time step if we are running backwards
    if(.not. RUN_BACKWARDS) then
      stf = source_time_function(dble(it-1)*deltat-hdur,hdur)
    else
      stf = source_time_function(dble(NSTEP-it)*deltat-hdur,hdur)
    endif
    accel(iglob_source) = accel(iglob_source) + stf*source_amp
  end subroutine addSourceAtGlobalLevel

  subroutine saveSeismograms()
    integer it

    open(unit=12,file=path//'seismogram',status='unknown')
    do it = 1,NSTEP
      write(12,*) sngl(dble(it-1)*deltat-hdur),sngl(seismogram(it))
    enddo
    close(12)
  end subroutine saveSeismograms

  subroutine saveLastFrame()
    ! Save last time frame of the forward run in order to be able to run backwards later if needed
    open(unit=12,file=path//'dump_last_snapshot',status='unknown',form='unformatted')
    write(12) displ
    write(12) veloc
    write(12) accel
    close(12)
  end subroutine saveLastFrame

  function func(phi,r,it)
    ! This function is used to calculate synthetics in axisymetric medium
    integer it
    double precision t,c,t2,func,phi,r
    c=sqrt(RIGIDITY/DENSITY)
    t=dble(it-1)*deltat
    t2=t-hdur
    func=source_time_function(t2-phi,hdur)/sqrt(phi**2-r**2/c**2)

    return
  end function func

  subroutine calculateSynthetics(it,syntheticDispl)
    integer i,j,it,ispec,iglob
    double precision r,c,t,t2
    double precision, dimension(NGLOB) :: syntheticDispl

    integer,parameter :: dbp = SELECTED_REAL_KIND (15,307)
    integer :: ngp = 40            ! # of Gauss Points
    real(dbp) :: xabsc(40), weig(40)

    c=sqrt(RIGIDITY/DENSITY)
    t=dble(it-1)*deltat
    do ispec = 1,NSPEC
      do i = 1,NGLL
        iglob = ibool(i,ispec)
        r=x(iglob)
        if ((ispec==1) .and. (i==1)) r=TINYVAL !To avoid r=0
        if (c*t/r <= 1.d0) then
          syntheticDispl(iglob) = ZERO
        else
        !  syntheticDispl(iglob)=(myQgauss(func, r/c, 1.01d0*r/c, 40,r,it)+myQgauss(func, 1.01d0*r/c, t, 40,r,it))
        !============Uncomment the 8 following lines and comment the line above if you have chosen a step source=====
          if (c*t/r <= c*hdur/r+1.d0) then
            syntheticDispl(iglob) = log(c*t/r+sqrt((c*t/r)**2-1)) !acosh(c*t/r)
          else if (c*t/r > c*hdur/r+1.d0) then
            syntheticDispl(iglob) = log(c*t/r+sqrt((c*t/r)**2-1))-log(c/r*(t-hdur)+sqrt((c/r*(t-hdur))**2-1))
          else
            print *,'hdur : ',hdur,' t : ',t,' r/c : ',r/c
            stop 'Error: in principle this stop statement should never happen... '
          endif
        endif
      enddo
    enddo
    syntheticDispl(:)=syntheticDispl(:)*source_amp*jacobian(1,1)*2.000d0/RIGIDITY

  end subroutine calculateSynthetics

  subroutine writeOutSyntheticSnapshots(it,displ,syntheticDispl)
    ! This subroutine has been used to validate AXISYM (to compare with synthetics)
    integer it,iglob
    double precision, dimension(NGLOB) :: syntheticDispl
    double precision, dimension(NGLOB) :: displ
    double precision :: L2norm1 = ZERO
    double precision :: L2norm2 = ZERO



    if(RUN_BACKWARDS) then
      write(moviefileSynth,"('synthetic_snapshot_backward2',i5.5)") it
    else
      write(moviefileSynth,"('synthetic_snapshot_forward_normal',i5.5)") it
    endif
    open(unit=15,file=path//moviefileSynth,status='unknown')
    do iglob = 1,NGLOB
        write(15,*) sngl(x(iglob)),sngl(syntheticDispl(iglob))
    enddo
    close(15)
  !  call calculateL2norms(syntheticDispl,L2norm1,L2norm2)
  !  call calculateH1norms(syntheticDispl,L2norm1,L2norm2)
    close(16)

  end subroutine writeOutSyntheticSnapshots

  subroutine calculateL2norms(syntheticDispl,L2norm1,L2norm2)
    ! This subroutine has been used to validate AXISYM (to compare with synthetics)
    integer iglob
    double precision, dimension(NGLOB) :: syntheticDispl
    double precision L2norm1,L2norm2
    do iglob = 1,floor(NGLOB/2.d0)
      L2norm1 = L2norm1 + syntheticDispl(iglob)**2
    enddo
    L2norm1 = sqrt(L2norm1)
    print *,'L2 norm 1 : ',L2norm1

    do iglob = ceiling(NGLOB/2.d0),NGLOB
      L2norm2 = L2norm2 + syntheticDispl(iglob)**2
    enddo
    L2norm2 = sqrt(L2norm2)
    print *,'L2 norm 2 : ',L2norm2
  end subroutine calculateL2norms

  subroutine calculateH1norms(syntheticDispl,L2norm1,L2norm2)
    ! This subroutine has been used to validate AXISYM (to compare with synthetics)
    integer iglob
    double precision, dimension(NGLOB) :: syntheticDispl
    double precision, dimension(NGLOB-1) :: syntheticDisplDeriv
    double precision L2norm1,L2norm2
    double precision :: H1norm1 = ZERO
    double precision :: H1norm2 = ZERO

    do iglob = 1,NGLOB-1
        syntheticDisplDeriv(iglob) = (syntheticDispl(iglob+1)-syntheticDispl(iglob))/(x(iglob+1)-x(iglob))
    enddo
    do iglob = 1,floor((NGLOB-1)/2.d0)
        H1norm1 = H1norm1 + abs(syntheticDisplDeriv(iglob))**2
    enddo
    print *,'H1 norm 1 : ',H1norm1

    do iglob = ceiling((NGLOB-1)/2.d0),NGLOB-1
        H1norm2 = H1norm2 + abs(syntheticDisplDeriv(iglob))**2
    enddo
    print *,'H1 norm 2 : ',H1norm2

  end subroutine calculateH1norms

end module otherVariables

!==============================================================================

module globalStiffMat
use grid
implicit none

! Formulation of stiffness matrix
! In general, when using SEM together with explicit newmark scheme(beta=0, gamma=1/2 with diagonal damping matrix C),
! we do not need to calculate element stiffness matrixs and then assemble them into global stiffness matrix.
! But in case of using SEM with implicit newmark scheme, global stiffness matrix is needed. Thus we include it.
  double precision :: jocobianl, xixl, B_matrix_left, B_matrix_right, element_stiffness_matrix_block
  double precision, dimension(NGLOB,NGLOB) :: global_stiffness_matrx

contains

  ! ====================== Implementation part ===============

  subroutine assembleGlobalStiffnessMatrix()
  ! Calculate the assembled global stiffness matrix

    integer :: i,j,ispec,i_interior,iglob_row,iglob_rol

    if(ASSEMBLE_GLOBAL_STIFFNESS_MATRIX)then
      global_stiffness_matrx = 0.
      do ispec = 1,NSPEC
        do i = 1,NGLL
          iglob_row = ibool(i,ispec)
          do j = 1,NGLL
            iglob_rol = ibool(j,ispec)
            element_stiffness_matrix_block = 0.d0
            do i_interior = 1,NGLL
              jocobianl = jacobian(i_interior,ispec)
              xixl = dxi_dx(i_interior,ispec)
              B_matrix_left = hprime(i_interior,i) * xixl
              B_matrix_right = hprime(i_interior,j) * xixl
              element_stiffness_matrix_block = element_stiffness_matrix_block + &
                                              B_matrix_left * mu(i_interior,ispec) * B_matrix_right * wgll(i_interior) * jocobianl
            enddo
            global_stiffness_matrx(iglob_row,iglob_rol) = global_stiffness_matrx(iglob_row,iglob_rol) + &
                                                          element_stiffness_matrix_block
          enddo
        enddo
      enddo
    endif

  end subroutine assembleGlobalStiffnessMatrix

end module globalStiffMat


