!     Last change:  TC   12 Sep 98    5:06 pm

!     Platform Information
!  1. Program Title: LF90
!     Command Line: lf90 <name> -nwrap -chk -vax -xref -wo -g
!     Working Directory: c:\lf90\bin

!     Files:
!     c:\lf90\bin\gaussm2.lib
!     c:\lf90\bin\gaussm2.mod
!     c:\wormyf90\pz1d\gaussm2.f90

!     VERSIONS: gaussm3
!   This is based on gaussm2 - but with major changes that WEIGHTS & ABSCISSAS are no longer common variables
!     They are obtained from gauleg as local variables. User must specify what order of integration (# Gaussian Pts)
!     Consequently, all Gaussian Quadrature routines are modified with CALL gauleg.


!
!     9/9/98 - All ALLOCATED variables (pointers & allocatables) have been deallocated except for
!   WEIGhts and ABSCissas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NOTES
! 1. This MODULE can be compiled on its own like header file.
!    Then a *.lib file is created. To use functions in module,
!     must have "USE module_name" statement.
! 2. FATAL -- Specification of assumed size may be done for last dimension only.
!        a(*,*) wrong      a(9,*) correct
! 3. But -- Specification of assumed shape can be done in ANY dimension
!        a(:,:)  O.K.
! 4. "Module" is different from "include". Non-module files which have functions
!    can be included when compiling main file by using "INCLUDE" inside the main
!    file.
! 5. Dummy arguments cannot be ALLOCATABLE. !!BUT!! Dummy arguments can have
!    assumed shapes such as a(:,:)
! 6. INTENT must be used with dummy argument only.
! 7. Array must not have both ALLOCATABLE and POINTER attributes, but
!    pointers can be ALLOCATED.
! 8. POINTER array (eg. ans) must be declared with a deferred shape specification
!    (see "Dynamic Arrays" in the Lahey Fortran 90 Language Reference).
!    ans(3,4)  wrong          ans(:,:) correct
! 9. When dummy arguments are pointers, their actual argument must be pointers
!    too; and vice versa.
! 10. Intrinsic function MAXLOC(find location of maximum element e.g. in a list)
!     returns an ARRAY.
! 11. Intrinsic function MATMUL for Matrix Multiplication can handle REAL of
!     selected_kind (eg. double precision)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  1. Topics: Gaussian elimination solution of linear equations & Gaussian Qudrature!
!   2. The variables in this module are GLOBALLY COMMON to all PROGRAM UNITS that "use"s
! this module and all INTERNAL subprograms within this module. The common variables are the
! Gauss-Legendre abscissas and weights for the Gaussian Qudrature.
!   3. The user must first call "gauleg" from the main program to set up all the abscissas
! and weights.

module gaussm3

                        ! Common Global variables within module !
    implicit none
   INTEGER, parameter :: dbp = SELECTED_REAL_KIND (15,307)
!   PRIVATE
   REAL (dbp) :: newv
   REAL(dbp)  :: EPS, M_PI
   PARAMETER (EPS=3.0d-15)        !EPS is the relative precision
   PARAMETER (M_PI=3.141592654d0)      ! Pi value

!   PUBLIC :: newv, EPS, M_PI, n, xabsc, weig, dbp, qgss2d

   INTERFACE

   END INTERFACE

   CONTAINS
!* This module has the following INTERNAL FUNCTIONS:
!* gauleg, qgauss, qgss3d, qgss2d, gsselm, identity_matrix
!* This module has the following INTERNAL SUBROUTINES:
!* linear_solver
!* They can call on each other without first specifying their type
!* NO INTERFACE nor EXTERNAL is required since they are INTERNAL functions

!********************************************************************************
!* Calculation of GAUSS-LEGENDRE abscissas and weights for Gaussian Quadrature
!* integration of polynomial functions.
!*      For normalized lower and upper limits of integration -1.0 & 1.0, and
!* given n, this routine calculates, arrays xabsc(1:n) and  weig(1:n) of length n,
!* containing the abscissas and weights of the Gauss-Legendre n-point quadrature
!* formula.  For detailed explanations finding weights & abscissas, see
!* "Numerical Recipes in Fortran */
!********************************************************************************
  SUBROUTINE  gauleg(ngp, xabsc, weig)

      implicit none
      INTEGER  i, j, m
      REAL(dbp)  p1, p2, p3, pp, z, z1
      INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
      REAL(dbp), INTENT(OUT) :: xabsc(ngp), weig(ngp)


     m = (ngp + 1) / 2
!* Roots are symmetric in the interval - so only need to find half of them  */

     do i = 1, m        ! Loop over the desired roots */

        z = cos( M_PI * (i-0.25d0) / (ngp+0.5d0) )
!*   Starting with the above approximation to the ith root,
!*          we enter the main loop of refinement by NEWTON'S method   */
100       p1 = 1.0d0
          p2 = 0.0d0
!*  Loop up the recurrence relation to get the Legendre
!*  polynomial evaluated at z                 */

          do j = 1, ngp
            p3 = p2
            p2 = p1
            p1 = ((2.0d0*j-1.0d0) * z * p2 - (j-1.0d0)*p3) / j
          enddo

!* p1 is now the desired Legendre polynomial. We next compute pp,
!* its derivative, by a standard relation involving also p2, the
!* polynomial of one lower order.      */
          pp = ngp*(z*p1-p2)/(z*z-1.0d0)
          z1 = z
          z = z1 - p1/pp             ! Newton's Method  */

          if (dabs(z-z1) > EPS) GOTO  100

        xabsc(i) =  - z                     ! Roots will be bewteen -1.0 & 1.0 */
        xabsc(ngp+1-i) =  + z                 ! and symmetric about the origin  */
        weig(i) = 2.0d0/((1.0d0-z*z)*pp*pp) ! Compute the weight and its       */
        weig(ngp+1-i) = weig(i)               ! symmetric counterpart         */

      enddo     ! i loop

   End subroutine gauleg

!********************************************************************************
!*     Returns the SINGLE integral of the function (of ONE VARIABLE) "func"
!* between x1 and x2 by N-point Gauss-Legendre integration. The function
!* is evaluated exactly N times at interior points in the range of
!* integration.       */
!********************************************************************************
   recursive function qgauss(func, x1, x2, ngp) RESULT(intgrl)
     implicit none
     REAL(dbp)  intgrl, x1, x2, func
     REAL(dbp)  xm, xl
     INTEGER j
     INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
     REAL(dbp) :: xabsc(ngp), weig(ngp)

     call gauleg(ngp, xabsc, weig)

     intgrl = 0.0d0
     xm = 0.5 * (x2 + x1)
     xl = 0.5 * (x2 - x1)
     do j = 1, ngp
       intgrl = intgrl + weig(j) * func( xm + xl*xabsc(j) )
     enddo
     intgrl = intgrl * xl;    !Scale the answer to the range of integration  */
   end function qgauss

    recursive function myQgauss(func, x1, x2, ngp,r,it) RESULT(intgrl)
      implicit none
      REAL(dbp)  intgrl, x1, x2, func
      REAL(dbp)  xm, xl
      INTEGER j
      INTEGER, INTENT(IN) :: ngp !# of Gauss Points
      INTEGER, INTENT(IN) :: it
      DOUBLE PRECISION, INTENT(IN) :: r
      REAL(dbp) :: xabsc(ngp), weig(ngp)

      call gauleg(ngp, xabsc, weig)

      intgrl = 0.0d0
      xm = 0.5 * (x2 + x1)
      xl = 0.5 * (x2 - x1)
      do j = 1, ngp
         intgrl = intgrl + weig(j) * func( xm + xl*xabsc(j) ,r,it)
      enddo
      intgrl = intgrl * xl;    !Scale the answer to the range of integration  */
    end function myQgauss



   recursive function qgaussmat1(func, x1, x2, ngp, frow, fcol) RESULT(intgrl)
      implicit none
      REAL(dbp), INTENT(IN) :: x1, x2
      INTEGER :: frow, fcol
    REAL(dbp) :: intgrl(frow, fcol), tmpm(frow,fcol)
      REAL(dbp) ::  func
      REAL(dbp) ::  xm, xl, arg
      INTEGER j

      INTEGER, INTENT(IN) :: ngp            ! # of Gauss Points
      REAL(dbp) :: xabsc(ngp), weig(ngp)

      call gauleg(ngp, xabsc, weig)

      intgrl(:,:) = 0.0d0
      tmpm(:,:) = 0.0d0
     xm = 0.5 * (x2 + x1)
    xl = 0.5 * (x2 - x1)
      do j = 1, ngp
         arg =  xm + xl*xabsc(j)
      PRINT *, 'szhgd ds'
          PRINT *,arg
         tmpm = func(arg)
         intgrl = intgrl + weig(j) * tmpm
      enddo

    intgrl = intgrl * xl;    !Scale the answer to the range of integration  */
   END function qgaussmat1

end module gaussm3


