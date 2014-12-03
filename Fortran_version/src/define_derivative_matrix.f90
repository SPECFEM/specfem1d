
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

  subroutine define_derivative_matrix(xigll,wgll,hprime)

  implicit none

  include "constants.h"

! Gauss-Lobatto-Legendre points of integration
  double precision, dimension(NGLL) :: xigll

! weights
  double precision, dimension(NGLL) :: wgll

! array with derivatives of Lagrange polynomials
  double precision, dimension(NGLL,NGLL) :: hprime

! function for calculating derivatives of Lagrange polynomials
  double precision, external :: lagrange_deriv_GLL

  integer i1,i2

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wgll,NGLL,GAUSSALPHA,GAUSSBETA)

! if number of points is odd, the middle abscissa is exactly zero
  if(mod(NGLL,2) /= 0) xigll((NGLL-1)/2+1) = 0.d0

! calculate derivatives of the Lagrange polynomials
! and precalculate some products in double precision
! hprime(i,j) = h'_j(xigll_i) by definition of the derivation matrix
  do i1=1,NGLL
    do i2=1,NGLL
      hprime(i2,i1) = lagrange_deriv_GLL(i1-1,i2-1,xigll,NGLL)
    enddo
  enddo

  end subroutine define_derivative_matrix

!____________________________________________________________________________________
!
!subroutine define_GLJ_derivation_matrix(xiglj,wxglj,hprimeBar,hprimeBarwglj)
! Calculate all that we need for the GLJ quadrature on axial elements :
! Weights, GLJ points and derivatives of polynomials at GLJ points.
!____________________________________________________________________________________
!

  subroutine define_GLJ_derivation_matrix(xiglj,wxglj,hprimeBar,hprimeBarwglj)

  implicit none

  include "constants.h"

! Gauss-Lobatto-Jacobi points of integration
  double precision, dimension(NGLJ) :: xiglj

! weights for the GLJ quadrature
  double precision, dimension(NGLJ) :: wxglj

! array with derivatives of GLJ polynomials
  double precision, dimension(NGLJ,NGLJ) :: hprimeBar,hprimeBarwglj

  double precision, parameter    :: alphaGLJ=0.d0,betaGLJ=1.d0

! function for calculating derivatives of GLJ polynomials
  double precision, external :: poly_deriv_GLJ

  integer i1,i2

! set up coordinates of the Gauss-Lobatto-Jacobi points
  call zwgljd(xiglj,wxglj,NGLJ,alphaGLJ,betaGLJ)

! calculate derivatives of the GLJ quadrature polynomials
! and precalculate some products in double precision
! hprimeBar(i,j) = hBar'_j(xiglj_i) by definition of the derivation matrix
  do i1=1,NGLJ
    do i2=1,NGLJ
      hprimeBar(i2,i1) = poly_deriv_GLJ(i1-1,i2-1,xiglj,NGLJ)
      hprimeBarwglj(i2,i1) = wxglj(i2) * hprimeBar(i2,i1)
    enddo
  enddo

  end subroutine define_GLJ_derivation_matrix
