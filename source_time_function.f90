
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

  double precision function source_time_function(t,hdur)

  implicit none

! source decay rate (also change in source spectrum if needed)
  double precision, parameter :: decay_rate = 2.628d0
  double precision, parameter :: PI = 3.141592653589793d0

  double precision t,hdur

  double precision alpha

  alpha = decay_rate/hdur

! Gaussian
! source_time_function = alpha*exp(-alpha*alpha*t*t)/sqrt(PI)

! Ricker wavelet
   source_time_function = -2.*(alpha**3)*t*exp(-alpha*alpha*t*t)/sqrt(PI)

  end function source_time_function

