
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

  use grid
  use otherVariables
  use globalStiffMat
  use GaussM3

  implicit none

  call defineDerivationMatrices()

  call makeGrid()

  call estimateTimeStep()

  call assembleGlobalMassMatrix()

  call setSourceReceivers()

  if(RUN_BACKWARDS) call readBack()

  if(ASSEMBLE_GLOBAL_STIFFNESS_MATRIX) call assembleGlobalStiffnessMatrix()

  call initializeTimeMarching()

  call mainTimeLoop()

  call saveSeismograms()

  if(.not. RUN_BACKWARDS) call saveLastFrame

end program wave

