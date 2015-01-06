#=====================================================================
#
#               S p e c f e m 1 D  V e r s i o n  1 . 0
#               ---------------------------------------
#
#                 Jeroen Tromp and Dimitri Komatitsch
#    Seismological Laboratory - California Institute of Technology
#         (c) California Institute of Technology November 2007
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#=====================================================================

O = obj

# Portland
#FC = pgf90
#FCFLAGS =-fast -Mnobounds -Minline -Mneginfo -Mdclchk -Knoieee -Minform=warn -fastsse -tp amd64e -Msmart

# Intel
#FC = ifort
#FCFLAGS = -O1 -vec-report0 -std95 -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage -check all -align sequence -assume byterecl -ftz -traceback -ftrapuv

# GNU gfortran
FC = gfortran
FCFLAGS = -std=gnu -fimplicit-none -frange-check -O2 -fmax-errors=10 -pedantic -pedantic-errors -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -fbounds-check

all: xwave xdiffusion

xwave: $O/wave.o $O/gll_library.o $O/lagrange_poly.o $O/source_time_function.o $O/define_derivative_matrix.o
	${FC} $(FCFLAGS) -o $@ $+

xdiffusion: $O/diffusion.o $O/gll_library.o $O/lagrange_poly.o $O/define_derivative_matrix.o
	${FC} $(FCFLAGS) -o $@ $+

$O/%.o: %.f90 constants.h
	${FC} $(FCFLAGS) -c -o $@ $<

.PHONY: clean purge
clean:
	rm -f $O/*.o xwave xdiffusion

purge: clean
	rm -f snapshot* seismogram
