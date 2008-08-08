
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
#F90 = pgf90
#FLAGS =-fast -Mnobounds -Minline -Mneginfo -Mdclchk -Knoieee -Minform=warn -fastsse -tp amd64e -Msmart

# Intel
#F90 = ifort
#FLAGS =-O3 -xP -vec-report0 -e03 -std03 -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage -check nobounds -align sequence -assume byterecl -fpe3 -ftz

# GNU gfortran
F90 = gfortran
FLAGS  = -std=f2003 -fimplicit-none -frange-check -O2 -Wunused-labels -Waliasing -Wampersand -Wsurprising -Wline-truncation -Wunderflow

wave: constants.h \
       $O/wave.o \
       $O/gll_library.o \
       $O/lagrange_poly.o \
       $O/source_time_function.o \
       $O/define_derivative_matrix.o
	${F90} $(FLAGS) -o xwave \
       $O/wave.o \
       $O/gll_library.o \
       $O/lagrange_poly.o \
       $O/source_time_function.o \
       $O/define_derivative_matrix.o

diffusion: constants.h \
       $O/diffusion.o \
       $O/gll_library.o \
       $O/lagrange_poly.o \
       $O/define_derivative_matrix.o
	${F90} $(FLAGS) -o xdiffusion \
       $O/diffusion.o \
       $O/gll_library.o \
       $O/lagrange_poly.o \
       $O/define_derivative_matrix.o

bak:
	cp *f90 *h Makefile bak

clean:
	rm -f $O/*.o *.o snapshot* seismogram xwave xdiffusion

all: clean wave diffusion

$O/wave.o: constants.h wave.f90
	${F90} $(FLAGS) -c -o $O/wave.o wave.f90

$O/diffusion.o: constants.h diffusion.f90
	${F90} $(FLAGS) -c -o $O/diffusion.o diffusion.f90

$O/define_derivative_matrix.o: constants.h define_derivative_matrix.f90
	${F90} $(FLAGS) -c -o $O/define_derivative_matrix.o define_derivative_matrix.f90

$O/gll_library.o: gll_library.f90
	${F90} $(FLAGS) -c -o $O/gll_library.o gll_library.f90

$O/lagrange_poly.o: lagrange_poly.f90
	${F90} $(FLAGS) -c -o $O/lagrange_poly.o lagrange_poly.f90

$O/source_time_function.o: source_time_function.f90
	${F90} $(FLAGS) -c -o $O/source_time_function.o source_time_function.f90

