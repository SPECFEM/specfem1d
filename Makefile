O = obj

# Dec Alpha
#F90 = f90
#FLAGS = -fast -warn truncated_source -warn argument_checking -warn unused -warn declarations -std95 -check nounderflow -check bounds

# Portland compiler
#F90 = pgf90
#FLAGS = -fast -Mbounds -Mneginfo -Mdclchk -Mstandard -Knoieee
#FLAGS = -fast -Mnobounds -Mneginfo -Mdclchk -Munroll=c:6 -Mstandard -Knoieee

# Absoft compiler
#F90 = f90
#FLAGS = -O2 -W132 -YEXT_NAMES=LCS  -s -B108 -YCFRL=1

# Intel Linux compiler
F90 = ifort
FLAGS = -O3 -e95 -implicitnone

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

