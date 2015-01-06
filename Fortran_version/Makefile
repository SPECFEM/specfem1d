##############################################################################
#
# Makefile for SPECFEM1D
# This makefile has been created by Alexis Bottero
#      (alexis DOT bottero aT gmail DOT com)
# make : create the executable bin/xspecfem1d
# make clean : remove the object files + module files
# make purge : # make clean + remove executable + remove OUTPUT_FILES
##############################################################################

#************ COMPILER OPTIONS **************#

DEBUG=no        # Debug mode or not
EXEC=xspecfem1d # Executable name

# Intel
#FC = ifort                    # Fortran Compiler
#FCFLAGS=#-W -Wall -Wextra -pedantic  # Compiler flags
#FCFLAGS_DEBUG = -O1 -vec-report0 -std95 -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage -check all -align sequence -assume byterecl -ftz -traceback -ftrapuv    # Compiler debug flags
#LDFLAGS = #Linker flags
#MOD_OPTION = -module #Specifies the directory <dir> where module (.mod) files should be placed when created and where they should be searched for

# GNU gfortran
FC = gfortran # Fortran Compiler
FCFLAGS=#-W -Wall -Wextra -pedantic # Compiler flags
FLAGS_DEBUG = -std=gnu -fimplicit-none -frange-check -O2 -fmax-errors=10 -pedantic -pedantic-errors -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -fbounds-check # Compiler debug flags
LDFLAGS = #Linker flags
MOD_OPTION = -J #Specifies the directory <dir> where module (.mod) files should be placed when created and where they should be searched for


#******* DO NOT CHANGE ANYTHING BELOW *******#

O := obj
S := src
B := bin

ifeq ($(DEBUG),yes)
	FCFLAGS := $(FCFLAGS_DEBUG)
endif

OBJ = \
	$O/modules.o \
	$O/define_derivative_matrix.o \
	$O/gaussm3.o \
	$O/gll_library.o \
	$O/lagrange_poly.o \
	$O/source_time_function.o \
	$O/specfem1d.o

MODS = \
	gaussm3 \
	globalstiffmat \
	grid \
	othervariables

# Default target
all: $B/$(EXEC)

# Module dependencies
$O/modules.o: $O/gaussm3.o
$O/define_derivative_matrix.o: $O/modules.o
$O/gll_library.o: $O/modules.o
$O/lagrange_poly.o: $O/modules.o
$O/source_time_function.o: $O/modules.o
$O/specfem1d.o: $O/modules.o

# To make the executable
$B/$(EXEC): $(OBJ)
	$(FC) -o $@ $(LDFLAGS) $+

$O/%.o: $S/%.f90
	$(FC) $(MOD_OPTION)$O -I$O -o $@ $(FCFLAGS) -c $<

# .PHONY is useful for example if a file called "clean" exists in the current directory -> "make clean" would not work
.PHONY: clean purge

# Clean remove all .o (if you don't want to see the messages use : @rm -rf *.o -> only the error messages will appear)
clean:
	rm -f $(OBJ) $(addprefix $O/, $(addsuffix .mod, $(MODS))) $B/$(EXEC)

# Clean and remove executable
purge: clean
	rm -f OUTPUT_FILES/* $S/*~ *~ _________plot_all_snapshots.gnu

