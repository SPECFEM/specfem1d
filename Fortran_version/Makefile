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

DEBUG=no                   # Debug mode or not
EXEC=xspecfem1d # Executable name

# Intel
#F90 = ifort                    # Fortran Compiler 
#FLAGS=#-W -Wall -Wextra -pedantic  # Compiler flags
#FLAGS_DEBUG = -O1 -vec-report0 -std95 -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments -warn ignore_loc -warn usage -check all -align sequence -assume byterecl -ftz -traceback -ftrapuv    # Compiler debug flags
#LDFLAGS = #Linker flags
#MOD_OPTION = -module #Specifies the directory <dir> where module (.mod) files should be placed when created and where they should be searched for

# GNU gfortran
F90 = gfortran # Fortran Compiler 
FLAGS=#-W -Wall -Wextra -pedantic # Compiler flags
FLAGS_DEBUG = -std=gnu -fimplicit-none -frange-check -O2 -fmax-errors=10 -pedantic -pedantic-errors -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -fbounds-check # Compiler debug flags
LDFLAGS = #Linker flags
MOD_OPTION = -J #Specifies the directory <dir> where module (.mod) files should be placed when created and where they should be searched for


#******* DO NOT CHANGE ANYTHING BELOW *******#

OBJ_DIR_NAME := obj
SRC_DIR_NAME := src
BIN_DIR_NAME := bin

PWD :=  $(shell pwd)
OBJ_DIR := $(PWD)/$(OBJ_DIR_NAME)
SRC_DIR := $(PWD)/$(SRC_DIR_NAME)
BIN_DIR := $(PWD)/$(BIN_DIR_NAME)

ifeq ($(DEBUG),yes)
	FLAGS := $(FLAGS_DEBUG)
endif

SRC = $(wildcard $(SRC_DIR)/*.f90)

# List of object files :
NAME := $(basename $(notdir $(SRC)))
OBJ := $(addprefix $(OBJ_DIR)/, $(addsuffix .o, $(NAME)))

# To make the executable
all: make_directory $(OBJ)
	$(F90) -o $(BIN_DIR)/$(EXEC) $(OBJ) $(LDFLAGS)

# Rule to make directory
make_directory: dirobj bin

$(OBJ_DIR)/gaussm3.o: $(SRC_DIR)/gaussm3.f90
	$(F90) $(MOD_OPTION)$(OBJ_DIR) -I$(OBJ_DIR) -c $< -o $@ $(FLAGS)
$(OBJ_DIR)/modules.o: $(SRC_DIR)/modules.f90
	$(F90) $(MOD_OPTION)$(OBJ_DIR) -I$(OBJ_DIR) -c $< -o $@ $(FLAGS)
$(OBJ_DIR)/define_derivative_matrix.o: $(SRC_DIR)/define_derivative_matrix.f90
	$(F90) $(MOD_OPTION)$(OBJ_DIR) -I$(OBJ_DIR) -c $< -o $@ $(FLAGS)
$(OBJ_DIR)/gll_library.o: $(SRC_DIR)/gll_library.f90
	$(F90) $(MOD_OPTION)$(OBJ_DIR) -I$(OBJ_DIR) -c $< -o $@ $(FLAGS)	
$(OBJ_DIR)/lagrange_poly.o: $(SRC_DIR)/lagrange_poly.f90
	$(F90) $(MOD_OPTION)$(OBJ_DIR) -I$(OBJ_DIR) -c $< -o $@ $(FLAGS)	
$(OBJ_DIR)/source_time_function.o: $(SRC_DIR)/source_time_function.f90
	$(F90) $(MOD_OPTION)$(OBJ_DIR) -I$(OBJ_DIR) -c $< -o $@ $(FLAGS)		
$(OBJ_DIR)/specfem1d.o: sleep_a_bit $(SRC_DIR)/specfem1d.f90
	$(F90) $(MOD_OPTION)$(OBJ_DIR) -I$(OBJ_DIR) -c $(SRC_DIR)/specfem1d.f90 -o $@ $(FLAGS)		

# This was added to solve a Makefile/Fortran bug... Indeed the creation of modules takes time and that create errors :
# use gaussm3
#       1
# Fatal Error: Can't open module file 'gaussm3.mod' for reading at (1): No such file or directory
# Even if at the end we see gaussm3.mod in obj/
# See : http://stackoverflow.com/questions/27268809/makefile-make-fails-but-make-d-and-make-n-hand-execution-works

sleep_a_bit:
	@sleep 0.5 # must be >0.005

# .PHONY is usefull for example if it exists a file named "clean" is the current directory -> "make clean" would not work
.PHONY: dirobj bin

# Create the directory $(OBJ_DIR) if needed :
ifeq ($(strip $( $(wildcard $(OBJ_DIR)) ) ), )
dirobj: 
	@mkdir -p $(OBJ_DIR)
else
dirobj: 
endif

# Create the directory $(BIN_DIR) if needed :
ifeq ($(strip $( $(wildcard $(BIN_DIR)) ) ), )
bin: 
	@mkdir -p $(BIN_DIR)
else
bin: 
endif 

.PHONY: clean purge 

# Clean remove all .o (if you don't want to see the messages use : @rm -rf *.o -> only the error messages will appear)
clean:
	rm -rf $(OBJ_DIR) $(SRC_DIR)/*~ *~ .plot_all_snapshots.gnu utils/*~ $(BIN_DIR)/*

# Clean + remove executable
purge: clean
	rm -rf $(BIN_DIR) OUTPUT_FILES
