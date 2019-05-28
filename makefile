#Set the compilers and linker
FC=mpif90
CC=mpicc
LD=mpif90

#Set different compiler options for when at home
home:	FC=gfortran
home: CC=gcc
home: LD=gfortran



#Set the objects
OBJS=																PETScVector.o															\
																		PETScMatrix.o															\
																		PETScNumMethods.o													\
																		PETScUseAll.o	     												\
																		PetTest.o


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fortran Compiler options
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MY_OPTIONS    =           -ffree-line-length-800 -o3
FC_FLAGS    =           -fcheck=bounds -ffree-line-length-800 -pg -g
FC_FLAGS_DEBUG = -fimplicit-none -fbounds-check -ffree-line-length-0 -fcheck=all -fbacktrace -g -DDEBUG
FC_FLAGS_HOME = -fbounds-check -ffree-line-length-0 -pg -DHOME

debug:	FC_FLAGS = $(FC_FLAGS_DEBUG)
home:		FC_FLAGS = $(FC_FLAGS_HOME)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# LINKER - FORTRAN LINKER (MPI w/ GNUCOMPILER gfortran)
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

LD_FLAGS = -O2

LD_FLAGS_DEBUG = -g -Og

debug:	LD_FLAGS = $(LD_FLAGS_DEBUG)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# PETSC (STATIC & DYNAMIC) LIBRARY
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

debug: PETSC_DIR = $(PETSC_DIR_DEBUG)

PETSC_INCLUDE_PATH = -I$(PETSC_DIR)/include/

PETSC_LINK_PATH = -L$(PETSC_DIR)/lib/

#--dev
PETSC_LOAD_PATH = -Wl,-rpath=$(PETSC_DIR)/lib/
PETSC_LOAD_PATH_HOME= -L$(PETSC_DIR)/lib/

# home:	PETSC_LOAD_PATH = $(PETSC_LOAD_PATH_HOME)

PETSC_LIB = -lpetsc

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# PATHS
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

INCLUDE_PATH = $(PETSC_INCLUDE_PATH)
LINK_PATH = $(PETSC_LINK_PATH)
LOAD_PATH = $(PETSC_LOAD_PATH)
LIBS = $(PETSC_LIB)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEFINE MODS: OBJS .o REPLACED BY .mod
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
MODS= $(OBJS:.o=.mod)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEFINE EXECUTABLE
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
EXEC=pet

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAKEFILE VARIABLE
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
DEFAULT=makefile

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAKE EXECUTABLE
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
all : $(EXEC)
debug : $(EXEC)
home: $(EXEC)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# LINK THE SOURCE (.o) INTO THE TARGET (EXEC) - explicit rule
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
$(EXEC): $(OBJS)
	$(LD) $(FC_FLAGS) $(LINK_PATH) -o $@ $^ $(LOAD_PATH) $(LIBS)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# COMPILE INFERRED SOURCE (.f90) INTO TARGET (.o) - inference rule
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%.o : %.f90
	$(FC) $(FC_FLAGS) -cpp -dM -c $< -o $@ $(INCLUDE_PATH) $(LINK_PATH) $(LIBS)

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# UPDATE OBJS
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%.f90 : $(DEFAULT)

#Clean the directory and any directories searched
.PHONY : clean
clean :
	rm -rf $(OBJS) $(EXEC) $(MODS) *.o *.mod
