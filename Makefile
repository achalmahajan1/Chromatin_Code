include ../../MakefileInc.mk

INCLUDE_DIRS = $(USERINCLUDE) -I../../
LIBRARY_DIRS =
LIBRARIES =
# Files

SRC =   ActiveForces.f90 BoundaryInt.f90 discrete.f90 factors.f90 gasdev.f90 gauss_trgl.f90 Geo_val.f90 Initconfig.f90 Initdip1.f90 Initdip.f90 Interp.f90 Leonard_jones.f90 SpringForce.f90 sdlp_3d_interp.f90 Ext_mod2a.f90 ../../Interface/FMMWrapper-f.f90 ../../FMM/FMMWrapper.cpp ../../Interface/FMMWrapper-c.cpp 
INC =   ../../FMM/FMMWrapper.h ../../Interface/FMMWrapper-c.h 

# Definitions
CXX = mpicxx
# enable link time optimization
CXXFLAGS += -flto -g # -DFMMTIMING #-check bounds -wrap-margin-

FC = mpif90
FFLAGS = $(subst -std=c++14, ,$(CXXFLAGS)) #-std03 # fortran 2003 in intel compiler

LINK := mpif90
LINKFLAGS := $(FFLAGS) $(LINKFLAGS) # -cxxlib # link with cxx runtime in intel compiler

EXE = Run1
OBJ = ../../Interface/FMMWrapper-f.o ActiveForces.o BoundaryInt.o discrete.o factors.o gasdev.o sdlp_3d_interp.o gauss_trgl.o Geo_val.o Initconfig.o Initdip1.o  Initdip.o Interp.o SpringForce.o Leonard_jones.o Ext_mod2a.o ../../FMM/FMMWrapper.o ../../Interface/FMMWrapper-c.o 

# Link rule
$(EXE):	$(OBJ) $(SRC)
	$(LINK) $(OBJ)  -o $(EXE) $(LINKFLAGS) $(SYSLIB) $(LIBRARY_DIRS) $(LIBRARIES) 

# Compilation rules
../../Interface/FMMWrapper-f.o: ../../Interface/FMMWrapper-f.f90
	$(FC) $(FFLAGS) $(INCLUDE_DIRS) -c ../../Interface/FMMWrapper-f.f90 -o ../../Interface/FMMWrapper-f.o
ActiveForces.o: ../../Interface/FMMWrapper-f.o ActiveForces.f90
	$(FC) $(FFLAGS) $(INCLUDE_DIRS) -c ActiveForces.f90 -o ActiveForces.o
BoundaryInt.o: ../../Interface/FMMWrapper-f.o BoundaryInt.f90
	$(FC) $(FFLAGS) $(INCLUDE_DIRS) -c BoundaryInt.f90 -o BoundaryInt.o	
discrete.o: ../../Interface/FMMWrapper-f.o discrete.f90
	$(FC) $(FFLAGS) $(INCLUDE_DIRS) -c discrete.f90 -o discrete.o	
factors.o: ../../Interface/FMMWrapper-f.o factors.f90
	$(FC) $(FFLAGS) $(INCLUDE_DIRS) -c factors.f90 -o factors.o
gasdev.o: ../../Interface/FMMWrapper-f.o gasdev.f90
	$(FC) $(FFLAGS) $(INCLUDE_DIRS) -c gasdev.f90 -o gasdev.o
gauss_trgl.o: ../../Interface/FMMWrapper-f.o gauss_trgl.f90
	$(FC) $(FFLAGS) $(INCLUDE_DIRS) -c gauss_trgl.f90 -o gauss_trgl.o
Geo_val.o: ../../Interface/FMMWrapper-f.o Geo_val.f90
	$(FC) $(FFLAGS) $(INCLUDE_DIRS) -c Geo_val.f90 -o Geo_val.o
Initconfig.o: ../../Interface/FMMWrapper-f.o Initconfig.f90
	$(FC) $(FFLAGS) $(INCLUDE_DIRS) -c Initconfig.f90 -o Initconfig.o
Initdip.o: ../../Interface/FMMWrapper-f.o Initdip.f90
	$(FC) $(FFLAGS) $(INCLUDE_DIRS) -c Initdip.f90 -o Initdip.o
Initdip1.o: ../../Interface/FMMWrapper-f.o Initdip1.f90
	$(FC) $(FFLAGS) $(INCLUDE_DIRS) -c Initdip1.f90 -o Initdip1.o
Interp.o: ../../Interface/FMMWrapper-f.o Interp.f90
	$(FC) $(FFLAGS) $(INCLUDE_DIRS) -c Interp.f90 -o Interp.o
Leonard_jones.o: ../../Interface/FMMWrapper-f.o Leonard_jones.f90
	$(FC) $(FFLAGS) $(INCLUDE_DIRS) -c Leonard_jones.f90 -o Leonard_jones.o
sdlp_3d_interp.o: ../../Interface/FMMWrapper-f.o sdlp_3d_interp.f90
	$(FC) $(FFLAGS) $(INCLUDE_DIRS) -c sdlp_3d_interp.f90 -o sdlp_3d_interp.o
SpringForce.o: ../../Interface/FMMWrapper-f.o SpringForce.f90
	$(FC) $(FFLAGS) $(INCLUDE_DIRS) -c SpringForce.f90 -o SpringForce.o	
Ext_mod2a.o: ../../Interface/FMMWrapper-f.o Ext_mod2a.f90
	$(FC) $(FFLAGS) $(INCLUDE_DIRS) -c Ext_mod2a.f90 -o Ext_mod2a.o												
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDE_DIRS) -c $*.cpp -o $*.o

# Individual dependencies
$(OBJ): $(INC)

all: $(EXE)

clean: 
	rm -f ./$(OBJ)
	rm -f ./*.mod
	rm -f ./$(EXE)
