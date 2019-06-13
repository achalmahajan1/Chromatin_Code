SRC = Variables.f90 ActiveForces.f90 discrete.f90 factors.f90 gasdev.f90 Geo_val.f90 gauss_trgl.f90 Greens_fs.f90 Initconfig.f90 Spring.f90 Initdip.f90 Interp.f90 sgf_3d_fs.f90 sdlp_3d_interp.f90 Leonard_jones.f90 BoundaryInt.f90 Ext_mod2a.f90

OBJ = Variables.o ActiveForces.o discrete.o factors.o gasdev.o Geo_val.o gauss_trgl.o Greens_fs.o Initconfig.o Spring.o Initdip.o Interp.o sgf_3d_fs.o sdlp_3d_interp.o Leonard_jones.o BoundaryInt.o Ext_mod2a.o

EXE = Run1
LINK = gfortran
LINKFLAGS = -O3 #Link rule

# Link rule
$(EXE):	$(OBJ) $(SRC)
	$(LINK) $(OBJ) -o $(EXE) $(LINKFLAGS)  

# Compilation rules	
Variables.o: 
	$(LINK) $(LINKFLAGS) -c Variables.f90 -o Variables.o
ActiveForces.o: 
	$(LINK) $(LINKFLAGS) -c ActiveForces.f90 -o ActiveForces.o
BoundaryInt.o: 
	$(LINK) $(LINKFLAGS) -c BoundaryInt.f90 -o BoundaryInt.o	
discrete.o: 
	$(LINK) $(LINKFLAGS) -c discrete.f90 -o discrete.o	
factors.o: 
	$(LINK) $(LINKFLAGS) -c factors.f90 -o factors.o
gasdev.o: 
	$(LINK) $(LINKFLAGS) -c gasdev.f90 -o gasdev.o
gauss_trgl.o: 
	$(LINK) $(LINKFLAGS) -c gauss_trgl.f90 -o gauss_trgl.o
Geo_val.o: 
	$(LINK) $(LINKFLAGS) -c Geo_val.f90 -o Geo_val.o
Greens_fs.o: 
	$(LINK) $(LINKFLAGS) -c Greens_fs.f90 -o Greens_fs.o
Initconfig.o: 
	$(LINK) $(LINKFLAGS) -c Initconfig.f90 -o Initconfig.o
Initdip.o:
	$(LINK) $(LINKFLAGS) -c Initdip.f90 -o Initdip.o
Interp.o: 
	$(LINK) $(LINKFLAGS) -c Interp.f90 -o Interp.o
Leonard_jones.o: 
	$(LINK) $(LINKFLAGS) -c Leonard_jones.f90 -o Leonard_jones.o
sdlp_3d_interp.o: 
	$(LINK) $(LINKFLAGS) -c sdlp_3d_interp.f90 -o sdlp_3d_interp.o
sgf_3d_fs.o: 
	$(LINK) $(LINKFLAGS) -c sgf_3d_fs.f90 -o sgf_3d_fs.o
Spring.o: 
	$(LINK) $(LINKFLAGS) -c Spring.f90 -o Spring.o
Ext_mod2a.o: 
	$(LINK) $(LINKFLAGS) -c Ext_mod2a.f90 -o Ext_mod2a.o												

	
clean: 
	rm ./$(OBJ)
	rm ./*.mod
	rm ./$(EXE)
	rm ./*.dat	
