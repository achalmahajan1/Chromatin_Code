gfortran -c Variables.f90 -o Variables.o
gfortran -c ActiveForces.f90 -o ActiveForces.o
gfortran -c discrete.f90 -o discrete.o
gfortran -c factors.f90 -o factors.o
gfortran -c gasdev.f90 -o gasdev.o
gfortran -c Geo_val.f90 -o Geo_val.o
gfortran -c gauss_trgl.f90 -o gauss_trgl.o
gfortran -c Greens_fs.f90 -o Greens_fs.o
gfortran -c Initconfig.f90 -o Initconfig.o
gfortran -c Spring.f90 -o Spring.o
gfortran -c Initdip.f90 -o Initdip.o
gfortran -c Interp.f90 -o Interp.o
gfortran -c sgf_3d_fs.f90 -o sgf_3d_fs.o
gfortran -c sdlp_3d_interp.f90 -o  sdlp_3d_interp.o
gfortran -c Leonard_jones.f90 -o Leonard_jones.o
gfortran -c BoundaryInt.f90 -o BoundaryInt.o
gfortran -c Ext_mod2a.f90 -o Ext_mod2a.o
gfortran -O3 -o Nu4 Variables.o ActiveForces.o discrete.o Spring.o factors.o gasdev.o gauss_trgl.o Geo_val.o Greens_fs.o Initconfig.o Initdip.o Interp.o Leonard_jones.o BoundaryInt.o sdlp_3d_interp.o sgf_3d_fs.o Ext_mod2a.o
