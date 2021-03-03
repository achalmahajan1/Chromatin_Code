PROGRAM Bead_Rod
        use fmmwrapper
        use iso_c_binding
        use mpi
        use omp_lib
	IMPLICIT NONE
	        INTEGER*4 :: rank,size,ierror,tag,status(MPI_STATUS_SIZE)
		INTEGER*8 :: N,i,j,Niter,k,m,l,npts,nn(5120,6),Nelm,i1,i2,i3,i4,i5,i6,Ntr,ind,filenum,&
		                idum,Nc,ierr,Ns,Nsc,j1,jj1,Nsb,Np,Np1  
		INTEGER*8, ALLOCATABLE, DIMENSION(:) :: jj
		INTEGER*8, ALLOCATABLE, DIMENSION(:,:) :: dip,dipth,tchk,patch
		REAL*8 :: del_t,T_end,pi,Rad_bead,tol,v1,v2,Rodl,Rs,dr,time,gasdev,xcm1,ycm1,zcm1,scale2,&
		          kon,koff,p(10242,3),start,finish,boa,coa,oot,scale1,x_axis,y_axis,z_axis,aoc,&
		          boc,Bpar,Hc,norm,Qbar,dr1,kton,ktoff,num,Fx,Fy,Fz,big
		REAL*8, ALLOCATABLE, DIMENSION(:) :: Fdp1,F_brown,noise,Bcho,x_temp,z_temp,y_temp,RM2,leng,&
		                                     Fdpx1,Fdpy1,Fdpz1,Fdpxa1,Fdpya1,Fdpza1,Fbx1,ux,uy,uz,&
		                                     Fljx1,Fljy1,Fljz1,Fby1,Fbz1,alpha,beta,gamma1,u_1
		REAL*8, ALLOCATABLE, DIMENSION(:,:) :: x,y,z,Fdpx,Fdpy,Fdpz,trel,tswitch,F_brownx,F_browny,F_brownz,&
		                                       Qx,Qy,Qz,Qo,Fbx,Fby,Fbz,Fdpxa,Fdpya,Fdpza,x_tr,tsth,&
		                                       Fljx,Fljy,Fljz,u_tr,Fdp,vna,RM1,x_temp1,y_temp1,z_temp1,&
						       trel1,xc,yc,zc,Fqx,Fqy,Fqz
                CHARACTER(LEN =30) :: filename                                       
		CHARACTER(LEN=38) :: path='/mnt/home/amahajan/ceph/PhaseSeg1b_18/'
                TYPE(c_ptr):: fmm                                                                            
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%VARIABLE ASSIGNMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		pi = 4.d0*atan(1.d0)
		!PRINT*,"Enter the number of beads to be simulated"
		!READ(*,*), N
		N = 1305  !No of beads in each chain
		Nc = 23   !No chains
		Ntr = 0  !No of tracer particles
                Nsb = 4!Number of patches for the crosslinks that can stay on
                Np = 90!Number of beads stays on in the patch (or size of a single patch)
                Np1 = N-5-Nsb*Np!Number of beads that cannot form Crosslinks (onlydipoles)
		ALLOCATE(x(N,Nc))
		ALLOCATE(y(N,Nc))
		ALLOCATE(z(N,Nc))
		ALLOCATE(xc(N,Nc))
                ALLOCATE(yc(N,Nc))
                ALLOCATE(zc(N,Nc))
		ALLOCATE(jj(N*Nc))
		ALLOCATE(x_temp1(N,Nc))
		ALLOCATE(y_temp1(N,Nc))
		ALLOCATE(z_temp1(N,Nc))
		ALLOCATE(x_temp(N*Nc+Ntr))
		ALLOCATE(y_temp(N*Nc+Ntr))
		ALLOCATE(z_temp(N*Nc+Ntr))
		ALLOCATE(Fljx(N,Nc))
		ALLOCATE(Fljy(N,Nc))
		ALLOCATE(Fljz(N,Nc))
		ALLOCATE(Fljx1(N*Nc))
		ALLOCATE(Fljy1(N*Nc))
		ALLOCATE(Fljz1(N*Nc))
		ALLOCATE(Fdpx(N,Nc))
		ALLOCATE(Fdpy(N,Nc))
		ALLOCATE(Fdpz(N,Nc))
		ALLOCATE(Fdpxa(N,Nc))
		ALLOCATE(Fdpya(N,Nc))
		ALLOCATE(Fdpza(N,Nc))
		ALLOCATE(Fdpx1(N*Nc))
		ALLOCATE(Fdpy1(N*Nc))
		ALLOCATE(Fdpz1(N*Nc))
		ALLOCATE(Fdpxa1(N*Nc))
		ALLOCATE(Fdpya1(N*Nc))
		ALLOCATE(Fdpza1(N*Nc))
		ALLOCATE(F_brownx(N,Nc))
		ALLOCATE(F_browny(N,Nc))
		ALLOCATE(F_brownz(N,Nc))
                ALLOCATE(Fqx(N,Nc))
                ALLOCATE(Fqy(N,Nc))
                ALLOCATE(Fqz(N,Nc))
		!ALLOCATE(F_brown(N,3))
		ALLOCATE(trel(N,Nc))
		ALLOCATE(trel1(N,Nc))
                ALLOCATE(tchk(N,Nc))
                ALLOCATE(tswitch(N,Nc))
                ALLOCATE(tsth(N,Nc))
                ALLOCATE(dip(N,Nc))
                ALLOCATE(dipth(N,Nc))
		ALLOCATE(x_tr(Ntr,3))
		ALLOCATE(u_tr(Ntr,3))
		ALLOCATE(Qx(N-1,Nc))
		ALLOCATE(Qy(N-1,Nc))
		ALLOCATE(Qz(N-1,Nc))
		ALLOCATE(Qo(N-1,Nc))
		ALLOCATE(Fbx1(N*Nc))
		ALLOCATE(Fby1(N*Nc))
		ALLOCATE(Fbz1(N*Nc))
		ALLOCATE(Fbx(N,Nc))
		ALLOCATE(Fby(N,Nc))
		ALLOCATE(Fbz(N,Nc))
		ALLOCATE(ux(N*Nc))
		ALLOCATE(uy(N*Nc))
		ALLOCATE(uz(N*Nc))
		ALLOCATE(Fdp(N*Nc,3))
		ALLOCATE(Fdp1(N*Nc*3))
		!ALLOCATE(u_1(N*Nc,3))
		ALLOCATE(patch(N,Nc))
                ALLOCATE(u_1(N*Nc*3))
		ALLOCATE(leng(Nc))
		call MPI_INIT(ierror)
                call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror)
                call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
                call MPI_BARRIER(MPI_COMM_WORLD, ierror); 
                fmm = create_fmm_wrapper(8, 5000, 0, 0, 0)
                call FMM_SetBox(fmm, -40.0d+0, 40.0d+0, -40.0d+0, 40.0d+0, -35.0d+0, 35.d+0);
		u_1 = 0.d0
		Rs = 28.d0 !(Effective radius of the Nucleus; see notes for scales)
		T_end = 200.d0
		dr = Rs*0.05d0
		time = 0.d0
		Fljx = 0.d0
		Fljy = 0.d0
		Fljz = 0.d0
		Fljx1 = 0.d0
		Fljy1 = 0.d0
		Fljz1 = 0.d0
		Fdpx = 0.d0
		Fdpy = 0.d0
		Fdpz = 0.d0
		Fbx1 = 0.d0
		Fby1 = 0.d0
		Fbz1 = 0.d0
		dr1 = 0.2d0
		Fbx = 0.d0
		Fby = 0.d0
		Fbz = 0.d0
		Hc = 3.d0
		Fdpxa = 0.d0
		Fdpya = 0.d0
		Fdpza = 0.d0
                Bpar = 6.d0
		Rodl = 1.d0 !(Length of the rod scaled with Sphere radius)
                kon = 2.d0 ! On/off rates for dipoles
                koff = 5.d0
		kton = 0.02d0
                ktoff = 0.1d0
		Rad_Bead = 0.1d0 !(Radius of the bead scaled with Sphere radius)
		ind = 1.d0
                ux = 0.d0
                uy = 0.d0
                uz = 0.d0
		idum = -140
		Nelm = 20
                !F_brown = 0.d0
		del_t = 0.000005d0
                big = 100000000.d0
		WRITE(*,*),"Radius of the bead:",Rad_Bead
                WRITE(*,*),"Number of bead per chains:",N
                WRITE(*,*),"Number of chains:",Nc
                WRITE(*,*),"Initial bond length:",Rodl
                WRITE(*,*),"Radius of the Spherical Cavity:",Rs
                WRITE(*,*),"Magnitude of FENE-p spring:",Hc                
                WRITE(*,*),"Off/On time for active forces chosen from Poisson Process:",kon,koff
                WRITE(*,*),"time step:",del_t
                WRITE(*,*),"End of the simulations:",T_end
		k = 0
		m = 0
		trel = 0.d0
		trel1 = 0.d0
                xc = 0.d0
                yc = 0.d0
                zc = 0.d0
                Fqx = 0.d0
                Fqy = 0.d0
                Fqz = 0.d0
		OPEN(UNIT=381, FILE='/mnt/home/amahajan/ceph/PhaseSeg1b_18/BEM.dat')
                DO l = 1,Nc
        		filenum = l
	        	WRITE(filename, '("x_chain",I2.2, ".dat")')filenum
	        	OPEN(FILE = path//filename,UNIT = 41 + l, FORM="unformatted")
	        	WRITE(filename, '("y_chain",I2.2, ".dat")')filenum
	        	OPEN(FILE = path//filename,UNIT = 71 + l, FORM="unformatted")
	        	WRITE(filename, '("z_chain",I2.2, ".dat")')filenum
	        	OPEN(FILE = path//filename,UNIT = 101 + l, FORM="unformatted")
	        	WRITE(filename, '("ForceTension_chain",I2.2, ".dat")')filenum
	        	OPEN(FILE = path//filename,UNIT = 131 + l, FORM="unformatted")
	        	WRITE(filename, '("ForceActive_chain",I2.2, ".dat")')filenum
	        	OPEN(FILE = path//filename,UNIT = 161 + l, FORM="unformatted")
	        	WRITE(filename, '("ForceLJ_chain",I2.2, ".dat")')filenum
	        	OPEN(FILE = path//filename,UNIT = 191 + l, FORM="unformatted")
		        WRITE(filename, '("ForceCross_chain",I2.2, ".dat")')filenum
                        OPEN(FILE = path//filename,UNIT = 221 + l, FORM="unformatted")
                ENDDO
                OPEN(UNIT=551, FILE='/mnt/home/amahajan/ceph/PhaseSeg1b_18/Velocity.dat', FORM="unformatted")
                OPEN(UNIT=24, FILE='/mnt/home/amahajan/ceph/PhaseSeg1b_18/Params.dat')
                WRITE(24,*),"Radius of the bead:",Rad_Bead
                WRITE(24,*),"Number of bead per chains:",N
                WRITE(24,*),"Number of chains:",Nc
                WRITE(24,*),"Initial bond length:",Rodl
                WRITE(24,*),"Radius of the Spherical Cavity:",Rs
                WRITE(24,*),"Magnitude of FENE-p spring:",Hc
                WRITE(24,*),"Off/On time for active forces chosen from Poisson Process:",kon,koff
                WRITE(24,*),"time step:",del_t
                WRITE(24,*),"End of the simulations:",T_end
                oot  = 1.0D0/3.0D0
                !boa = 1.d0 !For spheroid b/a = 1
                !coa = 0.8d0 !For oblate c < a and prolate c > a 
                aoc = 0.8d0
                boc = 0.8d0
                !scale1 = Rs/(boa*coa)**oot
                scale1 = Rs/(boc*aoc)**oot
                scale2 = (boc*aoc)**oot
                !x_axis = scale1        !Comment it for prolate 
                !y_axis = scale1*boa    !Comment it for prolate
                !z_axis = scale1*coa    !Comment it for prolate
                x_axis = scale1*aoc     !Comment it for oblate
                y_axis = scale1*boc     !Comment it for oblate
                z_axis = scale1         !Comment it for oblate
                CALL discrete(p,nn,npts,Nelm)
                DO i=1,npts                  ! scale
                        p(i,1) = x_axis*p(i,1)
                        p(i,2) = y_axis*p(i,2)
                        p(i,3) = z_axis*p(i,3)
                ENDDO

                ALLOCATE(RM1(Nelm*3,Nelm*3))
                ALLOCATE(RM2(9*Nelm*Nelm))
                ALLOCATE(alpha(Nelm))
                ALLOCATE(beta(Nelm))
                ALLOCATE(gamma1(Nelm))
                l = 0
                ierr = 1.d0
                !MAtrix for boundary integral evaluation \int G(x,x_0) dS saved
                !in a file and loaded before starting the simulations
                OPEN(32,FILE = "../SCTest/MatrixInverseG5120_Rs28_e08.txt", FORM="unformatted")
                !Need to read the initial position of the chains from the text
                !file. Chain is initialized as a random walk configuration. Each
                !chain is generated, and then placed manually inside the
                !elliposid domain. Eventually they are saved in
                !PosRWRs18_{01..23}.txt before being used in this code
                OPEN(261,FILE = "PosRWRs18_01.txt")
                OPEN(262,FILE = "PosRWRs18_02.txt")
                OPEN(263,FILE = "PosRWRs18_03.txt")
                OPEN(264,FILE = "PosRWRs18_04.txt")
                OPEN(265,FILE = "PosRWRs18_05.txt")
                OPEN(266,FILE = "PosRWRs18_06.txt")
                OPEN(267,FILE = "PosRWRs18_07.txt")
                OPEN(268,FILE = "PosRWRs18_08.txt")
                OPEN(269,FILE = "PosRWRs18_09.txt")
                OPEN(270,FILE = "PosRWRs18_10.txt")
                OPEN(271,FILE = "PosRWRs18_11.txt")
                OPEN(272,FILE = "PosRWRs18_12.txt")
                OPEN(273,FILE = "PosRWRs18_13.txt")
                OPEN(274,FILE = "PosRWRs18_14.txt")
                OPEN(275,FILE = "PosRWRs18_15.txt")
                OPEN(276,FILE = "PosRWRs18_16.txt")
                OPEN(277,FILE = "PosRWRs18_17.txt")
                OPEN(278,FILE = "PosRWRs18_18.txt")
                OPEN(279,FILE = "PosRWRs18_19.txt")
                OPEN(280,FILE = "PosRWRs18_20.txt")
                OPEN(281,FILE = "PosRWRs18_21.txt")
                OPEN(282,FILE = "PosRWRs18_22.txt")
                OPEN(283,FILE = "PosRWRs18_23.txt")
                CALL CPU_TIME(start)
                DO i = 1,Nelm*Nelm*9
                        READ(32,IOSTAT=ierr),RM2(i)! No need to read at every time step, just once at the beginning and can be saved in the RAM
                ENDDO
                DO i = 1,Nelm*3
                        DO j = 1,Nelm*3
                                l = l + 1
                                RM1(i,j) = RM2(l)
                        ENDDO
                ENDDO 
                DEALLOCATE(RM2)           
                DO i = 1,nelm
                        i1 = nn(i,1)
                        i2 = nn(i,2)
                        i3 = nn(i,3)
                        i4 = nn(i,4)
                        i5 = nn(i,5)
                        i6 = nn(i,6)
                        l = l + 1
                        WRITE(381,*),i,l,p(i1,1),p(i1,2),p(i1,3)
                        l = l + 1 
                        WRITE(381,*),i,l,p(i2,1),p(i2,2),p(i2,3)
                        l = l + 1
                        WRITE(381,*),i,l,p(i3,1),p(i3,2),p(i3,3) 
                        !PRINT*,i, p(i1,1),p(i1,2),p(i1,3)       
               ENDDO 
               !ALLOCATE(Bcho(3*N,3*N))
	       !ALLOCATE(noise(3*N,1))
               ALLOCATE(vna(npts,3))
! THIS gives the values of alpha, beta and gamma used           
                DO k=1,nelm                            
                        i1 = nn(k,1)
                        i2 = nn(k,2)
                        i3 = nn(k,3)
                        i4 = nn(k,4)
                        i5 = nn(k,5)
                        i6 = nn(k,6)
                        !PRINT*,"outside", p(i1,1),p(i1,2),p(i1,3)
                        CALL factors(p(i1,1),p(i1,2),p(i1,3),p(i2,1),p(i2,2),p(i2,3),p(i3,1),p(i3,2),&
                                p(i3,3),p(i4,1),p(i4,2),p(i4,3),p(i5,1),p(i5,2),p(i5,3),p(i6,1),&
                                p(i6,2),p(i6,3),alpha(k),beta(k),gamma1(k))
                                !PRINT*,"outside", DSQRT(p(i2,1)**2.d0 + p(i2,2)**2.d0 + p(i2,3)**2.d0)
                                !PAUSE
                ENDDO
                CALL Geo_val(nelm,npts,vna,nn,alpha,beta,gamma1,p)
                vna = -1.d0*vna
!******************************************INITIALIZE THE POSITION AND DIPOLES******************************************
			!CALL Initconfig(x,y,z,boa,coa,aoc,boc,Rs,Rodl,N,Nc)
                        
                        DO l = 1,Nc
			DO i = 1,N
			        READ(260+l,*),x(i,l),y(i,l),z(i,l)
			ENDDO
                        ENDDO
                        
                        CALL Initdip(trel,kon,koff,tswitch,dip,N,Nc)
                        CALL Initdip1(trel1,kton,ktoff,tsth,tswitch,dipth,dip,N-5,Nc,Nsb,Np,patch)
			DO j = 1,Nc
                                DO i = 1,N
                                        tchk(i,j) = 1.d0
                                        jj((j-1)*N+i) = (j-1)*N+i
                                ENDDO
                        ENDDO 

			DO l = 1,Nc
			        DO i=1,N
			                WRITE(41 + l),REAL(time),REAL(x(i,l)),INT(jj((l-1)*N+i))
			        ENDDO
			        DO i=1,N
			                WRITE(71 + l),REAL(time),REAL(y(i,l)),INT(jj((l-1)*N+i))
			        ENDDO
			        DO i=1,N
			                WRITE(101 + l),REAL(time),REAL(z(i,l)),INT(jj((l-1)*N+i))
			        ENDDO
			ENDDO
			DO j = 1,Nc
        		        leng(j) = 0.d0
		            	        DO i = 2,N
		                                leng(j) = leng(j) + DSQRT((x(i,j) - x(i-1,j))**2 + (y(i,j) &
		                                                - y(i-1,j))**2 + (z(i,j) - z(i-1,j))**2) 
		                        ENDDO 
		        ENDDO
		        PRINT*,"Length of",Nc,"chains",SUM(leng) 
			j = 1
!*********************BEGINNING OF HI (GREENS FUNTION FOR STOKESLET INSIDE A SPHERE- KIM and MAUL 1994 POF)*********		
			!DO i = 1,N
				!PRINT*, x(i), y(i), z(i), Rad_Bead
			!	x(i) = x(i)/Rs
			!	y(i) = y(i)/Rs
			!	z(i) = z(i)/Rs
			!ENDDO
			!Rad_Bead = Rad_Bead/Rs
			!CALL SphereGeo(x,y,z,Dxx, Dyy, Dzz, Dxy, Dyx, Dyz, Dzy, Dxz, Dzx, Rs,Rad_Bead,N)
			!Rad_Bead = Rad_Bead*Rs
			!DO i = 1,N
			!	x(i) = x(i)*Rs
			!	y(i) = y(i)*Rs
			!	z(i) = z(i)*Rs
				!PRINT*, i,Rs, x(i), y(i), z(i), Rad_Bead
			!ENDDO
!******************************END OF HI **********************************************		
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MAIN CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		DO WHILE(time .LE. T_end)
			m = m + 1
			time = del_t*m
	                x_temp = 0.d0
                        y_temp = 0.d0
                        z_temp = 0.d0
			
                        x_temp1 = x
                        y_temp1 = y
                        z_temp1 = z
			
                        DO i = 1,Nc
                                DO k = 1,N
					trel1(k,i) = trel1(k,i) + del_t
                                        trel(k,i) = trel(k,i) + del_t
                                !PRINT*,i,k,dipth(k,i),tsth(k,i),dip(k,i),tswitch(k,i)
                                ENDDO
                        ENDDO
                        
		        xcm1 = 0.d0
                        ycm1 = 0.d0
                        zcm1 = 0.d0  

			DO i = 1,Nc
                                DO j = 1,N
                                IF (trel1(j,i) .GT. tsth(j,i))THEN
                                        trel1(j,i) = 0.d0
                                         num = 0.d0
                                                IF(dipth(j,i) .EQ. 0)THEN
                                                        dipth(j,i) = 1.d0
                                                        DO WHILE(num .EQ. 0.d0)
                                                                num = rand(0)
                                                        ENDDO
                                                        tsth(j,i) = big!-log(num)/ktoff
                                                ELSE
                                                        DO WHILE(num .EQ. 0.d0)
                                                                num = rand(0)
                                                        ENDDO
                                                        jj(jj((i-1)*N+j)) = jj((i-1)*N+j)
                                                        dipth(j,i) = 0.d0
                                                        tsth(j,i) = -log(num)/kton
                                                        jj((i-1)*N+j) = (i-1)*N+j
                                               ENDIF
                                ENDIF
                                ENDDO
                        ENDDO


                        l = 0 
                        DO i = 1,Nc
                                DO j = 1,N
                                        l = l + 1
                                        x_temp(l) = x(j,i)
                                        y_temp(l) = y(j,i)
                                        z_temp(l) = z(j,i)
                                ENDDO
                        ENDDO                 			
                        l = 0
                        IF(time .GT. 10.d0)THEN
                        CALL ActiveForces(x,y,z,Fdpx,Fdpy,Fdpz,Fdpxa,Fdpya,Fdpza,N,Nc,kon,koff,dip,trel,tswitch) 			
                        ENDIF
                        CALL Leonard_jones(x,y,z,Fljx,Fljy,Fljz,jj,patch*dipth,N,Nc,Rad_Bead,Rodl,Rs,time)
			IF(MOD(m,1000) .EQ. 0.d0 .OR. m .EQ. 100)THEN !DON'T USE m .EQ.1 because if the beads are closer in initial config can lead to diverging
                        !Start with deterministic forces
                        DO i = 1,Nc
                                DO j = 1,N
                                l = l + 1
                                Fdp1(3*(l - 1) + 1) = Fdpx(j,i) + Fdpxa(j,i) + Fljx(j,i) + Fqx(j,i) + Fbx(j,i)
                                Fdp1(3*(l - 1) + 2) = Fdpy(j,i) + Fdpya(j,i) + Fljy(j,i) + Fqy(j,i) + Fby(j,i)
                                Fdp1(3*(l - 1) + 3) = Fdpz(j,i) + Fdpza(j,i) + Fljz(j,i) + Fqz(j,i) + Fbz(j,i)
                                ENDDO
                        ENDDO
                        !SUBROUTINE FOR INDUCED VELOCITY USING BOUNDARY INTEGRAL AND FAST MUTLIPOLE (Can be used as a black box)
                        u_1 = 0.d0
		        CALL BoundaryInt(Fdp1,x_temp,y_temp,z_temp,u_1,p,nn,Nelm,alpha,beta,gamma1,RM1,fmm,N*Nc)
			ENDIF 
			Fqx = 0.d0
                        Fqy = 0.d0
                        Fqz = 0.d0
!!!!%%%%%%%%%%%%%%%%%%%IN CASE TO MOVE WITH THE CENTRE OF MASS - UNCOMMENT THIS!!%%%%%%%%%%%%%%%%%%%%%%%%			      
			!IF(ANY(dip == 1))THEN
			!        DO i = 1,N	
        	        !	        xcm1 = x_temp(i) + xcm1
        		!                ycm1 = y_temp(i) + ycm1
        		!                zcm1 = z_temp(i) + zcm1
                        !        ENDDO
                        !        xcm1 = xcm1/N
                        !        ycm1 = ycm1/N
                        !        zcm1 = zcm1/N
                        !        DO i = 1,N
                        !                x_temp(i) = x_temp(i) - xcm1
                        !                y_temp(i) = y_temp(i) - ycm1
                        !                z_temp(i) = z_temp(i) - zcm1
                        !        ENDDO
                        !ENDIF	
!!!!%%%%%%%%%%%%%%%%%%%IN CASE TO MOVE WITH THE CENTRE OF MASS END!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
!!!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SWICTH ON OF THE CROSSLINKS!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           IF(time .GT. 5.d0)THEN
                           DO j = 1,Nc
                                DO i = 1,N
                                        IF(dipth(i,j) .EQ. 1.d0 .AND. jj((j-1)*N+i) .NE. (j-1)*N+i)THEN
                                                Ns = 1.d0
                                                Nsc = 1.d0
                                                jj1 = MOD(jj((j-1)*N+i),N)
                                                j1 = INT(jj((j-1)*N+i)/N)+1
                                                IF(jj1 .EQ. 0.d0)THEN
                                                        jj1 = N
                                                        j1 = j1-1
                                                ENDIF
                              CALL SpringForce(x_temp1(i,j),y_temp1(i,j),z_temp1(i,j),x_temp1(jj1,j1)&
                                                ,y_temp1(jj1,j1),z_temp1(jj1,j1),Fx,Fy,Fz,Ns,Nsc)
                                                !Fqx(i,j) = Fqx(i,j) + dipth(i,j)*Fx
                                                !Fqy(i,j) = Fqy(i,j) + dipth(i,j)*Fy
                                                !Fqz(i,j) = Fqz(i,j) + dipth(i,j)*Fz
                                                !Fqx(jj1,j1) = Fqx(jj1,j1) - dipth(i,j)*Fx
                                                !Fqy(jj1,j1) = Fqy(jj1,j1) - dipth(i,j)*Fy
                                                !Fqz(jj1,j1) = Fqz(jj1,j1) - dipth(i,j)*Fz
                Fqx(i,j) = Fqx(i,j) + patch(i,j)*patch(jj1,j1)*dipth(i,j)*Fx/(dipth(i,j)+dipth(jj1,j1))
                Fqy(i,j) = Fqy(i,j) + patch(i,j)*patch(jj1,j1)*dipth(i,j)*Fy/(dipth(i,j)+dipth(jj1,j1))
                Fqz(i,j) = Fqz(i,j) + patch(i,j)*patch(jj1,j1)*dipth(i,j)*Fz/(dipth(i,j)+dipth(jj1,j1))
                Fqx(jj1,j1) = Fqx(jj1,j1)-patch(i,j)*patch(jj1,j1)*dipth(i,j)*Fx/(dipth(i,j)+dipth(jj1,j1))
                Fqy(jj1,j1) = Fqy(jj1,j1)-patch(i,j)*patch(jj1,j1)*dipth(i,j)*Fy/(dipth(i,j)+dipth(jj1,j1))
                Fqz(jj1,j1) = Fqz(jj1,j1)-patch(i,j)*patch(jj1,j1)*dipth(i,j)*Fz/(dipth(i,j)+dipth(jj1,j1))
                                        ENDIF
                                ENDDO
                        ENDDO
                        ENDIF
                        l = 0
                        DO i = 1,Nc
                                DO j = 1,N
                                        l = l + 1
                                        x_temp1(j,i) = x_temp1(j,i)+del_t*(u_1(3*(l - 1) + 1)+Fljx(j,i)+Fqx(j,i))
                                        y_temp1(j,i) = y_temp1(j,i)+del_t*(u_1(3*(l - 1) + 2)+Fljy(j,i)+Fqy(j,i))
                                        z_temp1(j,i) = z_temp1(j,i)+del_t*(u_1(3*(l - 1) + 3)+Fljz(j,i)+Fqz(j,i))
                                ENDDO
                        ENDDO   
                        !PAUSE              
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%WITHOUT INTERACTION IN BROWNIAN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
			DO j = 1,Nc
			        DO i = 1,N
		         	        !F_brown = gasdev(idum) 
		                        F_browny(i,j) = 1.d0*DSQRT(2.d0*del_t)*gasdev(idum) 
		        	        y_temp1(i,j) = y_temp1(i,j) + F_browny(i,j)   
				        !F_brown = gasdev(idum)  
				        F_brownx(i,j) = 1.d0*DSQRT(2.d0*del_t)*gasdev(idum) 
				        x_temp1(i,j) = x_temp1(i,j)+F_brownx(i,j)
				        !F_brown = gasdev(idum) 
				        F_brownz(i,j) = 1.d0*DSQRT(2.d0*del_t)*gasdev(idum) 
				        z_temp1(i,j) = z_temp1(i,j)+F_brownz(i,j)
			        ENDDO
			ENDDO                   
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%WITHOUT INTERACTION IN BROWNIAN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        tchk = 0.d0
			DO j = 1,Nc
                        DO i = 1,N-1
                                Qx(i,j) = x_temp1(i,j) - x_temp1(i+1,j)
                	        Qy(i,j) = y_temp1(i,j) - y_temp1(i+1,j) 
        		        Qz(i,j) = z_temp1(i,j) - z_temp1(i+1,j)
			        Qo(i,j) = Qx(i,j)**2.d0 + Qy(i,j)**2.d0 + Qz(i,j)**2.d0
			ENDDO
		        Qbar = DSQRT(Qo(1,j))
        		Fbx(1,j) = -Hc*(Qx(1,j))/(1-(Qo(1,j)/Bpar))
        		Fby(1,j) = -Hc*(Qy(1,j))/(1-(Qo(1,j)/Bpar))
        		Fbz(1,j) = -Hc*(Qz(1,j))/(1-(Qo(1,j)/Bpar))
        		x(1,j) = x_temp1(1,j)+Fbx(1,j)*del_t
                	y(1,j) = y_temp1(1,j)+Fby(1,j)*del_t
        		z(1,j) = z_temp1(1,j)+Fbz(1,j)*del_t
        		norm = DSQRT((x(1,j)*scale2/aoc)**2.d0 + (y(1,j)*scale2/boc)**2.d0 + (z(1,j)*scale2)**2.d0)
					IF(norm .GT. (Rs-dr))THEN
				                x(1,j) = (x(1,j)/norm)*(Rs-dr)
				                y(1,j) = (y(1,j)/norm)*(Rs-dr)
				                z(1,j) = (z(1,j)/norm)*(Rs-dr)
				        ENDIF
				        
                        DO i = 2,N-1 
                        Fbx(i,j) = -Hc*(Qx(i,j))/(1-(Qo(i,j)/Bpar))-&
                                                        Hc*(-Qx(i-1,j))/(1-(Qo(i-1,j)/Bpar)) 
                        Fby(i,j) = -Hc*(Qy(i,j))/(1-(Qo(i,j)/Bpar))-&
                                                        Hc*(-Qy(i-1,j))/(1-(Qo(i-1,j)/Bpar)) 
                        Fbz(i,j) = -Hc*(Qz(i,j))/(1-(Qo(i,j)/Bpar))-&
                                                        Hc*(-Qz(i-1,j))/(1-(Qo(i-1,j)/Bpar))
			x(i,j) = x_temp1(i,j)+Fbx(i,j)*del_t
			y(i,j) = y_temp1(i,j)+Fby(i,j)*del_t
			z(i,j) = z_temp1(i,j)+Fbz(i,j)*del_t
        		norm = DSQRT((x(i,j)*scale2/aoc)**2.d0 + (y(i,j)*scale2/boc)**2.d0 + (z(i,j)*scale2)**2.d0)
					IF(norm .GT. (Rs-dr))THEN
					        !PRINT*,i,j,"YES"
				                x(i,j) = (x(i,j)/norm)*(Rs-dr)
				                y(i,j) = (y(i,j)/norm)*(Rs-dr)
			 	                z(i,j) = (z(i,j)/norm)*(Rs-dr)
				        ENDIF 	
			ENDDO
        		Qbar = DSQRT(Qo(N-1,j))
        		Fbx(N,j) = -Hc*(-Qx(N-1,j))/(1-(Qo(N-1,j)/Bpar)) 
                        Fby(N,j) = -Hc*(-Qy(N-1,j))/(1-(Qo(N-1,j)/Bpar)) 
                        Fbz(N,j) = -Hc*(-Qz(N-1,j))/(1-(Qo(N-1,j)/Bpar))
		        x(N,j) = x_temp1(N,j)+Fbx(N,j)*del_t
                	y(N,j) = y_temp1(N,j)+Fby(N,j)*del_t
        		z(N,j) = z_temp1(N,j)+Fbz(N,j)*del_t
        		norm = DSQRT((x(N,j)*scale2/aoc)**2.d0 + (y(N,j)*scale2/boc)**2.d0 + (z(N,j)*scale2)**2.d0)
					IF(norm .GT. (Rs-dr))THEN
				                x(N,j) = (x(N,j)/norm)*(Rs-dr)
				                y(N,j) = (y(N,j)/norm)*(Rs-dr)
				                z(N,j) = (z(N,j)/norm)*(Rs-dr)
				        ENDIF 
			ENDDO        				        	
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF SPRING FORCING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                DO i = 1,Ntr
                		        norm = DSQRT((x_tr(i,1)/aoc)**2.d0 + (x_tr(i,2)/boc)**2.d0 + x_tr(i,3)**2.d0)               
                		        IF(norm .GT. (Rs-dr))THEN
				                        x_tr(i,1) = (x_tr(i,1)/norm)*(Rs-dr)
				                        x_tr(i,2) = (x_tr(i,2)/norm)*(Rs-dr)
				                        x_tr(i,3) = (x_tr(i,3)/norm)*(Rs-dr)
				        ENDIF
                		ENDDO  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PRINT THE POSITIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                
				IF(MOD(m,2000) .EQ. 0.d0)THEN
				        DO j = 1,Nc
				                l = j
                                                DO i = 1,N
                			        WRITE(131+j),REAL(time),REAL(Fbx(i,j)),REAL(Fby(i,j)),REAL(Fbz(i,j))  
                			        WRITE(161+j),REAL(time),REAL(Fdpx(i,j)+Fdpxa(i,j)),&
                			                    REAL(Fdpy(i,j)+Fdpya(i,j)),REAL(Fdpz(i,j)+Fdpza(i,j)),&
                			                    INT(dip(i,j)),REAL(trel(i,j)),REAL(tswitch(i,j))        
                			        WRITE(191+j),REAL(time),REAL(Fljx(i,j)),REAL(Fljy(i,j)),REAL(Fljz(i,j))
                			        WRITE(221+j),REAL(time),REAL(Fqx(i,j)),REAL(Fqy(i,j)),REAL(Fqz(i,j)),&
                                                            INT(dipth(i,j)),REAL(trel1(i,j)),REAL(tsth(i,j))      
                                                ENDDO 
                		        ENDDO	
                		        DO i = 1,Nc*N
                                                WRITE(551),REAL(time),REAL(u_1(3*(i-1)+1)),REAL(u_1(3*(i-1)+2))&
                                                        ,REAL(u_1(3*(i-1)+3))
                		        ENDDO                 
	        		ENDIF 				               
                                !x = x_new
                                !y = y_new
                                !z = z_new
!%%%%%%%%%%%%%%%%%%%%%%MOBILITY TENSOR (GREENS FUNTION FOR STOKESLET INSIDE A SPHERE- KIM and MAUL 1994 POF)%%%%%%%%%%%%%%%
                                !IF(MOD(m,1000) .EQ. 0.d0)THEN
                                !        CALL SphereGeo(x,y,z,Dxx, Dyy, Dzz, Dxy, Dyx, Dyz, Dzy, Dxz, Dzx, Rs,Rad_Bead,N)
                                !ENDIF
				!DO i = 1,N
                                !	x(i) = x(i)*Rs
                                !	y(i) = y(i)*Rs
                                !	z(i) = z(i)*Rs
					!PRINT*, i, x(i), y(i), z(i)
                        	!ENDDO
				!Rad_Bead = Rad_Bead*Rs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PRINT VALUES IN A FILE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%				
				IF(MOD(m,1000) .EQ. 0.d0)THEN
        				DO j = 1,Nc
        				leng(j) = 0.d0
		                		DO i = 2,N
		                		        leng(j) = leng(j) + DSQRT((x(i,j) - x(i-1,j))**2 + (y(i,j) &
		                		                - y(i-1,j))**2 + (z(i,j) - z(i-1,j))**2) 
		                		ENDDO 
		                	ENDDO
		               ENDIF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PRINT ITERATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        			IF(MOD(m,2000) .EQ. 0.d0)THEN
				        DO l = 1,Nc
				        DO i=1,N
				                WRITE(41+l),REAL(time),REAL(x(i,l)),INT(jj((l-1)*N+i))
				        ENDDO
				        DO i=1,N
				                WRITE(71+l),REAL(time),REAL(y(i,l)),INT(jj((l-1)*N+i))
				        ENDDO
				        DO i=1,N
				                WRITE(101+l),REAL(time),REAL(z(i,l)),INT(jj((l-1)*N+i))
				        ENDDO
				        ENDDO
				ENDIF	
				
				IF(MOD(m,1000) .EQ. 0.d0)THEN
				        PRINT*,m,m*del_t, SUM(leng)
				ENDIF
				!PAUSE
		ENDDO	
		call delete_fmm_wrapper(fmm); 
                call MPI_BARRIER(MPI_COMM_WORLD, ierror); 
                call MPI_FINALIZE(ierror)				
		DEALLOCATE(x)
		DEALLOCATE(y)
		DEALLOCATE(z)
		DEALLOCATE(x_temp1)
		DEALLOCATE(y_temp1)
		DEALLOCATE(z_temp1)
		DEALLOCATE(x_temp)
		DEALLOCATE(y_temp)
		DEALLOCATE(z_temp)
		DEALLOCATE(Fljx)
		DEALLOCATE(Fljy)
		DEALLOCATE(Fljz)
		DEALLOCATE(Fljx1)
		DEALLOCATE(Fljy1)
		DEALLOCATE(Fljz1)
		DEALLOCATE(Fdpx)
		DEALLOCATE(Fdpy)
		DEALLOCATE(Fdpz)
		DEALLOCATE(Fdpxa)
		DEALLOCATE(Fdpya)
		DEALLOCATE(Fdpza)
		DEALLOCATE(Fdpx1)
		DEALLOCATE(Fdpy1)
		DEALLOCATE(Fdpz1)
		DEALLOCATE(Fdpxa1)
		DEALLOCATE(Fdpya1)
		DEALLOCATE(Fdpza1)
		DEALLOCATE(F_brownx)
		DEALLOCATE(F_browny)
		DEALLOCATE(F_brownz)
                !DEALLOCATE(F_brown)
		DEALLOCATE(trel)
		DEALLOCATE(tswitch)
		DEALLOCATE(dip)
		DEALLOCATE(x_tr)
		DEALLOCATE(u_tr)
		DEALLOCATE(Qx)
		DEALLOCATE(Qy)
		DEALLOCATE(Qz)
		DEALLOCATE(Qo)
		DEALLOCATE(Fbx1)
		DEALLOCATE(Fby1)
		DEALLOCATE(Fbz1)
		DEALLOCATE(Fbx)
		DEALLOCATE(Fby)
		DEALLOCATE(Fbz)
		DEALLOCATE(ux)
		DEALLOCATE(uy)
		DEALLOCATE(uz)
		DEALLOCATE(Fdp)
		DEALLOCATE(Fdp1)
		DEALLOCATE(u_1)
		DEALLOCATE(leng)
		DEALLOCATE(RM1)
                DEALLOCATE(alpha)
                DEALLOCATE(beta)
                DEALLOCATE(gamma1)
END PROGRAM Bead_Rod
!************************END OF THE MAIN PROGRAM************************************
!****************************END OF THE CODE****************************************
