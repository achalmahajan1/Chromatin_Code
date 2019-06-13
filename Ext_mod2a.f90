PROGRAM Bead_Rod
        USE Variables
	IMPLICIT NONE
		INTEGER*8 :: i,j,Niter,k,m,l,npts,Nelm,i1,i2,i3,i4,i5,i6,ind,filenum,&
		                idum  
		REAL*8 :: del_t,T_end,Rad_bead,tol,v1,v2,Rodl,dr,time,gasdev,xcm1,ycm1,zcm1,&
		          start,finish,boa,coa,oot,scale1,x_axis,y_axis,z_axis,norm
                REAL*8, ALLOCATABLE, DIMENSION(:) :: Fdp1,F_brown,noise,Bcho,x_temp,z_temp,y_temp,RM2,leng,&
		                                     ux,uy,uz
		REAL*8, ALLOCATABLE, DIMENSION(:,:) :: x,y,z,F_brownx,F_browny,F_brownz,&
		                                       x_tr,u_tr,Fdp,RM1
                CHARACTER(LEN =30) :: filename                                                                                                                   
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%VARIABLE ASSIGNMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ALLOCATE(x(N,Nc))
		ALLOCATE(y(N,Nc))
		ALLOCATE(z(N,Nc))
		ALLOCATE(x_temp1(N,Nc))
		ALLOCATE(y_temp1(N,Nc))
		ALLOCATE(z_temp1(N,Nc))
		ALLOCATE(x_temp(N*Nc+Ntr))
		ALLOCATE(y_temp(N*Nc+Ntr))
		ALLOCATE(z_temp(N*Nc+Ntr))
		ALLOCATE(Fljx(N,Nc))
		ALLOCATE(Fljy(N,Nc))
		ALLOCATE(Fljz(N,Nc))
		ALLOCATE(Fdpx(N,Nc))
		ALLOCATE(Fdpy(N,Nc))
		ALLOCATE(Fdpz(N,Nc))
		ALLOCATE(Fdpxa(N,Nc))
		ALLOCATE(Fdpya(N,Nc))
		ALLOCATE(Fdpza(N,Nc))
		ALLOCATE(F_brownx(N,Nc))
		ALLOCATE(F_browny(N,Nc))
		ALLOCATE(F_brownz(N,Nc))
		ALLOCATE(trel(N,Nc))
		ALLOCATE(tswitch(N,Nc))
		ALLOCATE(dip(N,Nc))
		ALLOCATE(x_tr(Ntr,3))
		ALLOCATE(u_tr(Ntr,3))
		ALLOCATE(Qx(N-1,Nc))
		ALLOCATE(Qy(N-1,Nc))
		ALLOCATE(Qz(N-1,Nc))
		ALLOCATE(Qo(N-1,Nc))
		ALLOCATE(Fbx(N,Nc))
		ALLOCATE(Fby(N,Nc))
		ALLOCATE(Fbz(N,Nc))
		ALLOCATE(ux(N*Nc))
		ALLOCATE(uy(N*Nc))
		ALLOCATE(uz(N*Nc))
		ALLOCATE(Fdp(N*Nc,3))
		ALLOCATE(Fdp1(N*Nc*3))
		ALLOCATE(u_1(N*Nc,3))
		ALLOCATE(leng(Nc))
		u_1 = 0.d0
		T_end = 40.d0
		dr = Rs*0.05d0
		time = 0.d0
		Fljx = 0.d0
		Fljy = 0.d0
		Fljz = 0.d0
		Fdpx = 0.d0
		Fdpy = 0.d0
		Fdpz = 0.d0
		Fbx = 0.d0
		Fby = 0.d0
		Fbz = 0.d0
		Fdpxa = 0.d0
		Fdpya = 0.d0
		Fdpza = 0.d0
		Rodl = 1.d0 !(Length of the rod scaled with Sphere radius)
		Rad_Bead = 0.1d0 !(Radius of the bead scaled with Sphere radius)
		ind = 1.d0
                ux = 0.d0
                uy = 0.d0
                uz = 0.d0
		idum = -140
		Nelm = 20
                !F_brown = 0.d0
		del_t = 0.000005d0
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
		OPEN(UNIT=181, FILE='BEM.dat')
                DO l = 1,Nc
        		filenum = l
	        	WRITE(filename, '("x_chain",I2.2, ".dat")')filenum
	        	OPEN(FILE = filename,UNIT = 41 + l, FORM="unformatted")
	        	WRITE(filename, '("y_chain",I2.2, ".dat")')filenum
	        	OPEN(FILE = filename,UNIT = 51 + l, FORM="unformatted")
	        	WRITE(filename, '("z_chain",I2.2, ".dat")')filenum
	        	OPEN(FILE = filename,UNIT = 61 + l, FORM="unformatted")
	        	WRITE(filename, '("ForceTension_chain",I2.2, ".dat")')filenum
	        	OPEN(FILE = filename,UNIT = 71 + l, FORM="unformatted")
	        	WRITE(filename, '("ForceActive_chain",I2.2, ".dat")')filenum
	        	OPEN(FILE = filename,UNIT = 81 + l, FORM="unformatted")
	        	WRITE(filename, '("ForceLJ_chain",I2.2, ".dat")')filenum
	        	OPEN(FILE = filename,UNIT = 91 + l, FORM="unformatted")
		ENDDO
                OPEN(UNIT=101, FILE='Velocity.dat', FORM="unformatted")
                OPEN(UNIT=24, FILE='Params.dat')
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
                
                !scale1 = Rs/(boa*coa)**oot
                scale1 = Rs/(boc*aoc)**oot
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
                OPEN(32,FILE = "../Extensile/MatrixInverseG1280.txt", FORM="unformatted")
                CALL CPU_TIME(start)
                DO i = 1,Nelm*Nelm*9
                        READ(32),RM2(i)! No need to read at every time step, just once at the beginning and can be saved in the RAM
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
                        WRITE(181,*),i,l,p(i1,1),p(i1,2),p(i1,3)
                        l = l + 1 
                        WRITE(181,*),i,l,p(i2,1),p(i2,2),p(i2,3)
                        l = l + 1
                        WRITE(181,*),i,l,p(i3,1),p(i3,2),p(i3,3) 
                        !PRINT*,i, p(i1,1),p(i1,2),p(i1,3)       
               ENDDO 
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
                                p(i6,2),p(i6,3))
                                !PRINT*,"outside", DSQRT(p(i2,1)**2.d0 + p(i2,2)**2.d0 + p(i2,3)**2.d0)
                                !PAUSE
                ENDDO
                CALL Geo_val(nelm,npts)
                vna = -1.d0*vna
!******************************************INITIALIZE THE POSITION AND DIPOLES******************************************
			CALL Initconfig(x,y,z,Rodl)
                        PRINT*,"Here"
                        CALL Initdip
                        PRINT*,"Here1"
                        DO l = 1,Nc
			        DO i=1,N
			                WRITE(41 + l),REAL(time),REAL(x(i,l))
			        ENDDO
			        DO i=1,N
			                WRITE(51 + l),REAL(time),REAL(y(i,l))
			        ENDDO
			        DO i=1,N
			                WRITE(61 + l),REAL(time),REAL(z(i,l))
			        ENDDO
			ENDDO
			DO j = 1,Nc
        		        leng(j) = 0.d0
		            	        DO i = 2,N
		                                leng(j) = leng(j) + DSQRT((x(i,j) - x(i-1,j))**2 + (y(i,j) &
		                                                - y(i-1,j))**2 + (z(i,j) - z(i-1,j))**2) 
		                        ENDDO 
		        ENDDO
		        PRINT*,"Length of",Nc,"chains", leng 
			j = 1	
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
                                        trel(k,i) = trel(k,i) + del_t
                                ENDDO
                        ENDDO
                        
		        xcm1 = 0.d0
                        ycm1 = 0.d0
                        zcm1 = 0.d0  
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
                        CALL ActiveForces(x,y,z) 		
                        CALL Leonard_jones(x,y,z,Rad_Bead,Rodl)
                        PRINT*,Fljx	
                        PAUSE
                        
			IF(MOD(m,1000) .EQ. 0.d0 .OR. m .EQ. 5)THEN !DON'T USE m .EQ.1 because if the beads are closer in initial config can lead to diverging
                        DO i = 1,Nc
                                DO j = 1,N
                                        l = l + 1
                                        Fdp(l,1) = Fdpx(j,i) + Fdpxa(j,i) + Fljx(j,i) + Fbx(j,i)
                                        Fdp(l,2) = Fdpy(j,i) + Fdpya(j,i) + Fljy(j,i) + Fby(j,i)
                                        Fdp(l,3) = Fdpz(j,i) + Fdpza(j,i) + Fljz(j,i) + Fbz(j,i)
                                ENDDO
                        ENDDO
                        u_1 = 0.d0
		        !CALL BoundaryIntFinal1(Fdp1,x_temp,y_temp,z_temp,u_1,p,nn,Nelm,alpha,beta,gamma1,RM1,fmm,N*Nc)
		        CALL BoundaryInt(Fdp,x_temp,y_temp,z_temp,Nelm,RM1)                         
			ENDIF 
			
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
                        l = 0
                        DO i = 1,Nc
                                DO j = 1,N
                                        l = l + 1
                                        x_temp1(j,i) = x_temp1(j,i)+del_t*(6.d0*Rad_bead*pi*u_1(l,1)+Fljx(j,i))
                                        y_temp1(j,i) = y_temp1(j,i)+del_t*(6.d0*Rad_bead*pi*u_1(l,2)+Fljy(j,i))
                                        z_temp1(j,i) = z_temp1(j,i)+del_t*(6.d0*Rad_bead*pi*u_1(l,3)+Fljz(j,i))
                                ENDDO
                        ENDDO   

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
                        CALL Spring
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF SPRING FORCING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%UPDATE THE POSITIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
                        DO j = 1,Nc	
                                DO i = 1,N		                         
          		                x(i,j) = x_temp1(i,j)+Fbx(i,j)*del_t
			                y(i,j) = y_temp1(i,j)+Fby(i,j)*del_t
			                z(i,j) = z_temp1(i,j)+Fbz(i,j)*del_t
			                norm = DSQRT((x(i,j)/aoc)**2.d0 + (y(i,j)/boc)**2.d0 + z(i,j)**2.d0)
					IF(norm .GT. (Rs-dr))THEN
				                x(i,j) = (x(i,j)/norm)*(Rs-dr)
				                y(i,j) = (y(i,j)/norm)*(Rs-dr)
				                z(i,j) = (z(i,j)/norm)*(Rs-dr)
				        ENDIF 
                                ENDDO	
                        ENDDO                                			        
                                DO i = 1,Ntr
                		        norm = DSQRT((x_tr(i,1)/aoc)**2.d0 + (x_tr(i,2)/boc)**2.d0 + x_tr(i,3)**2.d0)               
                		        IF(norm .GT. (Rs-dr))THEN
				                        x_tr(i,1) = (x_tr(i,1)/norm)*(Rs-dr)
				                        x_tr(i,2) = (x_tr(i,2)/norm)*(Rs-dr)
				                        x_tr(i,3) = (x_tr(i,3)/norm)*(Rs-dr)
				        ENDIF
                		ENDDO	        
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PRINT THE POSITIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                
				IF(MOD(m,5000) .EQ. 0.d0)THEN
				        DO j = 1,Nc
				                l = j
                                                DO i = 1,N
                			        WRITE(71+j),REAL(time),REAL(Fbx(i,j)),REAL(Fby(i,j)),REAL(Fbz(i,j))  
                			        WRITE(81+j),REAL(time),REAL(Fdpx(i,j)+Fdpxa(i,j)),&
                			                    REAL(Fdpy(i,j)+Fdpya(i,j)),REAL(Fdpz(i,j)+Fdpza(i,j)),&
                			                    REAL(dip(i,j)),REAL(trel(i,j)),REAL(tswitch(i,j))        
                			        WRITE(91+j),REAL(time),REAL(Fljx(i,j)),REAL(Fljy(i,j)),REAL(Fljz(i,j))
                			        !WRITE(101+j),time, F_brownx(1,j), F_browny(1,j), F_brownz(1,j)      
                                                ENDDO 
                		        ENDDO	
                		        DO i = 1,Nc*N
                                                WRITE(101),REAL(time),REAL(u_1(i,1)),REAL(u_1(i,2))&
                                                        ,REAL(u_1(i,3))
                		        ENDDO                 
	        		ENDIF 				               
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
        			IF(MOD(m,1000) .EQ. 0.d0)THEN
				        DO l = 1,Nc
				        DO i=1,N
				                WRITE(41+l),REAL(time),REAL(x(i,l))
				        ENDDO
				        DO i=1,N
				                WRITE(51+l),REAL(time),REAL(y(i,l))
				        ENDDO
				        DO i=1,N
				                WRITE(61+l),REAL(time),REAL(z(i,l))
				        ENDDO
				        ENDDO
				ENDIF	
				
				IF(MOD(m,1000) .EQ. 0.d0)THEN
				        PRINT*,m,m*del_t, leng
				ENDIF
				!PAUSE
		ENDDO			
		DEALLOCATE(x_temp)
		DEALLOCATE(y_temp)
		DEALLOCATE(z_temp)
		DEALLOCATE(Fljx)
		DEALLOCATE(Fljy)
		DEALLOCATE(Fljz)
		DEALLOCATE(Fdpx)
		DEALLOCATE(Fdpy)
		DEALLOCATE(Fdpz)
		DEALLOCATE(F_brownx)
		DEALLOCATE(F_browny)
		DEALLOCATE(F_brownz)
		DEALLOCATE(trel)
		DEALLOCATE(tswitch)
		DEALLOCATE(dip)
END PROGRAM Bead_Rod
!************************END OF THE MAIN PROGRAM************************************
!****************************END OF THE CODE****************************************
