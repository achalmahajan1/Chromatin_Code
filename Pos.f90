PROGRAM EXTRACT
	IMPLICIT NONE
	        TYPE Pos
			REAL*4 :: time
			REAL*4 :: xpos
		END TYPE Pos
		INTEGER*4 :: i,j,n,nobeads,Nframes,T_frame
		REAL*4 :: del_t, T_end
		TYPE(Pos), ALLOCATABLE, DIMENSION(:) :: X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, &
		                                      Z4,X5,Y5,Z5,X6,Y6,Z6,X7,Y7,Z7,X8,Y8, &
		                                      Z8
		nobeads = 750
		T_frame = 1000
		T_end =  10.d0
		del_t = 0.000005d0
		n = ((T_end/(del_t))/T_frame)*nobeads + 1
			!OPEN(23,FILE = "xpos.dat", FORM="unformatted")
			!OPEN(24,FILE = "ypos.dat", FORM="unformatted")
			!OPEN(25,FILE = "zpos.dat", FORM="unformatted")
			!	DO WHILE(i .GE. 0)	
			!		READ(23,IOSTAT = i)i
			!			n = n + 1
			!	ENDDO
				
			!PRINT*, "n", n
			ALLOCATE(X1(n))
			ALLOCATE(Y1(n))
			ALLOCATE(Z1(n))
			ALLOCATE(X2(n))
			ALLOCATE(Y2(n))
			ALLOCATE(Z2(n))
			ALLOCATE(X3(n))
			ALLOCATE(Y3(n))
			ALLOCATE(Z3(n))
			ALLOCATE(X4(n))
			ALLOCATE(Y4(n))
			ALLOCATE(Z4(n))
			!REWIND(23)
			!REWIND(24)
			!REWIND(25)
			OPEN(23,FILE = "x_chain01.dat", FORM="unformatted")
			OPEN(24,FILE = "y_chain01.dat", FORM="unformatted")
			OPEN(25,FILE = "z_chain01.dat", FORM="unformatted")
			OPEN(26,FILE = "x_chain02.dat", FORM="unformatted")
			OPEN(27,FILE = "y_chain02.dat", FORM="unformatted")
			OPEN(28,FILE = "z_chain02.dat", FORM="unformatted")
			OPEN(29,FILE = "x_chain03.dat", FORM="unformatted")
			OPEN(30,FILE = "y_chain03.dat", FORM="unformatted")
			OPEN(31,FILE = "z_chain03.dat", FORM="unformatted")
			OPEN(32,FILE = "x_chain04.dat", FORM="unformatted")
			OPEN(33,FILE = "y_chain04.dat", FORM="unformatted")
			OPEN(34,FILE = "z_chain04.dat", FORM="unformatted")
				DO j = 1,n
				        READ(23),X1(j)
				        READ(24),Y1(j)
				        READ(25),Z1(j)
				        READ(26),X2(j)
				        READ(27),Y2(j)
				        READ(28),Z2(j)
				        READ(29),X3(j)
				        READ(30),Y3(j)
				        READ(31),Z3(j)
				        READ(32),X4(j)
				        READ(33),Y4(j)
				        READ(34),Z4(j)
				        !IF(j .EQ. 1250)PRINT*,X2(j)%xpos, Y2(j)%xpos, Z2(j)%xpos, X4(j)%xpos
				        !PAUSE
				ENDDO
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                        
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MEAN SQUARE DISPLACEMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                       
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                  
                                        OPEN(41,FILE = "Pos1.dat")
                                        DO i = 1,N
                                                WRITE(41,*),X1(i)%time,X1(i)%xpos,Y1(i)%xpos,Z1(i)%xpos
                                        ENDDO
                                        CLOSE(41)
                                        OPEN(42,FILE = "Pos2.dat")
                                        DO i = 1,N
                                                WRITE(42,*),X2(i)%time,X2(i)%xpos,Y2(i)%xpos,Z2(i)%xpos
                                        ENDDO
                                        CLOSE(42)
                                        OPEN(43,FILE = "Pos3.dat")
                                        DO i = 1,N
                                                WRITE(43,*),X3(i)%time,X3(i)%xpos,Y3(i)%xpos,Z3(i)%xpos
                                        ENDDO
                                        CLOSE(43)
                                        OPEN(44,FILE = "Pos4.dat")
                                        DO i = 1,N
                                                WRITE(44,*),X4(i)%time,X4(i)%xpos,Y4(i)%xpos,Z4(i)%xpos
                                        ENDDO
                                        CLOSE(44)                   
END PROGRAM EXTRACT
