PROGRAM EXTRACT
	IMPLICIT NONE
	        TYPE Pos
			REAL*4 :: time
			REAL*4 :: xpos
		END TYPE Pos
		INTEGER*4 :: i, n, k, j, l, filenum, temp, s, nobeadsart, num, ip, ipp, bdl(2000000),&
					count1, T_frame, nobeads, Nc, i1, c, m, ic
		CHARACTER(LEN =30) :: filename
		REAL*4 :: del_t, T_end
		TYPE(Pos),ALLOCATABLE,DIMENSION(:) :: X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, &
		                                      Z4
		i = 0
		Nc = 4
		nobeads = 750
		T_frame = 1000
		T_end = 19.5d0
		del_t = 0.000005d0
		n = ((T_end/(del_t))/T_frame)*nobeads
		nobeadsart = nobeads*Nc
		count1 = 0
			!OPEN(23,FILE = "xpos.dat", FORM="unformatted")
			!OPEN(24,FILE = "ypos.dat", FORM="unformatted")
			!OPEN(25,FILE = "zpos.dat", FORM="unformatted")
			!	DO WHILE(i .LE. 0)	
			!		READ(23,IOSTAT = i)i
			!			n = n + 1
			!			PRINT*,i,n
			!	ENDDO
				
			PRINT*, n
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
				l = 0
				num = 1
				k = 0
				c = 0
				j = 0
				m = 0
				DO i = 1,n*Nc
						IF(mod(i,nobeadsart*l+1) .EQ. 0 .OR. i .EQ. 1)THEN
								l = l + 1
								filenum = l
								WRITE(filename, '("position",I4.4, ".pdb")' )filenum
								OPEN(FILE = filename,UNIT = 16)
                                                		WRITE(16,11)'HEADER    Rods positions'
        							WRITE(16,10)'COMPND    '
         							WRITE(16,10)'SOURCE    '
         							WRITE(16,10)'REMARK    '
                                                		num  = 1
                                                                k = 0
				                                c = 0
				                                j = 0
				                                m = 0
						ELSE
							    filenum = filenum
							    num = num + 1    

						ENDIF
						!PRINT*,(i-1)/nobeads,i
						!PRINT*,l,i,n*Nc,(i-1)/nobeads
						IF((i-1)/nobeads .EQ. 4*(l-1))THEN
                                                k = k + 1
						i1 = (l-1)*nobeads + k 
						!IF(l .EQ. 1)THEN
						!PRINT*,"YES",num,k,i1 
						!PAUSE
						!ENDIF
						!IF(l .EQ. 1)THEN
						!IF(DSQRT(X1(i1)%xpos**2.d0+Y1(i1)%xpos**2.d0+Z1(i1)%xpos**2.d0) .GT. (8.d0 - 8.d0*0.04999d0))THEN
						!PRINT*,"YES",i1,m,l,X2(i1)%xpos,Y2(i1)%xpos,Z2(i1)%xpos, DSQRT(X2(i1)%xpos**2.d0+&
						!        Y2(i1)%xpos**2.d0+Z2(i1)%xpos**2.d0)
						!PAUSE
						!ENDIF
						WRITE(16,12)'ATOM  ',1,'  CA  LYS ',k,'    ',X1(i1)%xpos*1.d0,Y1(i1)%xpos*1.d0,Z1(i1)%xpos*1.d0,1.0,1.0
						!PRINT*,DSQRT((X(i+1) - X(i))**2.d0+(Y(i+1) - Y(i))**2.d0 +(Z(i+1) - Z(i))**2.d0) 
						ELSEIF((i-1)/nobeads .EQ. 1 + 4*(l-1))THEN
						m = m + 1
						IF(m .EQ. 1)write(16,13)'TER'
						i1 = (l-1)*nobeads + m
						!PRINT*,"YES1",num,m,i1,l
						!PAUSE 
						!ENDIF
						!IF(DSQRT(X2(i1)%xpos**2.d0+Y2(i1)%xpos**2.d0+Z2(i1)%xpos**2.d0) .GT. (8.d0 - 8.d0*0.04999d0))THEN
						!PRINT*,"YES",i1,m,l,X2(i1)%xpos,Y2(i1)%xpos,Z2(i1)%xpos, DSQRT(X2(i1)%xpos**2.d0+&
						!        Y2(i1)%xpos**2.d0+Z2(i1)%xpos**2.d0)
						!PAUSE
						!ENDIF
						WRITE(16,12)'ATOM  ',1,'  CA  LEU ',m,'    ',X2(i1)%xpos*1.d0,Y2(i1)%xpos*1.d0,Z2(i1)%xpos*1.d0,1.0,1.0
						ELSEIF((i-1)/nobeads .EQ. 2 + 4*(l-1))THEN
						j = j + 1
						IF(j .EQ. 1)write(16,13)'TER'
						i1 = (l-1)*nobeads + j
						!IF(l .EQ. 2)THEN
						!PRINT*,"YES2",num,j,i1,l 
						!IF(DSQRT(X3(i1)%xpos**2.d0+Y3(i1)%xpos**2.d0+Z3(i1)%xpos**2.d0) .GT. (8.d0 - 8.d0*0.04999d0))THEN
						!PRINT*,"YES",i1,m,l,X2(i1)%xpos,Y2(i1)%xpos,Z2(i1)%xpos, DSQRT(X2(i1)%xpos**2.d0+&
						!        Y2(i1)%xpos**2.d0+Z2(i1)%xpos**2.d0)
						!PAUSE
						!ENDIF
						!PAUSE
						!ENDIF
						WRITE(16,12)'ATOM  ',1,'  CA  MET ',j,'    ',X3(i1)%xpos*1.d0,Y3(i1)%xpos*1.d0,Z3(i1)%xpos*1.d0,1.0,1.0
						ELSEIF((i-1)/nobeads .EQ. 3 + 4*(l-1))THEN
						c = c + 1
						IF(c .EQ. 1)write(16,13)'TER'
						i1 = (l-1)*nobeads + c
						!IF(DSQRT(X4(i1)%xpos**2.d0+Y4(i1)%xpos**2.d0+Z4(i1)%xpos**2.d0) .GT. (8.d0 - 8.d0*0.04999d0))THEN
						!PRINT*,"YES",i1,m,l,X2(i1)%xpos,Y2(i1)%xpos,Z2(i1)%xpos, DSQRT(X2(i1)%xpos**2.d0+&
						!        Y2(i1)%xpos**2.d0+Z2(i1)%xpos**2.d0)
						!PAUSE
						!ENDIF
						!IF(l .EQ. 2)THEN
						!PRINT*,"YES3",num,c,i1,l 
						!PAUSE
						!ENDIF
						WRITE(16,12)'ATOM  ',1,'  CA  HIS ',c,'    ',X4(i1)%xpos*1.d0,Y4(i1)%xpos*1.d0,Z4(i1)%xpos*1.d0,1.0,1.0
						ENDIF
				ENDDO
10      FORMAT(1a10)
11      FORMAT(1a24)
12      FORMAT(a6,i5,a10,i5,a4,3f8.3,2f6.2)
13      FORMAT(a3)
14      FORMAT(a64)	


      count1 = 0
      do ic = 1,Nc
         do ip = 1,nobeads-1
            count1 = count1+1
            bdl(count1) = ip+nobeads*(ic-1)
            count1 = count1+1
            bdl(count1) = ip+nobeads*(ic-1)+1
         enddo
      enddo
      open(unit=6,file='links.psf',status='unknown')
      write(6,20)'PSF'
      write(6,21)'       4 !NTITLE'
      write(6,21)' REMARKS BONDS  ' 
      write(6,*)
      write(6,21)'     3000 !NATOM '
 20   format(1a3)
 21   format(1a16)
 22   format(i8,a36,a26)
 23   format(1a27)
 24   format(8i8)
      
      do ip = 1,nobeads     
          
         write(6,22) ip,' MAIN 1    LYS  N    NH3    0.000000',&
             '       0.00000           0'            

      enddo
      do ip = 1,nobeads     
          
         write(6,22) ip+nobeads,' MAIN 1    LEU  N    NH3    0.000000',&
             '       0.00000           0'            

      enddo
      do ip = 1,nobeads     
          
         write(6,22) ip+2*nobeads,' MAIN 1    MET  N    NH3    0.000000',&
             '       0.00000           0'            

      enddo
      do ip = 1,nobeads     
          
         write(6,22) ip+3*nobeads,' MAIN 1    HIS  N    NH3    0.000000',&
             '       0.00000           0'            

      enddo



      write(6,*)
      write(6,23)'    2996 !NBOND: bonds     '
      
      do ip = 1,count1/8+1
         ipp=8*(ip-1)+1
         write(6,24) bdl(ipp),bdl(ipp+1),bdl(ipp+2),bdl(ipp+3),&
             bdl(ipp+4),bdl(ipp+5),bdl(ipp+6),bdl(ipp+7)

      enddo

      write(6,*)
      write(6,23)'       0 !NTHETA: angles   '

      write(6,*)
      write(6,23)'       0 !nobeadsHI: dihedrals  '

      write(6,*)
      write(6,23)'       0 !NIMPHI: impropers'

      write(6,*)
      write(6,23)'       0 !NDON: donors     '

      write(6,*)
      write(6,23)'       0 !NACC: acceptors  '

      write(6,*)
      write(6,23)'       0 !NNB              '

      write(6,*)
      write(6,23)'       0       0 !NGRP     '

      close(6)				
END PROGRAM EXTRACT
