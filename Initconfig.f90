SUBROUTINE Initconfig(x1,y1,z1,boa,coa,aoc,boc,Rs1,Rodl,N,Nc)    
!!This Subroutine initializes the position of beads with length of rod scaled to ellipsoidal (cavity) radius
!!Rod length (Cavity) radius is chosen to be the length scale in the problem
        IMPLICIT NONE
        INTEGER*8 :: i,N,Nc,j
        REAL*8 :: px,py,pz,pnorm,Rs,Rs1,dr,xnorm,x1(N,Nc),y1(N,Nc),z1(N,Nc),Rodl,boa,coa,aoc,boc,&
                        R_cx(Nc),R_cy(Nc),R_cz(Nc),num,dr1
        REAL*4 :: rand
	dr = Rs1*0.05d0
	num = 3.d0
	R_cx(1) = -Rs1/num
	R_cy(1) = 0.d0
	R_cz(1) = Rs1/num
	R_cx(2) = Rs1/num
	R_cy(2) = 0.d0
	R_cz(2) = -Rs1/num
	R_cx(3) = Rs1/num
	R_cy(3) = 0.d0
	R_cz(3) = Rs1/num
	R_cx(4) = -Rs1/num
	R_cy(4) = 0.d0
	R_cz(4) = -Rs1/num
	Rs = Rs1!/2.d0
	j = 1
	xnorm = 1000.d0	
	x1(1:N,j) = 1.d0
	y1(1:N,j) = 1000.d0
	z1(1:N,j) = -1.d0
	dr1 = 0.8d0
	DO WHILE(xnorm .GT. (Rs - dr) .OR. x1(1,j) .GT. -dr1 .OR. z1(1,j) .LT. dr1)
        	x1(1,j) = Rs*(rand(0) - 0.5d0)
                y1(1,j) = Rs*(rand(0) - 0.5d0)
                z1(1,j) = Rs*(rand(0) - 0.5d0) 
		xnorm = DSQRT((x1(1,j)/aoc-R_cx(j))**2.d0+(y1(1,j)/boc-R_cy(j))**2.d0+(z1(1,j)-R_cz(j))**2.d0)
	ENDDO
	
	DO i = 1,N-1
		xnorm = 1000.d0
		DO WHILE(xnorm .GT. (Rs - dr) .OR. x1(i+1,j) .GT. -dr1 .OR. z1(i+1,j) .LT. dr1)
			px = (rand(0) - 0.5d0)*2.d0
                        py = (rand(0) - 0.5d0)*2.d0
	        	pz = (rand(0) - 0.5d0)*2.d0
	        	pnorm = DSQRT(px**2.d0 + py**2.d0 + pz**2.d0)
	        	px = Rodl*px/pnorm
	        	py = Rodl*py/pnorm
	        	pz = Rodl*pz/pnorm
	        	x1(i+1,j) = x1(i,j) + px
			y1(i+1,j) = y1(i,j) + py
			z1(i+1,j) = z1(i,j) + pz
			xnorm = DSQRT((x1(i+1,j)/aoc)**2.d0+(y1(i+1,j)/boc)**2.d0+(z1(i+1,j))**2.d0)
			!PRINT*,i,j,x1(i+1,j),y1(i+1,j),z1(i+1,j),xnorm
		ENDDO
		!PRINT*,i+1,j,x1(i+1,j),y1(i+1,j),z1(i+1,j)
		!PAUSE
	ENDDO
	
	
	j = 2
	xnorm = 1000.d0	
	x1(1:N,j) = -1.d0
	y1(1:N,j) = 1000.d0
	z1(1:N,j) = 1.d0
	DO WHILE(xnorm .GT. (Rs - dr) .OR. x1(1,j) .LT. dr1 .OR. z1(1,j) .GT. -dr1)
        	x1(1,j) = Rs*(rand(0) - 0.5d0)
                y1(1,j) = Rs*(rand(0) - 0.5d0)
                z1(1,j) = Rs*(rand(0) - 0.5d0) 
		xnorm = DSQRT((x1(1,j)/aoc-R_cx(j))**2.d0+(y1(1,j)/boc-R_cy(j))**2.d0+(z1(1,j)-R_cz(j))**2.d0)
	ENDDO
	
	DO i = 1,N-1
		xnorm = 1000.d0
		DO WHILE(xnorm .GT. (Rs - dr) .OR. x1(i+1,j) .LT. dr1 .OR. z1(i+1,j) .GT. -dr1)
			px = (rand(0) - 0.5d0)*2.d0
                        py = (rand(0) - 0.5d0)*2.d0
	        	pz = (rand(0) - 0.5d0)*2.d0
	        	pnorm = DSQRT(px**2.d0 + py**2.d0 + pz**2.d0)
	        	px = Rodl*px/pnorm
	        	py = Rodl*py/pnorm
	        	pz = Rodl*pz/pnorm
	        	x1(i+1,j) = x1(i,j) + px
			y1(i+1,j) = y1(i,j) + py
			z1(i+1,j) = z1(i,j) + pz
			xnorm = DSQRT((x1(i+1,j)/aoc)**2.d0+(y1(i+1,j)/boc)**2.d0+(z1(i+1,j))**2.d0)
			!PRINT*,i,j,x1(i+1,j),y1(i+1,j),z1(i+1,j),xnorm
		ENDDO
		!PRINT*,i+1,j,x1(i+1,j),y1(i+1,j),z1(i+1,j)
		!PAUSE
	ENDDO
	
	
	j = 3
	xnorm = 1000.d0	
	x1(1:N,j) = -1.d0
	y1(1:N,j) = 1000.d0
	z1(1:N,j) = -1.d0
	DO WHILE(xnorm .GT. (Rs - dr) .OR. x1(1,j) .LT. dr1 .OR. z1(1,j) .LT. dr1)
        	x1(1,j) = Rs*(rand(0) - 0.5d0)
                y1(1,j) = Rs*(rand(0) - 0.5d0)
                z1(1,j) = Rs*(rand(0) - 0.5d0) 
		xnorm = DSQRT((x1(1,j)/aoc-R_cx(j))**2.d0+(y1(1,j)/boc-R_cy(j))**2.d0+(z1(1,j)-R_cz(j))**2.d0)
	ENDDO
	
	DO i = 1,N-1
		xnorm = 1000.d0
		DO WHILE(xnorm .GT. (Rs - dr) .OR. x1(i+1,j) .LT. dr1 .OR. z1(i+1,j) .LT. dr1)
			px = (rand(0) - 0.5d0)*2.d0
                        py = (rand(0) - 0.5d0)*2.d0
	        	pz = (rand(0) - 0.5d0)*2.d0
	        	pnorm = DSQRT(px**2.d0 + py**2.d0 + pz**2.d0)
	        	px = Rodl*px/pnorm
	        	py = Rodl*py/pnorm
	        	pz = Rodl*pz/pnorm
	        	x1(i+1,j) = x1(i,j) + px
			y1(i+1,j) = y1(i,j) + py
			z1(i+1,j) = z1(i,j) + pz
			xnorm = DSQRT((x1(i+1,j)/aoc)**2.d0+(y1(i+1,j)/boc)**2.d0+(z1(i+1,j))**2.d0)
			!PRINT*,i,j,x1(i+1,j),y1(i+1,j),z1(i+1,j),xnorm
		ENDDO
		!PRINT*,i+1,j,x1(i+1,j),y1(i+1,j),z1(i+1,j)
		!PAUSE
	ENDDO
	
	
	j = 4
	xnorm = 1000.d0	
	x1(1:N,j) = 1.d0
	y1(1:N,j) = 1000.d0
	z1(1:N,j) = 1.d0
	DO WHILE(xnorm .GT. (Rs - dr) .OR. x1(1,j) .GT. -dr1 .OR. z1(1,j) .GT. -dr1)
        	x1(1,j) = Rs*(rand(0) - 0.5d0)
                y1(1,j) = Rs*(rand(0) - 0.5d0)
                z1(1,j) = Rs*(rand(0) - 0.5d0) 
		xnorm = DSQRT((x1(1,j)/aoc-R_cx(j))**2.d0+(y1(1,j)/boc-R_cy(j))**2.d0+(z1(1,j)-R_cz(j))**2.d0)
	ENDDO
	
	DO i = 1,N-1
		xnorm = 1000.d0
		DO WHILE(xnorm .GT. (Rs - dr) .OR. x1(i+1,j) .GT. -dr1 .OR. z1(i+1,j) .GT. -dr1)
			px = (rand(0) - 0.5d0)*2.d0
                        py = (rand(0) - 0.5d0)*2.d0
	        	pz = (rand(0) - 0.5d0)*2.d0
	        	pnorm = DSQRT(px**2.d0 + py**2.d0 + pz**2.d0)
	        	px = Rodl*px/pnorm
	        	py = Rodl*py/pnorm
	        	pz = Rodl*pz/pnorm
	        	x1(i+1,j) = x1(i,j) + px
			y1(i+1,j) = y1(i,j) + py
			z1(i+1,j) = z1(i,j) + pz
			xnorm = DSQRT((x1(i+1,j)/aoc)**2.d0+(y1(i+1,j)/boc)**2.d0+(z1(i+1,j))**2.d0)
			!PRINT*,i,j,x1(i+1,j),y1(i+1,j),z1(i+1,j),xnorm
		ENDDO
		!PRINT*,i+1,j,x1(i+1,j),y1(i+1,j),z1(i+1,j)
		!PAUSE
	ENDDO
	
END SUBROUTINE Initconfig
