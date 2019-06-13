SUBROUTINE Initconfig(x1,y1,z1,Rodl)    
!!This Subroutine initializes the position of beads with length of rod scaled to ellipsoidal (cavity) radius
!!Rod length (Cavity) radius is chosen to be the length scale in the problem
        Use Variables
        IMPLICIT NONE
        INTEGER*8 :: i,j
        REAL*8 :: px,py,pz,pnorm,Rs1,dr,xnorm,x1(N,Nc),y1(N,Nc),z1(N,Nc),Rodl,&
                        R_cx(Nc),R_cy(Nc),R_cz(Nc),num
	dr = Rs*0.05d0
	num = 3.d0
	R_cx(1) = Rs/num
	R_cy(1) = -Rs/num
	R_cz(1) = Rs/num
	R_cx(2) = Rs/num
	R_cy(2) = -Rs/num
	R_cz(2) = -Rs/num
	!R_cx(3) = -Rs/num
	!R_cy(3) = Rs/num
	!R_cz(3) = Rs/num
	!R_cx(4) = -Rs/num
	!R_cy(4) = Rs/num
	!R_cz(4) = -Rs/num
	Rs1 = Rs/2.d0
	DO j = 1,Nc
	xnorm = 1000.d0	
	DO WHILE(xnorm .GT. (Rs1 - dr))
        	x1(1,j) = Rs1*(RAND(0) - 0.5d0)*2.d0
                y1(1,j) = Rs1*(RAND(0) - 0.5d0)*2.d0
                z1(1,j) = Rs1*(RAND(0) - 0.5d0)*2.d0 
		xnorm = DSQRT((x1(1,j)/aoc-R_cx(j))**2.d0+(y1(1,j)/boc-R_cy(j))**2.d0+(z1(1,j)-R_cz(j))**2.d0)
	ENDDO
	DO i = 1,N-1
		xnorm = 1000.d0
		DO WHILE(xnorm .GT. (Rs1 - dr))
			px = (RAND(0) - 0.5d0)*2.d0
                        py = (RAND(0) - 0.5d0)*2.d0
	        	pz = (RAND(0) - 0.5d0)*2.d0
	        	pnorm = DSQRT(px**2.d0 + py**2.d0 + pz**2.d0)
	        	px = Rodl*px/pnorm
	        	py = Rodl*py/pnorm
	        	pz = Rodl*pz/pnorm
	        	x1(i+1,j) = x1(i,j) + px
			y1(i+1,j) = y1(i,j) + py
			z1(i+1,j) = z1(i,j) + pz
			xnorm = DSQRT((x1(i+1,j)/aoc-R_cx(j))**2.d0+(y1(i+1,j)/boc-R_cy(j))**2.d0+(z1(i+1,j)-R_cz(j))**2.d0)
			!PRINT*,i,j,x1(i+1,j)
		ENDDO
	ENDDO
	ENDDO
END SUBROUTINE Initconfig
