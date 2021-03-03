SUBROUTINE SphereGeo(x1,y1,z1,Fdpx,Fdpy,Fdpz,ux,uy,uz,Rb,N,ind,Rs)    
        IMPLICIT NONE
        INTEGER*8 :: i, k, N, Nc, Ntr,l,ind
        REAL*8 :: x1(N), y1(N), z1(N), xo, yo, zo, Rx, Ry, Rz, Rbx, Rby, Rbz, dot, Rob, RRb, &
                        theta, co, si, Zt, C11, C21,uxij,uyij,uzij,uxyij,uzyij,uxzij,rij, &
            C22, C31, C32, XX, C33, ex, ey, ez, xob, yob, zob, RR, Dxx, Dyy, Dzz, Dxy, Dyx,Fdpx(N), &
            Dyz, Dzy, Dxz, Dzx, x, y, z, L13, L25, L15, pi, RRb3, RRb5, Ro, St_fac, alpha, Rb,Fdpy(N),Fdpz(N),&
            ux(N), uy(N), uz(N), Rs 
        pi = 4.d0*atan(1.d0)     
        !St_fac = 1.d0/(pi*8.d0)
        St_fac = (6.d0*Rb/8.d0)  !Factor to be multiplied due to scaling     
        ux = 0.d0
        uy = 0.d0
        uz = 0.d0
!%%%%%%%%%%%%%%%%%%Greens function inside the sphere (Singular Stokeslet)%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        IF(ind .EQ. 1)THEN
        DO i = 1,N
                DO k = 1,N
                        IF(i == k)THEN
                        Dxx = 1.d0
                        Dyy = 1.d0
                        Dzz = 1.d0
                        Dxy = 0.d0
                        Dyx = 0.d0
                        Dxz = 0.d0
                        Dzx = 0.d0
                        Dzy = 0.d0
                        Dyz = 0.d0
                        uy(i) = uy(i) 
                        ux(i) = ux(i) 
                        uz(i) = uz(i)
                        !uy(i) = uy(i) + (Dyx*Fdpx(k)+Dyy*Fdpy(k)+Dyz*Fdpz(k))
                        !ux(i) = ux(i) + (Dxx*Fdpx(k)+Dxy*Fdpy(k)+Dxz*Fdpz(k))
                        !uz(i) = uz(i) + (Dzx*Fdpx(k)+Dzy*Fdpy(k)+Dzz*Fdpz(k))
                        ELSE
                        x = x1(i)
                        y = y1(i)
                        z = z1(i)
                        XX = DSQRT(x**2.d0 + y**2.d0 + z**2.d0)! Evaluation point bead k
                        xo = x1(k)
                        yo = y1(k)
                        zo = z1(k)
                        Ro = DSQRT(xo**2.d0 + yo**2.d0 + zo**2.d0) !Position of the external forcing on bead i
                        ex = xo/Ro  
                        ey = yo/Ro
                        ez = zo/Ro
                        alpha = 1.d0**2.d0/Ro**2.d0 ! factor to calculate inverse point
                        xob = alpha*xo !xo is the position of the point force
                        yob = alpha*yo !xob is the image position of the point force
                        zob = alpha*zo
                        Rx = x - xo !Rx relative distance from the evaluation point to i
                        Ry = y - yo
                        Rz = z - zo
                        Rbx = x - xob !Rbx relative distance from the evaluation point to image i
                        Rby = y - yob
                        Rbz = z - zob
                        Rob = DSQRT(xob**2.d0 + yob**2.d0 + zob**2.d0)
                        dot = Rbx*ex + Rby*ey + Rbz*ez
                        RR = DSQRT(Rx**2.d0 + Ry**2.d0 + Rz**2.d0)
                        RRb = DSQRT(Rbx**2.d0 + Rby**2.d0 + Rbz**2.d0)
                        co = (ex*x + ey*y + ez*z)/XX
                        si = DSQRT(1.d0 - co**2.d0)
                        theta = ACOS(co)
                        Zt = Rob/XX - co
                        C11 = (1.d0 - 3.d0*Ro**2.d0)/(2.d0*Ro**3.d0)
                        C21 = (1.d0 - Ro**2.d0)/(Ro**4.d0)
                        C22 = ((1.d0 - Ro**2.d0)**2.d0)/(4.d0*Ro**5.d0)
                        C31 = 3.d0*(XX**2.d0 - 1.d0)*(1.d0 - Ro**2.d0)/Ro**2.d0
                        C32 = (3.d0*Ro**2.d0 - 5.d0)/(2.d0*Ro**3.d0)
                        C33 = 3.d0*(Ro**2.d0 - 1.d0)/Ro**3.d0
!***********************************************************************                        
!%%%%%%%%%%%%%%%%%%%%AXISYMMETRIC TERMS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
!***********************************************************************                      
!%%%%%%%%%%%%%%%%%%%%DIRECT STOKESLET TERM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        Dxx =  1.d0/RR + Rx**2.d0/RR**3.d0
                        Dyy =  1.d0/RR + Ry**2.d0/RR**3.d0
                        Dzz =  1.d0/RR + Rz**2.d0/RR**3.d0
                        Dxy =  Rx*Ry/RR**3.d0
                        Dxz =  Rx*Rz/RR**3.d0
                        Dzy =  Ry*Rz/RR**3.d0
                        Dyx =  Dxy
                        Dzx =  Dxz
                        Dyz =  Dzy               
!%%%%%%%%%%%%%%%%%%%%IMAGE STOKESLET TERM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
                        RRb3 = RRb**3.d0            
                        Dxx = Dxx + C11*(ex**2.d0/RRb + dot*Rbx*ex/RRb3)                  
                        Dyy = Dyy + C11*(ey**2.d0/RRb + dot*Rby*ey/RRb3)
                        Dzz = Dzz + C11*(ez**2.d0/RRb + dot*Rbz*ez/RRb3)                                  
		        Dxy = Dxy + C11*(ex*ey/RRb + dot*Rbx*ey/RRb3)
		        Dyx = Dyx + C11*(ex*ey/RRb + dot*Rby*ex/RRb3)
		        Dxz = Dxz + C11*(ex*ez/RRb + dot*Rbx*ez/RRb3)
		        Dzx = Dzx + C11*(ex*ez/RRb + dot*Rbz*ex/RRb3)
		        Dzy = Dzy + C11*(ey*ez/RRb + dot*Rbz*ey/RRb3)
		        Dyz = Dyz + C11*(ez*ey/RRb + dot*Rby*ez/RRb3)
!%%%%%%%%%%%%%%%%%%%%IMAGE STRESSLET TERMS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        RRb5 = RRb**5.d0
                        Dxx = Dxx - C21*((-dot*ex**2.d0 + Rbx*ex + dot*ex*ex)/RRb3 - 3.d0*dot**2.d0*Rbx*ex/RRb5)                  
                        Dyy = Dyy - C21*((-dot*ey**2.d0 + Rby*ey + dot*ey*ey)/RRb3 - 3.d0*dot**2.d0*Rby*ey/RRb5)
                        Dzz = Dzz- C21*((-dot*ez**2.d0 + Rbz*ez + dot*ez*ez)/RRb3 - 3.d0*dot**2.d0*Rbz*ez/RRb5)                                  
		        Dxy = Dxy - C21*((-dot*ex*ey + Rbx*ey + dot*ex*ey)/RRb3 - 3.d0*dot**2.d0*Rbx*ey/RRb5)
		        Dyx = Dyx - C21*((-dot*ex*ey + Rby*ex + dot*ex*ey)/RRb3 - 3.d0*dot**2.d0*Rby*ex/RRb5)
		        Dxz = Dxz - C21*((-dot*ex*ez + Rbx*ez + dot*ex*ez)/RRb3 - 3.d0*dot**2.d0*Rbx*ez/RRb5)
		        Dzx = Dzx - C21*((-dot*ez*ex + Rbz*ex + dot*ez*ex)/RRb3 - 3.d0*dot**2.d0*Rbz*ex/RRb5)
		        Dzy = Dzy - C21*((-dot*ez*ey + Rbz*ey + dot*ez*ey)/RRb3 - 3.d0*dot**2.d0*Rbz*ey/RRb5)
		        Dyz = Dyz - C21*((-dot*ez*ey + Rby*ez + dot*ez*ey)/RRb3 - 3.d0*dot**2.d0*Rby*ez/RRb5)
!%%%%%%%%%%%%%%%%%%%%%%IMAGE DEGENERATE QUADRUPOLE TERMS%%%%%%%%%%%%%%%%%%
                        Dxx = Dxx - C22*(2.d0*ex**2.d0/RRb3 - 6.d0*dot*Rbx*ex/RRb5)                  
                        Dyy = Dyy - C22*(2.d0*ey**2.d0/RRb3 - 6.d0*dot*Rby*ey/RRb5)
                        Dzz = Dzz - C22*(2.d0*ez**2.d0/RRb3 - 6.d0*dot*Rbz*ez/RRb5)                                  
		        Dxy = Dxy - C22*(2.d0*ex*ey/RRb3 - 6.d0*dot*Rbx*ey/RRb5)
		        Dyx = Dyx - C22*(2.d0*ex*ey/RRb3 - 6.d0*dot*Rby*ex/RRb5)
		        Dxz = Dxz - C22*(2.d0*ex*ez/RRb3 - 6.d0*dot*Rbx*ez/RRb5)
		        Dzx = Dzx - C22*(2.d0*ex*ez/RRb3 - 6.d0*dot*Rbz*ex/RRb5)
		        Dzy = Dzy - C22*(2.d0*ey*ez/RRb3 - 6.d0*dot*Rbz*ey/RRb5)
		        Dyz = Dyz - C22*(2.d0*ez*ey/RRb3 - 6.d0*dot*Rby*ez/RRb5)	
!***********************************************************************                        
!%%%%%%%%%%%%%%%%%%%%ASYMMETRIC TERMS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
!***********************************************************************            
!%%%%%%%%%%%%%%%%%%%%IMAGE STOKESLET TERM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  		                
		        Dxx = Dxx + C32*(1.d0/RRb -ex**2.d0/RRb - dot*Rbx*ex/RRb3 + Rbx**2.d0/RRb3)                  
                        Dyy = Dyy + C32*(1.d0/RRb -ey**2.d0/RRb - dot*Rby*ey/RRb3 + Rby**2.d0/RRb3)
                        Dzz = Dzz + C32*(1.d0/RRb -ez**2.d0/RRb - dot*Rbz*ez/RRb3 + Rbz**2.d0/RRb3)                                  
		        Dxy = Dxy + C32*(-ex*ey/RRb - dot*Rbx*ey/RRb3+ Rbx*Rby/RRb3)
		        Dyx = Dyx + C32*(-ex*ey/RRb - dot*Rby*ex/RRb3+ Rbx*Rby/RRb3)
		        Dxz = Dxz + C32*(-ex*ez/RRb - dot*Rbx*ez/RRb3+ Rbx*Rbz/RRb3)
		        Dzx = Dzx + C32*(-ex*ez/RRb - dot*Rbz*ex/RRb3+ Rbx*Rbz/RRb3)
		        Dzy = Dzy + C32*(-ey*ez/RRb - dot*Rbz*ey/RRb3+ Rby*Rbz/RRb3)
		        Dyz = Dyz + C32*(-ez*ey/RRb - dot*Rby*ez/RRb3+ Rbz*Rby/RRb3)
!%%%%%%%%%%%%%%%%%%%%IMAGE STRESSLET TERMS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        Dxx = Dxx + C21*(-Rbx*ex/RRb3 +(3.d0*dot**2.d0*Rbx*ex - 3.d0*dot*Rbx**2.d0)/RRb5 + dot/RRb3)                  
                        Dyy = Dyy + C21*(-Rby*ey/RRb3 +(3.d0*dot**2.d0*Rby*ey - 3.d0*dot*Rby**2.d0)/RRb5 + dot/RRb3)  
                        Dzz = Dzz + C21*(-Rbz*ez/RRb3 +(3.d0*dot**2.d0*Rbz*ez - 3.d0*dot*Rbz**2.d0)/RRb5 + dot/RRb3)                                
		        Dxy = Dxy + C21*((3.d0*dot**2.d0*Rbx*ey - 3.d0*dot*Rbx*Rby)/RRb5 -Rby*ex/RRb3)  
		        Dyx = Dyx + C21*((3.d0*dot**2.d0*Rby*ex - 3.d0*dot*Rby*Rbx)/RRb5 -Rbx*ey/RRb3) 
		        Dxz = Dxz + C21*((3.d0*dot**2.d0*Rbx*ez - 3.d0*dot*Rbx*Rbz)/RRb5 -Rbz*ex/RRb3)   
		        Dzx = Dzx + C21*((3.d0*dot**2.d0*Rbz*ex - 3.d0*dot*Rbz*Rbx)/RRb5 -Rbx*ez/RRb3)  
		        Dzy = Dzy + C21*((3.d0*dot**2.d0*Rbz*ey - 3.d0*dot*Rbz*Rby)/RRb5 -Rby*ez/RRb3) 
		        Dyz = Dyz + C21*((3.d0*dot**2.d0*Rby*ez - 3.d0*dot*Rby*Rbz)/RRb5 -Rbz*ey/RRb3)
!%%%%%%%%%%%%%%%%%%%%%%IMAGE DEGENERATE QUADRUPOLE TERMS%%%%%%%%%%%%%%%%%%
                        Dxx = (Dxx + C22*((-2.d0*ex**2.d0 + 2.d0)/RRb3 + (6.d0*dot*Rbx*ex - 6.d0*Rbx**2.d0)/RRb5))                  
                        Dyy = (Dyy + C22*((-2.d0*ey**2.d0 + 2.d0)/RRb3 + (6.d0*dot*Rby*ey - 6.d0*Rby**2.d0)/RRb5))
                        Dzz = (Dzz + C22*((-2.d0*ez**2.d0 + 2.d0)/RRb3 + (6.d0*dot*Rbz*ez - 6.d0*Rbz**2.d0)/RRb5))                                  
                        Dxy = (Dxy + C22*(-2.d0*ex*ey/RRb3 + (6.d0*dot*Rbx*ey- 6.d0*Rbx*Rby)/RRb5))
		        Dyx = (Dyx + C22*(-2.d0*ex*ey/RRb3 + (6.d0*dot*Rby*ex- 6.d0*Rbx*Rby)/RRb5))
		        Dxz = (Dxz + C22*(-2.d0*ex*ez/RRb3 + (6.d0*dot*Rbx*ez- 6.d0*Rbx*Rbz)/RRb5))
		        Dzx = (Dzx + C22*(-2.d0*ex*ez/RRb3 + (6.d0*dot*Rbz*ex- 6.d0*Rbx*Rbz)/RRb5))
		        Dzy = (Dzy + C22*(-2.d0*ey*ez/RRb3 + (6.d0*dot*Rbz*ey- 6.d0*Rbz*Rby)/RRb5))
		        Dyz = (Dyz + C22*(-2.d0*ez*ey/RRb3 + (6.d0*dot*Rby*ez- 6.d0*Rbz*Rby)/RRb5))
!%%%%%%%%%%%%%%%%%%%%%%IMAGE VECTOR HARMONI CS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        IF(THETA .EQ. 0.d0)THEN
                                L13 = 1.d0/RRb + XX/(2.d0*RRb**2.d0)
                                L15 = 1.d0/(3.d0*RRb**3.d0) + XX/(4.d0*RRb**4.d0)
                                L25 = 1.d0/(2.d0*RRb**2.d0) + (2.d0*XX)/(3.d0*RRb**3.d0) + (XX**2.d0)/(4.d0*RRb**4.d0)
                        ELSEIF(DABS(pi - THETA) .LE. 1.d-8)THEN
                               ! PRINT*, "YES"
                                L13 = 1.d0/RRb - XX/(2.d0*RRb**2.d0)
                                L15 = 1.d0/(3.d0*RRb**3.d0) - XX/(4.d0*RRb**4.d0)
                                L25 = 1.d0/(2.d0*RRb**2.d0) - (2.d0*XX)/(3.d0*RRb3) + (XX**2.d0)/(4.d0*RRb**4.d0)
                        ELSE
                                !PRINT*, "YES1"
                                L13 = 1.d0/RRb + (co/(XX*si**2d0))*(1.d0 - Zt*XX/RRb)
                                L15 = 1.d0/(3.d0*RRb3) - Zt*(co/(3.d0*si**2.d0*RRb3)) + (2.d0*co/(3.d0*&
                                                        XX**3.d0*si**4.d0))*(1.d0 - Zt*XX/RRb)
                                L25 = 2.d0*co*XX/(3.d0*RRb3) + XX*Zt/(3.d0*RRb3) - (1.d0/(XX**2.d0*si**2.d0))&
                                                *(1.d0/3.d0 + Zt*XX/(3.d0*RRb)+ co**2.d0*Zt*XX**3.d0/(3.d0*RRb3))&
                                                 + (2.d0/(3.d0*XX**2.d0*si**4.d0))*(1.d0 - co**2.d0*Zt*XX/RRb)
                        ENDIF
                        dot = x*ex + y*ey + z*ez
                        !PRINT*, C31,L13,pi, THETA
                        Dxx = Dxx + C31*((1.d0-ex*ex)*Ro*L13/2.d0 - (3.d0*Ro*(x - dot*ex)/2.d0)*(L15*x - ex*L25))                  
                        Dyy = Dyy + C31*((1.d0-ey*ey)*Ro*L13/2.d0 - (3.d0*Ro*(y - dot*ey)/2.d0)*(L15*y - ey*L25))
                        Dzz = Dzz + C31*((1.d0-ez*ez)*Ro*L13/2.d0 - (3.d0*Ro*(z - dot*ez)/2.d0)*(L15*z - ez*L25))                                  
		        Dxy = Dxy + C31*(-(ex*ey)*Ro*L13/2.d0 - (3.d0*Ro*(y - dot*ey)/2.d0)*(L15*x - ex*L25))
		        Dyx = Dyx + C31*(-(ex*ey)*Ro*L13/2.d0 - (3.d0*Ro*(x - dot*ex)/2.d0)*(L15*y - ey*L25))
		        Dxz = Dxz + C31*(-(ex*ez)*Ro*L13/2.d0 - (3.d0*Ro*(z - dot*ez)/2.d0)*(L15*x - ex*L25))
		        Dzx = Dzx + C31*(-(ex*ez)*Ro*L13/2.d0 - (3.d0*Ro*(x - dot*ex)/2.d0)*(L15*z - ez*L25))
		        Dyz = Dyz + C31*(-(ey*ez)*Ro*L13/2.d0 - (3.d0*Ro*(z - dot*ez)/2.d0)*(L15*y - ey*L25))
		        Dzy = Dzy + C31*(-(ey*ez)*Ro*L13/2.d0 - (3.d0*Ro*(y - dot*ey)/2.d0)*(L15*z - ez*L25))
!%%%%%%%%%%%%%%%%%%%%%%EXTRA CONSTANT TERM%%%%%%%%%%%%%%%%%%
                        Dxx = Dxx - C33*(1.d0 - ex**2.d0)/RRb                  
                        Dyy = Dyy - C33*(1.d0 - ey**2.d0)/RRb
                        Dzz = Dzz - C33*(1.d0 - ez**2.d0)/RRb                                 
                        Dxy = Dxy - C33*(-ex*ey)/RRb
		        Dyx = Dyx - C33*(-ex*ey)/RRb
		        Dxz = Dxz - C33*(-ex*ez)/RRb
		        Dzx = Dzx - C33*(-ex*ez)/RRb
		        Dzy = Dzy - C33*(-ez*ey)/RRb
		        Dyz = Dyz - C33*(-ez*ey)/RRb
		        Dxx = St_fac*Dxx                
                        Dyy = St_fac*Dyy
                        Dzz = St_fac*Dzz                            
                        Dxy = St_fac*Dxy
                	Dyx = St_fac*Dyx
	                Dxz = St_fac*Dxz
                	Dzx = St_fac*Dzx
	                Dzy = St_fac*Dzy
	                Dyz = St_fac*Dyz
		        uy(i) = uy(i) + (Dyx*Fdpx(k)+Dyy*Fdpy(k)+Dyz*Fdpz(k))/Rs
                        ux(i) = ux(i) + (Dxx*Fdpx(k)+Dxy*Fdpy(k)+Dxz*Fdpz(k))/Rs
                        uz(i) = uz(i) + (Dzx*Fdpx(k)+Dzy*Fdpy(k)+Dzz*Fdpz(k))/Rs      
		        ENDIF
		!PRINT*,i,k,ux(i),uy(i),uz(i)
		!PAUSE
		ENDDO
        ENDDO
!%%%%%%%%%%%%%%%%%%Free space Greens function(RPY form)%%%%%%%%%%%%%%%%%%%%%%%%%%%          
        ELSEIF(ind .EQ. 2)THEN
                DO i = 1,N
                        DO k = 1,N
                                uxij = x1(i) - x1(k)
                                uyij = y1(i) - y1(k)
                                uzij = z1(i) - z1(k)
                                uxyij = (x1(i) - x1(k))*(y1(i) - y1(k))
                                uxzij = (x1(i) - x1(k))*(z1(i) - z1(k))
                                uzyij = (z1(i) - z1(k))*(y1(i) - y1(k))
                                rij = DSQRT(uxij**2.d0 + uyij**2.d0 + uzij**2.d0)
                                IF(rij .GE. 2.d0*Rb)THEN
                                IF(i==k) THEN
                                                Dxx = 1.d0
                                                Dyy = 1.d0
                                                Dzz = 1.d0
                                                Dxz = 0.d0
                                                Dyz = 0.d0
                                                Dzy = 0.d0
                                                Dzx = 0.d0
                                                Dxy = 0.d0
                                                Dyx = 0.d0
                                uy(i) = uy(i) + (Dyx*Fdpx(k)+Dyy*Fdpy(k)+Dyz*Fdpz(k))
                                ux(i) = ux(i) + (Dxx*Fdpx(k)+Dxy*Fdpy(k)+Dxz*Fdpz(k))
                                uz(i) = uz(i) + (Dzx*Fdpx(k)+Dzy*Fdpy(k)+Dzz*Fdpz(k))                
                                ELSE        
                                Dxx = (3.d0*Rb)*((1.d0 + 2.d0*Rb**2/(3.d0*rij**2)) &
                                +(1.d0 - 2.d0*Rb**2/rij**2)*uxij**2/rij**2)/(4.d0 * rij)                    
                                Dxy = (3.d0*Rb)*((1.d0 - 2.d0*Rb**2/rij**2)*uxyij/rij**2)/(4.d0 * rij)                                 
                                Dyx = Dxy
                                Dyy = (3.d0*Rb)*((1.d0 + 2.d0*Rb**2/(3.d0*rij**2)) &
                                +(1.d0 - 2.d0*Rb**2/rij**2)*uyij**2/rij**2)/(4.d0 * rij)  
                                Dxz = (3.d0*Rb)*((1.d0 - 2.d0*Rb**2/rij**2)*uxzij/rij**2)/(4.d0 * rij)                
                                Dyz = (3.d0*Rb)*((1.d0 - 2.d0*Rb**2/rij**2)*uzyij/rij**2)/(4.d0 * rij)                               
                                Dzz = (3.d0*Rb)*((1.d0 + 2.d0*Rb**2/(3.d0*rij**2)) &
                                +(1.d0 - 2.d0*Rb**2/rij**2)*uzij**2/rij**2)/(4.d0 * rij) 	
                                Dzy = Dyz
                                Dzx = Dxz
                                uy(i) = uy(i) + (Dyx*Fdpx(k)+Dyy*Fdpy(k)+Dyz*Fdpz(k))
                                ux(i) = ux(i) + (Dxx*Fdpx(k)+Dxy*Fdpy(k)+Dxz*Fdpz(k))
                                uz(i) = uz(i) + (Dzx*Fdpx(k)+Dzy*Fdpy(k)+Dzz*Fdpz(k)) 			                              
                                ENDIF         
                                ELSEIF(rij .LT. 2.d0*Rb)THEN
                                IF(i==k) THEN
                                Dxx = 1.d0
                                Dyy = 1.d0
                                Dzz = 1.d0
                                Dxz = 0.d0
                                Dyz = 0.d0
                                Dzy = 0.d0
                                Dzx = 0.d0
                                Dxy = 0.d0
                                Dyx = 0.d0
                                uy(i) = uy(i) + (Dyx*Fdpx(k)+Dyy*Fdpy(k)+Dyz*Fdpz(k))
                                ux(i) = ux(i) + (Dxx*Fdpx(k)+Dxy*Fdpy(k)+Dxz*Fdpz(k))
                                uz(i) = uz(i) + (Dzx*Fdpx(k)+Dzy*Fdpy(k)+Dzz*Fdpz(k))
                                ELSE
                                Dxx = (3.d0*Rb)*(rij/(2.d0*Rb)*((8.d0/3.d0 - 3.d0*rij/(4.d0*Rb))&
                                + uxij**2/(4.d0*rij*Rb)))/(4.d0 * rij)                    
                                Dxy = (3.d0*Rb)*(uxyij/(4.d0*rij*Rb))/(4.d0*rij)                                                     
                                Dyz = (3.d0*Rb)*(uzyij/(4.d0*rij*Rb))/(4.d0*rij) 
                                Dxz = (3.d0*Rb)*(uxzij/(4.d0*rij*Rb))/(4.d0*rij) 
                                Dyx = Dxy
                                Dzy = Dyz
                                Dzx = Dxz
                                Dyy = (3.d0*Rb)*(rij/(2.d0*Rb)*((8.d0/3.d0 - 3.d0*rij/(4.d0*Rb))&
                                + uyij**2/(4.d0*rij*Rb)))/(4.d0 * rij) 
                                Dzz = (3.d0*Rb)*(rij/(2.d0*Rb)*((8.d0/3.d0 - 3.d0*rij/(4.d0*Rb))&
                                + uzij**2/(4.d0*rij*Rb)))/(4.d0 * rij) 
                                uy(i) = uy(i) + (Dyx*Fdpx(k)+Dyy*Fdpy(k)+Dyz*Fdpz(k))
                                ux(i) = ux(i) + (Dxx*Fdpx(k)+Dxy*Fdpy(k)+Dxz*Fdpz(k))
                                uz(i) = uz(i) + (Dzx*Fdpx(k)+Dzy*Fdpy(k)+Dzz*Fdpz(k))
                                ENDIF               
                                ENDIF                 
                        ENDDO
                ENDDO                        
        ENDIF				        
END SUBROUTINE SphereGeo

