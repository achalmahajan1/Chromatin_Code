SUBROUTINE Spring
        USE Variables
        IMPLICIT NONE
        INTEGER*8 :: j,i
        REAL*8 :: Qbar
                DO j = 1,Nc
                        DO i = 1,N-1
                                Qx(i,j) = x_temp1(i,j) - x_temp1(i+1,j)
                	        Qy(i,j) = y_temp1(i,j) - y_temp1(i+1,j) 
        		        Qz(i,j) = z_temp1(i,j) - z_temp1(i+1,j)
			        Qo(i,j) = Qx(i,j)**2.d0 + Qy(i,j)**2.d0 + Qz(i,j)**2.d0
			ENDDO
		        Qbar = DSQRT(Qo(1,j))
        		Fbx(1,j) = -Hc*(Qx(1,j)/Qbar)/(1-(Qo(1,j)/Bpar))
        		Fby(1,j) = -Hc*(Qy(1,j)/Qbar)/(1-(Qo(1,j)/Bpar))
        		Fbz(1,j) = -Hc*(Qz(1,j)/Qbar)/(1-(Qo(1,j)/Bpar))
                        DO i = 2,N-1 
                        Fbx(i,j) = -Hc*(Qx(i,j)/DSQRT(Qo(i,j)))/(1-(Qo(i,j)/Bpar))-&
                                                        Hc*(-Qx(i-1,j)/DSQRT(Qo(i-1,j)))/(1-(Qo(i-1,j)/Bpar)) 
                        Fby(i,j) = -Hc*(Qy(i,j)/DSQRT(Qo(i,j)))/(1-(Qo(i,j)/Bpar))-&
                                                        Hc*(-Qy(i-1,j)/DSQRT(Qo(i-1,j)))/(1-(Qo(i-1,j)/Bpar)) 
                        Fbz(i,j) = -Hc*(Qz(i,j)/DSQRT(Qo(i,j)))/(1-(Qo(i,j)/Bpar))-&
                                                        Hc*(-Qz(i-1,j)/DSQRT(Qo(i-1,j)))/(1-(Qo(i-1,j)/Bpar))		        
		        ENDDO
        		Qbar = DSQRT(Qo(N-1,j))
        		Fbx(N,j) = -Hc*(-Qx(N-1,j)/Qbar)/(1-(Qo(N-1,j)/Bpar)) 
                        Fby(N,j) = -Hc*(-Qy(N-1,j)/Qbar)/(1-(Qo(N-1,j)/Bpar)) 
                        Fbz(N,j) = -Hc*(-Qz(N-1,j)/Qbar)/(1-(Qo(N-1,j)/Bpar))
                ENDDO        				        	
END SUBROUTINE Spring              		
