SUBROUTINE ActiveForces(x,y,z,Fdpx,Fdpy,Fdpz,Fdpxa,Fdpya,Fdpza,N,Nc,kon,koff,dip,trel1,tswitch)
        IMPLICIT NONE
                INTEGER*8 :: N, i, j, k, Nc
                INTEGER*8 :: dip(N,Nc)
                REAL*4 :: rand
                REAL*8 :: kon,koff,Fm, num, leng
                REAL*8 :: trel1(N,Nc),tswitch(N,Nc)
                REAL*8 :: Fdpx(N,Nc), Fdpy(N,Nc), Fdpz(N,Nc), Fdpxa(N,Nc), Fdpya(N,Nc), &
                                Fdpza(N,Nc), x(N,Nc), y(N,Nc), z(N,Nc)

! Active force magnitude (negative for contractile)
                Fm = 0.d0
                leng = 0.d0
! Initialize force vector
                !DO i = 1,Nc
                                Fdpx = 0.d0
                                Fdpy = 0.d0
                                Fdpz = 0.d0
                                Fdpxa = 0.d0
                                Fdpya = 0.d0
                                Fdpza = 0.d0
                !ENDDO
                !j = 3
! Loop over links, update dipole states if needed
                DO i = 1,Nc
                        DO j = 1,N
                                IF (trel1(j,i) .GT. tswitch(j,i))THEN
                                        trel1(j,i) = 0.d0
                                                !PRINT*,"YES", trel1(j), j, dip(j)
                                                num = 0.d0
                                                IF (dip(j,i) .EQ. 0)THEN
                                                        dip(j,i) = 1.d0
                                                        DO WHILE(num .EQ. 0.d0)
                                                                num = rand(0)
                                                        ENDDO
                                                        tswitch(j,i) = -log(num)/koff
                                                ELSE
                                                        DO WHILE(num .EQ. 0.d0)
                                                                num = rand(0)
                                                        ENDDO
                                                        dip(j,i) = 0.d0
                                                        tswitch(j,i) = -log(num)/kon
                                               ENDIF
                                ENDIF
                        ENDDO
                ENDDO
! Update force density based on dipole states
                DO i = 1,Nc
                        DO j = 2,N-1
                               k = SIGN(1.d0,0.5d0 - rand(0))        
                               IF(k .EQ. -1)THEN
                                   leng = DSQRT((x(j+1,i)-x(j,i))**2.d0+(y(j+1,i)-y(j,i))**2.d0+(z(j+1,i)-z(j,i))**2.d0)
                                        Fdpx(j,i) = Fdpx(j,i)+k*Fm*dip(j,i)*(x(j+1,i) - x(j,i))/leng
                                        Fdpy(j,i) = Fdpy(j,i)+k*Fm*dip(j,i)*(y(j+1,i) - y(j,i))/leng
                                        Fdpz(j,i) = Fdpz(j,i)+k*Fm*dip(j,i)*(z(j+1,i) - z(j,i))/leng         
                                        Fdpxa(j+1,i) = Fdpxa(j+1,i)+Fm*dip(j,i)*(x(j+1,i) - x(j,i))/leng
                                        Fdpya(j+1,i) = Fdpya(j+1,i)+Fm*dip(j,i)*(y(j+1,i) - y(j,i))/leng
                                        Fdpza(j+1,i) = Fdpza(j+1,i)+Fm*dip(j,i)*(z(j+1,i) - z(j,i))/leng
                                ELSE
                                   leng = DSQRT((x(j,i)-x(j-1,i))**2.d0+(y(j,i)-y(j-1,i))**2.d0+(z(j,i)-z(j-1,i))**2.d0)     
                                        Fdpx(j,i) = Fdpx(j,i)+k*Fm*dip(j,i)*(x(j,i) - x(j-1,i))/leng
                                        Fdpy(j,i) = Fdpy(j,i)+k*Fm*dip(j,i)*(y(j,i) - y(j-1,i))/leng
                                        Fdpz(j,i) = Fdpz(j,i)+k*Fm*dip(j,i)*(z(j,i) - z(j-1,i))/leng
                                        Fdpxa(j-1,i) = Fdpxa(j-1,i)-Fm*dip(j,i)*(x(j,i) - x(j-1,i))/leng
                                        Fdpya(j-1,i) = Fdpya(j-1,i)-Fm*dip(j,i)*(y(j,i) - y(j-1,i))/leng
                                        Fdpza(j-1,i) = Fdpza(j-1,i)-Fm*dip(j,i)*(z(j,i) - z(j-1,i))/leng

                                ENDIF       
                        ENDDO
                ENDDO
END SUBROUTINE ActiveForces
