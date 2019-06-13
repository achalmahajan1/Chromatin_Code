SUBROUTINE ActiveForces(x,y,z) !IF
        USE Variables
        IMPLICIT NONE
                INTEGER*8 :: i, j, k
                REAL*8 ::num
                REAL*8 :: x(N,Nc), y(N,Nc), z(N,Nc)

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
                                IF (trel(j,i) .GT. tswitch(j,i))THEN
                                        trel(j,i) = 0.d0
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
                !dip = 1.d0
! Update force density based on dipole states
                DO i = 1,Nc
                        DO j = 2,N-1
                               k = SIGN(1.d0,0.5d0 - rand(0))        
                               IF(k .EQ. -1)THEN
                                        Fdpx(j,i) = Fdpx(j,i)+k*Fm*dip(j,i)*(x(j+1,i) - x(j,i))
                                        Fdpy(j,i) = Fdpy(j,i)+k*Fm*dip(j,i)*(y(j+1,i) - y(j,i))
                                        Fdpz(j,i) = Fdpz(j,i)+k*Fm*dip(j,i)*(z(j+1,i) - z(j,i))         
                                        Fdpxa(j+1,i) = Fdpxa(j+1,i)+Fm*dip(j,i)*(x(j+1,i) - x(j,i))
                                        Fdpya(j+1,i) = Fdpya(j+1,i)+Fm*dip(j,i)*(y(j+1,i) - y(j,i))
                                        Fdpza(j+1,i) = Fdpza(j+1,i)+Fm*dip(j,i)*(z(j+1,i) - z(j,i))
                                ELSE
                                        Fdpx(j,i) = Fdpx(j,i)+k*Fm*dip(j,i)*(x(j,i) - x(j-1,i))
                                        Fdpy(j,i) = Fdpy(j,i)+k*Fm*dip(j,i)*(y(j,i) - y(j-1,i))
                                        Fdpz(j,i) = Fdpz(j,i)+k*Fm*dip(j,i)*(z(j,i) - z(j-1,i))
                                        Fdpxa(j-1,i) = Fdpxa(j-1,i)-Fm*dip(j,i)*(x(j,i) - x(j-1,i))
                                        Fdpya(j-1,i) = Fdpya(j-1,i)-Fm*dip(j,i)*(y(j,i) - y(j-1,i))
                                        Fdpza(j-1,i) = Fdpza(j-1,i)-Fm*dip(j,i)*(z(j,i) - z(j-1,i))

                                ENDIF       
                        ENDDO
                ENDDO
END SUBROUTINE ActiveForces
