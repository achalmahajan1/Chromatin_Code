SUBROUTINE Initdip
        USE Variables
        IMPLICIT NONE
                INTEGER*8 :: i,j
                REAL*8 :: con

! Loop over chains
                DO i = 1,Nc

! Loop over links
                        DO j = 1,N

! Assign initial state randomly according to on/off rates
! Calculate time till next switch (Poisson process)
                                con = kon/(kon+koff)
                                        IF(rand(0) .LT. con)THEN
                                                dip(j,i) = 1
                                                trel(j,i) = 0.d0
                                                tswitch(j,i) = -log(rand(0))/koff
                                        ELSE
                                                dip(j,i) = 0
                                                trel(j,i) = 0.d0
                                                tswitch(j,i) = -log(rand(0))/kon
                                        ENDIF
                        ENDDO
                ENDDO
                PRINT*, dip, trel, tswitch
END SUBROUTINE Initdip
