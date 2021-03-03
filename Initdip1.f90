SUBROUTINE Initdip1(trel1,kon,koff,tswitch1,tswitch,dip1,dip,N,Nc,Nsb,Np,patch)
  IMPLICIT NONE
  INTEGER*8 :: i,j,N,Nc,dip1(N,Nc),dip(N,Nc),Nsb,Np,Np1,k,Nt,patch(N,Nc)
  REAL*8 :: con,trel1(N,Nc),tswitch(N,Nc),tswitch1(N,Nc),kon,koff,big
	Np1 = N-Nsb*Np
  Nt = Np+Np1/Nsb
  patch = 0.d0!0.d0 for normal beads, 1.d0 for crosslinks region
  ! Loop over chains
  big = 100000000.d0
  dip1 = 0.d0
  trel1 = 0.d0
  tswitch1 = big

  DO i = 1,Nc
    ! Loop over links
    DO j = 1,N/Nt
        con = kon/(kon+koff)
        patch(j*Nt,i) = 1.d0
        tswitch(j*Nt,i) = big
        dip(j*Nt,i) = 0

        IF(rand(0) .LT. con)THEN
          dip1(j*Nt,i) = 1
          trel1(j*Nt,i) = 0.d0
          tswitch1(j*Nt,i) = big!-log(rand(0))/koff
        ELSE
          dip1(j*Nt,i) = 0
          trel1(j*Nt,i) = 0.d0
          tswitch1(j*Nt,i) = -log(rand(0))/kon
        ENDIF

        DO k = 2,Np
        dip(j*Nt-k+1,i) = 0
        tswitch(j*Nt-k+1,i) = big
        patch(j*Nt-k+1,i) = 1.d0

        IF(rand(0) .LT. con)THEN
          dip1(j*Nt-k+1,i) = 1
          trel1(j*Nt-k+1,i) = 0.d0
          tswitch1(j*Nt-k+1,i) = big!-log(rand(0))/koff
        ELSE
          dip1(j*Nt-k+1,i) = 0
          trel1(j*Nt-k+1,i) = 0.d0
          tswitch1(j*Nt-k+1,i) = -log(rand(0))/kon
        ENDIF

        ENDDO
    ENDDO
  ENDDO
!  STOP
END SUBROUTINE Initdip1
