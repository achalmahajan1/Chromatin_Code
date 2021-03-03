!============================================================================
SUBROUTINE Leonard_jones(RX1,RY1,RZ1,fx1,fy1,fz1,jj,dipth,Ntot1,Nc,Rad,Rodl,Rs1,time)
  !----------------------------------------------------------------------------
  !    *******************************************************************
  !    ** ROUTINE TO SET UP LINKED LIST AND THE HEAD OF CHAIN ARRAYS    **
  !    ** TAKEN FROM ALLEN AND TILDESLEY (COMPUTER SIMULATION OF LIQUIDS) **
  !    ** PRINCIPAL VARIABLES:                                          **
  !    **                                                               **
  !    ** INTEGER N                  NUMBER OF ATOMS                    **
  !    ** INTEGER M                  NUMBER OF CELLS IN EACH DIRECTION  **
  !    ** INTEGER NCELL              TOTAL NUMBER OF CELLS (M**3)       **
  !    ** INTEGER LIST(N)            LINKED LIST OF ATOMS               **
  !    ** INTEGER HEAD(NCELL)        HEAD OF CHAIN FOR EACH CELL        **
  !    ** REAL    RX(N),RY(N),RZ(N)  POSITIONS                          **
  !    ** REAL    RCUT               THE CUTOFF DISTANCE FOR THE FORCE  **
  !    **                                                               **
  !    ** USAGE:                                                        **
  !    **                                                               **
  !    ** EACH ATOM IS SORTED INTO ONE OF THE M**3 SMALL CELLS.         **
  !    ** THE FIRST ATOM IN EACH CELL IS PLACED IN THE HEAD ARRAY.      **
  !    ** SUBSEQUENT ATOMS ARE PLACED IN THE LINKED LIST ARRAY.         **
  !    ** THE ROUTINE IS CALLED EVERY TIMESTEP BEFORE THE FORCE ROUTINE.**
  !    *******************************************************************
  IMPLICIT NONE
  INTEGER*8 :: N,NCELL,MAPSIZ,YCELL,IX, IY, IZ, IMAP, Ntot, iatm, jatm, count1, k, &
       count2, bead_ex, ICELL, I, vtemp, Ntot1, Nc, jj(Ntot1*Nc)
  REAL*8    :: CELLI, CELL, FX1(Ntot1,Nc), FY1(Ntot1,Nc), FZ1(Ntot1,Nc),Rs,Rs1
  INTEGER*8, PARAMETER :: M = 80
  INTEGER*8 :: JCELL0, JCELL,J, NABOR,dipth(Ntot1,Nc)
  REAL*8    :: RXI, RYI, RZI, FXIJ, FYIJ, FZIJ, FIJ, SIGSQ, FXI, FYI, FZI, SR2, SR6, &
       VIJ, WIJ,RX1(Ntot1,Nc), RY1(Ntot1,Nc), RZ1(Ntot1,Nc), RIJ, RIJ2, RXIJ, RYIJ, RZIJ, &
       Rcut2a,Rcut2,r2,dx,dy,dz,fr2,fr6,fpr,Rad,sigma,sigma2,Rodl,y_comp,x_comp,z_comp, &
       u,dist,time,phi,phia
  REAL*8, ALLOCATABLE, DIMENSION(:) :: x_temp, y_temp, z_temp, FX, FY, FZ
  INTEGER*8, ALLOCATABLE, DIMENSION(:) :: LIST,dipth1,jj1,MAP,HEAD,chain
  INTEGER*8 :: FIRST_CALL = 1
  SAVE :: MAP,HEAD

  !    *******************************************************************
  Rs = Rs1/0.8d0
  sigma = 10.d0*Rad
  sigma2 = sigma * sigma
  Rcut2 = sigma2
  bead_ex = 3.d0
  count1 = 0
  Ntot = Ntot1*Nc
  NCELL = M*M*M
  MAPSIZ = 13*NCELL
  ALLOCATE(x_temp((Ntot - 1)*bead_ex+Ntot))
  ALLOCATE(y_temp((Ntot - 1)*bead_ex+Ntot))
  ALLOCATE(z_temp((Ntot - 1)*bead_ex+Ntot))
  ALLOCATE(FX((Ntot - 1)*bead_ex+Ntot))
  ALLOCATE(chain((Ntot - 1)*bead_ex+Ntot))
  ALLOCATE(FY((Ntot - 1)*bead_ex+Ntot))
  ALLOCATE(FZ((Ntot - 1)*bead_ex+Ntot))
  ALLOCATE(jj1(Ntot))
  ALLOCATE(dipth1(Ntot))
  ALLOCATE(LIST((Ntot - 1)*bead_ex+Ntot))

  ! The map only needs to be allocated and calculated once, so let's do that
  IF (FIRST_CALL == 1) THEN
    ALLOCATE(HEAD(NCELL))
    ALLOCATE(MAP(MAPSIZ))
    DO IZ = 1,M
      DO IY = 1,M
        DO IX = 1, M
          IMAP = (YCELL(IX,IY,IZ)-1)*13
          MAP(IMAP+1) = YCELL(IX+1,IY,IZ)
          MAP(IMAP+2) = YCELL(IX+1,IY+1,IZ)
          MAP(IMAP+3) = YCELL(IX,IY+1,IZ)
          MAP(IMAP+4) = YCELL(IX-1,IY+1,IZ)
          MAP(IMAP+5) = YCELL(IX+1,IY,IZ-1)
          MAP(IMAP+6) = YCELL(IX+1,IY+1,IZ-1)
          MAP(IMAP+7) = YCELL(IX,IY+1,IZ-1)
          MAP(IMAP+8) = YCELL(IX-1,IY+1,IZ-1)
          MAP(IMAP+9) = YCELL(IX+1,IY,IZ+1)
          MAP(IMAP+10) = YCELL(IX+1,IY+1,IZ+1)
          MAP(IMAP+11) = YCELL(IX,IY+1,IZ+1)
          MAP(IMAP+12) = YCELL(IX-1,IY+1,IZ+1)
          MAP(IMAP+13) = YCELL(IX,IY,IZ+1)
        ENDDO
      ENDDO
    ENDDO
    FIRST_CALL = 0
  ENDIF

  count1 = 0
  k = 0
  DO i =1,Nc
    DO j = 1,Ntot1
      count1 = count1+1
      k = k + 1
      x_temp(count1) = RX1(j,i)
      y_temp(count1) = RY1(j,i)
      z_temp(count1) = RZ1(j,i)
      jj1(k) = jj(k)
      dipth1(k) = dipth(j,i)
      chain(count1) = j
    ENDDO
  ENDDO
  DO i =1,Nc
    DO j = 1,Ntot1-1
      count1 = count1+1
      x_temp(count1) = (RX1(j,i)*3.d0+RX1(j+1,i))/4.d0
      y_temp(count1) = (RY1(j,i)*3.d0+RY1(j+1,i))/4.d0
      z_temp(count1) = (RZ1(j,i)*3.d0+RZ1(j+1,i))/4.d0
      chain(count1) = j
    ENDDO
  ENDDO
  DO i =1,Nc
    DO j = 1,Ntot1-1
      count1 = count1+1
      x_temp(count1) = (RX1(j,i)+RX1(j+1,i))/2.d0
      y_temp(count1) = (RY1(j,i)+RY1(j+1,i))/2.d0
      z_temp(count1) = (RZ1(j,i)+RZ1(j+1,i))/2.d0
      chain(count1) = j
    ENDDO
  ENDDO
  DO i =1,Nc
    DO j = 1,Ntot1-1
      count1 = count1+1
      x_temp(count1) = (RX1(j,i)+RX1(j+1,i)*3.d0)/4.d0
      y_temp(count1) = (RY1(j,i)+RY1(j+1,i)*3.d0)/4.d0
      z_temp(count1) = (RZ1(j,i)+RZ1(j+1,i)*3.d0)/4.d0
      chain(count1) = j
    ENDDO
  ENDDO
  FX1 = 0.d0
  FY1 = 0.d0
  FZ1 = 0.d0

  !    ** ZERO HEAD OF CHAIN ARRAY **
  N = count1
  DO  ICELL = 1, NCELL
    HEAD(ICELL) = 0
  ENDDO
  CELLI = REAL ( M )/(2.d0*Rs)
  CELL  = 1.d0 / CELLI

  !    ** SORT ALL ATOMS **
  DO  I = 1,N
    ICELL = 1 + IDINT ( ( x_temp(I) + Rs) * CELLI ) + IDINT ( ( y_temp(I) + Rs)&
         * CELLI ) * M + IDINT ( ( z_temp(I) + Rs) * CELLI ) * M * M
    LIST(I)     = HEAD(ICELL)
    HEAD(ICELL) = I
  END DO

  !    *******************************************************************
  !    ** ROUTINE TO COMPUTE FORCES AND POTENTIAL USING A LINK LIST     **
  !    **                                                               **
  !    ** PRINCIPAL VARIABLES:                                          **
  !    **                                                               **
  !    ** INTEGER N                  NUMBER OF ATOMS                    **
  !    ** INTEGER M                  NUMBER OF CELLS IN EACH DIRECTION  **
  !    ** INTEGER NCELL              NUMBER OF SMALL CELLS (M**3)       **
  !    ** INTEGER MAPSIZ             SIZE OF CELL-CELL MAP              **
  !    ** INTEGER LIST(N)            THE LINKED LIST                    **
  !    ** INTEGER HEAD(NCELL)        THE HEAD OF CHAIN ARRAY            **
  !    ** INTEGER MAP(MAPSIZ)        LIST OF NEIGHBOURING CELLS         **
  !    ** REAL    RX(N),RY(N),RZ(N)  POSITIONS                          **
  !    ** REAL    FX(N),FY(N),FZ(N)  FORCES                             **
  !    ** REAL    SIGMA              THE LJ LENGTH PARAMETER            **
  !    ** REAL    RCUT               THE CUT-OFF DISTANCE               **
  !    ** REAL    V                  THE POTENTIAL ENERGY               **
  !    ** REAL    W                  THE VIRIAL                         **
  !    **                                                               **
  !    ** USAGE:                                                        **
  !    **                                                               **
  !    ** FORCE IS CALLED IN AN MD PROGRAM TO CALCULATE THE FORCE ON    **
  !    ** EACH ATOM. THE ROUTINE IS WRITTEN FOR A LIQUID OF LENNARD     **
  !    ** JONES ATOMS. SUBROUTINE FORCE REQUIRES A LINKED LIST, SET UP  **
  !    ** USING SUBROUTINE LINKS, AND THE MAP OF THE SMALL CELLS SET UP **
  !    ** USING SUBROUTINE MAPS.                                        **
  !    *******************************************************************
        IF(time .LT. 1.d0)THEN
                phi  = 1.d0
        ELSE
                phi = 5.d0
        ENDIF
        
  !    ** ZERO FORCES AND POTENTIAL **
  DO I = 1, N
    FX(I) = 0.d0
    FY(I) = 0.d0
    FZ(I) = 0.d0
  ENDDO
  !    ** LOOP OVER ALL CELLS **
  !
  !$OMP PARALLEL DO SCHEDULE(AUTO), &
  !$OMP PRIVATE(ICELL, I, J, RXI, RYI, RZI, FXI, FYI, FZI, RIJ, RIJ2, RXIJ, RYIJ, RZIJ, fpr, FXIJ, FYIJ, FZIJ, JCELL0, JCELL), &
  !$OMP SHARED(jj1),&
  !$OMP REDUCTION(+:FX, FY, FZ)
  DO ICELL = 1, NCELL
    I = HEAD(ICELL)
    !       ** LOOP OVER ALL MOLECULES IN THE CELL **
1000 IF ( I .GT. 0 ) THEN
      RXI = x_temp(I)
      RYI = y_temp(I)
      RZI = z_temp(I)
      FXI = FX(I)
      FYI = FY(I)
      FZI = FZ(I)
      !          ** LOOP OVER ALL MOLECULES BELOW I IN THE CURRENT CELL **
      J = LIST(I)
2000  IF ( J .GT. 0 ) THEN
        RXIJ  = RXI - x_temp(J)
        RYIJ  = RYI - y_temp(J)
        RZIJ  = RZI - z_temp(J)
        RIJ2 = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ
        RCUT2a = RCUT2
        phia = phi
        IF(i .LE. N .AND. j .LE. N)THEN
         IF(chain(i) .EQ. chain(j) .OR. chain(i) .EQ. chain(j)+1 .OR. chain(i) .EQ. chain(j)-1)THEN
                RCUT2a = 0.2d0*RCUT2
                phia = 0.01*phi
         ENDIF
        ELSEIF(i .LE. N  .AND. j .GT. N)THEN
         IF(chain(i) .EQ. chain(j) .OR. chain(i) .EQ. chain(j)+1)GO TO 1200
        ELSEIF(i .GT. N  .AND. j .LE. N)THEN
         IF(chain(i) .EQ. chain(j) .OR. chain(i) .EQ. chain(j)-1)GO TO 1200
        ELSEIF(i .GT. N  .AND. j .GT. N)THEN
         IF(chain(i) .EQ. chain(j))GO TO 1200
        ENDIF

        IF ( RIJ2 .LE. RCUT2a) THEN
          RIJ = DSQRT(RIJ2)
          RXIJ = RXIJ/RIJ
          RYIJ = RYIJ/RIJ
          RZIJ = RZIJ/RIJ

          fpr = phia/(RIJ*RIJ*RIJ+0.0001d0)
          FXIJ  = fpr * RXIJ
          FYIJ  = fpr * RYIJ
          FZIJ  = fpr * RZIJ
          FXI   = FXI + FXIJ
          FYI   = FYI + FYIJ
          FZI   = FZI + FZIJ
          FX(J) = FX(J) - FXIJ
          FY(J) = FY(J) - FYIJ
          FZ(J) = FZ(J) - FZIJ
        ENDIF
        IF(i .LT. Ntot .AND. j .LT. Ntot .AND. j .NE. i+1 .AND. time .GT. 5.d0 .AND. &
                                j .NE. i-1 .AND. RIJ2 .LE. 2.25d0*RCUT2) THEN
          IF(jj1(i) .EQ.  i .AND. jj1(j) .EQ. j .AND. dipth1(i) .EQ. 1.d0 .AND. dipth1(j) .EQ. 1.d0)THEN
            jj1(i) = j
            jj1(j) = i
          ENDIF
        ENDIF
1200        J = LIST(J)
        GO TO 2000
      ENDIF
      !          ** LOOP OVER NEIGHBOURING CELLS **
      JCELL0 = 13 * (ICELL - 1)
      DO NABOR = 1, 13
        JCELL = MAP ( JCELL0 + NABOR )
        !             ** LOOP OVER ALL MOLECULES IN NEIGHBOURING CELLS **
        J = HEAD(JCELL)
3000    IF ( J .NE. 0 ) THEN
          RXIJ  = RXI - x_temp(J)
          RYIJ  = RYI - y_temp(J)
          RZIJ  = RZI - z_temp(J)
          RIJ2 = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ
          phia = phi
          RCUT2a = RCUT2      
         IF(i .LE. N .AND. j .LE. N)THEN
          IF(chain(i) .EQ. chain(j) .OR. chain(i) .EQ. chain(j)+1 .OR. chain(i) .EQ. chain(j)-1)THEN
                RCUT2a = 0.2d0*RCUT2
                phia = 0.01*phi
          ENDIF
        ELSEIF(i .LE. N  .AND. j .GT. N)THEN
          IF(chain(i) .EQ. chain(j) .OR. chain(i) .EQ. chain(j)+1)GO TO 1300
        ELSEIF(i .GT. N  .AND. j .LE. N)THEN
          IF(chain(i) .EQ. chain(j) .OR. chain(i) .EQ. chain(j)-1)GO TO 1300
        ELSEIF(i .GT. N  .AND. j .GT. N)THEN
          IF(chain(i) .EQ. chain(j))GO TO 1300
        ENDIF
          IF ( RIJ2 .LE. RCUT2a) THEN
            RIJ = DSQRT(RIJ2)
            RXIJ = RXIJ/RIJ
            RYIJ = RYIJ/RIJ
            RZIJ = RZIJ/RIJ

            fpr = phia/(RIJ*RIJ*RIJ+0.0001d0)
            FXIJ  = fpr * RXIJ
            FYIJ  = fpr * RYIJ
            FZIJ  = fpr * RZIJ
            FXI   = FXI + FXIJ
            FYI   = FYI + FYIJ
            FZI   = FZI + FZIJ
            FX(J) = FX(J) - FXIJ
            FY(J) = FY(J) - FYIJ
            FZ(J) = FZ(J) - FZIJ
          ENDIF
         IF(i .LT. Ntot .AND. j .LT. Ntot .AND. j .NE. i+1 .AND. time .GT. 5.d0 &
                        .AND. j .NE. i-1 .AND. RIJ2 .LE. 2.25d0*RCUT2) THEN
   	IF(jj1(i) .EQ.  i .AND. jj1(j) .EQ. j .AND. dipth1(i) .EQ. 1.d0 .AND. dipth1(j) .EQ. 1.d0)THEN
              jj1(i) = j
              jj1(j) = i
   	ENDIF
         ENDIF
1300          J = LIST(J)
          GO TO 3000
        ENDIF
      ENDDO
      FX(I) = FXI
      FY(I) = FYI
      FZ(I) = FZI
      I = LIST(I)
      GO TO 1000
    ENDIF
  END DO
  !    ** INCORPORATE ENERGY FACTORS **
  count1 = 0
  k = 0
  DO i = 1,Nc
    DO j = 1,Ntot1
      count1 = count1+1
      k = k + 1
      FX1(j,i) =  FX1(j,i) + FX(count1)
      FY1(j,i) =  FY1(j,i) + FY(count1)
      FZ1(j,i) =  FZ1(j,i) + FZ(count1)
      jj(k) = jj1(k)
    ENDDO
  ENDDO
  DO i = 1,Nc
    DO j = 1,Ntot1-1
      count1 = count1+1
      FX1(j,i) = FX1(j,i) + 0.75d0*FX(count1)
      FY1(j,i) = FY1(j,i) + 0.75d0*FY(count1)
      FZ1(j,i) = FZ1(j,i) + 0.75d0*FZ(count1)
      FX1(j+1,i) = FX1(j+1,i) + 0.25d0*FX(count1)
      FY1(j+1,i) = FY1(j+1,i) + 0.25d0*FY(count1)
      FZ1(j+1,i) = FZ1(j+1,i) + 0.25d0*FZ(count1)
    ENDDO
  ENDDO
  DO i = 1,Nc
    DO j = 1,Ntot1-1
      count1 = count1+1
      FX1(j,i) = FX1(j,i) + 0.5d0*FX(count1)
      FY1(j,i) = FY1(j,i) + 0.5d0*FY(count1)
      FZ1(j,i) = FZ1(j,i) + 0.5d0*FZ(count1)
      FX1(j+1,i) = FX1(j+1,i) + 0.5d0*FX(count1)
      FY1(j+1,i) = FY1(j+1,i) + 0.5d0*FY(count1)
      FZ1(j+1,i) = FZ1(j+1,i) + 0.5d0*FZ(count1)
    ENDDO
  ENDDO
  DO i = 1,Nc
    DO j = 1,Ntot1-1
      count1 = count1+1
      FX1(j,i) = FX1(j,i) + 0.25d0*FX(count1)
      FY1(j,i) = FY1(j,i) + 0.25d0*FY(count1)
      FZ1(j,i) = FZ1(j,i) + 0.25d0*FZ(count1)
      FX1(j+1,i) = FX1(j+1,i) + 0.75d0*FX(count1)
      FY1(j+1,i) = FY1(j+1,i) + 0.75d0*FY(count1)
      FZ1(j+1,i) = FZ1(j+1,i) + 0.75d0*FZ(count1)
    ENDDO
  ENDDO

  DEALLOCATE(x_temp)
  DEALLOCATE(y_temp)
  DEALLOCATE(z_temp)
  DEALLOCATE(chain)
  DEALLOCATE(FX)
  DEALLOCATE(FY)
  DEALLOCATE(FZ)
  DEALLOCATE(jj1)
  DEALLOCATE(dipth1)
  DEALLOCATE(LIST)
END SUBROUTINE Leonard_jones

INTEGER*8 FUNCTION YCELL(IX,IY,IZ)
  IMPLICIT NONE
  INTEGER*8 :: IX, IY, IZ, M
  M = 80
  YCELL = 1 + MOD ( IX - 1 + M, M ) &
       + MOD ( IY - 1 + M, M ) * M &
       + MOD ( IZ - 1 + M, M ) * M * M
  RETURN
END FUNCTION YCELL
