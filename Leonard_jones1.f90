!============================================================================
SUBROUTINE Leonard_jones1(RX1,RY1,RZ1,fx1,fy1,fz1,Ntot,Nt,Rad,Rodl,Rs)
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
        INTEGER*8 :: N,NCELL,MAPSIZ,YCELL,IX,IY,IZ,IMAP, Ntot, iatm, jatm, count1, k, &
                        count2, bead_ex, sum1, ICELL,I,Nt,vtemp
        REAL*8    :: CELLI, RCUT, CELL, RX(Ntot), RY(Ntot), RZ(Ntot), FX1(Ntot), FY1(Ntot), FZ1(Ntot),Rs
        INTEGER*8, PARAMETER :: M = 60
        INTEGER*8 :: JCELL0, JCELL,J, NABOR  
        REAL*8    :: RXI, RYI, RZI, FXIJ, FYIJ, FZIJ, FIJ, RCUTSQ, SIGSQ, FXI, FYI, FZI, SR2, SR6, &
                        VIJ, WIJ,RX1(Ntot), RY1(Ntot), RZ1(Ntot), RIJSQ, RXIJ, RYIJ, RZIJ, Rcut2, &
                        r2,dx,dy,dz,fr2,fr6,fpr,Rad, sigma, sigma2, Rodl, y_comp, x_comp, z_comp, u , dist
        REAL*8, ALLOCATABLE, DIMENSION(:) :: x_temp, y_temp, z_temp, FX, FY, FZ
	INTEGER*8, ALLOCATABLE, DIMENSION(:) :: LIST,HEAD,MAP 
!    *******************************************************************
        sum1 = 0
        sigma = 10.d0*Rad
        sigma2 = sigma * sigma
        Rcut2 = sigma
        bead_ex = 3.d0
        count1 = 0
        NCELL = M*M*M
        MAPSIZ = 13*NCELL
        ALLOCATE(x_temp((Ntot - 1)*bead_ex+Ntot))
        ALLOCATE(y_temp((Ntot - 1)*bead_ex+Ntot))
        ALLOCATE(z_temp((Ntot - 1)*bead_ex+Ntot))
        ALLOCATE(FX((Ntot - 1)*bead_ex+Ntot))
        ALLOCATE(FY((Ntot - 1)*bead_ex+Ntot))
        ALLOCATE(FZ((Ntot - 1)*bead_ex+Ntot))
        ALLOCATE(LIST((Ntot - 1)*bead_ex+Ntot))
        ALLOCATE(HEAD(NCELL))
        ALLOCATE(MAP(MAPSIZ))
        DO j = 1,Ntot
                count1 = count1+1
                x_temp(count1) = RX1(j)
                y_temp(count1) = RY1(j)
                z_temp(count1) = RZ1(j)
        ENDDO      
        DO j = 1,Ntot-1
                IF(MOD(j,Nt) .EQ. 0)CYCLE
                count1 = count1+1
                x_temp(count1) = (RX1(j)*3.d0+RX1(j+1))/4.d0
                y_temp(count1) = (RY1(j)*3.d0+RY1(j+1))/4.d0
                z_temp(count1) = (RZ1(j)*3.d0+RZ1(j+1))/4.d0
        ENDDO
        DO j = 1,Ntot-1
                IF(MOD(j,Nt) .EQ. 0)CYCLE
                count1 = count1+1
                x_temp(count1) = (RX1(j)+RX1(j+1))/2.d0
                y_temp(count1) = (RY1(j)+RY1(j+1))/2.d0
                z_temp(count1) = (RZ1(j)+RZ1(j+1))/2.d0
        ENDDO
        DO j = 1,Ntot-1
                IF(MOD(j,Nt) .EQ. 0)CYCLE
                count1 = count1+1
                x_temp(count1) = (RX1(j)+RX1(j+1)*3.d0)/4.d0
                y_temp(count1) = (RY1(j)+RY1(j+1)*3.d0)/4.d0
                z_temp(count1) = (RZ1(j)+RZ1(j+1)*3.d0)/4.d0
        ENDDO
	RX = RX1
	RY = RY1
	RZ = RZ1
	FX1 = 0.d0
	FY1 = 0.d0
	FZ1 = 0.d0

                DO IZ = 1,M
                        DO IY = 1,M
                                DO IX = 1, M
                                        sum1 = sum1 + 1
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
			!PRINT*,ICELL, x_temp(I), y_temp(I), z_temp(I) 
                END DO
		!PAUSE
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

!    ** ZERO FORCES AND POTENTIAL **

                DO I = 1, N
                        FX(I) = 0.d0
                        FY(I) = 0.d0
                        FZ(I) = 0.d0
                ENDDO


!    ** LOOP OVER ALL CELLS **
                DO ICELL = 1, NCELL
                        I = HEAD(ICELL)
!       ** LOOP OVER ALL MOLECULES IN THE CELL **
1000                            IF ( I .GT. 0 ) THEN
                                        RXI = x_temp(I)
                                        RYI = y_temp(I)
                                        RZI = z_temp(I)
                                        FXI = FX(I)
                                        FYI = FY(I)
                                        FZI = FZ(I)
!          ** LOOP OVER ALL MOLECULES BELOW I IN THE CURRENT CELL **
                                        J = LIST(I)
2000                                              IF ( J .GT. 0 ) THEN
                                                        RXIJ  = RXI - x_temp(J)
                                                        RYIJ  = RYI - y_temp(J)
                                                        RZIJ  = RZI - z_temp(J)
                                                        RIJSQ = DSQRT(RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ) 
                                                        RXIJ = RXIJ/RIJSQ
                                                        RYIJ = RYIJ/RIJSQ
                                                        RZIJ = RZIJ/RIJSQ
                                                                IF ( RIJSQ .LE. RCUT2) THEN                                       
                                                                        fpr = 5.d0/(RIJSQ**3.d0+0.0001d0)
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
                                                        J = LIST(J)
                                                        GO TO 2000
                                                ENDIF
!          ** LOOP OVER NEIGHBOURING CELLS **
                                        JCELL0 = 13 * (ICELL - 1)
                                                DO NABOR = 1, 13
                                                        JCELL = MAP ( JCELL0 + NABOR )
!             ** LOOP OVER ALL MOLECULES IN NEIGHBOURING CELLS **
                                                        J = HEAD(JCELL)
3000                                                            IF ( J .NE. 0 ) THEN
                                                                        RXIJ  = RXI - x_temp(J)
                                                                        RYIJ  = RYI - y_temp(J)
                                                                        RZIJ  = RZI - z_temp(J)
                                                                        RIJSQ = DSQRT(RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ)
                                                                        RXIJ = RXIJ/RIJSQ
                                                                        RYIJ = RYIJ/RIJSQ
                                                                        RZIJ = RZIJ/RIJSQ
                                                                                IF ( RIJSQ .LE. RCUT2 ) THEN                                       
                                                                                        fpr = 5.d0/(RIJSQ**3.d0+0.0001d0)
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
                                                                        J = LIST(J)
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
                k = 0
!    ** INCORPORATE ENERGY FACTORS **
                count1 = 0
                DO j = 1,Ntot
                        count1 = count1+1
                        FX1(j) =  FX1(j) + FX(count1)
                        FY1(j) =  FY1(j) + FY(count1)
                        FZ1(j) =  FZ1(j) + FZ(count1)
                ENDDO      
                DO j = 1,Ntot-1
                IF(MOD(j,Nt) .EQ. 0)CYCLE
                        count1 = count1+1
                        FX1(j) = FX1(j) + 0.75d0*FX(count1)
                        FY1(j) = FY1(j) + 0.75d0*FY(count1)
                        FZ1(j) = FZ1(j) + 0.75d0*FZ(count1)
                        FX1(j+1) = FX1(j+1) + 0.25d0*FX(count1)
                        FY1(j+1) = FY1(j+1) + 0.25d0*FY(count1)
                        FZ1(j+1) = FZ1(j+1) + 0.25d0*FZ(count1)
                ENDDO
                DO j = 1,Ntot-1
                IF(MOD(j,Nt) .EQ. 0)CYCLE
                        count1 = count1+1
                        FX1(j) = FX1(j) + 0.5d0*FX(count1)
                        FY1(j) = FY1(j) + 0.5d0*FY(count1)
                        FZ1(j) = FZ1(j) + 0.5d0*FZ(count1)
                        FX1(j+1) = FX1(j+1) + 0.5d0*FX(count1)
                        FY1(j+1) = FY1(j+1) + 0.5d0*FY(count1)
                        FZ1(j+1) = FZ1(j+1) + 0.5d0*FZ(count1)
                ENDDO
                DO j = 1,Ntot-1
                IF(MOD(j,Nt) .EQ. 0)CYCLE
                        count1 = count1+1
                        FX1(j) = FX1(j) + 0.25d0*FX(count1)
                        FY1(j) = FY1(j) + 0.25d0*FY(count1)
                        FZ1(j) = FZ1(j) + 0.25d0*FZ(count1)
                        FX1(j+1) = FX1(j+1) + 0.75d0*FX(count1)
                        FY1(j+1) = FY1(j+1) + 0.75d0*FY(count1)
                        FZ1(j+1) = FZ1(j+1) + 0.75d0*FZ(count1)
                ENDDO
END SUBROUTINE Leonard_jones1      

INTEGER*8 FUNCTION YCELL(IX,IY,IZ)
    IMPLICIT NONE
    INTEGER*8 :: IX, IY, IZ, M
    M = 60
    YCELL = 1 + MOD ( IX - 1 + M, M ) &
                               + MOD ( IY - 1 + M, M ) * M &
                               + MOD ( IZ - 1 + M, M ) * M * M
    RETURN
END FUNCTION YCELL
