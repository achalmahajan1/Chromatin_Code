!============================================================================
SUBROUTINE Leonard_jones(RX1,RY1,RZ1,Rad,Rodl)
        USE Variables
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
        INTEGER*8 :: M, NCELL, MAPSIZ, YCELL, IX, IY, IZ, IMAP, Ntot, iatm, jatm, count1, k, &
                        count2, bead_ex, sum1, ICELL, I, vtemp,Ntot1
        REAL*8    :: CELLI, RCUT, CELL
        PARAMETER (M = 20, NCELL = M*M*M)
        PARAMETER (MAPSIZ = 13 * NCELL) 
        INTEGER*8 :: HEAD(NCELL), MAP(MAPSIZ), JCELL0, JCELL,J, NABOR  
        REAL*8    :: RXI, RYI, RZI, FXIJ, FYIJ, FZIJ, FIJ, RCUTSQ, SIGSQ, FXI, FYI, FZI, SR2, SR6, &
                        VIJ, WIJ,RX1(N,Nc), RY1(N,Nc), RZ1(N,Nc), RIJSQ, RXIJ, RYIJ, RZIJ, Rcut2, &
                        r2,dx,dy,dz,fr2,fr6,fpr,Rad, sigma, sigma2, Rodl, y_comp, x_comp, z_comp, u , dist
        REAL*8, ALLOCATABLE, DIMENSION(:) :: x_temp, y_temp, z_temp, FX, FY, FZ
	INTEGER*8, ALLOCATABLE, DIMENSION(:) :: LIST
!    *******************************************************************
        sum1 = 0
        sigma = 10.d0*Rad
        sigma2 = sigma * sigma
        Rcut2 = sigma
        bead_ex = 3.d0
        count1 = 0
        Ntot = N*Nc
        ALLOCATE(x_temp((Ntot - 1)*bead_ex+Ntot))
        ALLOCATE(y_temp((Ntot - 1)*bead_ex+Ntot))
        ALLOCATE(z_temp((Ntot - 1)*bead_ex+Ntot))
        ALLOCATE(FX((Ntot - 1)*bead_ex+Ntot))
        ALLOCATE(FY((Ntot - 1)*bead_ex+Ntot))
        ALLOCATE(FZ((Ntot - 1)*bead_ex+Ntot))
        ALLOCATE(LIST((Ntot - 1)*bead_ex+Ntot))
        count1 = 0
        DO i =1,Nc
        DO j = 1,N
                count1 = count1+1
                x_temp(count1) = RX1(j,i)
                y_temp(count1) = RY1(j,i)
                z_temp(count1) = RZ1(j,i)
        ENDDO      
        DO j = 1,N-1
                count1 = count1+1
                x_temp(count1) = (RX1(j,i)*3.d0+RX1(j+1,i))/4.d0
                y_temp(count1) = (RY1(j,i)*3.d0+RY1(j+1,i))/4.d0
                z_temp(count1) = (RZ1(j,i)*3.d0+RZ1(j+1,i))/4.d0
        ENDDO
        DO j = 1,N-1
                count1 = count1+1
                x_temp(count1) = (RX1(j,i)+RX1(j+1,i))/2.d0
                y_temp(count1) = (RY1(j,i)+RY1(j+1,i))/2.d0
                z_temp(count1) = (RZ1(j,i)+RZ1(j+1,i))/2.d0
        ENDDO
        DO j = 1,N-1
                count1 = count1+1
                x_temp(count1) = (RX1(j,i)+RX1(j+1,i)*3.d0)/4.d0
                y_temp(count1) = (RY1(j,i)+RY1(j+1,i)*3.d0)/4.d0
                z_temp(count1) = (RZ1(j,i)+RZ1(j+1,i)*3.d0)/4.d0
        ENDDO
        ENDDO
	Fljx = 0.d0
	Fljy = 0.d0
	Fljz = 0.d0
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
         Ntot1 = count1  
                DO  ICELL = 1, NCELL
                        HEAD(ICELL) = 0
                ENDDO
        CELLI = REAL ( M )/(2.d0*Rs)
        CELL  = 1.d0 / CELLI

!    ** SORT ALL ATOMS **
                DO  I = 1,Ntot1
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
                DO I = 1, Ntot1
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
                                                                        fpr = 2.d0/(RIJSQ**3.d0+0.0001d0)
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
                                                                                        fpr = 2.d0/(RIJSQ**3.d0+0.0001d0)
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
!    ** INCORPORATE ENERGY FACTORS **
        count1 = 0
        DO i = 1,Nc
                DO j = 1,N
                        count1 = count1+1
                        Fljx(j,i) =  Fljx(j,i) + FX(count1)
                        Fljy(j,i) =  Fljy(j,i) + FY(count1)
                        Fljz(j,i) =  Fljz(j,i) + FZ(count1)
                ENDDO      
                DO j = 1,N-1
                        count1 = count1+1
                        Fljx(j,i) = Fljx(j,i) + 0.75d0*FX(count1)
                        Fljy(j,i) = Fljy(j,i) + 0.75d0*FY(count1)
                        Fljz(j,i) = Fljz(j,i) + 0.75d0*FZ(count1)
                        Fljx(j+1,i) = Fljx(j+1,i) + 0.25d0*FX(count1)
                        Fljy(j+1,i) = Fljy(j+1,i) + 0.25d0*FY(count1)
                        Fljz(j+1,i) = Fljz(j+1,i) + 0.25d0*FZ(count1)
                ENDDO
                DO j = 1,N-1
                        count1 = count1+1
                        Fljx(j,i) = Fljx(j,i) + 0.5d0*FX(count1)
                        Fljy(j,i) = Fljy(j,i) + 0.5d0*FY(count1)
                        Fljz(j,i) = Fljz(j,i) + 0.5d0*FZ(count1)
                        Fljx(j+1,i) = Fljx(j+1,i) + 0.5d0*FX(count1)
                        Fljy(j+1,i) = Fljy(j+1,i) + 0.5d0*FY(count1)
                        Fljz(j+1,i) = Fljz(j+1,i) + 0.5d0*FZ(count1)
                ENDDO
                DO j = 1,N-1
                        count1 = count1+1
                        Fljx(j,i) = Fljx(j,i) + 0.25d0*FX(count1)
                        Fljy(j,i) = Fljy(j,i) + 0.25d0*FY(count1)
                        Fljz(j,i) = Fljz(j,i) + 0.25d0*FZ(count1)
                        Fljx(j+1,i) = Fljx(j+1,i) + 0.75d0*FX(count1)
                        Fljy(j+1,i) = Fljy(j+1,i) + 0.75d0*FY(count1)
                        Fljz(j+1,i) = Fljz(j+1,i) + 0.75d0*FZ(count1)
                ENDDO
        ENDDO
        PRINT*,Fljx
END SUBROUTINE Leonard_jones      

INTEGER*8 FUNCTION YCELL(IX,IY,IZ)
    IMPLICIT NONE
    INTEGER*8 :: IX, IY, IZ, M
    M = 20
    YCELL = 1 + MOD ( IX - 1 + M, M ) &
                               + MOD ( IY - 1 + M, M ) * M &
                               + MOD ( IZ - 1 + M, M ) * M * M
    RETURN
END FUNCTION YCELL
  
