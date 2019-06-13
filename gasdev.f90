! ***********************************************************************
! * GASDEV returns double precision figures from a Gaussian Distribution
! * with mean 0 and variance 1. 'IDUM' is set to a negative number to
! * initialize and then not changed. Reset IDUM every time, with a new
! * negative number (iteration number X -1, for instance) to avoid
! * repetition of saem sequence
! ***********************************************************************
    REAL(KIND=8) FUNCTION gasdev(idum)
		INTEGER, INTENT(INOUT):: idum
		INTEGER :: iset
		REAL(KIND=8) :: fac,gset,rsq,v1,v2
    REAL(KIND=8), EXTERNAL :: ran1

		!common/rndm/idum
		SAVE iset,gset
		DATA iset/0/
		IF(idum .LT. 0) iset=0
		IF(iset .EQ. 0)THEN
  1	 	        v1=2.d0*ran1(idum)-1.d0
		        v2=2.d0*ran1(idum)-1.d0
		        rsq=v1**2.d0+v2**2.d0
		                IF(rsq>=1.d0 .OR. rsq==0.d0)GOTO 1
		        fac=sqrt(-2.d0*log(rsq)/rsq)
		        gset=v1*fac
		        gasdev=v2*fac
		        iset=1
		ELSE
			gasdev=gset
			iset=0
		ENDIF
		RETURN
		END
    !--------------------------------------------------------------------------------------
		REAL(KIND=8)FUNCTION ran1(idum)
		!common/rndm/idum
		INTEGER, INTENT(INOUT):: idum
		INTEGER, PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB
                REAL(KIND=8), PARAMETER ::  AM=1.d0/IM,EPS=1.2d-7,RNMX=1.d0-EPS
		INTEGER :: j,k,iv(NTAB),iy
		SAVE iv,iy
		DATA iv /NTAB*0/, iy /0/
		IF(idum<=0 .OR. iy==0)THEN
			idum=MAX(-idum,1)
			DO j=NTAB+8,1,-1
				k=idum/IQ
				idum=IA*(idum-k*IQ)-IR*k
				        IF(idum<0) idum=idum+IM
				        IF(j<=NTAB) iv(j)=idum
			ENDDO
			iy=iv(1)
		ENDIF
		k=idum/IQ
		idum=IA*(idum-k*IQ)-IR*k
		IF(idum<0) idum=idum+IM
		j=1+iy/NDIV
		iy=iv(j)
		iv(j)=idum
		ran1=min(AM*iy,RNMX)
		RETURN
        END
