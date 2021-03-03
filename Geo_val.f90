SUBROUTINE Geo_val(nelm,npts,vna,n,alpha,beta,gamma1,p1)
        IMPLICIT NONE
                INTEGER*8 :: npts,i1,i2,i3,i4,i5,i6,n(5120,6),nelm,m,k,i,itally(npts)
                REAL*8 :: vna(npts,3),alpha(nelm),beta(nelm),gamma1(nelm),al,be,ga,alc,bec,gac,&
                           xxi(6), eet(6), xi, eta, par,p1(10242,3),vx(6),vy(6),vz(6),DxDxi,DyDxi,DzDxi,DxDet,&
                           DyDet,DzDet,hxi,het,hs,x,y,z
               
!============================================
! FDLIB, BEMLIB
!
! Copyright by C. Pozrikidis, 1999
! All rights reserved
!
! This program is to be used only under the
! stipulations of the licensing agreement
!============================================

!----------------------------------------
! Compute:
!
!  the surface area of the individual elements
!  x, y, and z moments over each element
!
!  the total particle surface area and volume
!
!  the mean curvature of each element
!
!  Also compute the average value of the
!  normal vector at each node.
!  This is done by computing the normal vector 
!  at the six nodes of each triangle, 
!  and then averaging the contributions
!----------------------------------------


!-----------
! initialize
!-----------

        vna(:,1) = 0.0D0
        vna(:,2) = 0.0D0
        vna(:,3) = 0.0D0
        itally = 0

!-----

        DO k=1,nelm
                i1 = n(k,1)
                i2 = n(k,2)
                i3 = n(k,3)
                i4 = n(k,4)
                i5 = n(k,5)
                i6 = n(k,6)
                al  = alpha(k)
                be  = beta(k)
                ga  = gamma1(k)
                alc = 1.0D0-al
                bec = 1.0D0-be
                gac = 1.0D0-ga
                xxi(1) = 0.0D0
                eet(1) = 0.0D0
                xxi(2) = 1.0D0
                eet(2) = 0.0D0
                xxi(3) = 0.0D0
                eet(3) = 1.0D0
                xxi(4) = al
                eet(4) = 0.0D0
                xxi(5) = ga
                eet(5) = gac
                xxi(6) = 0.0D0
                eet(6) = be
!------------------------------------------------------
! compute the average value of the normal vector
!
!         the mean curvature as a contour integral
!         using the nifty formula (4.2.10) 
!         of Pozrikidis (1997)
!------------------------------------------------------
                        DO i=1,6
                                xi  = xxi(i)
                                eta = eet(i)
                               CALL interp(p1(i1,1),p1(i1,2),p1(i1,3),p1(i2,1),p1(i2,2),p1(i2,3),p1(i3,1),&
                   p1(i3,2),p1(i3,3),p1(i4,1),p1(i4,2),p1(i4,3),p1(i5,1),p1(i5,2),p1(i5,3),p1(i6,1),p1(i6,2),p1(i6,3)&
                                ,al,be,ga,xi,eta,x,y,z,DxDxi,DyDxi,DzDxi,DxDet,DyDet,&
                                        DzDet,vx(i),vy(i),vz(i),hs)
                                m = n(k,i)
                                vna(m,1) = vna(m,1) + vx(i)
                                vna(m,2) = vna(m,2) + vy(i)
                                vna(m,3) = vna(m,3) + vz(i)
                                itally(m) = itally(m)+1
                                !PRINT*,p1(i1,1),p1(i2,1),xi,eta
                                !PAUSE
                        ENDDO
                        

        ENDDO
!-----------------------
! averaged normal vector
!-----------------------

      DO i=1,npts

        par = float(itally(i))

        vna(i,1) = vna(i,1)/par
        vna(i,2) = vna(i,2)/par
        vna(i,3) = vna(i,3)/par
        !PRINT*,m,par,vna(i,1),vna(i,2),vna(i,3)
        !PAUSE
        par = sqrt( vna(i,1)**2+vna(i,2)**2+vna(i,3)**2 )

        vna(i,1) = vna(i,1)/par
        vna(i,2) = vna(i,2)/par
        vna(i,3) = vna(i,3)/par
      END DO
        !PAUSE
!-----
! Done
!-----

        
END SUBROUTINE Geo_val
