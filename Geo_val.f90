SUBROUTINE Geo_val(nelm,npts)
        USE Variables
        IMPLICIT NONE
                INTEGER*8 :: npts,i1,i2,i3,i4,i5,i6,nelm,m,k,i,itally(npts)
                REAL*8 :: al,be,ga,alc,bec,gac,&
                           xxi(6), eet(6), xi, eta, par,vx(6),vy(6),vz(6),DxDxi,DyDxi,DzDxi,DxDet,&
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

        DO k = 1,nelm
                i1 = nn(k,1)
                i2 = nn(k,2)
                i3 = nn(k,3)
                i4 = nn(k,4)
                i5 = nn(k,5)
                i6 = nn(k,6)
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
                                CALL interp(p(i1,1),p(i1,2),p(i1,3),p(i2,1),p(i2,2),p(i2,3),p(i3,1),&
                   p(i3,2),p(i3,3),p(i4,1),p(i4,2),p(i4,3),p(i5,1),p(i5,2),p(i5,3),p(i6,1),p(i6,2),p(i6,3)&
                                ,al,be,ga,xi,eta,x,y,z,DxDxi,DyDxi,DzDxi,DxDet,DyDet,&
                                        DzDet,vx(i),vy(i),vz(i),hs)
                                m = nn(k,i)
                                vna(m,1) = vna(m,1) + vx(i)
                                vna(m,2) = vna(m,2) + vy(i)
                                vna(m,3) = vna(m,3) + vz(i)
                                itally(m) = itally(m)+1
                        ENDDO
                        

        ENDDO
!-----------------------
! averaged normal vector
!-----------------------

                DO i=1,npts

                        par = DFLOAT(itally(i))

                        vna(i,1) = vna(i,1)/par
                        vna(i,2) = vna(i,2)/par
                        vna(i,3) = vna(i,3)/par
                        par = DSQRT( vna(i,1)**2+vna(i,2)**2+vna(i,3)**2 )

                        vna(i,1) = vna(i,1)/par
                        vna(i,2) = vna(i,2)/par
                        vna(i,3) = vna(i,3)/par
                END DO
!-----
! Done
!-----

        
END SUBROUTINE Geo_val
