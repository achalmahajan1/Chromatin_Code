SUBROUTINE sgf_3d_fs(x,y,z,x0,y0,z0,Rb,Gxx,Gxy,Gxz,Gyx,Gyy,Gyz,Gzx,Gzy,Gzz,px,py,pz,Txxx,&
                        Txxy,Txxz,Tyxy,Tyxz,Tzxz,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz,Txzx,Txzy,&
                        Txzz,Tyzy,Tyzz,Tzzz)
        IMPLICIT NONE
        REAL*8 :: x,x0,y,y0,z,z0,dxx,dy,dz,dyy,dyz,dzz,r,r3,ri,ri3,Gxx,Gxy,Gxz,Gyy,Gyz,Gzz,pi,&
                        Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz,Txzx,Rb(2),&
                        Txzy,Txzz,Tyzy,Tyzz,Tzzz,dxy,dxz,Gyx,Gzx,Gzy,px,py,pz,cf,dx,Rad_Bead
!========================================
! FDLIB, BEMLIB
!
! Copyright by C. Pozrikidis, 1999
! All rights reserved
!
! This program is to be used only under the
! stipulations of the licensing agreement
!========================================

!---------------------------------------
! Free-space Green's function: Stokeslet
!
! Pozrikidis (1992, p. 23)
!
! Iopt =  1 generates only the Green's function
!      ne 1 generates the Green's function,
!           pressure, and stress
!---------------------------------------
      pi = 4.d0*atan(1.d0)
      dx = x-x0
      dy = y-y0
      dz = z-z0

      dxx = dx*dx
      dxy = dx*dy
      dxz = dx*dz
      dyy = dy*dy
      dyz = dy*dz
      dzz = dz*dz

      r = DSQRT(dxx+dyy+dzz)
      r3  = r*r*r
      ri  = 1.0D0/r
      ri3 = 1.0D0/r3
!%%%%%%%%%%%%OSEEN-BURGER%%%%%%%%%%%%%%%%
        !Gxx = ri + dxx*ri3
        !Gxy =      dxy*ri3
        !Gxz =      dxz*ri3
        !Gyy = ri + dyy*ri3
        !Gyz =      dyz*ri3
        !Gzz = ri + dzz*ri3

        !Gyx = Gxy
        !Gzx = Gxz
        !Gzy = Gyz 
!%%%%%%%%%%%%%%%%%%%ROTNE PARGER DIFFERENT SIZE PARTICLES%%%%%%%%%%%%%%%%%%%%%%
        IF(r .GE. (Rb(1)+Rb(2)))THEN    
                Gxx = ((1.d0+(Rb(1)+Rb(2))**2.d0/(3.d0*r**2)) &
                +(1.d0-(Rb(1)+Rb(2))**2.d0/r**2)*dxx/r**2)/(r)
                Gxy = ((1.d0-(Rb(1)+Rb(2))**2.d0/r**2)*dxy/r**2)/(r)
                Gxz = ((1.d0-(Rb(1)+Rb(2))**2.d0/r**2)*dxz/r**2)/(r)
                Gyy = ((1.d0+(Rb(1)+Rb(2))**2.d0/(3.d0*r**2)) &
                +(1.d0-(Rb(1)+Rb(2))**2.d0/r**2)*dyy/r**2)/(r)
                Gyz = ((1.d0-(Rb(1)+Rb(2))**2.d0/r**2)*dyz/r**2)/(r)
                Gzz = ((1.d0+(Rb(1)+Rb(2))**2.d0/(3.d0*r**2)) &
                +(1.d0-(Rb(1)+Rb(2))**2.d0/r**2)*dzz/r**2)/(r)
                Gyx = Gxy
                Gzx = Gxz
                Gzy = Gyz		                                       
        ELSEIF(r .LT. (Rb(1)+Rb(2)) .AND. r .GT. ABS(Rb(1)-Rb(2)))THEN         
                Gxx = (8.d0*pi/(6.d0*pi*Rb(1)*Rb(2)))*((16.d0*r**3*(Rb(1)+Rb(2))-&
                ((Rb(1)-Rb(2))**2+3.d0*r**2)**2)/(32.d0*r**3)+3.d0*((Rb(1)-Rb(2))**2-r**2)**2*dxx/(32.d0*r**5))                  
                Gxy = (8.d0*pi/(6.d0*pi*Rb(1)*Rb(2)))*(3.d0*((Rb(1)-Rb(2))**2-r**2)**2*dxy/(32.d0*r**5))                                                      
                Gyz = (8.d0*pi/(6.d0*pi*Rb(1)*Rb(2)))*(3.d0*((Rb(1)-Rb(2))**2-r**2)**2*dyz/(32.d0*r**5))   
                Gxz = (8.d0*pi/(6.d0*pi*Rb(1)*Rb(2)))*(3.d0*((Rb(1)-Rb(2))**2-r**2)**2*dxz/(32.d0*r**5))  
                Gyx = Gxy
                Gzy = Gyz
                Gzx = Gxz
                Gyy = (8.d0*pi/(6.d0*pi*Rb(1)*Rb(2)))*((16.d0*r**3*(Rb(1)+Rb(2))-&
                ((Rb(1)-Rb(2))**2+3.d0*r**2)**2)/(32.d0*r**3)+3.d0*((Rb(1)-Rb(2))**2-r**2)**2*dyy/(32.d0*r**5))      
                Gzz = (8.d0*pi/(6.d0*pi*Rb(1)*Rb(2)))*((16.d0*r**3*(Rb(1)+Rb(2))-&
                ((Rb(1)-Rb(2))**2+3.d0*r**2)**2)/(32.d0*r**3)+3.d0*((Rb(1)-Rb(2))**2-r**2)**2*dzz/(32.d0*r**5))
        ELSEIF(r .LT. ABS(Rb(1)-Rb(2)))THEN    
                Gxx = 8.d0*pi/(6.d0*pi*MAX(Rb(1),Rb(2)))
                Gyy = 8.d0*pi/(6.d0*pi*MAX(Rb(1),Rb(2)))
                Gzz = 8.d0*pi/(6.d0*pi*MAX(Rb(1),Rb(2)))
                Gxz = 0.d0
                Gzx = Gxz
                Gyz = 0.d0
                Gzy = Gyz
                Gxy = 0.d0
                Gyx = Gxy                     
        ENDIF          
      !PRINT*,"Greens inside",x,x0,y,y0

!---
! compute the stress tensor and the pressure
!---

      cf = -6.0D0/(r3*r*r)

      Txxx = dxx*dx * cf
      Txxy = dxy*dx * cf
      Txxz = dxz*dx * cf
      Tyxy = dyy*dx * cf
      Tyxz = dyz*dx * cf
      Tzxz = dzz*dx * cf

      Txyx = Txxy
      Txyy = Tyxy
      Txyz = Tyxz
      Tyyy = dyy*dy * cf
      Tyyz = dyz*dy * cf
      Tzyz = dzz*dy * cf

      Txzx = Txxz
      Txzy = Tyxz
      Txzz = Tzxz
      Tyzy = dyy*dz * cf
      Tyzz = dyz*dz * cf
      Tzzz = dzz*dz * cf

!---------
! pressure
!---------

      cf = 2.0D0*ri3
      px = dx * cf
      py = dy * cf
      pz = dz * cf


!-----
! done
!-----
END SUBROUTINE sgf_3d_fs

