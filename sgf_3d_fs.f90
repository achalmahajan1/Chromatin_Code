SUBROUTINE sgf_3d_fs(x,y,z,x0,y0,z0,Gxx,Gxy,Gxz,Gyx,Gyy,Gyz,Gzx,Gzy,Gzz,px,py,pz,Txxx,&
                        Txxy,Txxz,Tyxy,Tyxz,Tzxz,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz,Txzx,Txzy,&
                        Txzz,Tyzy,Tyzz,Tzzz)
        IMPLICIT NONE
        REAL*8 :: x,x0,y,y0,z,z0,dxx,dy,dz,dyy,dyz,dzz,r,r3,ri,ri3,Gxx,Gxy,Gxz,Gyy,Gyz,Gzz,&
                        Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz,Txzx,&
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
      Rad_Bead = 0.1d0
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
        Gxx = ri + dxx*ri3
        Gxy =      dxy*ri3
        Gxz =      dxz*ri3
        Gyy = ri + dyy*ri3
        Gyz =      dyz*ri3
        Gzz = ri + dzz*ri3

        Gyx = Gxy
        Gzx = Gxz
        Gzy = Gyz
!%%%%%%%%%%%%ROTNE-PRAGER%%%%%%%%%%%%%%%%
      !IF(r .GE. 2.d0*Rad_bead)THEN
      !  Gxx = ((1.d0 + 2.d0*Rad_Bead**2/(3.d0*r**2)) &
      !          +(1.d0 - 2.d0*Rad_Bead**2/r**2)*dxx/r**2)/(r)
      !  Gxy =      ((1.d0 - 2.d0*Rad_Bead**2/r**2)*dxy/r**2)/(r)
      !  Gxz =      ((1.d0 - 2.d0*Rad_Bead**2/r**2)*dxz/r**2)/(r)
      !  Gyy = ((1.d0 + 2.d0*Rad_Bead**2/(3.d0*r**2)) &
      !          +(1.d0 - 2.d0*Rad_Bead**2/r**2)*dyy/r**2)/(r)
      !  Gyz =      ((1.d0 - 2.d0*Rad_Bead**2/r**2)*dyz/r**2)/(r)
      !  Gzz = ((1.d0 + 2.d0*Rad_Bead**2/(3.d0*r**2)) &
      !          +(1.d0 - 2.d0*Rad_Bead**2/r**2)*dzz/r**2)/(r)
      !  Gyx = Gxy
      !  Gzx = Gxz
      !  Gzy = Gyz
     !ELSEIF(r .LT. 2.d0*Rad_bead)THEN
      !  Gxx = (r/(2.d0*Rad_Bead)*((8.d0/3.d0 - 3.d0*r/(4.d0*Rad_Bead))&
      !          + dxx/(4.d0*r*Rad_Bead)))/(r)
      !  Gxy = (dxy/(4.d0*r*Rad_Bead))/(r)
      !  Gxz = (dxz/(4.d0*r*Rad_Bead))/(r)
      !  Gyy = (r/(2.d0*Rad_Bead)*((8.d0/3.d0 - 3.d0*r/(4.d0*Rad_Bead))&
      !          + dyy/(4.d0*r*Rad_Bead)))/(r)
      !  Gyz = (dyz/(4.d0*r*Rad_Bead))/(r)
      !  Gzz = (r/(2.d0*Rad_Bead)*((8.d0/3.d0 - 3.d0*r/(4.d0*Rad_Bead))&
      !          + dzz/(4.d0*r*Rad_Bead)))/(r)
      !  Gyx = Gxy
      !  Gzx = Gxz
      !  Gzy = Gyz
      !ENDIF
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

