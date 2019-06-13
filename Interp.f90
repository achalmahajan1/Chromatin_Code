SUBROUTINE interp(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,al,be,ga,xi,eta&
                  ,x,y,z,DxDxi,DyDxi,DzDxi,DxDet,DyDet,DzDet,vnx,vny,vnz,hs)
        IMPLICIT NONE
                REAL*8 :: alc, bec, gac, alalc, bebec, gagac, al, be, ga, ph1, ph2, ph3, ph4, ph5, ph6, x,&
                          y, z, x1, y1, z1, x2, y2, z2, x5, y5, z5, x3, y3, z3, x4, y4, z4, x6, y6, z6, dph1,&
                          dph2, dph3, dph4, dph5, dph6, DxDxi, DyDxi, DzDxi, xi, eta, pph1, pph2, pph3, pph4,&
                          pph5, pph6, DxDet, DyDet, DzDet, vnx, vny, vnz, hs
! ===========================================
! FDLIB, BEMLIB
! 
! Copyright by C. Pozrikidis, 1999
! All rights reserved
!
! This program is to be used only under the
! stipulations of the li!ensing agreement
! ===========================================

!----------------------------------------------------
!  Interpolate over an element to compute geometrical 
!  variables including the following:
!
!  1) Position vector
!  2) Tangential ve!tors in the xi and eta dire!tions
!  3) Unit normal vector
!  4) Line and surface metrics
!
!  Set I!hoose = 1 to !ompute the position ve!tor only
!                2 for the position ve!tor
!                  and the rest of the variables
!-----------------------------------------------------------


!--------
! prepare
!--------

      alc = 1.0D0-al
      bec = 1.0D0-be
      gac = 1.0D0-ga

      alalc = al*alc
      bebec = be*bec
      gagac = ga*gac

!-------------------------
! evaluate basis fun!tions
!-------------------------

      ph2 = xi *(xi -al+eta*(al-ga)/gac)/alc
      ph3 = eta*(eta-be+xi *(be+ga-1.0D0)/ga)/bec
      ph4 = xi *(1.0D0-xi-eta)/alalc
      ph5 = xi*eta/gagac
      ph6 = eta*(1.0D0-xi-eta)/bebec
      ph1 = 1.0D0-ph2-ph3-ph4-ph5-ph6

!------------------------------------------
! interpolate the position vector (x, y, z)
!------------------------------------------

      x = x1*ph1 + x2*ph2 + x3*ph3 + x4*ph4 + x5*ph5 + x6*ph6
      y = y1*ph1 + y2*ph2 + y3*ph3 + y4*ph4 + y5*ph5 + y6*ph6
      z = z1*ph1 + z2*ph2 + z3*ph3 + z4*ph4 + z5*ph5 + z6*ph6

!------------------------------
      !if(Ichoose.eq.1) Go to 99
!------------------------------

!---
!  evaluate xi derivatives of basis fun!tions
!---

      dph2 =  (2.0D0*xi-al+eta*(al-ga)/gac)/alc
      dph3 =  eta*(be+ga-1.0D0)/(ga*bec)
      dph4 =  (1.0D0-2.0D0*xi-eta)/alalc
      dph5 =  eta/gagac
      dph6 = -eta/bebec
      dph1 = -dph2-dph3-dph4-dph5-dph6

!---
!  !ompute dx/dxi from xi derivatives of phi
!---

      DxDxi = x1*dph1 + x2*dph2 + x3*dph3 + x4*dph4 + x5*dph5 + x6*dph6
      DyDxi = y1*dph1 + y2*dph2 + y3*dph3 + y4*dph4 + y5*dph5 + y6*dph6
      DzDxi = z1*dph1 + z2*dph2 + z3*dph3 + z4*dph4 + z5*dph5 + z6*dph6

!---
!  evaluate eta derivatives of basis fun!tions
!---

      pph2 =  xi*(al-ga)/(alc*gac)
      pph3 =  (2.0D0*eta-be+xi*(be+ga-1.0D0)/ga)/bec
      pph4 =  -xi/alalc
      pph5 =   xi/gagac
      pph6 =  (1.0D0-xi-2.0D0*eta)/bebec
      pph1 = -pph2-pph3-pph4-pph5-pph6

!---
!  !ompute dx/deta from eta derivatives of phi 
!---

      DxDet = x1*pph1 + x2*pph2 + x3*pph3 + x4*pph4 + x5*pph5 + x6*pph6

      DyDet = y1*pph1 + y2*pph2 + y3*pph3 + y4*pph4 + y5*pph5 + y6*pph6

      DzDet = z1*pph1 + z2*pph2 + z3*pph3 + z4*pph4 + z5*pph5 + z6*pph6

!--
!  normal vector    vn = (DxDxi)x(DxDeta) 
!  surface metric   hs = norm(vn) 
!---

      vnx = DyDxi * DzDet - DyDet * DzDxi
      vny = DzDxi * DxDet - DzDet * DxDxi
      vnz = DxDxi * DyDet - DxDet * DyDxi

      hs = DSQRT(vnx*vnx + vny*vny + vnz*vnz)
!---
!  normalization
!---

      vnx = vnx/hs
      vny = vny/hs
      vnz = vnz/hs

!-----
! done
!-----


END SUBROUTINE interp
