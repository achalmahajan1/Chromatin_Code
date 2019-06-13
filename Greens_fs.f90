SUBROUTINE Greens_fs(x0,y0,z0,m,k,GExx,GExy,GExz,GEyx,GEyy,GEyz,GEzx,GEzy,GEzz,alpha1,beta1,&
                        gamma1,n,p1,u,vn,uxel,uyel,uzel,ind)
        IMPLICIT NONE
                INTEGER*8 :: k,i,i1,i2,i3,i4,i5,i6,mint,m,nn,n(6), ind
                REAL*8 :: pi, pi4, alpha1,beta1,gamma1, xi(7), wq(7), eta(7),p1(6,3),ux,uy,uz, &
                          x,x0,y,y0,z,z0,GExx,GExy,GExz,GEyx,GEyy,GEyz,GEzx,GEzy,GEzz,px,py,pz,F0(3),&
                          Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz,Txzx,&
                          Tyxx,Tyyx,Tzxx,Tzxy,Tyzx,Tzyx,Tzyy,Tzzx,Tzzy,Txzy,Txzz,Tyzy,Tyzz,Tzzz, DxDxi,&
                          DyDxi,DzDxi,DxDet,DyDet,DzDet,vnx,vny,vnz,hxi,het,hs,Gxx,Gxy,Gxz,Gyx,Gyy,Gyz,&
                          Gzx,Gzy,Gzz,fc,vn(6,3),cf,u(6,3), uxel, uyel, uzel
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
! Integrates the Green's function over a
! non-singular triangle numbered: k
!
! SYMBOLS:
! -------
! x0, y0, z0 : where the velocity is calculated, x_f, y_f, z_f : point forcing
! mint: order of triangle quadrature
!
! GE_ij: integrated ij component over the element
!----------------------------------------
!----------
! constants
!----------
                        mint = 7
                        nn = 7
                        pi  = 3.14159265358D0
                        pi4 = 4.D0*pi
                        cf = (1.d0/(2.d0*pi4))
!--------
! prepare
!--------

                        GExx = 0.D0
                        GExy = 0.D0
                        GExz = 0.D0
                        GEyx = 0.D0
                        GEyy = 0.D0
                        GEyz = 0.D0
                        GEzx = 0.D0
                        GEzy = 0.D0
                        GEzz = 0.D0
                        
                        uxel = 0.D0
                        uyel = 0.D0
                        uzel = 0.D0
                                IF(m .EQ. k)THEN
                                        nn = 6
                                        mint = 6
                                ENDIF  
        !CALL Gauss_leg(NQ,zz,ww)
                                CALL trgl_quad (nn,xi,eta,wq) !!This part of the code was changed from the normal code

                        i1 = n(1)
                        i2 = n(2)
                        i3 = n(3)
                        i4 = n(4)
                        i5 = n(5)
                        i6 = n(6)                     
                                DO i=1,mint
                CALL sdlp_3d_interp(p1(1,1),p1(1,2),p1(1,3),p1(2,1),p1(2,2),p1(2,3),p1(3,1),p1(3,2),p1(3,3)&
                                ,p1(4,1),p1(4,2),p1(4,3),p1(5,1),p1(5,2),p1(5,3),p1(6,1),p1(6,2),p1(6,3)&
                                ,vn(1,1),vn(1,2),vn(1,3),vn(2,1),vn(2,2),vn(2,3),vn(3,1),&
                                vn(3,2),vn(3,3),vn(4,1),vn(4,2),vn(4,3),vn(5,1),vn(5,2),&
                                vn(5,3),vn(6,1),vn(6,2),vn(6,3),u(1,1),u(1,2),u(1,3),u(2,1)&
                                ,u(2,2),u(2,3),u(3,1),u(3,2),u(3,3),u(4,1),u(4,2),u(4,3),u(5,1)&
                                ,u(5,2),u(5,3),u(6,1),u(6,2),u(6,3),alpha1,beta1,gamma1,xi(i),eta(i)&
                                ,x,y,z,vnx,vny,vnz,ux,uy,uz,hs)
                CALL sgf_3d_fs(x,y,z,x0,y0,z0,Gxx,Gxy,Gxz,Gyx,Gyy,Gyz,Gzx,Gzy,Gzz,px,py,pz,Txxx,Txxy,&
                                Txxz,Tyxy,Tyxz,Tzxz,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz)           
                                        fc = 0.5D0*hs*wq(i)
                                        GExx = GExx + gxx*fc
                                        GExy = GExy + gxy*fc
                                        GExz = GExz + gxz*fc
                                        GEyx = GEyx + gyx*fc
                                        GEyy = GEyy + gyy*fc
                                        GEyz = GEyz + gyz*fc
                                        GEzx = GEzx + gzx*fc
                                        GEzy = GEzy + gzy*fc
                                        GEzz = GEzz + gzz*fc
                
                                        Tyxx = Txxy
                                        Tzxx = Txxz
                                        Tzxy = Tyxz
                                        Tyyx = Txyy
                                        Tzyx = Txyz
                                        Tzyy = Tyyz
                                        Tyzx = Txzy
                                        Tzzx = Txzz
                                        Tzzy = Tyzz

                                        uxel = uxel +((ux*Txxx + uy*Tyxx + uz*Tzxx)*vnx+(ux*Txxy + uy*Tyxy + &
                                                        uz*Tzxy)*vny+(ux*Txxz + uy*Tyxz + uz*Tzxz)*vnz)*fc
                                        uyel = uyel +((ux*Txyx + uy*Tyyx + uz*Tzyx)*vnx+(ux*Txyy + uy*Tyyy + &
                                                        uz*Tzyy)*vny+(ux*Txyz + uy*Tyyz + uz*Tzyz)*vnz)*fc
                                        uzel = uzel +((ux*Txzx + uy*Tyzx + uz*Tzzx)*vnx+(ux*Txzy + uy*Tyzy + &
                                                        uz*Tzzy)*vny+(ux*Txzz + uy*Tyzz + uz*Tzzz)*vnz)*fc    
                                ENDDO
                                IF(ind .EQ. 1)THEN
                                        cf = (1.d0/(pi4))
                                ELSEIF(ind .EQ. 2)THEN
                                        cf = (1.d0/(2.d0*pi4))
                                ENDIF
                uxel = uxel*cf
                uyel = uyel*cf
                uzel = uzel*cf           
!-----
! Done
!-----
END SUBROUTINE Greens_fs
