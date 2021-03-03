SUBROUTINE Greens_fs(x0,y0,z0,m,k,GExx,GExy,GExz,GEyx,GEyy,GEyz,GEzx,GEzy,GEzz,alpha1,beta1,&
                        gamma1,n,p1,u,vn,fx,fy,fz)
        IMPLICIT NONE
                INTEGER*8 :: k,i,i1,i2,i3,i4,i5,i6,mint,m,nn,n(6), ind
                REAL*8 :: pi, pi4, alpha1,beta1,gamma1, xi(12), wq(12), eta(12),p1(6,3),ux,uy,uz, &
                          x,x0,y,y0,z,z0,GExx,GExy,GExz,GEyx,GEyy,GEyz,GEzx,GEzy,GEzz,px,py,pz,F0(3),&
                          Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz,Txzx,&
                          Tyxx,Tyyx,Tzxx,Tzxy,Tyzx,Tzyx,Tzyy,Tzzx,Tzzy,Txzy,Txzz,Tyzy,Tyzz,Tzzz, DxDxi,&
                          DyDxi,DzDxi,DxDet,DyDet,DzDet,vnx,vny,vnz,hxi,het,hs,Gxx,Gxy,Gxz,Gyx,Gyy,Gyz,&
                          Gzx,Gzy,Gzz,fc,vn(6,3),cf,u(6,3),fx,fy,fz, fc1
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
                        fc1 = 0.d0                 
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
                                !PRINT*,i,x,y,z,hs,wq(i)       
                                        fc = 0.5D0*hs*wq(i)
                                        !fc1 = fc1 + fc
                                        GExx = GExx + cf*gxx*fc*fx
                                        GExy = GExy + cf*gxy*fc*fy
                                        GExz = GExz + cf*gxz*fc*fz
                                        GEyx = GEyx + cf*gyx*fc*fx
                                        GEyy = GEyy + cf*gyy*fc*fy
                                        GEyz = GEyz + cf*gyz*fc*fz
                                        GEzx = GEzx + cf*gzx*fc*fx
                                        GEzy = GEzy + cf*gzy*fc*fy
                                        GEzz = GEzz + cf*gzz*fc*fz
                                        !GExx = GExx + cf*gxx
                                        !GExy = GExy + cf*gxy
                                        !GExz = GExz + cf*gxz
                                        !GEyx = GEyx + cf*gyx
                                        !GEyy = GEyy + cf*gyy
                                        !GEyz = GEyz + cf*gyz
                                        !GEzx = GEzx + cf*gzx
                                        !GEzy = GEzy + cf*gzy
                                        !GEzz = GEzz + cf*gzz                                                     
                                ENDDO
                                !fc1 = fc1/mint
                                !GExx = GExx*fc1*fx
                                !GExy = GExy*fc1*fy
                                !GExz = GExz*fc1*fz
                                !GEyx = GEyx*fc1*fx
                                !GEyy = GEyy*fc1*fy
                                !GEyz = GEyz*fc1*fz
                                !GEzx = GEzx*fc1*fx
                                !GEzy = GEzy*fc1*fy
                                !GEzz = GEzz*fc1*fz        
                                          
!-----
! Done
!-----
END SUBROUTINE Greens_fs
