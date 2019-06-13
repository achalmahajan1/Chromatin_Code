SUBROUTINE BoundaryInt(F_0,x_f,y_f,z_f,Nelm,RM1)
        USE Variables
        IMPLICIT NONE
                INTEGER*8 :: i,Nelm,j,k,l,i1,i2,i3,i4,i5,i6,ik2,k2,ind
                REAL*8 ::  GExx,GExy,GExz,GEyx,GEyy,GEyz,GEzx,GEzy,GEzz,cf,DxDet, DyDet,&
                           DzDet,DxDxi,DyDxi,DzDxi,xi,eta,x_f(N*Nc),y_f(N*Nc),z_f(N*Nc),Gxx,Gxy,&
                           Gxz,Gyx,Gyy,Gyz,Gzx,Gzy,Gzz,px,py,pz,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz,Txyx,Txyy,Txyz,&
                           Tyyy,Tyyz,Tzyz,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz,uxel,uyel,uzel,ux,uy,uz,&
                           ux_slp(N*Nc),uy_slp(N*Nc),uz_slp(N*Nc),ux_dlp(N*Nc),uy_dlp(N*Nc),&
                           uz_dlp(N*Nc),F_0(N*Nc,3),RM1(Nelm*3,Nelm*3),&
                           hs
               REAL*8, ALLOCATABLE, DIMENSION(:,:) :: u,vn,p1,ux1,uy1,uz1,u_0
               REAL*8, ALLOCATABLE, DIMENSION(:) :: x0,y0,z0,vnx0,vny0,vnz0,x_0,y_0,z_0,fx,&
                                                        fy,fz,RHS,F_tr         
                        !For n = Nelm, npts = 2562 & n = Nelm, npts = 2562           
                        !cf = 1.D0/(pi8*visc)                                  
                        xi  = 1.D0/3.D0
                        eta = 1.D0/3.D0
                        cf = (1.d0/(8.d0*pi))
                        ALLOCATE(ux1(Nelm,6))
                        ALLOCATE(uy1(Nelm,6))
                        ALLOCATE(uz1(Nelm,6))
                        ALLOCATE(u_0(Nelm,3))
                        ALLOCATE(p1(6,3))
                        ALLOCATE(u(6,3))
                        ALLOCATE(vn(6,3))
                        ALLOCATE(x0(Nelm))
                        ALLOCATE(y0(Nelm))
                        ALLOCATE(z0(Nelm))
                        ALLOCATE(vnx0(Nelm))
                        ALLOCATE(vny0(Nelm))
                        ALLOCATE(vnz0(Nelm))
                        ALLOCATE(x_0(Nelm))
                        ALLOCATE(y_0(Nelm))
                        ALLOCATE(z_0(Nelm))
                        ALLOCATE(fx(Nelm))
                        ALLOCATE(fy(Nelm))
                        ALLOCATE(fz(Nelm))
                        ALLOCATE(RHS(Nelm*3))
                        ALLOCATE(F_tr(Nelm*3))
                        !PRINT*,"YES",Nelm
                        u_0 = 0.d0
                        u = 0.d0
                        !OPEN(32,FILE = "MatrixInverseG1280.dat", FORM="unformatted")
                        ux1 = 0.d0
                        uy1 = 0.d0
                        uz1 = 0.d0
                        DO i=1,nelm
                                i1 = nn(i,1)
                                p1(1,:) = p(i1,:)
                                i2 = nn(i,2)
                                p1(2,:) = p(i2,:)
                                i3 = nn(i,3)
                                p1(3,:) = p(i3,:)
                                i4 = nn(i,4)
                                p1(4,:) = p(i4,:)
                                i5 = nn(i,5)
                                p1(5,:) = p(i5,:)
                                i6 = nn(i,6)
                                p1(6,:) = p(i6,:)
      CALL interp(p(i1,1),p(i1,2),p(i1,3),p(i2,1),p(i2,2),p(i2,3),p(i3,1),p(i3,2),p(i3,3),p(i4,1),p(i4,2),p(i4,3)&
                     ,p(i5,1),p(i5,2),p(i5,3),p(i6,1),p(i6,2),p(i6,3),alpha(i),beta(i),gamma1(i),xi,eta,x0(i),y0(i),&
                     z0(i),DxDxi,DyDxi,DzDxi,DxDet,DyDet,DzDet,vnx0(i),vny0(i),vnz0(i),hs)
                                x_0(i) = (p(i1,1) + p(i2,1) + p(i3,1))/3.d0
                                y_0(i) = (p(i1,2) + p(i2,2) + p(i3,2))/3.d0
                                z_0(i) = (p(i1,3) + p(i2,3) + p(i3,3))/3.d0   
                                DO j = 1,N*Nc
                 CALL sgf_3d_fs(x_f(j),y_f(j),z_f(j),x_0(i),y_0(i),z_0(i),Gxx,Gxy,Gxz,Gyx,Gyy,Gyz,Gzx,Gzy,Gzz,px,py,pz,Txxx,&
                        Txxy,Txxz,Tyxy,Tyxz,Tzxz,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz,Txzx,Txzy,&
                        Txzz,Tyzy,Tyzz,Tzzz)
                                        u_0(i,1) = u_0(i,1) + (Gxx*F_0(j,1) + Gxy*F_0(j,2) + Gxz*F_0(j,3))*cf  
                                        u_0(i,2) = u_0(i,2) + (Gyx*F_0(j,1) + Gyy*F_0(j,2) + Gyz*F_0(j,3))*cf
                                        u_0(i,3) = u_0(i,3) + (Gzx*F_0(j,1) + Gzy*F_0(j,2) + Gzz*F_0(j,3))*cf
                                                DO l = 1,6
       CALL sgf_3d_fs(p1(l,1),p1(l,2),p1(l,3),x_f(j),y_f(j),z_f(j),Gxx,Gxy,Gxz,Gyx,Gyy,Gyz,Gzx,Gzy,Gzz,px,py,pz,Txxx,&
                                Txxy,Txxz,Tyxy,Tyxz,Tzxz,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz)
                                                        ux1(i,l) = ux1(i,l) + cf*(Gxx*F_0(j,1) + Gxy*F_0(j,2) + Gxz*F_0(j,3))  
                                                        uy1(i,l) = uy1(i,l) + cf*(Gyx*F_0(j,1) + Gyy*F_0(j,2) + Gyz*F_0(j,3))
                                                        uz1(i,l) = uz1(i,l) + cf*(Gzx*F_0(j,1) + Gzy*F_0(j,2) + Gzz*F_0(j,3))
                                                ENDDO   
                                ENDDO
                        ENDDO   
                                                  
                DO i=1,nelm
                        ik2 = 3*i-2;
                        ux = 0.D0
                        uy = 0.D0
                        uz = 0.D0
                                DO j=1,nelm
                                        ind = 1
                                        k2 = 3*j-2;
                                        i1 = nn(j,1)
                                        p1(1,:) = p(i1,:)
                                        vn(1,:) = vna(i1,:)
                                        i2 = nn(j,2)
                                        p1(2,:) = p(i2,:)
                                        vn(2,:) = vna(i2,:)
                                        i3 = nn(j,3)
                                        p1(3,:) = p(i3,:)
                                        vn(3,:) = vna(i3,:)
                                        i4 = nn(j,4)
                                        p1(4,:) = p(i4,:)
                                        vn(4,:) = vna(i4,:)
                                        i5 = nn(j,5)
                                        p1(5,:) = p(i5,:)
                                        vn(5,:) = vna(i5,:)
                                        i6 = nn(j,6)
                                        p1(6,:) = p(i6,:)
                                        vn(6,:) = vna(i6,:)
                                                DO l = 1,6
                                                        u(l,1) = ux1(j,l)
                                                        u(l,2) = uy1(j,l)
                                                        u(l,3) = uz1(j,l)
                                                ENDDO
                                        CALL Greens_fs(x0(i),y0(i),z0(i),i,j,GExx,GExy,GExz,GEyx,&
                                         GEyy,GEyz,GEzx,GEzy,GEzz,alpha(j),beta(j),gamma1(j),nn(j,:),p1,u,vn,uxel,uyel,uzel,ind) 
                                        !ENDIF
                                        !RM(ik2,k2)       = cf*GExx
                                        !Rm(ik2,k2+1)     = cf*GExy
                                        !RM(ik2,k2+2)     = cf*GExz

                                        !RM(ik2+1,k2)     = cf*GEyx
                                        !RM(ik2+1,k2+1)   = cf*GEyy
                                        !RM(ik2+1,k2+2)   = cf*GEyz

                                        !RM(ik2+2,k2)     = cf*GEzx
                                        !RM(ik2+2,k2+1)   = cf*GEzy
                                        !RM(ik2+2,k2+2)   = cf*GEzz   
                                        ux  = ux + uxel
                                        uy  = uy + uyel
                                        uz  = uz + uzel                        
                                ENDDO
                                        RHS(ik2) = u_0(i,1) - ux
                                        RHS(ik2+1) = u_0(i,2) - uy
                                        RHS(ik2+2) = u_0(i,3) - uz
        ENDDO        

!%%%%%%%%%%%%CALCULATING THE INVERSE BY GAUSSIAN EMLIMINATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
               ! CALL Inverse(RM,RM1,3*nelm)
                !DO i = 1,Nelm*Nelm*9
                !        READ(32),RM2(i)
                !ENDDO     
                !l = 0
                !DO i = 1,Nelm*3
                !        DO j = 1,Nelm*3
                !                l = l + 1
                !                RM1(i,j) = RM2(l)
                !        ENDDO
                !ENDDO     
                l = 0 
                F_tr = 0.d0
                !DO i = 1,Nelm*3   
                !        temp = RM1(i,j)
                !        DO j = 1,Nelm*3
                !                F_tr(i) = F_tr(i) + temp*RHS(j)
                !        ENDDO
                !ENDDO                              
                F_tr = MATMUL(RM1,RHS)   
                DO i = 1,nelm
                        ik2 = 3*i-2;
                        fx(i) = F_tr(ik2)
                        fy(i) = F_tr(ik2+1)
                        fz(i) = F_tr(ik2+2)
                ENDDO         
                ux = 0.d0
                uy = 0.d0
                uz = 0.d0 
                ux_slp = 0.d0
                uy_slp = 0.d0
                uz_slp = 0.d0
                ux_dlp = 0.d0
                uy_dlp = 0.d0
                uz_dlp = 0.d0  
!%%%%%%%%%%%%CALCULATING THE VELOCITIES DUE TO THE TRACTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                        ind = 2
                        DO l = 1,N*Nc
                                GExx = 0.d0
                                GExy = 0.d0
                                GExz = 0.d0
                                GEyx = 0.d0
                                GEyy = 0.d0
                                GEyz = 0.d0
                                GEzx = 0.d0
                                GEzy = 0.d0
                                GEzz = 0.d0
                                DO j = 1,nelm
                                       ik2 = 3*j-2;
                                        !IF(i .NE. j)THEN     ! regular element !-----------------------
                                        i1 = nn(j,1)
                                        p1(1,:) = p(i1,:)
                                        vn(1,:) = vna(i1,:)
                                        i2 = nn(j,2)
                                        p1(2,:) = p(i2,:)
                                        vn(2,:) = vna(i2,:)
                                        i3 = nn(j,3)
                                        p1(3,:) = p(i3,:)
                                        vn(3,:) = vna(i3,:)
                                        i4 = nn(j,4)
                                        p1(4,:) = p(i4,:)
                                        vn(4,:) = vna(i4,:)
                                        i5 = nn(j,5)
                                        p1(5,:) = p(i5,:)
                                        vn(5,:) = vna(i5,:)
                                        i6 = nn(j,6)
                                        p1(6,:) = p(i6,:)
                                        vn(6,:) = vna(i6,:)
                                                DO i = 1,6
                                                        u(i,1) = ux1(j,i)
                                                        u(i,2) = uy1(j,i)
                                                        u(i,3) = uz1(j,i)
                                                ENDDO                                     
                                        CALL Greens_fs(x_f(l),y_f(l),z_f(l),i,j,GExx,GExy,GExz,GEyx,&
             GEyy,GEyz,GEzx,GEzy,GEzz,alpha(j),beta(j),gamma1(j),nn(j,:),p1,u,vn,uxel,uyel,uzel,ind)    
                                        ux_dlp(l) = ux_dlp(l) + uxel
                                        uy_dlp(l) = uy_dlp(l) + uyel
                                        uz_dlp(l) = uz_dlp(l) + uzel 
                                        ux = cf*GExx*fx(j) + cf*GExy*fy(j) + cf*GExz*fz(j) 
                                        uy = cf*GEyx*fx(j) + cf*GEyy*fy(j) + cf*GEyz*fz(j)
                                        uz = cf*GEzx*fx(j) + cf*GEzy*fy(j) + cf*GEzz*fz(j)
                                        ux_slp(l) = ux + ux_slp(l)
                                        uy_slp(l) = uy + uy_slp(l)
                                        uz_slp(l) = uz + uz_slp(l)
                                ENDDO
                        ENDDO
!---------------NOW TOTAL VELOCITY AT x_0------------------------------------ 
                        u_1 = 0.d0
                        DO k = 1,N*Nc
                        DO i = 1,N*Nc
                                IF(i .NE. k)THEN
             CALL sgf_3d_fs(x_f(i),y_f(i),z_f(i),x_f(k),y_f(k),z_f(k),Gxx,Gxy,Gxz,Gyx,Gyy,Gyz,Gzx,Gzy,Gzz,px,py,pz,Txxx,&
                                    Txxy,Txxz,Tyxy,Tyxz,Tzxz,Txyx,Txyy,Txyz,Tyyy,Tyyz,Tzyz,Txzx,Txzy,&
                                    Txzz,Tyzy,Tyzz,Tzzz)
                                   u_1(k,1) = u_1(k,1) + (Gxx*F_0(i,1) + Gxy*F_0(i,2) + Gxz*F_0(i,3))*cf 
                                   u_1(k,2) = u_1(k,2) + (Gyx*F_0(i,1) + Gyy*F_0(i,2) + Gyz*F_0(i,3))*cf 
                                   u_1(k,3) = u_1(k,3) + (Gzx*F_0(i,1) + Gzy*F_0(i,2) + Gzz*F_0(i,3))*cf 
                                ELSEIF(i .EQ. k)THEN
                                        u_1(k,1) = u_1(k,1) + F_0(i,1)*cf
                                        u_1(k,2) = u_1(k,2) + F_0(i,2)*cf
                                        u_1(k,3) = u_1(k,3) + F_0(i,3)*cf 
                                        !PRINT*,i,k,u_1(k,1),u_1(k,2),u_1(k,3)
                                ENDIF        
                        ENDDO
                        u_1(k,1) = u_1(k,1) - ux_slp(k) - ux_dlp(k) 
                        u_1(k,2) = u_1(k,2) - uy_slp(k) - uy_dlp(k)
                        u_1(k,3) = u_1(k,3) - uz_slp(k) - uz_dlp(k)
                        ENDDO
                        !CLOSE(32)
END SUBROUTINE BoundaryInt                             
