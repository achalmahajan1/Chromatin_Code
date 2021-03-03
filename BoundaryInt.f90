SUBROUTINE BoundaryInt(F_0,x_f,y_f,z_f,u_1,p,nn,Nelm,alpha,beta,gamma1,RM1,fmm,Npart)
        use fmmwrapper
        use iso_c_binding
        IMPLICIT NONE
!%%%This code solves the disturbance velocity inside the sphere with a point forcing inside the domain. The integral expression
!%%%involves evaluating a single layer potential. The disturbance velocity is sum of the imposed flow (point stokeslet in this case)
!%%%plus the contribution coming from the boundary to invoke no-slip at the wall. The code does not do well near the wall, and close
!%%%to the point forcing. The number of elements = 5120 (number of subdivisions = 4 of a regular icosahedron) which gives the optimum error in the velocity. Any increase
!%%%in the no of elements slows down the code. FMM is used to speed up the calculation of G_ij F_j with modified quadrature to account
!%%%%for multiple quadrature points.        
                INTEGER*4 :: nsrc, ntrg!, rank, size, ierror, tag, status(MPI_STATUS_SIZE), ntrg1
                INTEGER*8 :: i,N,Nelm,j,l,i1,k,i2,i3,i4,i5,i6,ik2,Npart,nn(Nelm,6),m,k1,nn1, mint
                REAL*8 ::  pi, pi4, cf,Gxx,Gxy,xi(12),wq(12),eta(12),x_f(Npart),y_f(Npart),z_f(Npart),Gxz,&
                           Gyx,Gyy,Gyz,Gzx,Gzy,Gzz,px,py,pz,Txxx,Txxy,Txxz,Tyxy,Tyxz,Tzxz,Txyx,Txyy,Txyz,&
                           Tyyy,Tyyz,Tzyz,Txzx,Txzy,Txzz,Tyzy,Tyzz,Tzzz,ux, uy, uz,gasdev,ux_slp(Npart),&
                           uy_slp(Npart),uz_slp(Npart),p(2*Nelm+2,3),fc,xi1,eta1,F_0(Npart*3),u_1(Npart*3),&
                           x,y,z,hs,RM1(Nelm*3,Nelm*3),vnx1,vny1,vnz1,gamma1(Nelm),alpha(Nelm),&
                           beta(Nelm),idum,Rad_Bead,start,finish,sumx,sumy,sumz,fac
               REAL*8, ALLOCATABLE, DIMENSION(:,:) :: u,RM,vn,p1
               REAL*8, ALLOCATABLE, DIMENSION(:) :: x_0,y_0,z_0,fx,srcValue1,trgValue1,fy,fz,RHS,F_tr,trgCoord,&
                                                    srcCoord,srcValue,trgValue,fc1,srcCoord1,trgCoord1
               TYPE(c_ptr):: fmm                                                
                        !For n = Nelm, npts = 2562 & n = Nelm, npts = 2562           
                        pi  = 3.14159265358D0
                        pi4 = 4.0D0*pi
                        !cf = 1.D0/(pi8*visc)                                  
                        xi1  = 1.D0/3.D0 !position of the centroid in the transformed domain
                        eta1 = 1.D0/3.D0 !position of the centroid in the transformed domain
                        cf = (1.d0/(8.d0*pi)) ! coefficient in front of the Stokeslet, scaled with the viscosity
                        Rad_Bead = 0.1d0 ! Radius of the spherical bead used in the friction coefficient
                        nn1 = 7 !Number of quadrature points to calculate the intgral
                        mint = 7
                        fac = 6.d0*Rad_Bead*pi ! Factor to be multiplied due to scaling
                        ALLOCATE(p1(6,3))
                        ALLOCATE(u(6,3))
                        ALLOCATE(vn(6,3))
                        ALLOCATE(x_0(Nelm))
                        ALLOCATE(y_0(Nelm))
                        ALLOCATE(z_0(Nelm))
                        ALLOCATE(fx(Nelm))
                        ALLOCATE(fy(Nelm))
                        ALLOCATE(fz(Nelm))
                        ALLOCATE(fc1(Nelm))
                        ALLOCATE(RHS(Nelm*3))
                        ALLOCATE(F_tr(Nelm*3))
                        ALLOCATE(trgCoord((Nelm+Npart)*3))
                        ALLOCATE(trgValue((Nelm+Npart)*3))
                        ALLOCATE(srcCoord(Npart*3))
                        ALLOCATE(srcValue(Npart*3))
                        ALLOCATE(trgCoord1(Npart*3))
                        ALLOCATE(trgValue1(Npart*3))
                        ALLOCATE(srcCoord1((Nelm*mint)*3))
                        ALLOCATE(srcValue1((Nelm*mint)*3))
                        !fmm = create_fmm_wrapper(8, 2000, 0, 0, 0) 
                        !call FMM_SetBox(fmm, -10.0d+0, 10.0d+0, -10.0d+0, 10.0d+0, -8.0d+0, 8.0d+0); 
                        idum = -1000    
                        m = 0
                        k1 = 0
                        u = 0.d0
                        fc1 = 0.d0
                        CALL trgl_quad (nn1,xi,eta,wq)
                        DO i=1,nelm !This part just needs to run once and then can be used at every time step
                                !m = m + 1
                                k1 = k1 + 1 
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
                                x_0(i) = (p(i1,1) + p(i2,1) + p(i3,1))/3.d0
                                y_0(i) = (p(i1,2) + p(i2,2) + p(i3,2))/3.d0
                                z_0(i) = (p(i1,3) + p(i2,3) + p(i3,3))/3.d0                                
                                trgCoord(3*(k1 - 1) + 1) = x_0(i)
                                trgCoord(3*(k1 - 1) + 2) = y_0(i)
                                trgCoord(3*(k1 - 1) + 3) = z_0(i)                                                               
                                DO j=1,mint
                                CALL sdlp_3d_interp(p1(1,1),p1(1,2),p1(1,3),p1(2,1),p1(2,2),p1(2,3),p1(3,1),p1(3,2),p1(3,3)&
                                ,p1(4,1),p1(4,2),p1(4,3),p1(5,1),p1(5,2),p1(5,3),p1(6,1),p1(6,2),p1(6,3)&
                                ,vn(1,1),vn(1,2),vn(1,3),vn(2,1),vn(2,2),vn(2,3),vn(3,1),&
                                vn(3,2),vn(3,3),vn(4,1),vn(4,2),vn(4,3),vn(5,1),vn(5,2),&
                                vn(5,3),vn(6,1),vn(6,2),vn(6,3),u(1,1),u(1,2),u(1,3),u(2,1)&
                                ,u(2,2),u(2,3),u(3,1),u(3,2),u(3,3),u(4,1),u(4,2),u(4,3),u(5,1)&
                                ,u(5,2),u(5,3),u(6,1),u(6,2),u(6,3),alpha(i),beta(i),gamma1(i),xi(j),eta(j)&
                                ,x,y,z,vnx1,vny1,vnz1,ux,uy,uz,hs)
                                m = m + 1
                                srcCoord1(3*(m - 1) + 1) = x
                                srcCoord1(3*(m - 1) + 2) = y
                                srcCoord1(3*(m - 1) + 3) = z      
                                fc = 0.5D0*hs*wq(j)
                                fc1(i) = fc1(i) + fc                                                                                          
                                ENDDO
                                fc1(i) = fc1(i)/mint
                        ENDDO  
                        !PAUSE
                        trgValue1 = 0.d0
                        trgValue = 0.d0 
                        DO j = 1,Npart
                                k1 = k1 + 1
                                srcCoord(3*(j - 1) + 1) = x_f(j)
                                srcCoord(3*(j - 1) + 2) = y_f(j)
                                srcCoord(3*(j - 1) + 3) = z_f(j)
                                trgCoord(3*(k1 - 1) + 1) = x_f(j)!+0.01d0
                                trgCoord(3*(k1 - 1) + 2) = y_f(j)!+0.01d0
                                trgCoord(3*(k1 - 1) + 3) = z_f(j)!+0.01d0
                                srcValue(3*(j - 1) + 1) = F_0(3*(j - 1) + 1)
                                srcValue(3*(j - 1) + 2) = F_0(3*(j - 1) + 2)
                                srcValue(3*(j - 1) + 3) = F_0(3*(j - 1) + 3) 
                                trgCoord1(3*(j - 1) + 1) = x_f(j)!+0.01d0
                                trgCoord1(3*(j - 1) + 2) = y_f(j)!+0.01d0
                                trgCoord1(3*(j - 1) + 3) = z_f(j)!+0.01d0
                        ENDDO   
                        ntrg = INT(Nelm+Npart)
                        nsrc = INT(Npart)
                        !CALL CPU_TIME(start)
                        call FMM_UpdateTree(fmm, trgCoord, srcCoord, ntrg, nsrc); 
                        call FMM_Evaluate(fmm, trgValue, srcValue, ntrg, nsrc); 
                        !CALL CPU_TIME(finish)
                        !WRITE(*,*),"Run time 1= ", finish-start
                        m = 0
                        sumx = 0.d0
                        sumy = 0.d0
                        sumz = 0.d0
                        DO i=1,nelm
!%%%%%%%%%%%%Giving back the vectors from the main list to individual list%%%%%%%%%%%%%%%%%%%%%%                                
                                ik2 = 3*i-2;
                                RHS(ik2) = 2.d0*trgValue(3*(i - 1) + 1)
                                RHS(ik2+1) = 2.d0*trgValue(3*(i - 1) + 2)
                                RHS(ik2+2) = 2.d0*trgValue(3*(i - 1) + 3)
                                sumx = sumx + trgValue(3*(i - 1) + 1)
                                sumy = sumy + trgValue(3*(i - 1) + 2)
                                sumz = sumz + trgValue(3*(i - 1) + 3)
                        ENDDO                                              
!%%%%%%%%%%%%CALCULATING THE TRACTION FORCES ON THE BOUNDARY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
                F_tr = 0.d0                            
                !CALL CPU_TIME(start)!Solve the linear problem for the traction by setting up the right hand side
                F_tr = MATMUL(RM1,RHS) !Inverse is already saved in a file and read at the start of the simulations and stored
                m = 0                  !The same inverted matrix is used again and again at every step to calculate the
                DO i = 1,nelm
                        ik2 = 3*i-2;
                        fx(i) = F_tr(ik2)
                        fy(i) = F_tr(ik2+1)
                        fz(i) = F_tr(ik2+2)
                        DO j=1,mint
                        m = m + 1
                        srcValue1(3*(m - 1) + 1) = fc1(i)*fx(i)
                        srcValue1(3*(m - 1) + 2) = fc1(i)*fy(i)
                        srcValue1(3*(m - 1) + 3) = fc1(i)*fz(i)
                        ENDDO
                ENDDO 
                !CALL CPU_TIME(finish)
                !WRITE(*,*),"Run time = 2",finish-start
                m = 0
                ux_slp = 0.d0
                uy_slp = 0.d0
                uz_slp = 0.d0 
!%%%%%%%%%%%%CALCULATING THE VELOCITIES INSIDE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                        !CALL CPU_TIME(start)
                        ntrg = INT(Npart)
                        nsrc = INT(Nelm*mint)
                        call FMM_UpdateTree(fmm, trgCoord1, srcCoord1, ntrg, nsrc); !CALL Fast Multipole
                        call FMM_Evaluate(fmm, trgValue1, srcValue1, ntrg, nsrc);   !CALL Fast Multipole
                        DO l = 1,Npart                               
                                        ux_slp(l) = trgValue1(3*(l - 1) + 1)
                                        uy_slp(l) = trgValue1(3*(l - 1) + 2)
                                        uz_slp(l) = trgValue1(3*(l - 1) + 3)
                        ENDDO
                        !CALL CPU_TIME(finish)
		        !WRITE(*,*),"Run time = 3",finish-start
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF FMM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                               
                        !CALL CPU_TIME(start)
                        k1 = Nelm
                        DO k = 1,Npart
                        k1 = k1 + 1
                        u_1(3*(k - 1) + 1) = (trgValue(3*(k1 - 1) + 1) - ux_slp(k))*fac !+ F_0(3*(k - 1) + 1)
                        u_1(3*(k - 1) + 2) = (trgValue(3*(k1 - 1) + 2) - uy_slp(k))*fac !+ F_0(3*(k - 1) + 2)
                        u_1(3*(k - 1) + 3) = (trgValue(3*(k1 - 1) + 3) - uz_slp(k))*fac !+ F_0(3*(k - 1) + 3)
                        !u_1(3*(k - 1) + 1) = (- ux_slp(k) + F_0(3*(k - 1) + 1)*cf)*fac
                        !u_1(3*(k - 1) + 2) = (- uy_slp(k) + F_0(3*(k - 1) + 2)*cf)*fac
                        !u_1(3*(k - 1) + 3) = (- uz_slp(k) + F_0(3*(k - 1) + 3)*cf)*fac
                        ENDDO 
                        !u_1(1) = F_0(1)
                        !u_1(2) = F_0(2)
                        !u_1(3) = F_0(3)                     
                        !CALL CPU_TIME(finish)
		        !WRITE(*,*),"Run time = 4", finish-start 
                        CLOSE(32)
                        CLOSE(21)
                        CALL FMM_DataClear(fmm);                      
                        DEALLOCATE(srcCoord)
                        DEALLOCATE(srcValue)
                        DEALLOCATE(trgCoord)
                        DEALLOCATE(trgValue)
                        DEALLOCATE(p1)
                        DEALLOCATE(u)
                        DEALLOCATE(vn)
                        DEALLOCATE(x_0)
                        DEALLOCATE(y_0)
                        DEALLOCATE(z_0)
                        DEALLOCATE(fx)
                        DEALLOCATE(fy)
                        DEALLOCATE(fz)
                        DEALLOCATE(fc1)
                        DEALLOCATE(RHS)
                        DEALLOCATE(F_tr)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END OF THE CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END SUBROUTINE BoundaryInt                           
