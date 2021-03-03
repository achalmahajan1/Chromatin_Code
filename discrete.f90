SUBROUTINE discrete(p,n,Npts,nelm)
        IMPLICIT NONE
                INTEGER*8 :: i, nelm, j, k, l, Ndiv, nm, Npts, Iflag, n(5120,6)
                REAL*8 :: rm, ru, c0, c1, c2, c3, c4, c5, c6, c7, rad, eps=0.00000001d0
                REAL*8 :: VX(12), VY(12), VZ(12), x(5120,6), y(5120,6), z(5120,6), &
                                xn(5120,6), yn(5120,6), zn(5120,6), p(10242,3)                
        !Dividing over the icosahedron
                !N = 12
                Ndiv = 4
                ru = 0.25D0*dsqrt(10.0D0+2.0D0*dsqrt(5.0D0))
                rm = 0.25D0*(1.0D0+dsqrt(5.0D0))   ! twice the golden ratio
                c0 = 0.D0
               !c1 = 0.9512 D0
                c1 = ru
               !c2 = 0.8507 D0
                c2 = 2.0*ru/dsqrt(5.0D0)
               !c3 = 0.8090 D0
                c3 = rm
               !c4 = 0.4253 D0
                c4 = ru/dsqrt(5.0D0)
               !c5 = 0.2629 D0
                c5 = dsqrt( ru**2-c3**2-c4**2)
                c6 = 0.5D0
               !c7 = 0.6882 D0
                c7 = dsqrt( ru**2-c4**2-c6**2)
               ! PRINT*,ru, rm, c0, c1, c2, c3, c4, c5, c6, c7
                VX(1) =  c0
                VY(1) =  c0
                VZ(1) =  c1

                VX(2) =  c0
                VY(2) =  c2
                VZ(2) =  c4

                VX(3) =  c3
                VY(3) =  c5
                VZ(3) =  c4

                VX(4) =  c6
                VY(4) = -c7
                VZ(4) =  c4

                VX(5) = -c6
                VY(5) = -c7
                VZ(5) =  c4

                VX(6) = -c3
                VY(6) =  c5
                VZ(6) =  c4

                VX(7) = -c6
                VY(7) =  c7
                VZ(7) = -c4

                VX(8) =  c6
                VY(8) =  c7
                VZ(8) = -c4

                VX(9) =  c3
                VY(9) = -c5
                VZ(9) = -c4

                VX(10) =  c0
                VY(10) = -c2
                VZ(10) = -c4

      VX(11) = -c3
      VY(11) = -c5
      VZ(11) = -c4

      VX(12) =  c0
      VY(12) =  c0
      VZ(12) = -c1
      !OPEN(23,FILE = "pos.dat")
      !DO i = 1,N
      !          WRITE(23,*), VX(i), VY(i), VZ(i)
      !ENDDO
      x(1,1) = VX(1)   ! first element
      y(1,1) = VY(1)
      z(1,1) = VZ(1)
      x(1,2) = VX(3)
      y(1,2) = VY(3)
      z(1,2) = VZ(3)
      x(1,3) = VX(2)
      y(1,3) = VY(2)
      z(1,3) = VZ(2)
        !---
      x(2,1) = VX(1)
      y(2,1) = VY(1)
      z(2,1) = VZ(1)
      x(2,2) = VX(4)
      y(2,2) = VY(4)
      z(2,2) = VZ(4)
      x(2,3) = VX(3)
      y(2,3) = VY(3)
      z(2,3) = VZ(3)
                                    !---
      x(3,1) = VX(1)
      y(3,1) = VY(1)
      z(3,1) = VZ(1)
      x(3,2) = VX(5)
      y(3,2) = VY(5)
      z(3,2) = VZ(5)
      x(3,3) = VX(4)
      y(3,3) = VY(4)
      z(3,3) = VZ(4)
                                    !---
      x(4,1) = VX(1)
      y(4,1) = VY(1)
      z(4,1) = VZ(1)
      x(4,2) = VX(6)
      y(4,2) = VY(6)
      z(4,2) = VZ(6)
      x(4,3) = VX(5)
      y(4,3) = VY(5)
      z(4,3) = VZ(5)
                                    !---
      x(5,1) = VX(1)
      y(5,1) = VY(1)
      z(5,1) = VZ(1)
      x(5,2) = VX(2)
      y(5,2) = VY(2)
      z(5,2) = VZ(2)
      x(5,3) = VX(6)
      y(5,3) = VY(6)
      z(5,3) = VZ(6)
                                    !---
      x(6,1) = VX(2)
      y(6,1) = VY(2)
      z(6,1) = VZ(2)
      x(6,2) = VX(3)
      y(6,2) = VY(3)
      z(6,2) = VZ(3)
      x(6,3) = VX(8)
      y(6,3) = VY(8)
      z(6,3) = VZ(8)
                                    !---
      x(7,1) = VX(3)
      y(7,1) = VY(3)
      z(7,1) = VZ(3)
      x(7,2) = VX(4)
      y(7,2) = VY(4)
      z(7,2) = VZ(4)
      x(7,3) = VX(9)
      y(7,3) = VY(9)
      z(7,3) = VZ(9)
                                    !---
      x(8,1) = VX(4)
      y(8,1) = VY(4)
      z(8,1) = VZ(4)
      x(8,2) = VX(5)
      y(8,2) = VY(5)
      z(8,2) = VZ(5)
      x(8,3) = VX(10)
      y(8,3) = VY(10)
      z(8,3) = VZ(10)
                                    !---
      x(9,1) = VX(5)
      y(9,1) = VY(5)
      z(9,1) = VZ(5)
      x(9,2) = VX(6)
      y(9,2) = VY(6)
      z(9,2) = VZ(6)
      x(9,3) = VX(11)
      y(9,3) = VY(11)
      z(9,3) = VZ(11)
                                    !---
      x(10,1) = VX(6)
      y(10,1) = VY(6)
      z(10,1) = VZ(6)
      x(10,2) = VX(2)
      y(10,2) = VY(2)
      z(10,2) = VZ(2)
      x(10,3) = VX(7)
      y(10,3) = VY(7)
      z(10,3) = VZ(7)
                                    !---
      x(11,1) = VX(2)
      y(11,1) = VY(2)
      z(11,1) = VZ(2)
      x(11,2) = VX(8)
      y(11,2) = VY(8)
      z(11,2) = VZ(8)
      x(11,3) = VX(7)
      y(11,3) = VY(7)
      z(11,3) = VZ(7)
                                    !---
      x(12,1) = VX(3)
      y(12,1) = VY(3)
      z(12,1) = VZ(3)
      x(12,2) = VX(9)
      y(12,2) = VY(9)
      z(12,2) = VZ(9)
      x(12,3) = VX(8)
      y(12,3) = VY(8)
      z(12,3) = VZ(8)
                                    !---
      x(13,1) = VX(4)
      y(13,1) = VY(4)
      z(13,1) = VZ(4)
      x(13,2) = VX(10)
      y(13,2) = VY(10)
      z(13,2) = VZ(10)
      x(13,3) = VX(9)
      y(13,3) = VY(9)
      z(13,3) = VZ(9)
                                    !---
      x(14,1) = VX(5)
      y(14,1) = VY(5)
      z(14,1) = VZ(5)
      x(14,2) = VX(11)
      y(14,2) = VY(11)
      z(14,2) = VZ(11)
      x(14,3) = VX(10)
      y(14,3) = VY(10)
      z(14,3) = VZ(10)
                                    !---
      x(15,1) = VX(6)
      y(15,1) = VY(6)
      z(15,1) = VZ(6)
      x(15,2) = VX(7)
      y(15,2) = VY(7)
      z(15,2) = VZ(7)
      x(15,3) = VX(11)
      y(15,3) = VY(11)
      z(15,3) = VZ(11)
                                    !---
      x(16,1) = VX(7)
      y(16,1) = VY(7)
      z(16,1) = VZ(7)
      x(16,2) = VX(8)
      y(16,2) = VY(8)
      z(16,2) = VZ(8)
      x(16,3) = VX(12)
      y(16,3) = VY(12)
      z(16,3) = VZ(12)
                                    !---
      x(17,1) = VX(8)
      y(17,1) = VY(8)
      z(17,1) = VZ(8)
      x(17,2) = VX(9)
      y(17,2) = VY(9)
      z(17,2) = VZ(9)
      x(17,3) = VX(12)
      y(17,3) = VY(12)
      z(17,3) = VZ(12)
                                    !--- 
      x(18,1) = VX(9)
      y(18,1) = VY(9)
      z(18,1) = VZ(9)
      x(18,2) = VX(10)
      y(18,2) = VY(10)
      z(18,2) = VZ(10)
      x(18,3) = VX(12)
      y(18,3) = VY(12)
      z(18,3) = VZ(12)
                                    !---
      x(19,1) = VX(10)
      y(19,1) = VY(10)
      z(19,1) = VZ(10)
      x(19,2) = VX(11)
      y(19,2) = VY(11)
      z(19,2) = VZ(11)
      x(19,3) = VX(12)
      y(19,3) = VY(12)
      z(19,3) = VZ(12)
                                    !---
      x(20,1) = VX(11)
      y(20,1) = VY(11)
      z(20,1) = VZ(11)
      x(20,2) = VX(7)
      y(20,2) = VY(7)
      z(20,2) = VZ(7)
      x(20,3) = VX(12)
      y(20,3) = VY(12)
      z(20,3) = VZ(12)
        !------------------------------------------
        ! compute the mid-points of the three sides
        ! of the 20 first-generation elements
        !
        ! midpoints are numbered 4, 5, 6
        !------------------------------------------

      Do i=1,Nelm

       x(i,4) = 0.5D0*(x(i,1)+x(i,2))
       y(i,4) = 0.5D0*(y(i,1)+y(i,2))
       z(i,4) = 0.5D0*(z(i,1)+z(i,2))

       x(i,5) = 0.5D0*(x(i,2)+x(i,3))
       y(i,5) = 0.5D0*(y(i,2)+y(i,3))
       z(i,5) = 0.5D0*(z(i,2)+z(i,3))

       x(i,6) = 0.5D0*(x(i,3)+x(i,1))
       y(i,6) = 0.5D0*(y(i,3)+y(i,1))
       z(i,6) = 0.5D0*(z(i,3)+z(i,1))

      END DO         
      !OPEN(23,FILE = "pos.dat")
      !DO i = 1,Nelm
      !          DO j= 1,3
      !                  WRITE(23,*), X(i,j), Y(i,j), Z(i,j)
      !          ENDDO
      !ENDDO
      !OPEN(24,FILE = "pos1.dat")
      !DO i = 1,Nelm
      !          DO j= 4,6
      !                  WRITE(24,*), X(i,j), Y(i,j), Z(i,j)
      !          ENDDO
      !ENDDO  
      
                !---
                ! project the nodes onto the unit sphere
                !---

        DO k=1,Nelm
                DO l=1,6
                        rad = DSQRT(x(k,l)**2+y(k,l)**2+z(k,l)**2) 
                        x(k,l) = x(k,l)/rad
                        y(k,l) = y(k,l)/rad
                        z(k,l) = z(k,l)/rad
                ENDDO
        ENDDO
        !OPEN(23,FILE = "pos.dat")
        !DO i = 1,Nelm
        !        DO j= 1,3
        !                WRITE(23,*), X(i,j), Y(i,j), Z(i,j)
        !        ENDDO
        !ENDDO
        !OPEN(24,FILE = "pos1.dat")
        !DO i = 1,Nelm
        !        DO j= 4,6
        !                WRITE(24,*), X(i,j), Y(i,j), Z(i,j)
        !        ENDDO
        !ENDDO
        !-------------------------------------------
        ! compute the local element node coordinates
        ! for discretization levels 1 through Ndiv
        !-------------------------------------------

      DO i=1,Ndiv

       nm = 0      ! count the new elements arising by sub-division
                   ! four element will be generated during each pass
       Do j=1,Nelm ! over old elements

        !---
        ! assign corner points to sub-elements
        ! these will become the "new" elements
        !---

        nm = nm+1

        xn(nm,1) = x(j,1)                  !  first sub-element
        yn(nm,1) = y(j,1)
        zn(nm,1) = z(j,1)

        xn(nm,2) = x(j,4)
        yn(nm,2) = y(j,4) 
        zn(nm,2) = z(j,4)

        xn(nm,3) = x(j,6)
        yn(nm,3) = y(j,6)
        zn(nm,3) = z(j,6)

        xn(nm,4) = 0.5D0*(xn(nm,1)+xn(nm,2))
        yn(nm,4) = 0.5D0*(yn(nm,1)+yn(nm,2))
        zn(nm,4) = 0.5D0*(zn(nm,1)+zn(nm,2))

        xn(nm,5) = 0.5D0*(xn(nm,2)+xn(nm,3))
        yn(nm,5) = 0.5D0*(yn(nm,2)+yn(nm,3))
        zn(nm,5) = 0.5D0*(zn(nm,2)+zn(nm,3))

        xn(nm,6) = 0.5D0*(xn(nm,3)+xn(nm,1))
        yn(nm,6) = 0.5D0*(yn(nm,3)+yn(nm,1))
        zn(nm,6) = 0.5D0*(zn(nm,3)+zn(nm,1))

        nm = nm+1

        xn(nm,1) = x(j,4)                !  second sub-element
        yn(nm,1) = y(j,4)
        zn(nm,1) = z(j,4)

        xn(nm,2) = x(j,2)
        yn(nm,2) = y(j,2)
        zn(nm,2) = z(j,2)

        xn(nm,3) = x(j,5)
        yn(nm,3) = y(j,5)
        zn(nm,3) = z(j,5)

        xn(nm,4) = 0.5D0*(xn(nm,1)+xn(nm,2))
        yn(nm,4) = 0.5D0*(yn(nm,1)+yn(nm,2))
        zn(nm,4) = 0.5D0*(zn(nm,1)+zn(nm,2))

        xn(nm,5) = 0.5D0*(xn(nm,2)+xn(nm,3))
        yn(nm,5) = 0.5D0*(yn(nm,2)+yn(nm,3))
        zn(nm,5) = 0.5D0*(zn(nm,2)+zn(nm,3))

        xn(nm,6) = 0.5D0*(xn(nm,3)+xn(nm,1))
        yn(nm,6) = 0.5D0*(yn(nm,3)+yn(nm,1))
        zn(nm,6) = 0.5D0*(zn(nm,3)+zn(nm,1))

        nm = nm+1

        xn(nm,1) = x(j,6)                !  third sub-element
        yn(nm,1) = y(j,6)
        zn(nm,1) = z(j,6)
 
        xn(nm,2) = x(j,5)
        yn(nm,2) = y(j,5)
        zn(nm,2) = z(j,5)

        xn(nm,3) = x(j,3)
        yn(nm,3) = y(j,3)
        zn(nm,3) = z(j,3)

        xn(nm,4) = 0.5D0*(xn(nm,1)+xn(nm,2))
        yn(nm,4) = 0.5D0*(yn(nm,1)+yn(nm,2))
        zn(nm,4) = 0.5D0*(zn(nm,1)+zn(nm,2))

        xn(nm,5) = 0.5D0*(xn(nm,2)+xn(nm,3))
        yn(nm,5) = 0.5D0*(yn(nm,2)+yn(nm,3))
        zn(nm,5) = 0.5D0*(zn(nm,2)+zn(nm,3))

        xn(nm,6) = 0.5D0*(xn(nm,3)+xn(nm,1))
        yn(nm,6) = 0.5D0*(yn(nm,3)+yn(nm,1))
        zn(nm,6) = 0.5D0*(zn(nm,3)+zn(nm,1))

        nm = nm+1

        xn(nm,1) = x(j,4)                !  fourth sub-element
        yn(nm,1) = y(j,4)
        zn(nm,1) = z(j,4)

        xn(nm,2) = x(j,5)
        yn(nm,2) = y(j,5)
        zn(nm,2) = z(j,5)

        xn(nm,3) = x(j,6)
        yn(nm,3) = y(j,6)
        zn(nm,3) = z(j,6)

        xn(nm,4) = 0.5D0*(xn(nm,1)+xn(nm,2))    ! mid points
        yn(nm,4) = 0.5D0*(yn(nm,1)+yn(nm,2))
        zn(nm,4) = 0.5D0*(zn(nm,1)+zn(nm,2))

        xn(nm,5) = 0.5D0*(xn(nm,2)+xn(nm,3))
        yn(nm,5) = 0.5D0*(yn(nm,2)+yn(nm,3))
        zn(nm,5) = 0.5D0*(zn(nm,2)+zn(nm,3))

        xn(nm,6) = 0.5D0*(xn(nm,3)+xn(nm,1))
        yn(nm,6) = 0.5D0*(yn(nm,3)+yn(nm,1))
        zn(nm,6) = 0.5D0*(zn(nm,3)+zn(nm,1))

       End Do                      !  end of old-element loop

!--------------------------------------
! number of elements has been increased
! by a factor of four
!--------------------------------------

       Nelm = 4*Nelm

!---
! relabel the new points
! and place them in the master list
!---
        DO k=1,Nelm
                DO l=1,6

                x(k,l) = xn(k,l)
                y(k,l) = yn(k,l)
                z(k,l) = zn(k,l)

        !--- project onto the unit sphere

                rad = dsqrt(x(k,l)**2+y(k,l)**2+z(k,l)**2)
                x(k,l) = x(k,l)/rad
                y(k,l) = y(k,l)/rad
                z(k,l) = z(k,l)/rad

                xn(k,l) = 0.0D0   ! zero just in case
                yn(k,l) = 0.0D0
                zn(k,l) = 0.0D0
                ENDDO
        ENDDO
       ! !----------

        ENDDO            !  end of discretization- level loop
        !-----------------------------------------
        !PRINT*,"inside",Nelm, x(5120,1),x(32,1)

!-----------------------------------
! Generate a list of global nodes by looping 
! over all elementsand adding nodes not found
! in the current list.
!
! Fill in the connectivity table n(i,j) 
! containing node numbers of element points 1-6
!-----------------------------------

!---
! six nodes of the first element are
! entered mannualy
!---
        !OPEN(26,FILE = "pos.dat")
        !j = 0
        !DO k = 1,Nelm
        !        DO l= 1,3
        !                j = j + 1
        !                WRITE(26,*), X(k,l), Y(k,l), Z(k,l)
        !                !PRINT*,j
        !        ENDDO
       !ENDDO
       !PRINT*,"inside",Nelm,j
       !OPEN(27,FILE = "pos1.dat")
       !DO k = 1,Nelm
       !         DO l= 4,6
       !                WRITE(27,*), X(k,l), Y(k,l), Z(k,l)
       !         ENDDO
       !ENDDO
       
      p(1,1) = x(1,1)
      p(1,2) = y(1,1)
      p(1,3) = z(1,1)

      p(2,1) = x(1,2)
      p(2,2) = y(1,2)
      p(2,3) = z(1,2)

      p(3,1) = x(1,3)
      p(3,2) = y(1,3)
      p(3,3) = z(1,3)

      p(4,1) = x(1,4)
      p(4,2) = y(1,4)
      p(4,3) = z(1,4)

      p(5,1) = x(1,5)
      p(5,2) = y(1,5)
      p(5,3) = z(1,5)

      p(6,1) = x(1,6)
      p(6,2) = y(1,6)
      p(6,3) = z(1,6)

      n(1,1) = 1  ! first  node of first element is global node 1
      n(1,2) = 2  ! second node of first element is global node 2
      n(1,3) = 3  ! third  node of first element is global node 3
      n(1,4) = 4
      n(1,5) = 5
      n(1,6) = 6  ! sixth  node of first element is global node 6

      Npts = 6

!---
! loop over further elements
!
! Iflag=0 will signal a new global node
!---

      Do i=2,Nelm        ! loop over elements
       Do j=1,6          ! loop over element nodes

        Iflag=0

         Do k=1,Npts
          If(abs(x(i,j)-p(k,1)).le.eps) then
           If(abs(y(i,j)-p(k,2)).le.eps) then
            If(abs(z(i,j)-p(k,3)).le.eps) then

             Iflag = 1    ! the node has been recorded previously
             n(i,j) = k   ! the jth local node of element i 
                          ! is the kth global node 
             !PRINT*,"1","i",i,"j",j,"npts",npts,'nij',n(i,j),"k",k             
            !PRINT*,"YES", npts, j, k , Iflag, n(i,j), p(i,j), x(i,j)             
            End If
           End If
          End If
         End Do
        
         If(Iflag.eq.0) then     ! record the node

          Npts = Npts+1          ! one more global node
          p(Npts,1) = x(i,j)
          p(Npts,2) = y(i,j)
          p(Npts,3) = z(i,j)
          rad = dsqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2)
          n(i,j) = Npts   ! the jth local node of element i
           !PRINT*,"2","i",i,"j",j,"npts",npts,'nij',n(i,j)
          !PAUSE                ! is the new global node 
         !PRINT*,"YES1", npts, j, k , Iflag, n(i,j), Npts, p(Npts,1), x(i,j)        
         End If
               !PAUSE
       End Do
      End Do                      !  end of loop over elements
      
            !PRINT*,"inside npoints", npts, nelm  
      Do i=1,Npts

       rad = dsqrt(p(i,1)**2+p(i,2)**2+p(i,3)**2)

       p(i,1) = p(i,1)/rad
       p(i,2) = p(i,2)/rad
       p(i,3) = p(i,3)/rad

      End Do
END SUBROUTINE discrete                  
