SUBROUTINE Cholesky(A,n)
        IMPLICIT NONE
                INTEGER*4 :: n,j      
                REAL*8  :: A(n,n) 

                DO j = 1,n
                        A(j,j) = DSQRT(A(j,j) - DOT_PRODUCT(A(j,1:j-1),A(j,1:j-1)))
                                IF(j < n) A(j+1:n,j) = (A(j+1:n,j) - MATMUL(A(j+1:n,1:j-1),A(j,1:j-1)))/A(j,j)
                        A(j,j+1:n)=0.d0
                ENDDO
END SUBROUTINE Cholesky


!SUBROUTINE Cholesky(A,n)
!        IMPLICIT NONE
!                INTEGER*4 :: n,j      
!                REAL*8  :: A(n,n) 

!                DO j = 1,n
!                        A(j+1:n,j+1:n) = (A(j+1:n,j+1:n) - MATMUL(A(j+1:n,j),A(j+1:n,j)))/A(j,j)
!                        A(j,j) = DSQRT(A(j,j))
!                        A(j+1:n,j)=A(j+1:n,j)/A(j,j)
!                        A(j,j+1:n)=0.d0
!                ENDDO
!END SUBROUTINE Cholesky


!SUBROUTINE Cholesky(A,n)
!        DO i=1:n
!                A(i+1:n,i+1:n)=A(i+1:n,i+1:n)-A(i+1:n,i)*A(i+1:n,i)/A(i,i)
!                A(i,i)=DSQRT(A(i,i))
!                A(i+1:n,i)=A(i+1:n,i)/A(i,i)
!                A(i,i+1:n)=0.d0
!        ENDDO
!END SUBROUTINE Cholesky
