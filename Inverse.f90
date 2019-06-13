SUBROUTINE inverse(a,c,n)
        IMPLICIT NONE 
                INTEGER*4 n, i, j ,k 
                REAL*8 :: a(n,n), c(n,n), L(n,n), U(n,n), b(n), d(n), x(n), coeff
                
! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
                L=0.d0
                U=0.d0
                b=0.d0

! step 1: forward elimination
                DO k=1, n-1
                        DO i=k+1,n
                                coeff=a(i,k)/a(k,k)
                                L(i,k) = coeff
                                        DO j=k+1,n
                                                a(i,j) = a(i,j)-coeff*a(k,j)
                                        END DO
                        END DO
                END DO

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
                DO i=1,n
                        L(i,i) = 1.d0
                ENDDO
! U matrix is the upper triangular part of A
                DO j=1,n
                        DO i=1,j
                                U(i,j) = a(i,j)
                        END DO
                END DO

! Step 3: compute columns of the inverse matrix C
                DO k=1,n
                        b(k)=1.d0
                        d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
                                DO i=2,n
                                        d(i)=b(i)
                                                DO j=1,i-1
                                                        d(i) = d(i) - L(i,j)*d(j)
                                                END DO
                                END DO
! Step 3b: Solve Ux=d using the back substitution
                x(n)=d(n)/U(n,n)
                                DO i = n-1,1,-1
                                        x(i) = d(i)
                                                DO j=n,i+1,-1
                                                        x(i)=x(i)-U(i,j)*x(j)
                                                END DO
                                        x(i) = x(i)/u(i,i)
                                END DO
! Step 3c: fill the solutions x(n) into column k of C
                                DO i=1,n
                                        c(i,k) = x(i)
                                END DO
                                b(k)=0.d0
               END DO
END SUBROUTINE inverse
