SUBROUTINE SpringForce(x,y,z,xc,yc,zc,Fqx,Fqy,Fqz,N,Nc)
        IMPLICIT NONE
                INTEGER*8 :: N, i, j, Nc
                REAL*8 :: H, xc, yc, zc, Q, Qo
                REAL*8 :: Fqx, Fqy, Fqz, x, y, z, dx, dy, dz
                H = 15.d0
                Qo = 8.d0
                Fqx = 0.d0
                Fqy = 0.d0
                Fqz = 0.d0
                !DO i = 1,Nc
                !        DO j = 1,N
                                dx = xc - x
                                dy = yc - y
                                dz = zc - z
                                Q = DSQRT(dx**2.d0 + dy**2.d0 + dz**2.d0)
                                Fqx = H*dx/(1-Q**2.d0/Qo)
                                Fqy = H*dy/(1-Q**2.d0/Qo)
                                Fqz = H*dz/(1-Q**2.d0/Qo)
                 !       ENDDO
                 !ENDDO   
END SUBROUTINE SpringForce                     
