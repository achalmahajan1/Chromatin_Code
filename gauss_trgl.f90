SUBROUTINE trgl_quad (n,xi,eta,w)
                IMPLICIT NONE
                
        !-----------------------------------------
        ! Copyright by C. Pozrikidis, 1999
        ! All rights reserved.
        !
        ! This program is to be used only under the
        ! stipulations of the licensing agreement.
        !----------------------------------------

        !------------------------------------
        ! Abscissas and Weight Factors for 
        ! Gaussian Integration over a triangle
        !------------------------------------
        REAL*8 :: alpha, beta, gamma1, delta, omega1, omega2, al, be, ga, de, rh, qa,&
                  ru, o1, o2, o3      
        INTEGER*8 :: n
        REAL*8 :: xi(n),eta(n),w(n)
!--------------------
      If(n.eq.1) then
!--------------------

        xi(1) = 1.0/3.0
        eta(1) = 1.0/3.0
        w(1) = 1.0        
!--------------------------
      Else If(n.eq.12) then
!--------------------------

      al = 0.873821971016996
      be = 0.249286745170910
      ga = 0.501426509658179
      de = 0.063089014491502
      rh = 0.636502499121399
      qa = 0.310352451033785
      ru = 0.053145049844816
      o1 = 0.050844906370207
      o2 = 0.116786275726379
      o3 = 0.082851075618374

      xi(1)  = de
      xi(2)  = al
      xi(3)  = de
      xi(4)  = be
      xi(5)  = ga
      xi(6)  = be
      xi(7)  = qa
      xi(8)  = ru
      xi(9)  = rh
      xi(10) = qa
      xi(11) = ru
      xi(12) = rh

      eta(1)  = de
      eta(2)  = de
      eta(3)  = al
      eta(4)  = be
      eta(5)  = be
      eta(6)  = ga
      eta(7)  = ru
      eta(8)  = qa
      eta(9)  = qa
      eta(10) = rh
      eta(11) = rh
      eta(12) = ru

      w(1)  = o1
      w(2)  = o1
      w(3)  = o1
      w(4)  = o2
      w(5)  = o2
      w(6)  = o2
      w(7)  = o3
      w(8)  = o3
      w(9)  = o3
      w(10) = o3
      w(11) = o3
      w(12) = o3     
!-----------
      ELSEIF(n .EQ. 7)THEN
!-----------
         alpha  = 0.797426958353087
         beta   = 0.470142064105115
         gamma1  = 0.059715871789770
         delta  = 0.101286507323456
         omega1 = 0.125939180544827
         omega2 = 0.132394152788506

         xi(1)  = delta 
         xi(2)  = alpha
         xi(3)  = delta
         xi(4)  = beta
         xi(5)  = gamma1
         xi(6)  = beta
         xi(7)  = 1.d0/3.d0
         eta(1) = delta
         eta(2) = delta
         eta(3) = alpha
         eta(4) = beta
         eta(5) = beta
         eta(6) = gamma1
         eta(7) = 1.d0/3.d0
         w(1)   = omega1
         w(2)   = omega1
         w(3)   = omega1
         w(4)   = omega2
         w(5)   = omega2
         w(6)   = omega2
         w(7)   = 0.225d0
         
         
!-----------
      ELSEIF(n .EQ. 6)THEN
!-----------

         alpha  = 0.816847572980459
         beta   = 0.445948490915965
         gamma1  = 0.108103018168070
         delta  = 0.091576213509771
         omega1 = 0.109951743655322
         omega2 = 0.223381589678011

         xi(1)  = delta
         xi(2)  = alpha
         xi(3)  = delta
         xi(4)  = beta
         xi(5)  = gamma1
         xi(6)  = beta
         eta(1) = delta
         eta(2) = delta
         eta(3) = alpha
         eta(4) = beta
         eta(5) = beta
         eta(6) = gamma1
         w(1)   = omega1
         w(2)   = omega1
         w(3)   = omega1
         w(4)   = omega2
         w(5)   = omega2
         w(6)   = omega2
        ENDIF         
END SUBROUTINE trgl_quad
