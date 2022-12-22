       program pytest 
         implicit none 
         double precision :: q, rr, den, eps 
         double precision :: peryev 
                                                                        
         integer :: i 
                                                                        
         rr   = 1 
         den  = 0.01 
         eps  = 1d-7 
                                                                        
        do i=1,200 
         q = i * 0.02 
         write(6,*)i,q,peryev(q,rr,den,eps) 
        enddo 
                                                                        
                                                                        
       end program pytest 
                                                                        
                                                                        
!  percus-yevick s(q)                                                   
! -------------------                                                   
! lit. m.s. wertheim, prl 10 (1963) 321                                 
! -------------------------------------                                 
!                                                                       
       function peryev(q,rr,den,eps) 
!      -----------------------------                                    
!      q     = impulsuebertrag                                          
!      rr    = teilchenradius                                           
!      den   = teilchendichte                                           
!      eps   = realteil ---> 0                                          
!          
       implicit none                                                             
       double precision, parameter :: pi=3.141592654d0
       double complex              :: ct, cst, clt, cgt 
     
       
       double precision            :: eta, dimag, deps, dr, dq, q, r, rr, den, eps 
       double precision            :: peryev 
!                                                                       
       r   = 2*rr 
! --- r ist der durchmesser ! der teilchen !                            
       dr  = r 
       dq  = q 
       deps= eps 
       eta = pi*den*r*r*r/6.d0 
!                                                                       
       ct  = cmplx ( deps , -dq*dr ) 
       cst = (((1.d0-eta)**2*ct+6.d0*eta*(1.d0-eta))*ct                 &
     &      +18.d0*eta*eta)*ct                                          &
     &      -12.d0*eta*(1.d0+2.d0*eta)                                  
       clt = 12.d0*eta*((1.d0+0.5d0*eta)*ct+(1.d0+2.d0*eta)) 
!                                                                       
       cgt = (ct*clt) / ( 12.d0*eta*(clt+cst*cdexp(ct)) ) 
!                                                                       
       peryev =  1.d0 + den * (4.d0*pi/dq) * dr*dr * imag( cgt ) 
!                                                                       
       return 
      END                                           
