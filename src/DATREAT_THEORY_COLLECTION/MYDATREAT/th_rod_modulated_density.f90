      function th_rodmod(x,pa,thnam,parnam,npar,iadda,ini) 
!     ===================================================               
!                                                                       
! -------> density modulated rod rho(z) ~ (cos(2Pi/rep*z)**2 + c) randomly orien in space 
!                                                                       
!                                                                       
       implicit none 
                                                                        
       character*8 thnam,parnam(20) 
       real*4 pa(20) 
       integer ini, npar, nparx 
                                                                        
       real*8 adapint, a, b 
       real*8 epsilon, erraccu, eps 
       integer maxit 
       parameter (maxit=3000) 
                                                                        
       real*8 pi 
       parameter (pi=3.141592654d0) 
       real*8 Navog 
       parameter (Navog=6.022045d23) 
                                                                        
                                                                        
                                                                        
       real*8 amplitu, d_cos 
       real*4 rhosolv4, x, th_rodmod 
                                                                        
! Entry of the additional parameters !                                  
       integer Nlay 
       double precision q,r,rep,const 
       common/crirke/ q,r,rep,const,Nlay 
                                                                        
                                                                        
       double precision y 
       integer i 
                                                                        
       double precision rkernel 
       external rkernel 
                                                                        
                                                                        
! ----- initialisation -----                                            
       if(ini.eq.0) then 
         thnam = 'modrod  ' 
         nparx = 6 
         if(npar.lt.nparx) then 
           write(6,1)thnam,nparx,npar 
    1      format(' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
           th_rodmod  = 0 
           return 
         endif 
         npar = nparx 
!        --------------> set the number of parameters                   
                                  ! eventually concentration factor     
         parnam(1) = 'amplitu' 
                                  ! radius of rod                       
         parnam(2) = 'radius ' 
                                  ! modulation period along z                  
         parnam(3) = 'd_cos  ' 
                                  ! constant c                          
         parnam(4) = 'const  ' 
                                  ! number of oszillations              
         parnam(5) = 'nrep   ' 
                                  ! epsilon for integration             
         parnam(6) = 'epsilon' 
                                                                        
!                                                                       
         th_rodmod = 0 
         return 
       endif 
!                                                                       
! ---- calculate theory here length unit is cm -----                    
!                                                                       
       q = x 
                                                                        
       amplitu = pa(1) 
       r       = abs(pa(2)) 
       d_cos   = abs(pa(3)) 
       const   = pa(4) 
       Nlay    = NINT(pa(5)/4) 
       epsilon = pa(6) 
       if(epsilon.le.0.0d0) epsilon = 1d-9 
                                                                        
       rep     = d_cos*2 
!                                                                       
! ACHTUNG: wegen der Definition der Modulation als cos**2 ist die Period
!          ebenso ist -Nlay...Nlay equivalent zu 4*Nlay oszillationen   
                                                                        
                                                                        
! hole Streulaengendichte des Loesungsmittels aus der Parameterliste    
!       call parget('bsolv   ',rhosolv4,iadda,ier)                      
!       rhosolv = rhosolv4                                              
                                                                        
!!!!!!!!!!!!! the computation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
       erraccu = 0 
                                                                        
       y = 0 
       a = 0 
       b = Pi/32 
       y = y+adapint(rkernel,a,b,epsilon,maxit,erraccu) 
       a = b 
       b = Pi/2-Pi/32 
       y = y+adapint(rkernel,a,b,epsilon,maxit,erraccu) 
       a = b 
       b = Pi/2 
       y = y+adapint(rkernel,a,b,epsilon,maxit,erraccu) 
                                                                        
       y = y/q 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
                                                                        
       th_rodmod = amplitu * y 
                                                                        
                                                                        
       return 
      END                                           
                                                                        
                                                                        
                                                                        
                                                                        
       double precision function modcos(q,c,r,N) 
!      -----------------------------------------                        
!      (cosine**2 modulation along the rod)                             
!                                                                       
       implicit none 
                                                                        
       double precision Pi 
       parameter(Pi=3.141592654d0) 
                                                                        
                                  !   q                                 
       double precision q 
                                  !   added constant                    
       double precision c 
                                  !   repetition distance, periode      
       double precision r 
                                  !   length of half of rod in terms of 
       integer          N 
                                                                        
!! The fouriertramsform squared of:                                     
!!             y = cos(2*Pi/r * z)**2 + const                           
!! modcos =  | int(y(z)*exp(i*q*z),-N*rep..+N*rep) |**2                 
!!                                                                      
        double precision qh,L,Norm 
        L = 2*N*r 
        Norm = 4*N**2*r**2*c**2+4*N**2*r**2*c+N**2*r**2 
                                                                        
                                                                        
        qh = q**2*r**2-16*Pi**2 
        if(abs(qh).lt.1d-20) then 
          modcos = (N*r/2)**2 / Norm 
          return 
        endif 
                                                                        
!        modcos = -4*(-1+cos(q*N*r)**2)*                                
!     *          (q**2*r**2+c*q**2*r**2-c*16*Pi**2-8*Pi**2)**2/         
!     *          ((qh**2) * q**2)                                       
                                                                        
        modcos = -4*(-1+cos(q*N*r)**2)*                                 &
     &          (q**2*r**2+c*qh-8*Pi**2)**2/                            &
     &          ((qh**2) * q**2)                                        
                                                                        
        modcos = modcos / Norm 
                                                                        
       return 
      END                                           
                                                                        
                                                                        
       double precision function rodq(q,r) 
!      -----------------------------------                              
!      Rod |F(q_perp)|**2                                               
!                                                                       
       implicit none 
       double precision Pi 
       parameter(Pi=3.141592654d0) 
                                                                        
                                  !   q component perpendicular to rod  
       double precision q 
                                  !   radius of the rod                 
       double precision r 
                                                                        
                                                                        
!      rodq = Pi*r*r* exp(-(q*r/2)**2)                                  
                                        ! normalized                    
       rodq =         exp(-(q*r/2)**2) 
                                                                        
       return 
      END                                           
                                                                        
                                                                        
       double precision function rkernel(phi) 
!      --------------------------------------                           
!      Integrand function                                               
!                                                                       
       implicit none 
                                  ! integration variable                
       double precision phi 
       double precision qlong 
       double precision qperp 
       double precision rodq, modcos 
                                                                        
! Entry of the additional parameters !                                  
       integer Nlay 
       double precision q,r,rep,const 
       common/crirke/ q,r,rep,const,Nlay 
                                                                        
                                                                        
       qlong = q*cos(phi) 
       qperp = q*sin(phi) 
                                                                        
       rkernel = rodq(qperp,r)*modcos(qlong,const,rep,Nlay)*qperp 
                                                                        
       return 
      END                                           
