                                                                        
       subroutine NzimPnaX(q,t,temp,Dr,etasolv,N,R,pmax,pwidth,nue,alpha,Sq,Sqt)                  
!      =========================================================================                  
!                                                                       
! Rouse expression for a chain of finite length:                        
! Input parameters:                                                     
!    q     ----> momentum transfer in A**-1                             
!    t     ----> time in nano-sec                                       
!    temp  ----> temperature in K                                       
!    Dr    ----> center of mass diffusion constant in A**2/ns, if 0 <-- 
!    etasolv --> solvent viscosity units?                               
!    N     ----> number of chain segments                               
!    R     ----> end-to-end distance of the polymer molecule            
!    pmax  ----> maximum flexible mode nr. (transition)                 
!    pwidth----> width of modenr. transition                            
!    nue   ----> chain expansion exponent (Gaussian=1/2)                
!    alpha ----> stiffness-parameter (only for relaxations, Norbornenes.
! Output parameters:                                                    
!    Sq    <--- S(Q)                                                    
!    Sqt   <--- S(Q,t)                                                  
! ------------------------------------------------------------          
!                                                                       
       implicit none 
                                                                        
       double precision kb, pi 
       parameter(kb=1.380662d-23) 
       parameter(pi=3.141592654d0) 
                                                                        
       double precision q,t,temp,Dr,etasolv,R,Sq,Sqt 
       double precision pmax, pwidth, nue, alpha 
       integer N, nn,mm,p 
                                                                        
       double precision l, tau_p, kbt, Sq0, arg1, arg2 
                                                           ! ,arg20     
       double precision a0,e0, ff2, ffc,    arg10 
       double precision aa1 , aa2, rsi 
       double precision rate_p 
       double precision fdf 
    
       double precision :: cosarray(N,N), ewfac(N)

                                                                    
                                                                        
 !      integer iout                                                    
                                                                        
       if(N.le.0) then 
         Sq = 999 
         Sqt= 999 
         write(6,*)'Error Number of chain segments is <= 0!',N 
         return 
       endif 
                                                                        
! ---- determine the segment length l ----                              
       l = R/(dfloat(N)**nue) 
                                                                        
                                ! in Joule = kg*m**2/s**2               
       kbt = temp*kb 
                                ! Re in meter                           
       Rsi = R*1d-10 
                                                                        
! ---- set the diffusion constant if input is zero --- !                
       if(Dr.eq.0.0d0) then 
         ! do a linear interpolation between values for Gaussian and goo
         ! limting factors give in Doi-Edwards                          
         fdf= 0.196d0 + (0.203d0-0.196d0)*(nue-0.5d0)/(0.6d0-0.5d0) 
                                         ! in m**2/s                    
         Dr = fdf*kbt/(Rsi*etasolv) 
                                         ! in A**2/ns                   
         Dr = Dr*1d20/1d9 
       endif 
                    


!$OMP PARALLEL DO     
       do nn=1,N
        do p=1,N
         cosarray(nn,p) = cos((pi*p*nn)/dfloat(N)) / sqrt((dfloat(p)**(2*nue+1.0d0)+alpha*dfloat(p)**(4))) 
        enddo
       enddo
!$OMP END PARALLEL DO   


!$OMP PARALLEL DO PRIVATE(p,tau_p,rate_p,a0,e0)   
       do p=1,N
! rouse:    tau_p = 2*W*(1-cos((pi*p)/dfloat(N)))                      !
            tau_p = 0.325*etasolv*Rsi**3/kbt /                          &
     &              (dfloat(p)**(3*nue)+alpha*dfloat(p)**(4-nue))       
                                                                       !
!                                 --- since we use R-end-to-end durectly
!                                     here the chain statistics is impli
!                   ----- this is the Gaussian chain factor (hope that i
                                                                       !
            tau_p = tau_p * 1d9 
                                                                        
            rate_p = transf_p(dfloat(p),pmax,pwidth)/tau_p 
                                                                        
            a0    = -t*rate_p 
            e0    = 1.0d0-exp(a0) 
 
            ewfac(p) = (1d0-exp(a0))
       enddo
!$OMP END PARALLEL DO    



                                                    
! ---- init sums ----                                                   
       Sq0 = 0 
       Sq  = 0 
       Sqt = 0 

       ff2  = -2*(R*q)**2/(3*pi**2) 

                
!$OMP PARALLEL DO PRIVATE(arg1,arg10) REDUCTION(+:Sq,Sqt,arg2)
                                                        
! ---- Do the sums -----                                                
       do nn = 1,N 
        do mm = 1,N 

          arg2 = sum(cosarray(mm,:)*cosarray(nn,:)* ewfac(:) )
                                                  
          Sq  = Sq  + exp(-(q**2)*(       (abs(nn-mm)**(2*nue))*(l**2)/6.0d0) ) 
          Sqt = Sqt + exp(-(q**2)*(Dr*t + (abs(nn-mm)**(2*nue))*(l**2)/6.0d0) +arg2*ff2)                                                                         
        enddo 
       enddo 
!$OMP END PARALLEL DO

                                                                        
       Sq  = Sq /N 
       Sqt = Sqt/N 
                                                                        
!       if(iout().gt.0)write(6,'(1x,5E14.6)')q,t,Sq,Sqt, Sqt/Sq         
                                                                        
       return 
CONTAINS                                           
                                                                        
                                                                        
       double precision function transf_p(x,p,w) 
!      -----------------------------------------                        
       implicit none 
       double precision x, p, w, f, arg 
                                                                        
       if(w.le.0.0d0) then 
         if(x.gt.p) then 
          f = 0.0d0 
         else 
          f = 1.0d0 
         endif 
       else 
          arg = (x-p)/w 
          if(arg.lt.-200d0) arg = -200d0 
          if(arg.gt. 200d0) arg =  200d0 
          f   = 1.0d0/(1.0d0+exp(arg)) 
       endif 
                                                                        
       transf_p = f 
                                                                        
       return 
      END   

END                                        
