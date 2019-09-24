      function th_nzimmaz (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!      use theory_description 
      implicit none 
      real    :: th_ndebye
      character(len=8) :: thnam, parnam (*) 
      real    :: pa (*) 
      real    :: x , xh
      integer :: mbuf, nparx, ier, ini, npar, iadda
      integer, intent(inout) :: nopar       
      character(len=80), intent(inout) :: napar(mbuf) 
      real, intent(inout) :: params(mbuf) 
     
      double precision, parameter :: Pi = 4*atan(1d0)
      integer                     :: actual_record_address

      double precision  ::     temp, qz, tau, eta, yz, SQ_rouse, a,b,etasolv 
      double precision  ::     R, Sq0, Sq, Sqt 
      double precision  ::     a0, sumq, sumnorm, q_width, dqw, qzz, fn 
      double precision  ::     diff, dr 
      real              ::     qget, tget 
      real              ::     th_nzimmaz 

      integer           ::     i, nqw 
      integer           ::     N 
      double precision  ::     pmax, pwidth, nue, alpha 
      double precision  ::     tau_zif
!                                                                       
! ----- initialisation -----                                            
       if(ini.eq.0) then 
         thnam = 'nzimmaz ' 
         nparx = 12 
         if(npar.lt.nparx) then 
           write(6,1)thnam,nparx,npar 
    1      format(' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
           th_nzimmaz = 0 
           return 
         endif 
         npar = nparx 
!        --------------> set the number of parameters                   
                                   ! Intensität tau=0                  
         parnam(1) = 'intensit' 
                                   ! Viscosität                        
         parnam(2) = 'etasolv' 
                                   !number Chain elements               
         parnam(3) = 'n_segmnt' 
                                   ! : re**2 = N * b**2 in Angstroem    
         parnam(4) = 're      ' 
                                   ! temp                               
         parnam(5) = 'temp    ' 
                                   ! Diffusion coefficient              
         parnam(6) = 'com_diff' 
                                   ! q width of detector (width for an a
         parnam(7) = 'q_width ' 
                                   ! : maximum mode that may move > 0   
         parnam(8) = 'p_max   ' 
                                   ! : width of transition regime f=1/(1
         parnam(9) = 'p_width ' 
                                   ! : chain expansion parameter (Gaussi
         parnam(10)= 'nue     ' 
                                   ! : chain stiffness descriptor alpha*
         parnam(11)= 'alpha   ' 
                                   ! : zif tau in ns                   
         parnam(12)= 'tau_zif ' 
                                                                        
!                                                                       
         th_nzimmaz = 0 
         return 
       endif 

       iadda = actual_record_address()
!                                                                       
! ---- calculate theory here -----                                      
       tau      = x 
       a0       = pa(1) 
       etasolv  = abs(pa(2)) 
       N        = nint(pa(3)) 
       R        = abs(pa(4)) 
       temp     = pa(5) 
                              ! in cm**2/sec                            
       diff     = abs(pa(6)) 
       q_width  = pa(7) 
       pmax     = abs(pa(8)) 
       pwidth   = abs(pa(9)) 
       nue      = abs(pa(10)) 
       alpha    = abs(pa(11)) 
       tau_zif  = abs(pa(12)) 
                                                                        
                                       ! Default Gaussian               
       if(nue.eq.0.0d0) nue = 0.5D0 
                                                                        
                                       ! in A**2/ns                     
       diff     = diff * 1d-9 / 1d-16 
                                                                        
                                                                        
       qget = 0.01 
       call        parget('q       ',qget,iadda,ier) 
       qz   = qget 
       if(ier.ne.0) write(6,*)'Warning q not found' 
       if(temp.eq.0.0d0) then 
         tget = 300.0 
         call        parget('temp    ',tget,iadda,ier) 
         temp = tget 
       endif
       if(temp <= 1d0) then
         write(*,*)"WARNING: temp is not properly defined! temp= ",temp
       endif 
                                                                        
                                                                        
      th_nzimmaz     = 0 
      sumq     = 0 
      sumnorm = 0 
      nqw     = 15 
      dqw     = 4*q_width/nqw 
      if(q_width.eq.0) then 
        nqw = 0 
        q_width = 1.0d0 
      endif 
                                                                        
      do i=-nqw,nqw 
       qzz = qz + i*dqw 
       if(qzz.gt.0) then 
        fn   = dexp( -(i*dqw/q_width)**2) 
                                                                        
! --- get center of mass diffusion only effectve if diff .ne. 0  ---    
         dr = diff 
                                                                        
                                                                        
         call NzimPnaY(qzz,tau,temp,dr,etasolv,N,R,                     &
     &                pmax,pwidth,nue,alpha, tau_zif, Sq,Sqt)                     
                                                                       
                                                                        
         sumq = sumq + fn*Sqt 
         sumnorm = sumnorm + Sq*fn 
                                                                        
                                                                        
        endif 
       enddo 
                                                                        
       if(sumnorm.gt.0.0d0) then 
          th_nzimmaz = a0*sumq/sumnorm 
       else 
          th_nzimmaz = 0 
       endif 
                                                                        
                                        ! in cm**2/s                    
       dr        = dr /( 1d-9 / 1d-16 ) 
       call        parset('diff    ',sngl(dr),iadda,ier) 
!                                                                       
       return 

CONTAINS                                          

                                                                        
                                                                        
       subroutine NzimPnaY(q,t,temp,Dr,etasolv,N,R,pmax,pwidth,nue,alpha,tzif,Sq,Sqt)                  
!      ==============================================================================                 
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
!    tzif  ----> ZIF model tau_i IN NanoSeconds !.
! Output parameters:                                                    
!    Sq    <--- S(Q)                                                    
!    Sqt   <--- S(Q,t)                                                  
! ------------------------------------------------------------          
!                                                                       
       implicit none 
                                                                        
       double precision kb, pi 
       parameter(kb=1.380662d-23) 
       parameter(pi=3.141592654d0) 
                                                                        
       double precision, intent(in)     :: q,t,temp
       double precision, intent(inout)  :: Dr
       double precision, intent(in)     :: etasolv,R
       double precision, intent(in)     :: pmax, pwidth, nue, alpha, tzif 
       double precision, intent(out)    :: Sq,Sqt 


       integer :: N, nn,mm,p 
                                                                        
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
                    (dfloat(p)**(3*nue)+alpha*dfloat(p)**(4-nue))       
                                                                       !
!                                 --- since we use R-end-to-end durectly
!                                     here the chain statistics is impli
!                   ----- this is the Gaussian chain factor (hope that i
                                                                       !
            tau_p = tau_p * 1d9  + tzif 
                                                                        
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
 
!      write(6,'(1x,5E14.6)')q,t,Sq,Sqt, Sqt/Sq                                                                         

       return 
 
 end subroutine NzimPnaY                                      
                                                                        
                                                                        
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
      END function transf_p   


END  function th_nzimmaz                                        
