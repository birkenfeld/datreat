 FUNCTION th_repmode(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  Mode modified Rouse modified according to simulation analysis by Grest et al.
!  J.T. Kalathi, S.K. Kumar, M. Rubinstein and G.S. Grest; Macromolecules (2014) 47, 6925-6931
      use theory_description 
      implicit none 
      real    :: th_repmode
      character(len=8) :: thnam, parnam (*) 
      real    :: pa (*) 
      real    :: x , xh
      integer :: mbuf, nparx, ier, ini, npar, iadda
      integer, intent(inout) :: nopar       
      character(len=80), intent(inout) :: napar(mbuf) 
      real, intent(inout) :: params(mbuf) 
     
      double precision, parameter :: Pi = 4*atan(1d0)
      integer                     :: actual_record_address
     
! the internal parameter representation 
     double precision :: ampli      ! prefactor                                                                       
     double precision :: wl4        ! reference scale for the Rouse rate                                              
     double precision :: re         ! end-end radius of chains                                                        
     double precision :: d          ! tube diameter (ree of subchain)                                                 
     double precision :: n          ! effective number of segments in summation                                       
     double precision :: wscalinf   ! high mode limit factor for w                                                    
     double precision :: wscal1     ! lowmode scale factor for w                                                      
     double precision :: wnu        ! sharpness factor for transition function of W                                   
     double precision :: betainf    ! high mode limit for beta                                                        
     double precision :: betamin    ! minimum value of beta                                                           
     double precision :: sigmabet   ! width parameter for beta dip                                                    
     double precision :: db2        ! parameter for the limit at low mode numbers                                     
     double precision :: ar         ! amplitude prefactor for mode amplitude expression (propto C-infinity)           
     double precision :: c          ! coefficient c in mode amplitude expression <x**2>=ar*(1-c/sqrt(N/p))/(4*sin(p*pi
! the recin parameter representation 
     double precision :: q          ! q-value                                                                         
! the reout parameter representation 
     double precision :: l          ! effective segment length                                                        
     double precision :: ne         ! number of segments of length l in entanglement                                  
     double precision :: w          ! basic rate                                                                      
     double precision :: w1         ! p=1 rate                                                                        
     double precision :: wn         ! p=n rate                                                                        
 
     double precision   :: t
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'repmode'
       nparx =       14
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_repmode = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " Mode modified Rouse modified according to simulation analysis by Grest et al."
       th_citation(idesc)     = " J.T. Kalathi, S.K. Kumar, M. Rubinstein and G.S. Grest; Macromolecules (2014) 47, 6925-6931"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'wl4     '  ! reference scale for the Rouse rate                                              
        parnam ( 3) = 're      '  ! end-end radius of chains                                                        
        parnam ( 4) = 'd       '  ! tube diameter (ree of subchain)                                                 
        parnam ( 5) = 'n       '  ! effective number of segments in summation                                       
        parnam ( 6) = 'wscalinf'  ! high mode limit factor for w                                                    
        parnam ( 7) = 'wscal1  '  ! lowmode scale factor for w                                                      
        parnam ( 8) = 'wnu     '  ! sharpness factor for transition function of W                                   
        parnam ( 9) = 'betainf '  ! high mode limit for beta                                                        
        parnam (10) = 'betamin '  ! minimum value of beta                                                           
        parnam (11) = 'sigmabet'  ! width parameter for beta dip                                                    
        parnam (12) = 'db2     '  ! parameter for the limit at low mode numbers                                     
        parnam (13) = 'ar      '  ! amplitude prefactor for mode amplitude expression (propto C-infinity)           
        parnam (14) = 'c       '  ! coefficient c in mode amplitude expression <x**2>=ar*(1-c/sqrt(N/p))/(4*sin(p*pi
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "reference scale for the Rouse rate" !//cr//parspace//&
        th_param_desc( 3,idesc) = "end-end radius of chains" !//cr//parspace//&
        th_param_desc( 4,idesc) = "tube diameter (ree of subchain)" !//cr//parspace//&
        th_param_desc( 5,idesc) = "effective number of segments in summation" !//cr//parspace//&
        th_param_desc( 6,idesc) = "high mode limit factor for w" !//cr//parspace//&
        th_param_desc( 7,idesc) = "lowmode scale factor for w" !//cr//parspace//&
        th_param_desc( 8,idesc) = "sharpness factor for transition function of W" !//cr//parspace//&
        th_param_desc( 9,idesc) = "high mode limit for beta" !//cr//parspace//&
        th_param_desc(10,idesc) = "minimum value of beta" !//cr//parspace//&
        th_param_desc(11,idesc) = "width parameter for beta dip" !//cr//parspace//&
        th_param_desc(12,idesc) = "parameter for the limit at low mode numbers" !//cr//parspace//&
        th_param_desc(13,idesc) = "amplitude prefactor for mode amplitude expression (propto C-infinity)" !//cr//parspace//&
        th_param_desc(14,idesc) = "coefficient c in mode amplitude expression <x**2>=ar*(1-c/sqrt(N/p))/(4*sin(p*pi/2/N)**2)" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q-value"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "l        > effective segment length"
        th_out_param(  2,idesc) = "ne       > number of segments of length l in entanglement"
        th_out_param(  3,idesc) = "w        > basic rate"
        th_out_param(  4,idesc) = "w1       > p=1 rate"
        th_out_param(  5,idesc) = "wn       > p=n rate"
! 
        th_repmode = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      wl4      =      pa( 2)
      re       =      pa( 3)
      d        =      pa( 4)
      n        =      pa( 5)
      wscalinf =      pa( 6)
      wscal1   =      pa( 7)
      wnu      =      pa( 8)
      betainf  =      pa( 9)
      betamin  =      pa(10)
      sigmabet =      pa(11)
      db2      =      pa(12)
      ar       =      pa(13)
      c        =      pa(14)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: q-value
      xh =  0.1
      call parget('q       ',xh,iadda,ier)
      q        = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     t   = x
     th_repmode = ! INSERT RESULTING VALUE OF TH EVALUATION HERE
 
! ---- writing computed parameters to the record >>>  
      call parset('l       ',sngl(l),iadda,ier)
      call parset('ne      ',sngl(ne),iadda,ier)
      call parset('w       ',sngl(w),iadda,ier)
      call parset('w1      ',sngl(w1),iadda,ier)
      call parset('wn      ',sngl(wn),iadda,ier)
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 



       subroutine NrouseY(q,t,temp,Dr,wl4,N,R, W, l,ifix, Sq,Sqt)
!      ========================================================
!
! Rouse expression for a chain of finite length:
! Input parameters:
!    q     ----> momentum transfer in A**-1
!    t     ----> time in nano-sec
!    temp  ----> temperature in K
!    Dr    ----> center of mass diffusion constant in A**2/ns, if 0 <-- Rouse-expectation
!    wl4   ----> friction coefficient in A**4/ns
!    N     ----> number of chain segments
!    R     ----> end-to-end distance of the polymer molecule
!    ifix  ----> if 0: normal Rouse, else segment(N) fix
! Output parameters:
!    W     <--- "Rouse factor" 3kT/(xi*l**2); R**2=N*l**2
!    l     <--- "Segment length l"
!    Sq    <--- S(Q)
!    Sqt   <--- S(Q,t)
! ------------------------------------------------------------
!
       implicit none

       double precision kb, pi
       parameter(kb=1.380662d-23)
       parameter(pi=3.141592654d0)

       double precision q,t,temp,Dr,xi,R, W,Sq,Sqt, wl4
       integer N, nn,mm,ifix,ip

       double precision l, tau_p, kbt, Sq0, arg1, arg2
       double precision a0,e0, ff2, ffc,    arg10,arg20
       double precision aa1 , aa2
       double precision p, p0fix, pfac

       integer iout

       if(N.le.0) then
         W  = 999
         Sq = 999
         Sqt= 999
         write(6,*)'Error Number of chain segments is <= 0!',N
         return
       endif

! ---- determine the segment length l ----
       l = sqrt(R**2/N)

! ---- and the Rousefactor ----
       kbt = temp*kb            ! in Joule = kg*m**2/s**2
       kbt = kbt * 100          ! in         kg*A**2/ns**2
       xi  = 3*kbt*l**2 / wl4
       W   = 3*kbt/(xi*(l**2))  ! in 1/ns

! ---- set the diffusion constant if input is zero --- !
       if(Dr.eq.0.0d0) then
         Dr = kbt/(N*xi)
       endif


       if(ifix.eq.0) then
         p0fix = 0
         pfac  = 1
       else
         p0fix = -0.5d0
         pfac  = 0.5d0   ! this is wrong ? !
         pfac=1         !! I think this has to be like this !!
       endif

! ---- init sums ----
       Sq0 = 0
       Sq  = 0
       Sqt = 0

! ---- Do the sums -----
       do nn = 1,N
        do mm = 1,N
          arg1 = -(q**2)*(Dr*t + abs(nn-mm)*(l**2)/6.0d0)
          arg10= -(q**2)*(       abs(nn-mm)*(l**2)/6.0d0)
          ff2  = -2*N*(l*q)**2/(3*pi**2)

          arg2 = 0
          arg20= 0
          do ip=1,N
            p = ip+p0fix
            tau_p = 2*W*(1-cos((pi*p)/dfloat(N)))
            a0    = -t*tau_p
            if(a0.lt.-200.0d0) a0 = -200.0d0
            e0    = 1.0d0-exp(a0)

            ffc   = cos((pi*p*nn)/dfloat(N)) * cos((pi*p*mm)/dfloat(N))
            ffc   = ffc / (p**2)
            ffc   = ffc*pfac

            arg2  = arg2  + ffc*e0
            arg20 = arg20 + ffc

          enddo
          arg2  = arg2  * ff2
          arg20 = arg20 * ff2

          aa1 = arg10
          if(aa1.lt.-300.0d0) aa1 = -300.0d0
          if(aa1.gt. 300.0d0) aa1 =  300.0d0

          aa2 = arg1+arg2
          if(aa2.lt.-300.0d0) aa2 = -300.0d0
          if(aa2.gt. 300.0d0) aa2 =  300.0d0


          Sq  = Sq  + exp(aa1)
          Sqt = Sqt + exp(aa2)

        enddo
       enddo

       Sq  = Sq /N
       Sqt = Sqt/N

       if(iout().gt.1)write(6,'(1x,a,6E14.6)')
     *        'q,t,Sq,Sqt, Sqt/Sq, w=', q,t,Sq,Sqt, Sqt/Sq, w

       return
       end

 end function th_repmode
