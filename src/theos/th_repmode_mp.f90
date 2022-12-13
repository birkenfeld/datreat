 FUNCTION th_repmode_mp(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  Mode modified Rouse modified according to simulation analysis by Grest et al.
!  J.T. Kalathi, S.K. Kumar, M. Rubinstein and G.S. Grest; Macromolecules (2014) 47, 6925-6931
      use theory_description 
      implicit none 
      real    :: th_repmode_mp
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
     integer          :: n          ! effective number of segments in summation                                       
     double precision :: wscalinf   ! high mode limit factor for w                                                    
     double precision :: wscal1     ! lowmode scale factor for w                                                      
     double precision :: wnu        ! sharpness factor for transition function of W                                   
     double precision :: betainf    ! high mode limit for beta                                                        
     double precision :: betamin    ! minimum value of beta                                                           
     double precision :: sigmabet   ! width parameter for beta dip                                                    
     double precision :: db2        ! parameter for the limit at low mode numbers                                     
     double precision :: ar         ! amplitude prefactor for mode amplitude expression (propto C-infinity)           
     double precision :: c          ! coefficient c in mode amplitude expression <x**2>=ar*(1-c/sqrt(N/p))/(4*sin(p*pi
     double precision :: nescalw    ! scaling factor for Ne in Wscale expression
     double precision :: nescalb    ! scaling factor for Ne in beta   expression
! the recin parameter representation 
     double precision :: q          ! q-value                                                                         
! the reout parameter representation 
     double precision :: l          ! effective segment length                                                        
     double precision :: ne         ! number of segments of length l in entanglement                                  
     double precision :: w          ! basic rate                                                                      
     double precision :: w1         ! p=1 rate                                                                        
     double precision :: wn         ! p=n rate                                                                        
 
     double precision   :: t

     double precision   :: sqt, sq
     double precision   :: temp
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'repmode2'
       nparx =       16
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_repmode_mp = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " Mode modified Rouse modified according to simulation analysis by Grest et al." //cr//parspace//&
                                " ACCELERATED (OMP+IMPROVED ALG.) VERSION!"
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
        parnam (14) = 'cc      '  ! coefficient c in mode amplitude expression <x**2>=ar*(1-c/sqrt(N/p))/(4*sin(p*pi
        parnam (15) = 'nescalw '  ! Ne scaling in Wscal
        parnam (16) = 'nescalb '  ! Ne scaling in beta 
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
        th_param_desc(15,idesc) = "scaling factor for Ne in Wscal expression" !//cr//parspace//&
        th_param_desc(16,idesc) = "scaling factor for Ne in beta expression " !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q-value (default 0.1)"
        th_file_param(  2,idesc) = "temp     > T-value (default 500K)"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "l        > effective segment length"
        th_out_param(  2,idesc) = "ne       > number of segments of length l in entanglement"
        th_out_param(  3,idesc) = "w        > basic rate"
        th_out_param(  4,idesc) = "w1       > p=1 rate"
        th_out_param(  5,idesc) = "wn       > p=n rate"
! 
        th_repmode_mp = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      abs( pa( 1) )
      wl4      =      abs( pa( 2) )
      re       =      abs( pa( 3) )
      d        =      abs( pa( 4) )
      N        =     nint( pa( 5) )
      wscalinf =      abs( pa( 6) )
      wscal1   =      abs( pa( 7) )
      wnu      =      abs( pa( 8) )
      betainf  =      abs( pa( 9) )
      betamin  =      abs( pa(10) )
      sigmabet =      abs( pa(11) )
      db2      =      abs( pa(12) )
      ar       =      abs( pa(13) )
      c        =      abs( pa(14) )
      nescalw  =      abs( pa(15) )
      nescalb  =      abs( pa(16) )
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: q-value
      xh =  0.1
      call parget('q       ',xh,iadda,ier)
      q        = xh
! >>> extract: temp-value
      xh =  500d0
      call parget('temp    ',xh,iadda,ier)
      temp     = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     t   = x
     call nrepmod(q,t,N,Sq,Sqt)
     th_repmode_mp = Sqt/Sq
 
! ---- writing computed parameters to the record >>>  
      call parset('l       ',sngl(l),iadda)
      call parset('ne      ',sngl(ne),iadda)
      call parset('w       ',sngl(w),iadda)
      call parset('w1      ',sngl(w1),iadda)
      call parset('wn      ',sngl(wn),iadda)
 
 

 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 



       subroutine nrepmod(q,t,N,Sq,Sqt)
!      ==============================
!
! Rouse expression for a chain of finite length:
! Input parameters:
!    q     ----> momentum transfer in A**-1
!    t     ----> time in nano-sec
!    n     ----> length of chain
!    Sq    <--- S(Q)
!    Sqt   <--- S(Q,t)
! ------------------------------------------------------------
!
       implicit none

       double precision kb, pi
       parameter(kb=1.380662d-23)
       parameter(pi=3.141592654d0)

       double precision, intent(in)   :: q,t
       integer         , intent(in)   :: N
       double precision, intent(out)  :: sq, sqt
  
       integer                        :: nn ,mm, p

       double precision ::  tau_p, kbt, Sq0, arg1, arg2
       double precision ::  a0,e0, ff2, ffc
       double precision ::  aa1 , aa2
       double precision ::  xi, beta, wscale
 
 !      integer :: iout

       integer, parameter :: maxp = 1000
       double precision, save   :: beta_vec(maxp)
       double precision, save   :: wscal_vec(maxp)
       double precision, save   :: amod_vec(maxp)

 
! save previsou variable for recomputation of beta, wscal vectors
       integer         , save   :: n0          = 0
       double precision, save   :: ne0         = 0
       double precision, save   :: wscalinf0   = 0
       double precision, save   :: wscal10     = 0 
       double precision, save   :: wnu0        = 0 

       double precision, save   :: betainf0    = 0
       double precision, save   :: betamin0    = 0
       double precision, save   :: db20        = 0
       double precision, save   :: sigmabet0   = 0

       double precision, save   :: ar0         = 0
       double precision, save   :: c0          = 0

       double precision, save   :: nescalw0    = 0
       double precision, save   :: nescalb0    = 0

    
       double precision :: cosarray(N,N), ewfac(N)
       integer          :: mp


 

       if(N.le.0) then
         W  = 999
         Sq = 999
         Sqt= 999
         write(6,*)'Error Number of chain segments is <= 0!',N
         return
       endif

! ---- determine the segment length l and the number of effective segments per entanglement----
       l  = sqrt(re**2/N)
       ne = (d/l)**2
       
! ---- and the Rousefactor ----
       kbt = temp*kb            ! in Joule = kg*m**2/s**2
       kbt = kbt * 100          ! in         kg*A**2/ns**2
       xi  = 3*kbt*l**2 / wl4
       W   = 3*kbt/(xi*(l**2))  ! in 1/ns

! prepare the mode dependent scaling factor as vectors (to speed up and create info files) !
iv:    if(   ( n        .ne.  n0        ) .or.  &       
             ( ne       .ne.  ne0       ) .or.  &       
             ( wscalinf .ne.  wscalinf0 ) .or.  &       
             ( wscal1   .ne.  wscal10   ) .or.  &       
             ( wnu      .ne.  wnu0      ) .or.  &       
             ( betainf  .ne.  betainf0  ) .or.  &       
             ( betamin  .ne.  betamin0  ) .or.  &       
             ( db2      .ne.  db20      ) .or.  &       
             ( sigmabet .ne.  sigmabet0 ) .or.  &            
             ( ar       .ne.  ar0       ) .or.  &            
             ( c        .ne.  c0        ) .or.  &    
             ( nescalw  .ne.  nescalw0  ) .or.  &    
             ( nescalb  .ne.  nescalb0  )       &    
             ) then    

               n0        =  n        
               ne0       =  ne       
               wscalinf0 =  wscalinf 
               wscal10   =  wscal1   
               wnu0      =  wnu      
               betainf0  =  betainf  
               betamin0  =  betamin  
               db20      =  db2       
               sigmabet0 =  sigmabet
               ar0       =  ar 
               c0        =  c

               mp        = min(N,maxp)

!$OMP PARALLEL DO     
             do p=1,N
               wscal_vec(p)  = Wscalfun( dfloat(n)/dfloat(p), ne * nescalw )
               beta_vec(p)   = betafun(  dfloat(n)/dfloat(p), ne * nescalb )
               amod_vec(p)   = amplimod(  dfloat(n)/dfloat(p) )
             enddo
!$OMP END PARALLEL DO 

             open(51,file="repmode_scales.dat")
               write(51,'(a,i7)   ') "n        := ",n       
               write(51,'(a,e14.7)') "ne       := ",ne      
               write(51,'(a,e14.7)') "wscalinf := ",wscalinf
               write(51,'(a,e14.7)') "wscal1   := ",wscal1  
               write(51,'(a,e14.7)') "wnu      := ",wnu     
               write(51,'(a,e14.7)') "betainf  := ",betainf 
               write(51,'(a,e14.7)') "betamin  := ",betamin 
               write(51,'(a,e14.7)') "db2      := ",db2      
               write(51,'(a,e14.7)') "sigmabet := ",sigmabet
               write(51,'(a,e14.7)') "ar       := ",ar
               write(51,'(a,e14.7)') "c        := ",c
               write(51,'(a)'      ) "  p    N      N/p     W_scale(p)      beta(p)  "
               write(51,'(i7,i7,4e15.7)')(p,n,dfloat(n)/dfloat(p),  wscal_vec(p),  beta_vec(p), amod_vec(p), p=1,min(N,maxp))
             close(51)
   
        endif iv  
 




!$OMP PARALLEL DO     
       do nn=1,N
        do p=1,N
         cosarray(nn,p) = cos((pi*p*nn)/dfloat(N)) * sqrt(amod_vec(p)) / ((2*N/pi)*sin(p*pi/(2*dfloat(N))))
        enddo
       enddo
!$OMP END PARALLEL DO   


!$OMP PARALLEL DO   
       do p=1,N
          ewfac(p) = 1.0d0-exp(-(t*2*W*Wscal_vec(p)*(1-cos((pi*p)/dfloat(N))))**beta_vec(p))
       enddo
!$OMP END PARALLEL DO    


! ---- init sums ----
       Sq0 = 0
       Sq  = 0
       Sqt = 0
       ff2  = -2*N*(l*q)**2/(3*pi**2)

!$OMP PARALLEL DO REDUCTION(+:Sq,Sqt) PRIVATE(arg2)
! ---- Do the sums -----
       do nn = 1,N
        do mm = 1,N

          arg2 = sum(cosarray(nn,1:mp) *cosarray(mm,1:mp) *  ewfac(1:mp) ) * ff2
          Sq  = Sq  + exp(-(q**2)*( abs(nn-mm)*(l**2)/6.0d0))
          Sqt = Sqt + exp(-(q**2)*( abs(nn-mm)*(l**2)/6.0d0) + arg2)

        enddo
       enddo
!$OMP END PARALLEL DO

       Sq  = Sq /N
       Sqt = Sqt/N

!       if(iout().gt.1)write(6,'(1x,a,6E14.6)')'q,t,Sq,Sqt, Sqt/Sq, w=', q,t,Sq,Sqt, Sqt/Sq, w

       return
       end

    double precision function Fermi(x)
      double precision, intent(in)    :: x
      Fermi =  1d0/(exp(x)+1d0)
    end function Fermi


    double precision function wscalfun( np, nent )
      implicit none
      double precision, intent(in)    :: np      ! N/p
      double precision, intent(in)    :: nent    ! Ne

      double precision :: alpha

      alpha = log(wscalinf/wscal1)
    
      wscalfun=wscalinf*exp(-alpha)*exp(alpha*Fermi(log(Np/Nent)*wnu))

    end function wscalfun


    double precision function betafun( np, nent )
      implicit none
      double precision, intent(in)    :: np      ! N/p
      double precision, intent(in)    :: nent    ! Ne


      double precision :: b0, db, mu

      b0 = betainf
      db = betainf-betamin
      mu = 1d0

      betafun=erf(betainf-db*exp(-log(Np/Nent)**2/sigmabet**2)+db2*log(Np/Nent)**2*Fermi(-log(Np/Nent)*mu))

    end function betafun

    double precision function amplimod( np )
      implicit none
      double precision, intent(in)    :: np      ! N/p
 
      
      amplimod = ( ar * (1d0-c/sqrt(np)) )

    end function amplimod



 end function th_repmode_mp
