      FUNCTION th_nrouconfy (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ================================================================================               
!                                                                       
! -------> th_nrouconfy <--------
! experimental intro if tau_e 
! Totally exploratory, must be refiend or abandoned
! just a tool to tray some funny ideas       
!                                                                       
!                
      implicit none  
      double precision, parameter :: Pi                  =     3.141592653589793238462643d0    
                                                 
      CHARACTER(len=8), intent(inout) :: thnam, parnam (20) 
      real, intent(in)                :: x,  pa (20)
      integer, intent(in)             :: ini, mbuf
      integer, intent(inout)          :: npar
      integer, intent(inout)          :: nopar            ! Anzahl der Parameter data
      character*80, intent(inout)     :: napar(mbuf)      ! name des parameters n
      real, intent(inout)             :: params(mbuf)     ! value des parameters n
      real                            :: th_nrouconfy

      integer          :: ier, nparx

      double precision :: temp, tau, q, wl4, l 
      double precision :: a0 , qz, R, dr, W, a, b
      double precision :: diff , tau_e
      double precision :: sq, sqt
      REAL             :: qget, tget 
      integer          :: n

      integer          :: iadda 
      common/thiadd/iadda

 
 ! ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'nrouconf'
         nparx = 9
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,' exceeds current max. = ',i8)
           th_nrouconfy = 0
           return
         endif
         npar = nparx
!         --------------> set the number of parameters
         parnam(1) = 'intensit'
         parnam(2) = 'wl4     '   ! rouse rate
         parnam(3) = 'n       '   ! segments 
         parnam(4) = 're      '   ! re of one block
         parnam(5) = 'diff    '   ! diffusion 
         parnam(6) = 'tau_e   '   ! confinement time
         parnam(7) = 'temp    '   ! temperature (maybe overwritten)
         parnam(8) = 'a       '   ! amplitude mod
         parnam(9) = 'b       '   ! rate mod

   
         th_nrouconfy = 0
         return
       endif
 
!  ---- calculate theory here -----
       tau      = x
       a0       = pa(1)
       Wl4      = abs(pa(2))
       N        = nint(pa(3))
       R        = abs(pa(4))
       diff     = abs(pa(5)) 
       tau_e    = abs(pa(6))
       temp     = pa(7)
       a        = 2*atan(abs(pa(8)))/Pi
       b        = abs(pa(9))

       qget = 0.01
       call        parget('q       ',qget,iadda,ier)
       if(ier.eq.0) qz   = qget
       call        parget('temp    ',tget,iadda,ier)
       if(ier.eq.0) temp = tget
             

       call NrouseC(qz, tau,  temp, diff, wl4, N, R, tau_e, a,b,    W, l, Sq,Sqt)

       th_nrouconfy =  a0 * Sqt/Sq

       call        parset('l       ',sngl(l),iadda,ier)      ! in ns A units
       dr        = dr /( 1d-9 / 1d-16 ) ! in cm**2/s
       call        parset('diff    ',sngl(dr),iadda,ier)
       call        parset('w       ',sngl(w),iadda,ier)
 
     CONTAINS
       


 subroutine NrouseC(q,t,temp,Dr,wl4,N,R,tau_e,A,B, W, l, Sq,Sqt)
!===========================================================
       implicit none
!
! Rouse expression for a chain of finite length with time controlled confinenment (experimental):
! Input parameters:
 double precision, intent(in)  :: q             !    q     ----> momentum transfer in A**-1
 double precision, intent(in)  :: t             !    t     ----> time in nano-sec
 double precision, intent(in)  :: temp          !    temp  ----> temperature in K
 double precision, intent(inout)  :: Dr         !    Dr    ----> center of mass diffusion constant in A**2/ns, if 0 <-- Rouse-expectation
 double precision, intent(in)  :: wl4           !    wl4   ----> friction coefficient in A**4/ns
 integer         , intent(in)  :: N             !    N     ----> number of chain segments
 double precision, intent(in)  :: R             !    R     ----> end-to-end distance of the polymer molecule
 double precision, intent(in)  :: tau_e         !    tau_e ----> control time
 double precision, intent(in)  :: A             !    amplitude distribution
 double precision, intent(in)  :: B             !    rate modifier local rep
! Output parameters:
 double precision, intent(out) :: W             !    W     <--- "Rouse factor" 3kT/(xi*l**2); R**2=N*l**2
 double precision, intent(out) :: l             !    l     <--- "Segment length l"
 double precision, intent(out) :: Sq            !    Sq    <--- S(Q)
 double precision, intent(out) :: Sqt           !    Sqt   <--- S(Q,t)
! ------------------------------------------------------------
!


       double precision, parameter :: kb=1.380662d-23
       double precision, parameter :: pi= 4d0*atan(1d0)

       double precision  :: xi
       integer           :: nn,mm,ip

       double precision  :: tau_p, kbt, Sq0,  arg1, arg2
       double precision  :: a0,e0, ff2, ffc,  arg10,arg20
       double precision  :: aa1 , aa2
       double precision  :: p
       double precision  :: tc

!       integer, external :: iout
       
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

    
!       p0fix = 0
!       pfac  = 1
 

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
            p = ip
            tau_p = 2*W*(1-cos((pi*p)/dfloat(N)))

            tc = 0.5d0 * sqrt(Pi) * tau_e * ERF( t / tau_e )   ! hier findet die Begrenzung statt !
!            tc = tau_e * atan( t / tau_e )   ! hier findet die Begrenzung statt !

!%            a0    = -tc * tau_p
!%            if(a0.lt.-200.0d0) a0 = -200.0d0
!%            e0    = 1.0d0-exp(a0)

            e0 = A*(1d0-exp(-t*tau_p * B)) + (1d0-A)*(1d0-exp(-tc*tau_p))

            
            ffc   = cos((pi*p*nn)/dfloat(N)) * cos((pi*p*mm)/dfloat(N))
            ffc   = ffc / (p**2)
            ffc   = ffc

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

!       if(iout().gt.1)write(6,'(1x,a,6E14.6)') &
!             'q,t,Sq,Sqt, Sqt/Sq, w=', q,t,Sq,Sqt, Sqt/Sq, w 

       
       end subroutine NrouseC


       
     end FUNCTION th_nrouconfy
