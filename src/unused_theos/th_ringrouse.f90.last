      FUNCTION th_ringrouse (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ================================================================================               
!                                                                       
! -------> nrouse_sh4 <--------
! selfhealing multimodal distribution mix of Rouse        
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
      real                            :: th_ringrouse

      integer          :: ier, nparx

      double precision :: temp, tau, q, wl4, l 
      double precision :: a0 , qz, R, dr
      double precision :: sq, sqt
      double precision :: diff, r02, tc, a_cross, nu_subdiff, rr
      double precision :: amod(1000)    ! quick an dirty upper limit for mode modifiers
      
      REAL             :: qget, tget 
      integer          :: n
      logical          :: sqof0

      integer          :: iadda 
      common/thiadd/iadda

 
 ! ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'ringrous'
 !##!        thnam = 'nrouse_h'
         nparx = 18
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,' exceeds current max. = ',i8)
           th_ringrouse = 0
           return
         endif
         npar = nparx
!         --------------> set the number of parameters
         parnam(1) = 'intensit'
         parnam(2) = 'wl4     '
         parnam(3) = 'n_segmnt'   ! segments in ring
         parnam(4) = 're_lin  '   ! re of one ring as linear chain
         parnam(5) = 'diff    '   ! long tiome cm-diffusion A**2/ns
         parnam(6) = 'r02     '   ! r**2(tc) in A**2 (sublinear diffusion)
         parnam(7) = 'tc      '   ! time scale of sublinear diffusion (ns)
         parnam(8) = 'nu_subd '   ! exponent of sublinear diffusion
         parnam(9) = 'a_cross '   ! crossover exponent for sublin--> linear diffusion
         parnam(10)= 'a1_mod  '   ! mode modifying factor
         parnam(11)= 'a2_mod  '   ! mode modifying factor
         parnam(12)= 'a3_mod  '   ! mode modifying factor
         parnam(13)= 'a4_mod  '   ! mode modifying factor
         parnam(14)= 'a5_mod  '   ! mode modifying factor
         parnam(15)= 'a6_mod  '   ! mode modifying factor
         parnam(16)= 'a7_mod  '   ! mode modifying factor
         parnam(17)= 'a8_mod  '   ! mode modifying factor
         parnam(18)= 'a9_mod  '   ! mode modifying factor

         th_ringrouse = 0
         return
       endif
 
       sqof0 = .false.
!  ---- calculate theory here -----
       tau        = x
       a0         = pa(1)
       Wl4        = abs(pa(2))
       N          = nint(pa(3))
       R          = abs(pa(4))
       diff       = abs(pa(5)) 
       r02        = abs(pa(6))
       tc         = abs(pa(7))
       nu_subdiff = abs(pa(8))
       a_cross    = abs(pa(9)) ! 1=>smooth crossover .... 100=>sharp crossover
       amod       = 1d0
       amod(1:9)  = abs(pa(10:18))

       qget = 0.01
       call        parget('q       ',qget,iadda,ier)
       qz   = qget
       if(ier.ne.0) then 
         write(6,*)'Assume Q as independend variable!'
         qz   = x
         tau = 0.0d0
         call parget('tau     ',qget,iadda,ier)
         if(ier.eq.0) then
          tau = qget
         endif
         write(6,*)'tau = ',tau 
         sqof0 = .true.
       endif
             
       call Nrouse_ring(qz,tau,Wl4,N,R,amod,l, Sq,Sqt)
 
       if(sqof0) then
          th_ringrouse =  a0 * Sqt
       else
          th_ringrouse =  a0 * Sqt/Sq
       endif

!!! Multiply the subdiff expression !!!
       rr = ( (r02*(tau/tc)**nu_subdiff)**a_cross + (6*diff*tau)**a_cross)**(1d0/a_cross)
       th_ringrouse = th_ringrouse * exp( -qz*qz*rr/6d0 )
         

 


       call        parset('l       ',sngl(l),iadda,ier)      ! in ns A units
       call        parset('wl4     ',sngl(wl4),iadda,ier)
 
       return
       end




       subroutine Nrouse_ring(q,t,Wl4,Nb,R,amod,l,SqAve,SqtAve)
!      ======================================================
!
! Rouse expression for a ring polymer
! with Nb segments each, the ideal Ree of the linear version of the polymer is given by R.
! amod contains mode modifiers
! Input parameters:
!    q     ----> momentum transfer in A**-1
!    t     ----> time in nano-sec
!    Wl4   ----> Rouse rate in A**4/ns
!    N     ----> number of chain segments in one block
!    R     ----> end-to-end distance of the polymer molecule
!    amodm ----> mode modiufiers
! Output parameters:
!    l     <--- segment length
!    Sq    <--- S(Q)
!    Sqt   <--- S(Q,t)
! ------------------------------------------------------------
!
       implicit none

       double precision kb, pi
       parameter(kb=1.380662d-23)
       parameter(pi=3.141592654d0)

       double precision, intent(in)     ::  q,t,Wl4,R
       double precision, intent(in)     ::  amod(1000)
       double precision, intent(out)    ::  l, SqAve, SqtAve
       integer,          intent(in)     ::  Nb

       integer                          ::  nn,mm,p

       double precision :: W, rate_p, kbt, Sq0, arg1, arg2, Dr
       double precision :: a0,e0, ff2, ffc,    arg10,arg20
       double precision :: aa1 , aa2, weight, Sqt, Sq, weightsum, wmode


 

       integer          :: N, Mblocks, n1, n2, mmax
       double precision :: pcindex

       integer          :: iout, itrans
        
       if(Nb.le.0) then
         W  = 999
         Sq = 999
         Sqt= 999
         write(6,*)'Error Number of chain segments is <= 0!',Nb
         return
       endif

! ---- determine the segment length l ----
       l = sqrt(R**2/Nb)       
       
! ---- and the Rousefactor ----

       W   = Wl4 / l**4  

! ---- init sums ----
       Sq0 = 0
       SqAve  = 0
       SqtAve = 0
       weightsum = 0


           Sq0 = 0
           Sq  = 0
           Sqt = 0

           N      =  Nb 
           n1     =  1
           n2     =  N

! ---- Do the sums -----
           do nn = n1,n2
            do mm = n1,n2
              arg1 = -(q**2)*(       abs(nn-mm)*(N-abs(nn-mm))*(l**2)/(6.0d0*N))
              arg10= -(q**2)*(       abs(nn-mm)*(N-abs(nn-mm))*(l**2)/(6.0d0*N))
              ff2  = -2*N*(l*q)**2/(3*pi**2)
        
              arg2 = 0
              arg20= 0
              do p = 2,N,2           
                rate_p = 2*W*(1-cos((pi*p)/dfloat(N)))   ! only even modes for ring
                 
                wmode = amod(p/2)
                   
                a0    = -t * rate_p

                e0    = 1.0d0-exp(a0)
                
                ffc   = cos((pi*p*(nn-mm))/dfloat(N)) 
                ffc   = ffc / (p**2)
                ffc   = ffc * wmode                    !####new
    
                arg2  = arg2  + ffc*e0  
                arg20 = arg20 + ffc
    
              enddo   
              arg2  = arg2  * ff2
              arg20 = arg20 * ff2
    
              aa1 = arg10
              aa2 = arg1+arg2
    
    
              Sq  = Sq  + exp(aa1)
              Sqt = Sqt + exp(aa2)
    
            enddo
           enddo
    
           Sq  = Sq /Nb
           Sqt = Sqt/Nb
    

       SqAve  = Sq
       SqtAve = Sqt

       return
       end

