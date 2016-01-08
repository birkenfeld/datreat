      FUNCTION th_nrouse_sh (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ================================================================================               
!                                                                       
! -------> nrouse_sh <--------
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
      real                            :: th_nrouse_sh

      integer          :: ier, nparx

      double precision :: temp, tau, q, wl4, l 
      double precision :: a0 , qz, R, dr
      double precision :: m_aver, wcut
      double precision :: diffscal 
      double precision :: sq, sqt
      REAL             :: qget, tget 
      integer          :: n
      logical          :: sqof0

      integer          :: iadda 
      common/thiadd/iadda

 
 ! ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'nrouse_s'
         nparx = 7
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,' exceeds current max. = ',i8)
           th_nrouse_sh = 0
           return
         endif
         npar = nparx
!         --------------> set the number of parameters
         parnam(1) = 'intensit'
         parnam(2) = 'wl4     '
         parnam(3) = 'n_segmnt'   ! segments in one block
         parnam(4) = 're      '   ! re of one block
         parnam(5) = 'diffscal'   ! diffusion scaling
         parnam(6) = 'm_aver  '   ! average aggregation number
         parnam(7) = 'wcut    '   ! weight cutoff
   
         th_nrouse_sh = 0
         return
       endif
 
       sqof0 = .false.
!  ---- calculate theory here -----
       tau      = x
       a0       = pa(1)
       Wl4      = abs(pa(2))
       N        = nint(pa(3))
       R        = abs(pa(4))
       diffscal = abs(pa(5)) 
       m_aver   = pa(6)
       wcut     = pa(7)

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

       dr = diffscal
             
       call Nrouse_sh(qz,tau,dr,Wl4,N,R, m_aver, wcut, l,Sq,Sqt)
 
       if(sqof0) then
          th_nrouse_sh =  a0 * Sqt
       else
          th_nrouse_sh =  a0 * Sqt/Sq
       endif

 


       call        parset('l       ',sngl(l),iadda,ier)      ! in ns A units
       dr        = dr /( 1d-9 / 1d-16 ) ! in cm**2/s
       call        parset('diffav  ',sngl(dr),iadda,ier)
 
       return
       end




       subroutine Nrouse_sh(q,t,DrScal,Wl4,Nb,R,m_aver,wcut,l,Daver, SqAve,SqtAve)
!      ===========================================================================
!
! Rouse expression for a chain of finite length:
! Input parameters:
!    q     ----> momentum transfer in A**-1
!    t     ----> time in nano-sec
!    Dr    ----> scaling factor for Rouse diffusuion
!    Wl4   ----> Rouse rate in A**4/ns
!    N     ----> number of chain segments in one block
!    R     ----> end-to-end distance of the polymer molecule
!    wcut  ----> weight cutoff for mblocks summation
!    m_aver----> average number of blocks
! Output parameters:
!    Daver <--- averaged effective diffusuion
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
       double precision, intent(in)     ::  m_aver
       double precision, intent(in)     ::  wcut
       double precision, intent(in)     ::  DrScal
       double precision, intent(out)    ::  Daver, l, SqAve, SqtAve
       integer,          intent(in)     ::  Nb

       integer                          ::  nn,mm,p

       double precision :: W, tau_p, kbt, Sq0, arg1, arg2, Dr
       double precision :: a0,e0, ff2, ffc,    arg10,arg20
       double precision :: aa1 , aa2, weight, Sqt, Sq, weightsum

 

       integer          :: N, Mblocks, n1, n2, mmax
       double precision :: pcindex

       integer          :: iout, itrans
        
       if(N.le.0) then
         W  = 999
         Sq = 999
         Sqt= 999
         write(6,*)'Error Number of chain segments is <= 0!',N
         return
       endif

 

       pcindex = 1d0-1d0/m_aver
       mmax    = log(wcut)/log(pcindex)


! ---- determine the segment length l ----
       l = sqrt(R**2/Nb)       
       
! ---- and the Rousefactor ----

       W   = Wl4 / l**4


! ---- init sums ----
       Sq0 = 0
       SqAve  = 0
       SqtAve = 0
       weightsum = 0
       Daver     = 0


mloop: do mblocks=1,mmax

           Sq0 = 0
           Sq  = 0
           Sqt = 0

           weight    =  mblocks * pcindex ** mblocks
           weightsum = weightsum + weight
           N      =  Nb * mblocks
           n1     =  max(1,N/2-Nb/2)
           n2     =  min(N,n1+Nb)

           Dr    = DrScal * Wl4 / ( 3*Nb*mblocks *l**2) 
           Daver = Daver + Dr * weight




! ---- Do the sums -----
           do nn = n1,n2
            do mm = n1,n2
              arg1 = -(q**2)*(Dr*t + abs(nn-mm)*(l**2)/6.0d0)
              arg10= -(q**2)*(       abs(nn-mm)*(l**2)/6.0d0)
              ff2  = -2*N*(l*q)**2/(3*pi**2)
        
              arg2 = 0
              arg20= 0
              do p = 1,N
                tau_p = 2*W*(1-cos((pi*p)/dfloat(N))) 
                a0    = -t*tau_p
                if(a0.lt.-200.0d0) a0 = -200.0d0
                e0    = 1.0d0-exp(a0)
                
                ffc   = cos((pi*p*nn)/dfloat(N)) * cos((pi*p*mm)/dfloat(N))
                ffc   = ffc / (p**2)
    
                arg2  = arg2  + ffc*e0
                arg20 = arg20 + ffc
    
              enddo   
              arg2  = arg2  * ff2
              arg20 = arg20 * ff2
    
              aa1 = arg10
!              if(aa1.lt.-300.0d0) aa1 = -300.0d0
!              if(aa1.gt. 300.0d0) aa1 =  300.0d0
    
              aa2 = arg1+arg2
!              if(aa2.lt.-300.0d0) aa2 = -300.0d0
!              if(aa2.gt. 300.0d0) aa2 =  300.0d0
    
    
              Sq  = Sq  + exp(aa1)
              Sqt = Sqt + exp(aa2)
    
            enddo
           enddo
    
           Sq  = Sq /Nb
           Sqt = Sqt/Nb
    
           SqAve  = SqAve  + Sq * weight
           SqtAve = SqtAve + Sqt* weight
    
       enddo mloop

       SqAve  = SqAve  + Sq / weightsum
       SqtAve = SqtAve + Sqt/ weightsum
       Daver  = Daver       / weightsum


       if(iout().gt.0)write(6,'(1x,6E14.6)')q,t,SqAve,SqtAve, SqtAve/SqAve, w 

       return
       end

