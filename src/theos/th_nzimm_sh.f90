      FUNCTION th_nzimm_sh (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
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
      real                            :: th_nzimm_sh

      integer          :: ier, nparx

      double precision :: temp, tau, q, l , taulim
      double precision :: etasolv, nue, alpha
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
         thnam = 'nzimm_s'
         nparx = 10
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,' exceeds current max. = ',i8)
           th_nzimm_sh = 0
           return
         endif
         npar = nparx
!         --------------> set the number of parameters
         parnam(1) = 'intensit'
         parnam(2) = 'etasolv '
         parnam(3) = 'n_segmnt'   ! segments in one block
         parnam(4) = 're      '   ! re of one block
         parnam(5) = 'diffscal'   ! diffusion scaling
         parnam(6) = 'm_aver  '   ! average aggregation number
         parnam(7) = 'wcut    '   ! weight cutoff in terms of mmax = wcut * (n at max of w)
         parnam(8) = 'nue     '   ! exponent for R=l*N**nue 
         parnam(9) = 'taulim  '   ! mode cutoff by weight exp(-tau(p)/taulim)
         parnam(10)= 'alpha   '   ! stiffness paramneter 
  
   
         th_nzimm_sh = 0
         return
       endif
 
       sqof0 = .false.
!  ---- calculate theory here -----
       tau      = x
       a0       = pa(1)
       etasolv  = abs(pa(2))
       N        = nint(pa(3))
       R        = abs(pa(4))
       diffscal = abs(pa(5)) 
       m_aver   = abs(pa(6))
       wcut     = abs(pa(7))
       nue      = abs(pa(8))
       taulim   = abs(pa(9))
       alpha    = pa(10)



       tget = 300.0
       call        parget('temp    ',tget,iadda,ier)
       temp = tget

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
             
       call Nzimm_sh(qz,tau,diffscal,etasolv,nue,alpha,temp,N,R, m_aver, wcut, l,dr,taulim, Sq,Sqt)
 
       if(sqof0) then
          th_nzimm_sh =  a0 * Sqt
       else
          th_nzimm_sh =  a0 * Sqt/Sq
       endif

 


       call        parset('l       ',sngl(l),iadda)      ! in ns A units
       dr        = dr /( 1d-9 / 1d-16 ) ! in cm**2/s
       call        parset('diffav  ',sngl(dr),iadda)
  
       return
       end




       subroutine Nzimm_sh(q,t,DrScal,etasolv,nue,alpha,temp,Nb,Rb,m_aver,wcut,l,Daver,taulim, SqAve,SqtAve)
!      ===========================================================================
!
! Rouse expression for an ensemble of polydisperse chains made up of pieces
! with Nb segments each, the ideal Ree of the Nb-wide pieces is given by R.
! The system contains of a polydispers mixture of 1, 2, 3 .... combinations of these
! chains, the Nb segments are either deuterated or protonated. The assumption is
! a random placement of h or d subchains, as expected as a result of a system with
! dynamical links.
! The size distribution is assume to be that of a polycondensation model with
! an average number of subchains: m_aver. wcut controlls the cutoff of the 
! summation over number of subchains which is approximately at m_aver * wcut.
!           (more precise -wcut/ln(p), with -1/ln(p) the n of the peak of the distribution.
! The diffusions are taken from Rousemodel with a scaling factor DrScal
! The marke block is inserted at relative position fpos in the long aggreate chain.  
! 
! Input parameters:
!    q     ----> momentum transfer in A**-1
!    t     ----> time in nano-sec
!    Dr    ----> scaling factor for Rouse diffusuion
!    etasolv----> solvent viscosity in SI
!    nue   ----> chain expansion coefficient (Gaussian = 0.5)
!    alpha ----> stiffness coeff
!    N     ----> number of chain segments in one block
!    Rb    ----> end-to-end distance of the polymer molecule
!    m_aver----> average number of blocks
!    wcut  ----> weight cutoff for mblocks summation in terms of wcut * location of max of distribution
! Output parameters:
!    Daver <--- averaged effective diffusuion
!    taulim<--- mode cut off beyond that time (chain breaking)         
!    l     <--- segment length
!    Sq    <--- S(Q)
!    Sqt   <--- S(Q,t)
! ------------------------------------------------------------
!
       implicit none

       double precision kb, pi
       parameter(kb=1.380662d-23)
       parameter(pi=3.141592654d0)

       double precision, intent(in)     ::  q,t,etasolv, nue, alpha ,Rb
       double precision, intent(in)     ::  m_aver
       double precision, intent(in)     ::  wcut
       double precision, intent(in)     ::  temp
       double precision, intent(in)     ::  DrScal
       double precision, intent(in)     ::  taulim
       double precision, intent(out)    ::  Daver, l, SqAve, SqtAve
       integer,          intent(in)     ::  Nb

       integer                          ::  nn,mm,p

       double precision :: tau_p, rate_p, kbt, Sq0, arg1, arg2, Dr, R, Rsi
       double precision :: a0,e0, ff2, ffc,    arg10,arg20, fdf
       double precision :: aa1 , aa2, weight, Sqt, Sq, weightsum, wmode

 

       integer          :: N, Mblocks, n1, n2, mmax
       double precision :: pcindex

       integer          :: iout, itrans
        
       if(Nb.le.0) then
         Sq = 999
         Sqt= 999
         write(6,*)'Error Number of chain segments is <= 0!',Nb
         return
       endif

 

       pcindex = 1d0-1d0/m_aver
       mmax    = nint(-1d0/log(pcindex) * wcut)


! ---- determine the segment length l ----
       l   = Rb / dfloat(Nb)**nue       
       kbt = temp*kb            ! in Joule = kg*m**2/s**2

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

           weight    =  mblocks * (pcindex ** (mblocks-1)) * (1-pcindex)**2 
           weightsum = weightsum + weight
           N      =  Nb * mblocks
           n1     =  1
           n2     =  N

           R      =  l * dfloat(N)**nue 
           Rsi    =  R * 1d-10            ! Re in meter                   

         ! do a linear interpolation between values for Gaussian and good solvent chain
         ! limting factors give in Doi-Edwards
           fdf= 0.196d0 + (0.203d0-0.196d0)*(nue-0.5d0)/(0.6d0-0.5d0)
           Dr = fdf*kbt/(Rsi*etasolv)      ! in m**2/s
           Dr = Dr*1d20/1d9                ! in A**2/ns

           Dr    = DrScal * Dr
           Daver = Daver  + Dr * weight

! ---- Do the sums -----
           do nn = n1,n2
            do mm = n1,n2
              arg1 = -(q**2)*(Dr*t + (abs(nn-mm)**(2*nue))*(l**2)/6.0d0)
              arg10= -(q**2)*(       (abs(nn-mm)**(2*nue))*(l**2)/6.0d0)
              ff2  = -2*(R*q)**2/(3*pi**2)
        
              arg2 = 0
              arg20= 0
              do p = 1,N
                tau_p = 0.325*etasolv*Rsi**3/kbt / (dfloat(p)**(3*nue)+alpha*dfloat(p)**(4-nue))   ! in s (?? 4-nue ok? fuer nue=1/2 ja)
!                                 --- since we use R-end-to-end durectly
!                                     here the chain statistics is implicit
!                   ----- this is the Gaussian chain factor (hope that it still is ok!
                tau_p  = tau_p * 1d9                                        ! in ns
 
                rate_p = 1d0 / tau_p 

                wmode = exp(-1d0/(rate_p*taulim))        !####new
                    
                a0    = -t * rate_p
                if(a0.lt.-200.0d0) a0 = -200.0d0
                e0    = 1.0d0-exp(a0)
                
                ffc   = cos((pi*p*nn)/dfloat(N)) * cos((pi*p*mm)/dfloat(N))
                ffc   = ffc/(dfloat(p)**(2*nue+1.0d0)+alpha*dfloat(p)**(4))

                ffc   = ffc * wmode                    !####new
    
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

       SqAve  = SqAve  + SqAve / weightsum
       SqtAve = SqtAve + SqtAve/ weightsum
       Daver  = Daver          / weightsum



!! checking
!!        write(6,*)
!!        write(6,*)"DrScal  = ",DrScal
!!        write(6,*)"etasolv = ",etasolv
!!        write(6,*)"nue     = ",nue
!!        write(6,*)"alpha   = ",alpha
!!        write(6,*)"temp    = ",temp
!!        write(6,*)"Nb      = ",Nb
!!        write(6,*)"Rb      = ",Rb
!!        write(6,*)"m_aver  = ",m_aver
!!        write(6,*)"wcut    = ",wcut
!!        write(6,*)"l       = ",l
!!        write(6,*)"taulim  = ",taulim
!!
!!        write(6,'(1x,7E14.6)')q,t,SqAve,SqtAve, SqtAve/SqAve, weightsum, daver   ! checking 
!!
       return
       end

