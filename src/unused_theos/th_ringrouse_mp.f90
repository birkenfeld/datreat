      FUNCTION th_ringrouse_mp (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ================================================================================               
!                                                                       
! -------> nrouse_sh4 <--------
! selfhealing multimodal distribution mix of Rouse        
!                                                                       
!                
      implicit none  
      double precision, parameter :: Pi                  =     3.141592653589793238462643d0    
                                                 
      CHARACTER(len=8), intent(inout) :: thnam, parnam (*) 
      real, intent(in)                :: x,  pa (*)
      integer, intent(in)             :: ini, mbuf
      integer, intent(inout)          :: npar
      integer, intent(inout)          :: nopar            ! Anzahl der Parameter data
      character(len=80), intent(inout):: napar(mbuf)      ! name des parameters n
      real, intent(inout)             :: params(mbuf)     ! value des parameters n
      real                            :: th_ringrouse_mp

      integer          :: ier, nparx

      double precision :: temp, tau, q, wl4, l 
      double precision :: a0 , qz, R, dr
      double precision :: sq, sqt
      double precision :: diff, r02, tc, a_cross, nu_subdiff, rr
      double precision :: ampmod(1000)    ! quick an dirty upper limit for mode modifiers
      double precision :: pmin, pwidth
      
      REAL             :: qget, tget 
      integer          :: n
      logical          :: sqof0

      integer          :: iadda 
      common/thiadd/iadda

 
 ! ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'ringrous'
 !##!        thnam = 'nrouse_h'
         nparx = 20
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,' exceeds current max. = ',i8)
           th_ringrouse_mp = 0
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
         parnam(19)= 'pmin    '   ! additional mode spectrum multiplication factor (Fermi) a(p) = 1 /(1+exp(-(p-pmin)/pwidth))
         parnam(20)= 'pwidth  '   ! additional mode spectrum width factor

         th_ringrouse_mp = 0
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
       ampmod     = 1d0
       ampmod(1:9)= abs(pa(10:18))
       pmin       = pa(19)
       pwidth     = abs(pa(20))
       
       if(pwidth == 0d0) then
          pwidth = 10
          pmin   = -100
       endif

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
             
       call Nrouse_ring(qz,tau,Wl4,N,R,ampmod,pmin,pwidth,l, Sq,Sqt)
 
       if(sqof0) then
          th_ringrouse_mp =  a0 * Sqt
       else
          th_ringrouse_mp =  a0 * Sqt/Sq
       endif

!!! Multiply the subdiff expression !!!
       rr = ( (r02*(tau/tc)**nu_subdiff)**a_cross + (6*diff*tau)**a_cross)**(1d0/a_cross)
       th_ringrouse_mp = th_ringrouse_mp * exp( -qz*qz*rr/6d0 )
         

 


       call        parset('l       ',sngl(l),iadda,ier)      ! in ns A units
       call        parset('wl4     ',sngl(wl4),iadda,ier)
 
CONTAINS


       subroutine Nrouse_ring(q,t,Wl4,Nb,R,ampmod,pmin,pwidth,l,SqAve,SqtAve)
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
       double precision, intent(in)     ::  ampmod(1000)
       double precision, intent(in)     ::  pmin         ! pmin in ring mode counting (i.e p--> p/2)
       double precision, intent(in)     ::  pwidth 
       double precision, intent(out)    ::  l, SqAve, SqtAve
       integer,          intent(in)     ::  Nb

       integer                          ::  nn,mm,p

       double precision :: W, rate_p, kbt, Sq0, arg1, arg2, Dr
       double precision :: a0,e0, ff2, ffc,    arg10,arg20
       double precision :: aa1 , aa2, weight, Sqt, Sq, weightsum, wmode

       
       double precision :: cosarray(0:Nb,Nb/2), ewfac(Nb/2)


 

       integer          :: N, ip, N2

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
       N2 = N/2

!$OMP PARALLEL DO     
       do nn=0,N
        do ip=1,N2
 
         cosarray(nn,ip) = cos((pi*2*ip*(nn))/dfloat(N))  &
                          *  ampmod(ip) / (1d0+exp(-(ip-pmin)/pwidth)) &
                          /  (2*ip)**2
        enddo
       enddo
!$OMP END PARALLEL DO   


!$OMP PARALLEL DO    
       do ip=1,N/2
         ewfac(ip) = 1.0d0-exp(-t * 2*W*(1-cos((pi*ip*2)/dfloat(N))) )
       enddo
!$OMP END PARALLEL DO    



       Sq0 = 0
       SqAve  = 0
       SqtAve = 0
       weightsum = 0


           Sq0 = 0
           Sq  = 0
           Sqt = 0

           N      =  Nb 
           ff2  = -2*N*(l*q)**2/(3*pi**2)

! ---- Do the sums -----
!$OMP PARALLEL DO REDUCTION(+:Sq,Sqt)
           do nn = 1,N
            do mm = 1,N

    
              Sq  = Sq  + exp(-(q**2)*( abs(nn-mm)*(N-abs(nn-mm))*(l**2)/(6.0d0*N)))
              Sqt = Sqt + exp(-(q**2)*( abs(nn-mm)*(N-abs(nn-mm))*(l**2)/(6.0d0*N)) &
                              + ff2*sum(cosarray(abs(nn-mm),1:N2)*ewfac(1:N2)) )
    
            enddo
           enddo
!$OMP END PARALLEL DO
    
           Sq  = Sq /Nb
           Sqt = Sqt/Nb
    

       SqAve  = Sq
       SqtAve = Sqt

       return
       end

end function th_ringrouse_mp
