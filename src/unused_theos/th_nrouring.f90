 FUNCTION th_nrouring(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
! -------> nrouse+ring <--------
!
!
      use theory_description 
      implicit none 
      real    :: th_nrouring
      character(len=8) :: thnam, parnam (*) 
      real    :: pa (*) 
      real    :: x , xh
      integer :: mbuf, nparx, ier, ini, npar, iadda
      integer, intent(inout) :: nopar       
      character(len=80), intent(inout) :: napar(mbuf) 
      real, intent(inout) :: params(mbuf) 
     
      double precision, parameter :: Pi = 4*atan(1d0)
 



       real*8    temp, qz, tau, eta, yz, SQ_rouse, a,b,xi
       real*8    R, W, l, Sq0, Sq, Sqt, wl4
       real*8    Rn, Rr
       real*8    a0, sumsqt, sumnorm, q_width, qzz
       real*8    diff, dr, drr, diffring, phirouse
       real*4    qget, tget
       integer   n, nmin, Nn, Nb, isring1, isring2
       logical   sqof0

 
       common/thiadd/iadda

!
! ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'nrouserr'
         nparx = 11
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,&
           ' exceeds current max. = ',i8)
           th_nrouring = 0
           return
         endif
         npar = nparx

! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = "Combination (additive) of two linear or ring Rouse chains f*S1(q,t)+(f-1)*S2(S(Q,t)"
       th_citation(idesc)     = ""

!        --------------> set the number of parameters
         parnam(1)  = 'intensit'
         parnam(2)  = 'wl4'
         parnam(3)  = 'n_segmnt'
         parnam(4)  = 're1      '   ! : re**2 = N * b**2
         parnam(5)  = 'temp    '
         parnam(6)  = 'diff1   '
         parnam(7)  = 'diff2   '
         parnam(8)  = 'phi12   '   ! : 0 normal Rouse, 1 fixed segment (N)
         parnam(9)  = 're2     '   ! : re (thought as linear strech of ring)
         parnam(10) = 'isring1 '   ! : re (thought as linear strech of ring)
         parnam(11) = 'isring2 '   ! : re (thought as linear strech of ring)

!
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "rouse rate common to both components" !//cr//parspace//&
        th_param_desc( 3,idesc) = "effective segment number in summation" !//cr//parspace//&
        th_param_desc( 4,idesc) = "end-to-end radius of component 1 if seen as linear Gaussian chain !" !//cr//parspace//&
        th_param_desc( 5,idesc) = "temperature" !//cr//parspace//&
        th_param_desc( 6,idesc) = "diffusion (cm**1/s) of component 1 (if=0, Rouse default)" !//cr//parspace//&
        th_param_desc( 7,idesc) = "diffusion (cm**1/s) of component 2 (if=0, Rouse default)" !//cr//parspace//&

        th_param_desc( 8,idesc) = "volume fraction f of component1 (internal limitation 0..1)" !//cr//parspace//&
        th_param_desc( 9,idesc) = "end-to-end radius of component 2 if seen as linear Gaussian chain !" !//cr//parspace//&
        th_param_desc(10,idesc) = "switch linear=0, ring > 0 for component 1" !//cr//parspace//&
        th_param_desc(11,idesc) = "switch linear=0, ring > 0 for component 2" !//cr//parspace//&

! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > scattering wavevector"
        th_file_param(  2,idesc) = "temp     > temperature"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(1,idesc)  = "l  effective segment length "
        th_out_param(2,idesc)  = "w  rouse rate "
        th_out_param(3,idesc)  = "wl4  rouse rate wl**4 "
        th_out_param(4,idesc)  = "diff  diffusion component1 "
        th_out_param(5,idesc)  = "diffring  diffusion component 2"

! 



         th_nrouring = 0
         return
       endif
!
       sqof0 = .false.
! ---- calculate theory here -----
       tau      = x
       a0       = pa(1)
       wl4      = abs(pa(2))  ! in A**4/ns
       N        = nint(pa(3))
       Rn       = abs(pa(4))  ! in A
       temp     = pa(5)
       diff     = abs(pa(6))  ! in cm**2/sec
       diffring = abs(pa(7))
       phirouse = min(abs(pa(8)),1.0)
       Rr       = abs(pa(9))
       isring1  = nint(pa(10))
       isring2  = nint(pa(11))


       diff     = diff * 1d-9 / 1d-16  ! in A**2/ns

!       if(epsilon.eq.0.0d0) epsilon = 1.0d-3

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
       if(temp.eq.0.0d0) then
         tget = 300.0 
         call        parget('temp    ',tget,iadda,ier)
         temp = tget
       endif


      th_nrouring     = 0
      sumsqt  = 0
      sumnorm = 0

      qzz  = qz        
      if(qzz.gt.0) then 

! --- include center of mass diffusion ---
!         a = fn * a0 * dexp(-qzz*qzz*diff*tau)
 
          dr = diff
          Nn = N
          if(isring1 == 0) then
            call NrouseY(qzz,tau,temp,dr,wl4,Nn,Rn, W, l, Sq,Sqt)
          else
            call Nrouse_ring(qzz,tau,dr,wl4,Nn,Rn, l, Sq,Sqt)
          endif
          sumsqt  = sumsqt  + phirouse * Sqt       
          sumnorm = sumnorm + phirouse * Sq

          Nb = N * min(3,nint((Rr/Rn)**2))
          drr = diffring
          if(isring2 == 0) then
            call NrouseY(qzz,tau,temp,drr,wl4,Nb,Rr, W, l, Sq,Sqt)
          else
            call Nrouse_ring(qzz,tau,drr,Wl4,Nb,Rr,l,Sq,Sqt)
          endif
          sumsqt  = sumsqt  + (1d0-phirouse) * Sqt       
          sumnorm = sumnorm + (1d0-phirouse) * Sq

        endif

       if(sumnorm.gt.0.0d0) then
          if(sqof0) then
            th_nrouring = sumsqt
          else
            th_nrouring = sumsqt/sumnorm
          endif
       else
          th_nrouring = 0
       endif 
       th_nrouring = th_nrouring*a0

       call        parset('l       ',sngl(l),iadda,ier)      ! in ns A units
       call        parset('w       ',sngl(W),iadda,ier)      !     "
       call        parset('wl4     ',sngl(W*l**4),iadda,ier) !     "
       dr        = dr /( 1d-9 / 1d-16 ) ! in cm**2/s
       call        parset('diff    ',sngl(dr),iadda,ier)
       drr       = drr /( 1d-9 / 1d-16 ) ! in cm**2/s
       call        parset('diffring',sngl(drr),iadda,ier)
!

CONTAINS

       subroutine NrouseY(q,t,temp,Dr,wl4,N,R, W, l, Sq,Sqt)
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
       integer N, nn,mm,ip

       double precision l, tau_p, kbt, Sq0, arg1, arg2
       double precision a0,e0, ff2, ffc,    arg10,arg20
       double precision aa1 , aa2
       double precision p, p0fix, pfac

!       integer iout
       
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

    
!       if(ifix.eq.0) then
         p0fix = 0
         pfac  = 1
!       else
!         p0fix = -0.5d0
!         pfac  = 0.5d0   ! this is wrong ? !
!         pfac=1         !! I think this has to be like this !!
!       endif

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

 !      if(iout().gt.1)write(6,'(1x,a,6E14.6)')
 !    *        'q,t,Sq,Sqt, Sqt/Sq, w=', q,t,Sq,Sqt, Sqt/Sq, w 

       return
       end subroutine NrouseY



       subroutine Nrouse_ring(q,t,Dr,Wl4,Nb,R,l,SqAve,SqtAve)
!      ====================================================
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
       double precision, intent(inout)  ::  dr
       double precision, intent(out)    ::  l, SqAve, SqtAve
       integer,          intent(in)     ::  Nb

       integer                          ::  nn,mm,p

       double precision :: W, rate_p, kbt, Sq0, arg1, arg2, xi
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
! ---- and the Rousefactor ----
       kbt = temp*kb            ! in Joule = kg*m**2/s**2
       kbt = kbt * 100          ! in         kg*A**2/ns**2
       xi  = 3*kbt*l**2 / wl4
       W   = 3*kbt/(xi*(l**2))  ! in 1/ns




! ---- set the diffusion constant if input is zero --- !
       if(Dr.eq.0.0d0) then
         Dr = kbt/(N*xi)
       endif


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
                 
!                wmode = amod(p/2)
!                wmode = wmode  / (1d0+exp(-(p/2-pmin)/pwidth))
                wmode = 1d0
                   
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
       SqtAve = Sqt * exp(-q*q*dr*t)

       return
       end subroutine Nrouse_ring







      
       end function th_nrouring


