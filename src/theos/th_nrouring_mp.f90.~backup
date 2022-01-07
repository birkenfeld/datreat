 FUNCTION th_nrouring_mp(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
! -------> nrouse+ring <--------
!
!
      use theory_description 
      implicit none 
      real    :: th_nrouring_mp
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
           th_nrouring_mp = 0
           return
         endif
         npar = nparx

! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = "Combination (additive) of two linear or ring Rouse chains f*S1(q,t)+(f-1)*S2(S(Q,t)"&
                                //cr//parspace//"ACCELERATED VERSION (OMP+ALG.)"

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



         th_nrouring_mp = 0
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


      th_nrouring_mp     = 0
      sumsqt  = 0
      sumnorm = 0

      qzz  = qz        
      if(qzz.gt.0) then 

! --- include center of mass diffusion ---
!         a = fn * a0 * dexp(-qzz*qzz*diff*tau)
 
          dr = diff
          Nn = N
          if(isring1 == 0) then
            call NrousePZ(qzz,tau,temp,dr,wl4,Nn,Rn, W, l, 1d0, dble(Nn), Sq,Sqt)
          else
! write(*,*)"A:",dr,wl4,Nn,Rn,l
            call Nrouse_ring(qzz,tau,temp,dr,wl4,Nn,Rn, l, Sq,Sqt)
          endif
          sumsqt  = sumsqt  + phirouse * Sqt       
          sumnorm = sumnorm + phirouse * Sq

          Nb = N * min(3,nint((Rr/Rn)**2))
          drr = diffring
          if(isring2 == 0) then
            call NrousePZ(qzz,tau,temp,drr,wl4,Nb,Rr, W, l, 1d0, dble(Nb), Sq,Sqt)
          else
! write(*,*)"B:",drr,wl4,Nb,Rr,l
            call Nrouse_ring(qzz,tau,temp,drr,Wl4,Nb,Rr,l,Sq,Sqt)
          endif
          sumsqt  = sumsqt  + (1d0-phirouse) * Sqt       
          sumnorm = sumnorm + (1d0-phirouse) * Sq

        endif

       if(sumnorm.gt.0.0d0) then
          if(sqof0) then
            th_nrouring_mp = sumsqt
          else
            th_nrouring_mp = sumsqt/sumnorm
          endif
       else
          th_nrouring_mp = 0
       endif 
       th_nrouring_mp = th_nrouring_mp*a0

       call        parset('l       ',sngl(l),iadda,ier)      ! in ns A units
       call        parset('w       ',sngl(W),iadda,ier)      !     "
       call        parset('wl4     ',sngl(W*l**4),iadda,ier) !     "
       dr        = dr /( 1d-9 / 1d-16 ) ! in cm**2/s
       call        parset('diff    ',sngl(dr),iadda,ier)
       drr       = drr /( 1d-9 / 1d-16 ) ! in cm**2/s
       call        parset('diffring',sngl(drr),iadda,ier)
!

CONTAINS


       subroutine Nrouse_ring(q,t,temp,dr,Wl4,Nb,R,l,SqAve,SqtAve)
!      ===================================================
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
       double precision, intent(in)     ::  temp
       double precision, intent(inout)  ::  dr
       double precision, intent(out)    ::  l, SqAve, SqtAve
       integer,          intent(in)     ::  Nb

       integer                          ::  nn,mm,p

       double precision :: W, rate_p, kbt, Sq0, xi
       double precision :: ff2
       double precision :: weight, Sqt, Sq, weightsum

       double precision :: ewfac2(Nb/2)
       double precision :: cosarray2(0:Nb,Nb/2)


 

       integer          :: N, ip, N2
       integer          :: iout, itrans
        
       if(Nb.le.0) then
         W  = 999
         Sq = 999
         Sqt= 999
         write(6,*)'Error Number of chain segments is <= 0!',Nb
         return
       endif
     
       N = Nb

! ---- determine the segment length l ----
       l = sqrt(R**2/Nb)       
       
! ---- and the Rousefactor ----

       W   = Wl4 / l**4  

       kbt = temp*kb            ! in Joule = kg*m**2/s**2
       kbt = kbt * 100          ! in         kg*A**2/ns**2
       xi  = 3*kbt*l**2 / wl4
!       W   = 3*kbt/(xi*(l**2))  ! in 1/ns

! ---- set the diffusion constant if input is zero --- !
       if(Dr.eq.0.0d0) then
         Dr = kbt/(N*xi)
       endif

! ---- init sums ----
       N2 = N/2

! write(*,*)"test2:",Nb,N,N2,R

!$OMP PARALLEL DO     
       do nn=0,N
        do ip=1,N2
 
         cosarray2(nn,ip) = cos((pi*2*ip*(nn))/dfloat(N))  &
                          /  (2*ip)**2
        enddo
       enddo
!$OMP END PARALLEL DO   

! write(*,*)"test2a:"

!$OMP PARALLEL DO    
       do ip=1,N/2
         ewfac2(ip) = 1.0d0-exp(-t * 2*W*(1-cos((pi*ip*2)/dfloat(N))) )
       enddo
!$OMP END PARALLEL DO    
! write(*,*)"test2b:"



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
                              + ff2*sum(cosarray2(abs(nn-mm),1:N2)*ewfac2(1:N2)) )
    
            enddo
           enddo
!$OMP END PARALLEL DO
    
           Sq  = Sq /Nb
           Sqt = Sqt/Nb

! write(*,*)"test2c:", Nb, Sq, Sqt
    

       SqAve  = Sq
       SqtAve = Sqt * exp(-q*q*dr*t)

       return
       end

       subroutine NrousePZ(q,t,temp,Dr,wl4,N,R, W, l,pmin,pmax, Sq,Sqt)
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
!    pmin  ----> minimum p
!    pmax  ----> maximum p
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

       double precision q,t,temp,Dr,xi,R, W,Sq,Sqt, wl4, pmin, pmax
       integer, intent(in) ::  N
       integer ::  nn,mm,ifix,ip

       double precision l, tau_p, kbt, Sq0, arg1, arg2
       double precision a0,e0, ff2, ffc,    arg10,arg20
       double precision aa1 , aa2
       double precision p, p0fix, pfac

       double precision :: cosarray(N,N), ewfac(N)

       integer :: ipmin, ipmax, i

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

!$OMP PARALLEL DO     
       do nn=1,N
        do ip=1,N
         cosarray(nn,ip) = cos((pi*ip*nn)/dfloat(N)) / ip 
        enddo
       enddo
!$OMP END PARALLEL DO   

!$OMP PARALLEL DO    
       do i=1,N
         ewfac(i) = (1d0-exp(-2*W*(1-cos((pi*i)/dfloat(N)))*t)) 
       enddo
!$OMP END PARALLEL DO    

       ipmin = max(1,nint(pmin))
       ipmax = min(N,nint(pmax))
!write(*,*)"TEST: ",ipmin, ipmax, N
! ---- init sums ----
       Sq0 = 0
       Sq  = 0
       Sqt = 0
       ff2  = -2*N*(l*q)**2/(3*pi**2)

! ---- Do the sums -----

!$OMP PARALLEL DO REDUCTION(+:Sq,Sqt)
       do nn = 1,N
        do mm = 1,N

          Sq  = Sq  + exp(-(q**2)*(       abs(nn-mm)*(l**2)/6.0d0))
          Sqt = Sqt + exp(-(q**2)*(Dr*t + abs(nn-mm)*(l**2)/6.0d0) + &
                ff2* sum(cosarray(nn,ipmin:ipmax) * cosarray(mm,ipmin:ipmax) *  ewfac(ipmin:ipmax) ))

        enddo
       enddo
!$OMP END PARALLEL DO

       Sq  = Sq /N
       Sqt = Sqt/N

!!       write(6,'(1x,a,6E14.6)')'q,t,Sq,Sqt, Sqt/Sq, w=', q,t,Sq,Sqt, Sqt/Sq, w 

       return
       end 
 

end function th_nrouring_mp
