 FUNCTION th_dolgushev(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
! -------> nrouse+ring <--------
!
!
      use theory_description 
      implicit none 
      real    :: th_dolgushev
      character(len=8) :: thnam, parnam (*) 
      real    :: pa (*) 
      real    :: x , xh
      integer :: mbuf, nparx, ier, ini, npar, iadda
      integer, intent(inout) :: nopar       
      character(len=80), intent(inout) :: napar(mbuf) 
      real, intent(inout) :: params(mbuf) 
     
      double precision, parameter :: Pi = 4*atan(1d0)
 


       double precision ::    tau_e
       double precision ::    temp, qz, tau, eta, yz, SQ_rouse, a,b,xi
       double precision ::    R, W, l, Sq0, Sq, Sqt, wl4
       double precision ::    Rn
       double precision ::    a0
       double precision ::    diff12, diff3, dtube 
       real             ::    qget, tget
       integer          ::    N, nskip, pmax
       logical          ::    sofq0

 
       common/thiadd/iadda

!
! ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'dolgushe'
         nparx = 10
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,&
           ' exceeds current max. = ',i8)
           th_dolgushev = 0
           return
         endif
         npar = nparx

! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = "S(Q,t) [eq(52) of cit.] of a Rouse chain in a cylindrical potential"
       th_citation(idesc)     = "Maxim Dolgushev, Margarita Krutyeva, Macromol. Theory and Simul. 2012, 21, 565 "

!        --------------> set the number of parameters
         parnam(1)  = 'intensit'
         parnam(2)  = 'wl4'
         parnam(3)  = 'n_segmnt'
         parnam(4)  = 're      '   ! : re**2 = N * b**2
         parnam(5)  = 'temp    '
         parnam(6)  = 'd_orho  '
         parnam(7)  = 'd_para  '
         parnam(8)  = 'dtube   '   ! : lateral width ('tube')
         parnam(9)  = 'nskip   '   ! : speed up params
         parnam(10) = 'pmax    '   ! :  "
!
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "rouse rate common to both components" !//cr//parspace//&
        th_param_desc( 3,idesc) = "effective segment number in summation" !//cr//parspace//&
        th_param_desc( 4,idesc) = "end-to-end radius of component 1 if seen as linear Gaussian chain !" !//cr//parspace//&
        th_param_desc( 5,idesc) = "temperature" !//cr//parspace//&
        th_param_desc( 6,idesc) = "diffusion (cm**2/s) of component alpha 1,2 (if 0==> STD" !//cr//parspace//&
        th_param_desc( 7,idesc) = "diffusion (cm**2/s) of component 3 along tube (if 0==> STD)" !//cr//parspace//&

        th_param_desc( 8,idesc) = "tube diameter (lateral potential width)" !//cr//parspace//&
        th_param_desc( 9,idesc) = " step width n,m sums to speed up" !//cr//parspace//&
        th_param_desc(10,idesc) = " maximum mode number ... to speed up" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > scattering wavevector"
        th_file_param(  2,idesc) = "temp     > temperature"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(1,idesc)  = "l  effective segment length in A"
        th_out_param(2,idesc)  = "w  rouse rate in ns "
        th_out_param(3,idesc)  = "wl4  rouse rate wl**4 in A**4/ns"
        th_out_param(4,idesc)  = "diff12  diffusion component ortho in cm**2/s "
        th_out_param(5,idesc)  = "diff3   diffusion component parallel in cm**2/s"
        th_out_param(6,idesc)  = "'entanglement time' in ns "

! 
         th_dolgushev = 0
         return
       endif
!
       sofq0 = .false.
! ---- calculate theory here -----
       tau      = max(1d-6,x)
       a0       = pa(1)
       wl4      = abs(pa(2))  ! in A**4/ns
       N        = nint(pa(3))
       Rn       = abs(pa(4))  ! in A
       temp     = pa(5)
       diff12   = abs(pa(6))  ! in cm**2/sec  
       diff3    = abs(pa(7))
       dtube    = abs(pa(8))
       nskip    = max(1,nint(abs(pa(9 ))))
       pmax     = nint(abs(pa(10)))
       if(pmax == 0) pmax = N-1

!??       diff     = diff * 1d-9 / 1d-16  ! in A**2/ns

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
         sofq0 = .true.
       endif

       if(temp.eq.0.0d0) then
         tget = 300.0 
         call        parget('temp    ',tget,iadda,ier)
         temp = tget
       endif

       diff12 = diff12 * 1d7
       diff3  = diff3  * 1d7 

       th_dolgushev     = 0       
       call  N_dolgu(qz, tau, temp, diff12, diff3, wl4, N, nskip, pmax,  Rn, dtube, W, l, tau_e, Sq,Sqt)

       if(sofq0) then
          th_dolgushev = Sq * a0
!!? write(6,*)"sq ", Sq, a0
       else
          th_dolgushev = Sqt/Sq * a0
!!? write(6,*)"sqt ",sq, Sqt, a0
       endif

       call        parset('l       ',sngl(l),iadda,ier)      ! in ns A units
       call        parset('w       ',sngl(W),iadda,ier)      !     "
       call        parset('wl4     ',sngl(W*l**4),iadda,ier) !     "
       call        parset('diff12  ',sngl(diff12*1d-7),iadda,ier) ! in cm**2/s
       call        parset('diff3   ',sngl(diff3*1d-7),iadda,ier)
       call        parset('tau_e   ',sngl(tau_e),iadda,ier)
!

CONTAINS

       subroutine N_dolgu(q, t, temp, D_ortho0, D_para, wl4, N, nskip, pmax, R, dtube, W, l, tau_e, Sq,Sqt)
!      ====================================================================================================
!
! Dolgushev+Kruteva Model of a Rouse chain in a scylindrical confining potential (linear tube):
! with angular averaging.
! Input parameters:
!    q     ----> momentum transfer in A**-1
!    t     ----> time in nano-sec
!    temp  ----> temperature in K
!    D_ortho    ----> center of mass diffusion factor across tube,  if 0 <-- Rouse-expectation
!    D_para     ----> center of mass diffusion constant along tube, if 0 <-- Rouse-expectation
!    wl4   ----> friction coefficient in A**4/ns
!    N     ----> number of chain segments
!    R     ----> end-to-end distance of the polymer molecule
!    dtube ----> lateral tube dimension
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

       double precision, parameter       :: eps=1d-6

       double precision, intent(in)      :: q
       double precision, intent(in)      :: t
       double precision, intent(in)      :: temp
       double precision, intent(inout)   :: D_ortho0
       double precision, intent(inout)   :: D_para
       double precision, intent(in)      :: wl4
       integer         , intent(in)      :: N
       integer         , intent(in)      :: nskip
       integer         , intent(inout)   :: pmax
       double precision, intent(in)      :: R
       double precision, intent(in)      :: dtube
       double precision, intent(out)     :: w
       double precision, intent(out)     :: l
       double precision, intent(out)     :: tau_e
       double precision, intent(out)     :: Sq
       double precision, intent(out)     :: Sqt
 


       integer            ::  nn,mm,ip

       double precision   ::  xi 
       double precision   ::  tau_12p, tau_3p, kbt
       double precision   ::  X12, X3, e12, e3, D_ortho
       double precision   ::  A, A0
       double precision   ::  B, B0

       double precision   ::  p, cosn, cosm, dcnm, mcnm
!       integer iout
       
       if(N.le.0) then
         W  = 999
         Sq = 999
         Sqt= 999
         write(6,*)'Error Number of chain segments is <= 0!',N
         return
       endif

       if(pmax >= N) pmax = N-1
       if(pmax < 1 ) pmax = 1
!!? write(6,*)"th_dolgu:"
!!? write(6,*)"q       ", q       
!!? write(6,*)"t       ", t      
!!? write(6,*)"temp    ", temp   
!!? write(6,*)"D_ortho ", D_ortho
!!? write(6,*)"D_para  ", D_para 
!!? write(6,*)"wl4     ", wl4    
!!? write(6,*)"N       ", N      
!!? write(6,*)"R       ", R      
!!? write(6,*)"dtube   ", dtube    

! ---- determine the segment length l ----
       l = sqrt(R**2/N)       
       
! ---- and the Rousefactor ----
       kbt   = temp*kb            ! in Joule = kg*m**2/s**2
       kbt   = kbt * 100          ! in         kg*A**2/ns**2
       xi    = 3*kbt*l**2 / wl4
       W     = 3*kbt/(xi*(l**2))  ! in 1/ns
       tau_e = (3*dtube**2)/(8*W*l**2)

! ---- set the diffusion constant if input is zero --- !
       if(D_ortho0.eq.0.0d0) then
         D_ortho0 = (dtube**2) / (8*N*t) 
       endif
       D_ortho = D_ortho0 * (1d0-exp(-t/tau_e)) 

       if(D_para.eq.0.0d0) then
         D_para = kbt/(N*xi)
       endif

!!? write(6,*)"xi      ", xi       
!!? write(6,*)"W       ", W      
!!? write(6,*)"kbt     ", kbt   
!!? write(6,*)"D_ortho ", D_ortho
!!? write(6,*)"D_para  ", D_para 
!!? write(6,*)"tau_e   ", tau_e    



! ---- init sums ----
       Sq  = 0
       Sqt = 0

! ---- Do the sums -----
dn:    do nn = 1,N,nskip
dm:     do mm = 1,N,nskip
          A   = D_ortho * t
          A0  = 0
          B   = D_para  * t
          B0  = 0
   
dp:       do ip=1,pmax
            p = ip
            tau_3p   = 1d0 / (4*W * sin((pi*p)/dfloat(2*N))**2)
            tau_12p  = 1d0  /( 1d0/tau_3p + 1d0/tau_e )

            e12    = 1.0d0-exp(-t/tau_12p)
            e3     = 1.0d0-exp(-t/tau_3p)

            X3    = 1d0 / ( (24d0*N/l**2) * sin((pi*p)/dfloat(2*N))**2 ) 
            X12   = 1d0 / ( (1d0/X3)      + (16d0*N/dtube**2) )

            cosn  =  cos((pi*p*nn)/dfloat(N)) 
            cosm  =  cos((pi*p*mm)/dfloat(N)) 
            dcnm  = (cosn-cosm)**2
            mcnm  = 2*cosn*cosm

            A0    = A0 + 2 *  X12 *   dcnm 
            A     = A  + 2 *  X12 * ( dcnm  + mcnm * e12 )

            B0    = B0 + 2 *  X3  *   dcnm 
            B     = B  + 2 *  X3  * ( dcnm  + mcnm * e3  )

          enddo dp

!!? write(6,*)nn,mm,B-A,B0-A0

          if(B0-A0 > eps) then
            Sq  = Sq  + erf(q*sqrt(B0-A0))   * exp(-A0*q**2)  / (q*sqrt(B0-A0))
          else
            Sq =  Sq  + 2/sqrt(Pi) * exp(-A0*q**2) 
          endif
          Sqt = Sqt + erf(q*sqrt(B-A))     * exp(-A*q**2)   / (q*sqrt(B-A))

        enddo dm
       enddo dn

       Sq  = Sq  * sqrt(Pi)/( 2d0*N) * nskip**2
       Sqt = Sqt * sqrt(Pi)/( 2d0*N) * nskip**2



!!? write(6,'(1x,a,6E14.6)') 'q,t,Sq,Sqt, Sqt/Sq, w=', q,t,Sq,Sqt, Sqt/Sq, w 

       return
       end subroutine N_dolgu


      
       end function th_dolgushev


