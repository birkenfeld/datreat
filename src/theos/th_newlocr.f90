 FUNCTION th_newlocr(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  generalized local reptation expression along the lines of deGennes, 
!  but with finite summation of integrals and lenthscale, 
!  timescale and fluctuation ratio as parameters
!  Journal de Physique, 1981, 42 (5), pp.735-740. <10.1051/jphys:01981004205073500>
      use theory_description 
      implicit none 
      real    :: th_newlocr
      character(len=8) :: thnam, parnam (*) 
      real    :: pa (*) 
      real    :: x , xh
      integer :: mbuf, nparx, ier, ini, npar, iadda
      integer, intent(inout) :: nopar       
      character(len=80), intent(inout) :: napar(mbuf) 
      real, intent(inout) :: params(mbuf) 
     
      double precision, parameter :: Pi = 4*atan(1d0)
      integer                     :: actual_record_address
     
   
                                                                    
!     double precision :: sqdeb0, sqdebt
      double precision :: ampli     ! prefactor       
      double precision :: d         ! length scale coarse chain 
      double precision :: ne        ! entanglemen ne length        
      integer          :: n         ! n chain-length      
      double precision :: wl4       ! rouse rate    
      double precision :: b         ! effective fluctuation   
      double precision :: lseg      ! segment length 
      double precision :: wmod      ! tube W modifcation factor
      double precision :: lsegmod   ! blob lmodification factor
      integer          :: inc       ! if=1 inc für locrep ?!

      integer          :: Z         ! n/ne

!
 
! the recin parameter representation 
     double precision :: q          ! q-value                                                                         

     double precision :: W0, temp, Re, lseg_nr
     integer          :: nr
     double precision :: sqt, sqt0
     double precision :: sqdebt, sqdeb0
     double precision, parameter   :: teps = 1d-6
 
     double precision   :: t
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'newlocr'
       nparx =        10
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_newlocr = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " generalized local reptation expression along the lines of deGennes, " //cr//parspace//&
                                " but lumped to finite sums," //cr//parspace//&
                                " . "
       th_citation(idesc)     = " Journal de Physique, 1981, 42 (5), pp.735-740. <10.1051/jphys:01981004205073500>"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor       
        parnam ( 2) = 'd       '  ! length scale coarse chain 
        parnam ( 3) = 'ne      '  ! entanglemen ne length        
        parnam ( 4) = 'n       '  ! n chain-length      
        parnam ( 5) = 'wl4     '  ! rouse rate    
        parnam ( 6) = 'b       '  ! effective fluctuation   
        parnam ( 7) = 'lseg    '  ! segment length 
        parnam ( 8) = 'wrmod   '  ! tube W modification factor 
        parnam ( 9) = 'lsegmod '  ! tube W modification factor 
        parnam (10) = 'inc     '  ! 
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "coarse tube d" !//cr//parspace//&
        th_param_desc( 3,idesc) = "entanglement ne" !//cr//parspace//&
        th_param_desc( 4,idesc) = "total n" !//cr//parspace//&
        th_param_desc( 5,idesc) = "rouse rate" !//cr//parspace//&
        th_param_desc( 6,idesc) = "fluctuation intensity of locrep" !//cr//parspace//&
        th_param_desc( 7,idesc) = "segment length" !//cr//parspace//&
        th_param_desc( 8,idesc) = "tube W modification factor " !//cr//parspace//&
        th_param_desc( 9,idesc) = "lseg Rouseblob modification factor " !//cr//parspace//&
        th_param_desc(10,idesc) = "if =1 pseudo inc locrep " !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
 ! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q-value"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)   = " "
!
                                             
! >>>>> describe parameters >>>>>>> 
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > momentum transfer"
        th_file_param(  2,idesc) = "temp     > temperature"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(1,idesc)  = "Re end to end of blob "



 
        th_newlocr = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----



      ampli    =        pa( 1)   
      d        =   abs( pa( 2) )
      ne       =   abs( pa( 3) ) 
      n        =  nint( pa( 4) )
      wl4      =   abs( pa( 5))
      b        =   abs( pa( 6))
      lseg     =   abs( pa( 7))
      wmod     =   abs( pa( 8))
      lsegmod  =   abs( pa( 9))
      inc      =  nint( pa(10))

      if(wmod    == 0d0) wmod    = 1d0
      if(lsegmod == 0d0) lsegmod = 1d0

     if(ne > 500d0) ne = 500

! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: q-value
      xh =     0.1d0
      call parget('q       ',xh,iadda,ier)
      q        = xh

      xh =     400d0
      call parget('temp    ',xh,iadda,ier)
      temp     = xh



      t        = x
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 

     Re       = sqrt(ne)*lseg    !! modification to enable Ne fittting (fractional ne)
     Nr       = nint(ne)
     lseg_nr  = Re/sqrt(dble(Nr)) * lsegmod

!  the Rosue-blob
     call NrousePlain(q,t,wl4,lseg_nr,Nr, Sqt0,Sqt)

     W0        = wl4/lseg**4 * wmod
     Z         = nint(n/ne)
     Sqdeb0 = nndeblr(q, teps, d, b ,Z, nint(Ne), lseg, W0 ) 
     Sqdebt = nndeblr(q, t,    d, b ,Z, nint(Ne), lseg, W0 )

     th_newlocr = ampli * sqt/sqt0  * Sqdebt/Sqdeb0
!    ===============================================

     call parset('re      ',sngl(re),iadda)
     call parset('lseg_nr ',sngl(lseg_nr),iadda)
     call parset('wtube   ',sngl(W0),iadda)

                   


CONTAINS

  function wlocrep_diffusion_of_n(t,n,d,lseg,Ne,W) result(val)
    implicit none
    double precision, intent(in) :: t    !! time
    integer         , intent(in) :: n    !! distance in number os steps of caorse brained tube
    double precision, intent(in) :: d    !! tube step-length
    double precision, intent(in) :: lseg !! chain segemnt length (single bond)
    integer         , intent(in) :: Ne   !! number of segments per blob
    double precision, intent(in) :: W    !! Rouse rate

    double precision             :: val

    double precision :: delta_t

    if(inc==1 .and. n .ne. 0) then; val=0; return; endif

    delta_t = t *  W*lseg**2 / Ne

    val = 0.5d0*(  erf( d*(2*n+1) / (4d0*sqrt(delta_t) )) &
                  -erf( d*(2*n-1) / (4d0*sqrt(delta_t) )) )

  end function wlocrep_diffusion_of_n 



function nndeblr(q, t, d, b ,Z, Ne, lseg, W ) result(val)
   implicit none
   double precision ,intent(in) :: q           ! Q-value
   double precision ,intent(in) :: t           ! time
   double precision ,intent(in) :: d           ! effective (tube) segment length
   double precision ,intent(in) :: b           ! relative fluctuation amplitued (1/3 ??)
   integer          ,intent(in) :: Z           ! number of steps (coarse grainded tube)
   integer          ,intent(in) :: Ne          ! number of segments per blob
   double precision ,intent(in) :: lseg        ! chain segment length
   double precision ,intent(in) :: W           ! rouse rate
   double precision             :: val

   integer          :: i, j
   double precision :: eterms(0:n-1)
   double precision :: qd2, mu

   qd2 = (1d0/6d0) * (q*d)**2
   
!  Rg  = sqrt((l*l*dble(n)**mu) / 6d0) )
!$OMP PARALLEL DO
   do i=0,Z-1
     eterms(i) = exp(-qd2*dble(i)**mu) * ( 1d0 + b * wlocrep_diffusion_of_n(t,i,d,lseg,Ne,W)  )
   enddo
!$OMP END PARALLEL DO
   val = 0
!$OMP PARALLEL DO REDUCTION(+:val)
   do i=1,Z
    do j=1,Z
      val = val + eterms(abs(i-j))
    enddo
   enddo
!$OMP END PARALLEL DO

   val = val / (Z*Z)


end function nndeblr




       subroutine NrousePlain(q,t,wl4,l,N, Sq,Sqt)
!      ==================================================
!
! Rouse expression for a chain of finite length:
! Input parameters:
!    q     ----> momentum transfer in A**-1
!    t     ----> time in nano-sec
!    wl4   ----> friction coefficient in A**4/ns
!    l     ----> segment length
!    N     ----> number of chain segments
! ------------------------------------------------------------
!
       implicit none

       double precision kb, pi
       parameter(kb=1.380662d-23)
       parameter(pi=3.141592654d0)

       double precision, intent(in) :: q
       double precision, intent(in) :: t
       double precision, intent(in) :: wl4
       double precision, intent(in) :: l
       integer,          intent(in) :: N 
       double precision, intent(out):: Sq,Sqt

       integer nn,mm,ip

       double precision :: Sq0
       double precision :: ff2
       double precision :: W

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
      
! ---- and the Rousefactor ----
       W   = wl4/l**4  ! in 1/ns


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

       ipmin = 1
       ipmax = N

! ---- init sums ----
       Sq0 = 0
       Sq  = 0
       Sqt = 0
       ff2  = -2*N*(l*q)**2/(3*pi**2)

! ---- Do the sums -----

!$OMP PARALLEL DO REDUCTION(+:Sq,Sqt)
       do nn = 1,N
        do mm = 1,N

          Sq  = Sq  + exp(-(q**2)*(  abs(nn-mm)*(l**2)/6.0d0))
          Sqt = Sqt + exp(-(q**2)*(  abs(nn-mm)*(l**2)/6.0d0) + &
                ff2* sum(cosarray(nn,ipmin:ipmax) * cosarray(mm,ipmin:ipmax) *  ewfac(ipmin:ipmax) ))

        enddo
       enddo
!$OMP END PARALLEL DO

       Sq  = Sq /N
       Sqt = Sqt/N
       return
       end subroutine NrousePlain
 






 end function th_newlocr

 
