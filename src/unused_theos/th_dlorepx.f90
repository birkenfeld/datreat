 FUNCTION th_dlocrepx(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  discrete local reptation Ansatz
!  ??
      use theory_description 
      implicit none 
      real    :: th_dlocrepx
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
     double precision :: nn         ! totla number of segements/chain                             
     double precision :: ne         ! entanglement N, integer: do not fit                                             
     double precision :: l          ! segment length                                                                  
     double precision :: b          ! b-factor local reptation                                                        
     double precision :: w          ! rouse rate                                                                      
     double precision :: fw1        ! modification factor w-locrecp                                                   
     double precision :: fw2        ! modification factor Rouse blob                                                  
     double precision :: fa         ! modification factor step length                                                 
! the recin parameter representation 
     double precision :: q          ! q-value                                                                         
! the reout parameter representation 
     double precision :: re         ! Ree Rouse                                                              
 
     double precision :: th
 
   double precision :: dne
   double precision :: t
   double precision :: sqtloc, sqtloc0, sqtrous, sqtrous0
   double precision :: sqtloc_w, sqtloc0_w, sqtrous_w, sqtrous0_w
   double precision :: weight, wn
   integer          :: nz, nr
   integer          :: mode
   integer          :: iw
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'dlocrepx'
       nparx =       10
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_dlocrepx = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " discrete local reptation Ansatz:"//cr//parspace//&
                                "    scales     Wrous=W*fw2"//cr//parspace//&
                                "               Wloc =(W/ne^2)*fw1"//cr//parspace//&
                                "               B    =b*(1/3)"//cr//parspace//&
                                "               a    =l*sqrt(ne)*fa"//cr//parspace//&
                                "    Version with +-5% ne averaging (--> ne can be fitted)"
                                
       th_citation(idesc)     = " ??"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'nn      '  ! number of total segements/chain                              
        parnam ( 3) = 'ne      '  ! entanglement N, integer: do not fit                                             
        parnam ( 4) = 'l       '  ! segment length                                                                  
        parnam ( 5) = 'b       '  ! b-factor local reptation                                                        
        parnam ( 6) = 'w       '  ! rouse rate                                                                      
        parnam ( 7) = 'fw1     '  ! modification factor w-locrecp                                                   
        parnam ( 8) = 'fw2     '  ! modification factor Rouse blob                                                  
        parnam ( 9) = 'fa      '  ! modification factor step length                                                 
        parnam (10) = 'mode    '  ! show what ..                                                 
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "total number of segments pre chain" !//cr//parspace//&
        th_param_desc( 3,idesc) = "entanglement Ne " !//cr//parspace//&
        th_param_desc( 4,idesc) = "segment length" !//cr//parspace//&
        th_param_desc( 5,idesc) = "b-factor local reptation" !//cr//parspace//&
        th_param_desc( 6,idesc) = "rouse rate" !//cr//parspace//&
        th_param_desc( 7,idesc) = "modification factor w-locrecp" !//cr//parspace//&
        th_param_desc( 8,idesc) = "modification factor Rouse blob" !//cr//parspace//&
        th_param_desc( 9,idesc) = "modification factor step length" !//cr//parspace//&
        th_param_desc(10,idesc) = "selector: 0=full, 1=locrep only, 2=rouse only" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q-value"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "re       > Ree rouse"
! 
        th_dlocrepx = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      nn       =  abs(pa( 2))
      ne       =  abs(pa( 3))
      l        =  abs(pa( 4))
      b        =  abs(pa( 5))
      w        =  abs(pa( 6))
      fw1      =  abs(pa( 7))
      fw2      =  abs(pa( 8))
      fa       =  abs(pa( 9))
      mode     = nint(pa(10))
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: q-value
      xh =      0.1
      call parget('q       ',xh,iadda,ier)
      q        = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     t   = x
     nz  = nint(nn/ne)
     nr  = nint(ne)
     dne = 0.05d0 * ne

     sqtrous0 = 0
     sqtrous  = 0
     sqtloc0  = 0
     sqtloc   = 0
     wn       = 0


     do iw = nr-nint(dne), nr+nint(dne)
        nz  = nint(nn/dble(nr))
        weight = exp(- (((iw-ne)/dne)*3)**2 ) 

        call NrousePb(q,t, nr, W*fw2, l,1,nr, sqtrous0_w,sqtrous_w)
        call Nlocrep (q,t, nz, W*fw1/ne**2, l*fa*sqrt(ne), B*(1d0/3d0),   Sqtloc0_w,  Sqtloc_w)
!!!!                                
        sqtrous0 = sqtrous0 + weight *  sqtrous0_w 
        sqtrous  = sqtrous  + weight *  sqtrous_w  
        sqtloc0  = sqtloc0  + weight *  sqtloc0_w  
        sqtloc   = sqtloc   + weight *  sqtloc_w  
        wn       = wn + weight 

     enddo


     

!!--------------------------------^ ==> check ob hier ne oder ne**2 ??!!

     select case(mode)
     case(0)
       th  = ampli * sqtrous*Sqtloc/(sqtrous0*sqtloc0)
     case(1)
       th  = ampli * Sqtloc/(sqtloc0)
     case(2)
       th  = ampli * sqtrous/(sqtrous0)
     case(3) 
       th = Sqtloc / wn
     case(4)
       th = Sqtrous / wn
     case default
       th  = ampli * sqtrous*Sqtloc/(sqtrous0*sqtloc0)
     end select

! write(*,*)"T1: ", ampli, q, t, nr, w, fw2, l, sqtrous0,sqtrous
! write(*,*)"T2: ", nz, fw1, Sqtloc0,  Sqtloc, th


     th_dlocrepx = th
 
! ---- writing computed parameters to the record >>>  
      re = l * sqrt(ne)
      call parset('re      ',sngl(re),iadda)
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 


       subroutine NrousePb(q,t, N, W, l,pmin,pmax, Sq,Sqt)
!      ==================================================
!
! Rouse expression for a chain of finite length (without com diffusion):
!
       implicit none
       double precision, intent(in) :: q              ! momentum transfer
       double precision, intent(in) :: t              ! time
       integer         , intent(in) :: n              ! number of segments
       double precision, intent(in) :: w              ! Rouse rate
       double precision, intent(in) :: l              ! segment length
       integer         , intent(in) :: pmin           ! lowest mode nr
       integer         , intent(in) :: pmax           ! highest mode nr
       double precision, intent(out):: Sq             ! Sq(t=0)
       double precision, intent(out):: Sqt            ! Sq(t)




       double precision kb, pi
       parameter(pi=3.141592654d0)

       integer          :: nn,mm,ip

       double precision :: tau_p,  Sq0, arg1, arg2
       double precision :: a0,e0, ff2, ffc,    arg10,arg20
       double precision :: aa1 , aa2
       double precision :: p, p0fix, pfac

       double precision :: cosarray(N,N), ewfac(N)

       integer :: ipmin, ipmax, i

       if(N.le.0) then
         Sq = 999
         Sqt= 999
         write(6,*)'Error Number of chain segments is <= 0!',N
         return
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

       ipmin = max(1,pmin)
       ipmax = min(N,pmax)

! ---- init sums ----
       Sq  = 0
       Sqt = 0
       ff2  = -2*N*(l*q)**2/(3*pi**2)

! ---- Do the sums -----

!$OMP PARALLEL DO REDUCTION(+:Sq,Sqt)
       do nn = 1,N
        do mm = 1,N

          Sq  = Sq  + exp(-(q**2)*(       abs(nn-mm)*(l**2)/6.0d0))
          Sqt = Sqt + exp(-(q**2)*(       abs(nn-mm)*(l**2)/6.0d0) + &
                ff2* sum(cosarray(nn,ipmin:ipmax) * cosarray(mm,ipmin:ipmax) *  ewfac(ipmin:ipmax) ))

        enddo
       enddo
!$OMP END PARALLEL DO

       Sq  = Sq /N
       Sqt = Sqt/N

       end subroutine NrousePb




       subroutine Nlocrep(q,t, N, W, l, B, Sq,Sqt)
!      ===============================================
!
! Local repatation as finite sum
!
       implicit none
       double precision, intent(in) :: q              ! momentum transfer
       double precision, intent(in) :: t              ! time
       integer         , intent(in) :: n              ! number of segments (=Z=steps)
!       integer         , intent(in) :: ne             ! number of segments (ne)
       double precision, intent(in) :: w              ! Rouse rate
       double precision, intent(in) :: l              ! segment length
       double precision, intent(in) :: B              ! b-factor
       double precision, intent(out):: Sq             ! Sq(t=0)
       double precision, intent(out):: Sqt            ! Sq(t)




       double precision kb, pi
       parameter(pi=3.141592654d0)

       integer          ::  nn,mm

       double precision :: tau_p, arg1, arg2
       double precision :: a0,e0, ff2, ffc,    arg10,arg20
       double precision :: aa1 , aa2
       double precision :: p, p0fix, pfac

       double precision :: aq, sig, s1, s2



       double precision            :: distar(0:N), distar0(0:N), ss
       double precision, parameter :: Delta = 0.5d0

       integer :: ipmin, ipmax, i


       aq  = l*l * q*q 
!!??       sig = 4*w*l**2
!!TBD: hier l oder a=l*sqrt(ne) ????
       distar(0)  =  2*Delta + B*erf(Delta/(2*sqrt(W*t)))
       distar0(0) =  2*Delta + B

!$OMP PARALLEL DO PRIVATE(ss)
       do nn=1,N
!         ss = nn
!        distar(nn)  = (1d0 + B*exp(-ss**2/(sig*t))/(sqrt(Pi*sig*t))) * exp(-aq*ss/l)
         distar(nn)  = (2*Delta+(1d0/2)*B*erf((Delta+nn)/(2*sqrt(W*t))) &
                              +(1/2d0)*B*erf((Delta-nn)/(2*sqrt(W*t)))) * exp(-aq*nn/(6d0))

         distar0(nn) = exp(-aq*nn/(6d0))



!! Alternate full integral (numerically problematic erf diff at large aq, get fix.. !
!!        distar(nn)=exp(-2*aq*nn)*(B*exp(aq*(W*aq*t+2*nn))*erf((-2*W*aq*t+Delta-nn)*W**(-1d0/2)*t**(-1d0/2)/2)*aq+ &
!!                                  B*erf((2*W*aq*t+Delta+nn)*W**(-1d0/2)*t**(-1d0/2)/2)*exp(aq*(W*aq*t+2*nn))*aq+  &
!!                                  2*l*(exp(aq*(Delta+nn))-exp(aq*(-Delta+nn))))/aq/l/2



       enddo
!$OMP END PARALLEL DO

! ---- init sums ----
       Sq  = 0
       Sqt = 0
! ---- Do the sums -----

!$OMP PARALLEL DO REDUCTION(+:Sq,Sqt)
       do nn = 1,N
        do mm = 1,N
          Sqt   = Sqt   + distar(abs(nn-mm))
          Sq    = Sq    + distar0(abs(nn-mm))
        enddo
       enddo
!$OMP END PARALLEL DO

       Sq  = Sq /N
       Sqt = Sqt/N

       end subroutine Nlocrep



 end function th_dlocrepx
