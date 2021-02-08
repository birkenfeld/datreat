 FUNCTION th_dlocrep4(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  discrete local reptation Ansatz
!  ??
      use theory_description 
      implicit none 
      real    :: th_dlocrep4
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
     double precision :: nn         ! total number of segments                             
     double precision :: ne         ! entanglement N, integer: do not fit                                             
     double precision :: l          ! segment length                                                                  
     double precision :: b          ! b-factor local reptation                                                        
     double precision :: w          ! rouse rate                                                                      
     double precision :: fw1        ! modification factor w-locrecp                                                   
     double precision :: fw2        ! modification factor Rouse blob                                                  
     double precision :: fa         ! modification factor step length   

     double precision :: eps_yy                                               

  
     double precision :: alpha0     ! prefactor f alpha function
     double precision :: talphamax  ! max tau in log
     double precision :: talphawd   ! width in taulog
     double precision :: aoffset    ! offset





! the recin parameter representation 
     double precision :: q          ! q-value                                                                         
! the reout parameter representation 
     double precision :: re         ! Ree Rouse                                                              
 
     double precision :: th
 




   double precision :: t
   double precision :: sqtloc, sqtloc0, sqtrous, sqtrous0, wout
   integer          :: nz, nr
   integer          :: mode
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'dlocrep4'
       nparx =       15
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_dlocrep4 = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " discrete local reptation Ansatz"
       th_citation(idesc)     = " ??"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'n       '  ! total number of segments                               
        parnam ( 3) = 'ne      '  ! entanglement N, integer: do not fit                                             
        parnam ( 4) = 'l       '  ! segment length                                                                  
        parnam ( 5) = 'b       '  ! b-factor local reptation                                                        
        parnam ( 6) = 'w       '  ! rouse rate                                                                      
        parnam ( 7) = 'fw1     '  ! modification factor w-locrecp                                                   
        parnam ( 8) = 'fw2     '  ! modification factor Rouse blob                                                  
        parnam ( 9) = 'fa      '  ! modification factor step length                                                 
        parnam (10) = 'mode    '  ! show what ..                                                 
        parnam (11) = 'eps_yy  '  ! show what ..   
        parnam (12) = 'alpha0  '  !                                                     
        parnam (13) = 'talpmax '  !                                                     
        parnam (14) = 'talpwd  '  !                                                    
        parnam (15) = 'alpoffs '  !                                                
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "total number of segments: do not fit!" !//cr//parspace//&
        th_param_desc( 3,idesc) = "entanglement N, integer: do not fit" !//cr//parspace//&
        th_param_desc( 4,idesc) = "segment length" !//cr//parspace//&
        th_param_desc( 5,idesc) = "b-factor local reptation" !//cr//parspace//&
        th_param_desc( 6,idesc) = "rouse rate" !//cr//parspace//&
        th_param_desc( 7,idesc) = "modification factor w-locrecp" !//cr//parspace//&
        th_param_desc( 8,idesc) = "modification factor Rouse blob" !//cr//parspace//&
        th_param_desc( 9,idesc) = "modification factor step length" !//cr//parspace//&
        th_param_desc(10,idesc) = "selector: 0 use Snm(2D), 1 use Snm(1D), 2 use Snm(0D) 3 Snm(2D)chk" !//cr//parspace//&
        th_param_desc(11,idesc) = "epsilon default 1/4 (for mode 2 only)" !//cr//parspace//&
        th_param_desc(12,idesc) = "prefactor of non Gaussianity alpha(t) expression" !//cr//parspace//&
        th_param_desc(13,idesc) = "t of max in log distribution alpha(t)" !//cr//parspace//&
        th_param_desc(14,idesc) = "width of Gaussion max. in log t of alpha(t)" !//cr//parspace//&
        th_param_desc(15,idesc) = "offset in alpha(t)" !//cr//parspace//&

! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q-value"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "re       > Ree rouse"
! 
        th_dlocrep4 = 0.0
 
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
      eps_yy   =  abs(pa(11))

      alpha0   =      pa(12)
      talphamax=  abs(pa(13))
      talphawd =  abs(pa(14))
      aoffset  =      pa(15)
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
!     call NrousePb  (q,t, nr, W*fw2,               l,1,nr, sqtrous0,sqtrous)
     re = l * sqrt(ne)
     call nrousalpha(q,t,509d0,1d-20,W*l**4*fw2, Nr, Re, Wout, l, 1d0, ne, sqtrous0,sqtrous)
     call Nlocrep (q,t, nz, W*fw1/ne**2, l*fa*sqrt(ne), B,   Sqtloc0,  Sqtloc)
!!--------------------------------^ ==> check ob hier ne oder ne**2 ??!!

     th  = ampli * sqtrous*Sqtloc/(sqtrous0*sqtloc0)
! write(*,*)"T1: ", ampli, q, t, nr, w, fw2, l, sqtrous0,sqtrous
! write(*,*)"T2: ", nz, fw1, Sqtloc0,  Sqtloc, th


     th_dlocrep4 = th
 
! ---- writing computed parameters to the record >>>  
      call parset('re      ',sngl(re),iadda,ier)
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 
 

       subroutine nrousalpha(q,t,temp,Dr,wl4,N,R, W, l,pmin,pmax, Sq,Sqt)
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
       integer N, nn,mm,ifix,ip

       double precision l, tau_p, kbt, Sq0, arg1, arg2
       double precision a0,e0, ff2, ffc,    arg10,arg20
       double precision aa1 , aa2
       double precision p, p0fix, pfac

       double precision :: cosarray(N,N), ewfac(N)
       double precision :: rmm, fqq, fqq0

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

! ---- init sums ----
       Sq0 = 0
       Sq  = 0
       Sqt = 0
       ff2  = -2*N*(l*q)**2/(3*pi**2)

! ---- Do the sums -----

       rmm = 0
!$OMP PARALLEL DO REDUCTION(+:rmm)
       do mm = 1,N
             rmm = rmm + 4d0*N*l**2/(pi**2) * &
                   sum(cosarray(mm,ipmin:ipmax) * cosarray(mm,ipmin:ipmax) *  ewfac(ipmin:ipmax) )
       enddo
!$OMP END PARALLEL DO
       rmm = rmm/N

       fqq  = 1d0 -q**2 * rmm/12d0 * alpha(t)       !! see Guenza PHYSICAL REVIEW E 89, 052603 (2014) Eq(4)
       fqq0 = 1d0 -q**2 * rmm/12d0 * alpha(0d0)     !! see Guenza PHYSICAL REVIEW E 89, 052603 (2014) Eq(4)


!$OMP PARALLEL DO REDUCTION(+:Sq,Sqt)
       do nn = 1,N
        do mm = 1,N

          Sq  = Sq  + exp(- fqq0*(q**2)*(       abs(nn-mm)*(l**2)/6.0d0))
          Sqt = Sqt + exp(- fqq* (q**2)*(Dr*t + abs(nn-mm)*(l**2)/6.0d0) + &
                fqq*ff2* sum(cosarray(nn,ipmin:ipmax) * cosarray(mm,ipmin:ipmax) *  ewfac(ipmin:ipmax) ))

        enddo
       enddo
!$OMP END PARALLEL DO

       Sq  = Sq /N
       Sqt = Sqt/N

!!       write(6,'(1x,a,6E14.6)')'q,t,Sq,Sqt, Sqt/Sq, w=', q,t,Sq,Sqt, Sqt/Sq, w 

       return
       end
 

function alpha(t) result(a) !! see Guenza PHYSICAL REVIEW E 89, 052603 (2014) Eq(4)
   double precision, intent(in) :: t
   double precision             :: a

! since "contained!  in th_nrosueaplha, the decribing parameters if not declared explictitly here
!                                        are shared (common) with those of the th-function
! the model
   a = alpha0 * exp(-(((log(t+1d-3)-log(talphamax))/talphawd)**2) / 2d0) + aoffset
   
end function alpha






 
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



       double precision            :: distar(0:N), distar0(0:N), ss, ss0
       double precision, parameter :: Delta = 0.5d0

       integer :: ipmin, ipmax, i


       aq  = l*l * q*q 

!$OMP PARALLEL DO PRIVATE(ss,ss0)
       do nn=0,N
!         ss = nn
!        distar(nn)  = (1d0 + B*exp(-ss**2/(sig*t))/(sqrt(Pi*sig*t))) * exp(-aq*ss/l)
         select case(mode)
         case(0)
           ss  = Snm(w*t,B,nn)
           ss0 = Snm(0d0,B,nn)
         case(1)
           ss  = Snm_s2(w*t,B,nn)
           ss0 = Snm_s2(0d0,B,nn)
         case(2)
           ss  = Snm_yy(w*t,B,nn)
           ss0 = Snm_yy(0d0,B,nn)
         case(3)
           ss  = Snm_s20(w*t,B,nn)
           ss0 = Snm_s20(0d0,B,nn)
         case default
           ss  = 0d0
           ss0 = 1d0
         end select
         distar(nn)  = ss  * exp(-aq*nn/(6d0))
         distar0(nn) = ss0 * exp(-aq*nn/(6d0))
!                      wt=0 wird in Snm abgefangen, ggf. hier?

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



function Snm(wt0, b, nm) result(t49)
 implicit none
 double precision, intent(in) :: wt0       ! ==> W/Ne**2 * t
 double precision, intent(in) :: b         ! ==> modifier for B (default=1)
 integer         , intent(in) :: nm        ! |n-m|


! y:=1+beta * exp(-(nm+u-v)**2/(4*wt));
! Snm:=int(y,u=0..1,v=0..1)
! with
! beta=b/sqrt(36*Pi*wt)

! Maple expression:
! Snm:=
! sqrt(wt)*b*sqrt(exp(nm/wt))/(3*sqrt(Pi)*(exp(nm^2/wt))^(1/4)*(exp(1/wt))^(1/4))+
! sqrt(wt)*b/(3*sqrt(Pi)*(exp(nm^2/wt))^(1/4)*sqrt(exp(nm/wt))*(exp(1/wt))^(1/4))-
! 2*sqrt(wt)*b/(3*sqrt(Pi)*(exp(nm^2/wt))^(1/4))-(1/6)*b*erf(-1/(2*sqrt(wt))+
! nm/(2*sqrt(wt)))+(1/6)*b*erf(-1/(2*sqrt(wt))+nm/(2*sqrt(wt)))*nm+
! (1/6)*b*erf(1/(2*sqrt(wt))+nm/(2*sqrt(wt)))+(1/6)*b*erf(1/(2*sqrt(wt))+
! nm/(2*sqrt(wt)))*nm-(1/3)*erf(nm/(2*sqrt(wt)))*b*nm+1
! 


 double precision :: t49
 real(kind=16)    :: wt, t1,t2,t4,t5,t6,t7,t9,t10,t11,t13,t14,t16,t17,t18,t30,t31,t33,t34,t39,t40,t45

      if(wt0 < 1d-8) then
        wt = 1d-8
      else
        wt = wt0
      endif 

      t1 = sqrt(wt)
      t2 = sqrt(0.3141592654D1)
      t4 = 0.1D1 / t2 * t1
      t5 = b * t4
      t6 = 0.1D1 / wt
      t7 = nm ** 2
      t9 = exp(t7 * t6)
      t10 = t9 ** (0.1D1 / 0.4D1)
      t11 = 0.1D1 / t10
      t13 = exp(nm * t6)
      t14 = sqrt(t13)
      t16 = exp(t6)
      t17 = t16 ** (0.1D1 / 0.4D1)
      t18 = 0.1D1 / t17
      t30 = 0.1D1 / t1
      t31 = nm * t30
      t33 = erf(-t30 / 0.2D1 + t31 / 0.2D1)
      t34 = t33 * b
      t39 = erf(t30 / 0.2D1 + t31 / 0.2D1)
      t40 = t39 * b
      t45 = erf(t31 / 0.2D1)
      t49 = t18 * t14 * t11 * t5 / 0.3D1 + t18 / t14 * t11 * t5 / 0.3D1 &
          - 0.2D1 / 0.3D1 * t11 * b * t4 - t34 / 0.6D1 + nm * t34 / 0.6D1 + &
            t40 / 0.6D1 + nm * t40 / 0.6D1 - nm * b * t45 / 0.3D1 + 0.1D1
      if(isnan(t49)) t49=1

end function Snm



function Snm_yy(wt0, b, nm) result(yyy)
 implicit none
 double precision, intent(in) :: wt0       ! ==> W/Ne**2 * t
 double precision, intent(in) :: b         ! ==> modifier for B (default=1)
 integer         , intent(in) :: nm        ! |n-m|

 double precision :: yyy
 double precision :: eps=1d0/4d0, wt
 double precision, parameter :: Pi=4d0*atan(1d0)

 wt=wt0
 if(wt<1d-6) wt=1d-6

 eps = eps_yy


 yyy = 1d0 + b*exp(-(eps+nm)**2/(4*wt0))/(6*sqrt(Pi*wt))

end function Snm_yy


function Snm_s2(wt0, b, nm) result(s2)
 implicit none
 double precision, intent(in) :: wt0       ! ==> W/Ne**2 * t
 double precision, intent(in) :: b         ! ==> modifier for B (default=1)
 integer         , intent(in) :: nm        ! |n-m|

 double precision :: s2
 double precision :: eps=1d0/4d0, wt
 double precision, parameter :: Pi=4d0*atan(1d0)

 wt=wt0
 if(wt<1d-6) wt=1d-6

 s2 = 1d0 + b*(erf((2d0*nm+1d0)/(4d0*sqrt(wt)))-erf((2d0*nm-1d0)/(4d0*sqrt(wt))))/6d0

end function Snm_s2


function Snm_s20(wt0, b, nm) result(s2)
 implicit none
 double precision, intent(in) :: wt0       ! ==> W/Ne**2 * t
 double precision, intent(in) :: b         ! ==> modifier for B (default=1)
 integer         , intent(in) :: nm        ! |n-m|

 double precision :: s2
 double precision :: eps=1d0/4d0, wt, swt, acc
 double precision, parameter :: Pi=4d0*atan(1d0)
 double precision, parameter :: sPi=sqrt(4d0*atan(1d0))
 double precision            :: nmm, nmp

 wt=wt0
 if(wt<1d-6) wt=1d-6
 swt = sqrt(wt)
 nmm =nm-1
 nmp =nm+1

 acc =        wt*exp(-(nmm)**2/(4*wt))
 acc = acc +  wt*exp(-(nmp)**2/(4*wt))
 acc = acc -2*wt*exp(-(dble(nm)**2)/(4*wt))
 acc = acc +  sPi*swt*(-nmm)*erf((-nmm)/(2*swt))/2
 acc = acc +  sPi*swt*( nmp)*erf(( nmp)/(2*swt))/2
 acc = acc -  sPi*swt*2*erf(nm/(2*swt))*nm/2

 s2 = 1d0 + b * acc / (3*sPi*swt)


end function Snm_s20

 end function th_dlocrep4
