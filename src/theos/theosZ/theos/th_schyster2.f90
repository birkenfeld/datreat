 MODULE hysteresis_empirical2
implicit none

SAVE 

integer,parameter :: max_comp = 6
integer           :: n_comp   = 4
double precision  :: h_comp(max_comp) =   [ &
                        2.686883d0        , &
                        0.2304476d0       , &    
                        0.2658331d0       , &    
                       -0.7871898d0       , 0d0, 0d0 ]

double precision  :: h_offset(max_comp) = [ & 
                        0.7532015d-03     , &   
                       -0.3018093d0       , &   
                       -0.4798850d0       , &   
                        0.6171697d-01     , 0d0, 0d0  ]   

double precision  :: h_hsharp(max_comp) = [ &
                       11.76612d0         , &    
                       6.250469d0         , &    
                       6.184682d0         , &    
                       16.06952d0         , 0d0, 0d0  ]      

double precision  :: h_beta(max_comp)   = [ &
                        0.5394006d0       , &   
                        0.2124426d0       , &   
                        0.7496565d-01     , &   
                        0.5801371d0       , 0d0, 0d0 ]  

double precision :: h_const =  0.5495992d-02 
double precision :: h_slope =  -0.3783910d-01

double precision :: mb_slope0 = 100d0/2     !! may be required to make slopüe h-dependent
double precision :: mb_slope1 = 100d0/2     !! may be required to make slopüe h-dependent
double precision :: mb_slope2 = 1d0         !! slope of slope with h
double precision :: mb_slope3 = 1d0         !! slope of slope with h

double precision :: h_last = 0
double precision :: h_turn = 1

integer          :: h_sense = 1             !! hysteresis loop sense (1 ==> "counter clockwise")
logical          :: increasing_last
logical          :: h_pristine = .true.     !! if true we start unmagnetized



contains

     function m_sc_upper(hb) result(m)
       implicit none
       double precision, intent(in) :: hb
       double precision             :: m
       integer :: i
        m = h_const + h_slope * hb
        do i=1,n_comp
          m = m + h_comp(i) * ( 1d0+(h_hsharp(i) * (hb-h_offset(i)))**2 ) ** (-h_beta(i))
        enddo

     end function m_sc_upper


     function m_sc_lower(hb) result(m)
       implicit none
       double precision, intent(in) :: hb
       double precision             :: m
       m = -m_sc_upper(-hb)
     end function m_sc_lower


     function m_transition_lowup(hb, hturn) result(m)
       implicit none
       double precision, intent(in) :: hb
       double precision, intent(in) :: hturn
       double precision             :: m
       double precision :: A, Ba, Bb
       A   = m_sc_upper(hb)
       if(abs(hb) > abs(hturn) .and. hb*hturn > 0) then
        m   = m_sc_lower(hb)
        return
       endif
       Bb  = m_sc_lower(hturn)
!       if(h_pristine) A=0
!       if(h_pristine) Bb=0
       m   = A + (Bb-A) * exp(-mb_slope0*max(0d0, -(hb-hturn)) )
     end function m_transition_lowup


     function m_transition_uplow(hb, hturn) result(m)
       implicit none
       double precision, intent(in) :: hb
       double precision, intent(in) :: hturn
       double precision             :: m

       double precision :: A, Ba, Bb
       A   = m_sc_lower(hb)
       Bb  = m_sc_upper(hturn)
!       if(h_pristine) A=0
       if(h_pristine) Bb=0
       m   = A + (Bb-A) * exp(-mb_slope0*max(0d0, (hb-hturn)) )
     end function m_transition_uplow

     function m_transition_lowup_R(hb, hturn) result(m)
       implicit none
       double precision, intent(in) :: hb
       double precision, intent(in) :: hturn
       double precision             :: m
       double precision :: A, Ba, Bb
       A   = m_sc_upper(hb)   
       Bb  = m_sc_lower(hturn)
       m   = A + (Bb-A) * exp(-mb_slope0*max(0d0, (hb-hturn)) )
     end function m_transition_lowup_R


     function m_transition_uplow_R(hb, hturn) result(m)
       implicit none
       double precision, intent(in) :: hb
       double precision, intent(in) :: hturn
       double precision             :: m

       double precision :: A, Ba, Bb
       A   = m_sc_lower(hb)
       Bb  = m_sc_upper(hturn)
       m   = A + (Bb-A) * exp(-mb_slope0*max(0d0, -(hb-hturn)) )
     end function m_transition_uplow_R


     function mhyster(h) result(m)   !! this is for one direction of "looping"
       implicit none
       double precision, intent(in) :: h
       double precision             :: m
       logical                      :: increasing

!       increasing =  (h > h_last) 
       increasing =  (h > h_last+1d-5) 
       if( (increasing .and. .not. increasing_last) .or. &
           (.not. increasing .and. increasing_last)) then
         h_turn = h_last
         h_pristine = .false.
       endif
           
       h_last = h
       increasing_last = increasing

!! tentative:
       mb_slope0 = mb_slope1 + h**2 * mb_slope2 + h**4 * mb_slope3


  
       if(h_sense==1) then                        !! depending on the sense with which the
                                                  !! hysteresiscurves are probed
          if(increasing) then
             m = m_transition_uplow(h, h_turn)
          else
             m = m_transition_lowup(h, h_turn) 
          endif
       else
          if(.not.increasing) then
             m = m_transition_uplow_R(h, h_turn)
          else
             m = m_transition_lowup_R(h, h_turn) 
          endif
       endif
       
     end function mhyster

END MODULE hysteresis_empirical2


FUNCTION th_schyste2(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  phenomenological approach to sc-hysteresis-curves, with history and B reversal
!  mm
      use theory_description 
      use hysteresis_empirical2

      implicit none 
      real    :: th_schyste2
      character(len=8) :: thnam, parnam (*) 
      real    :: pa (*) 
      real    :: x , xh
      integer :: mbuf, nparx, ier, ini, npar, iadda
      integer, intent(inout) :: nopar       
      character(len=80), intent(inout) :: napar(mbuf) 
      real, intent(inout) :: params(mbuf) 
     
      double precision, parameter :: Pi = 4*atan(1d0)
      integer                     :: actual_record_address
      integer                     :: actual_point
     
! the internal parameter representation 
     double precision :: ampli      ! prefactor   
 
     double precision :: amp_mod(4), off_mod(4)                                                                    
! the reout parameter representation 
 
     double precision :: th
 
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'schyste2'
       nparx =        5+16+2
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_schyste2 = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " phenomenological approach to sc-hysteresis-curves, with history and B reversal"
       th_citation(idesc)     = " mm"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'h_slope0 '  ! slope of reversal approach (=field expulsion)   
        parnam ( 3) = 'h_slope1 '  ! slope of reversal approach (=field expulsion) 
        parnam ( 4) = 'h_slope2 '  ! slope of reversal approach (=field expulsion) 
        parnam ( 5) = 'h_sense  '  ! sense of running through hysteresis loops 
        parnam ( 6) = 'amp_mod1' ! sense of running through hysteresis loops                                       
        parnam ( 7) = 'amp_mod2' ! sense of running through hysteresis loops                                       
        parnam ( 8) = 'amp_mod3' ! sense of running through hysteresis loops                                       
        parnam ( 9) = 'amp_mod4' ! sense of running through hysteresis loops                                       
        parnam (10) = 'off_mod1'  ! sense of running through hysteresis loops                                       
        parnam (11) = 'off_mod2'  ! sense of running through hysteresis loops                                       
        parnam (12) = 'off_mod3'  ! sense of running through hysteresis loops                                       
        parnam (13) = 'off_mod4'  ! sense of running through hysteresis loops                                       
        parnam (14) = 'hsh_mod1'  ! sense of running through hysteresis loops                                       
        parnam (15) = 'hsh_mod2'  ! sense of running through hysteresis loops                                       
        parnam (16) = 'hsh_mod3'  ! sense of running through hysteresis loops                                       
        parnam (17) = 'hsh_mod4'  ! sense of running through hysteresis loops                                       
        parnam (18) = 'bet_mod1'  ! sense of running through hysteresis loops                                       
        parnam (19) = 'bet_mod2'  ! sense of running through hysteresis loops                                       
        parnam (20) = 'bet_mod3'  ! sense of running through hysteresis loops                                       
        parnam (21) = 'bet_mod4'  ! sense of running through hysteresis loops                                       
        parnam (22) = 'bgr_mod'  ! sense of running through hysteresis loops                                       
        parnam (23) = 'slo_mod'  ! sense of running through hysteresis loops                                       
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "slope of reversal approach (=field expulsion)" !//cr//parspace//&
        th_param_desc( 3,idesc) = "slope1" !//cr//parspace//&
        th_param_desc( 4,idesc) = "slope2" !//cr//parspace//&
        th_param_desc( 5,idesc) = "sense of running through hysteresis loops" !//cr//parspace//&
        th_param_desc( 6:21,idesc) = "new h_comp, h_off, h_sharp, h_beta parameters " !//cr//parspace//&
        th_param_desc( 22:23,idesc) = "bgr and slope for hysteresis " !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "         > "
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
! 
        th_schyste2 = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli     =      pa( 1)
      mb_slope1 =      pa( 2)
      mb_slope2 =      pa( 3)
      mb_slope3 =      pa( 4)
      h_sense   = nint(pa( 5))
      h_comp(1) =   pa( 6) ! 'amp_mod1' 
      h_comp(2) =   pa( 7) ! 'amp_mod2' 
      h_comp(3) =   pa( 8) ! 'amp_mod3' 
      h_comp(4) =   pa( 9) ! 'amp_mod4' 
      h_offset(1) =   pa(10) ! 'off_mod1' 
      h_offset(2) =   pa(11) ! 'off_mod2' 
      h_offset(3) =   pa(12) ! 'off_mod3' 
      h_offset(4) =   pa(13) ! 'off_mod4' 
      h_hsharp(1) =   pa(14) ! 'hsh_mod1' 
      h_hsharp(2) =   pa(15) ! 'hsh_mod2' 
      h_hsharp(3) =   pa(16) ! 'hsh_mod3' 
      h_hsharp(4) =   pa(17) ! 'hsh_mod4' 
      h_beta(1) =   pa(18) ! 'bet_mod1' 
      h_beta(2) =   pa(19) ! 'bet_mod2' 
      h_beta(3) =   pa(20) ! 'bet_mod3' 
      h_beta(4) =   pa(21) ! 'bet_mod4' 
      h_const   =   pa(22) ! 'bgr_mod'  
      h_slope   =   pa(23) ! 'slo_mod'  


! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
     
      if(actual_point() == 10)  h_pristine = .true. 
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     th  = ampli * mhyster(dble(x))
     th_schyste2 = th
 
! ---- writing computed parameters to the record >>>  
 
! subroutines and functions entered here are private to this theory and share its variables 
 
 end function th_schyste2
