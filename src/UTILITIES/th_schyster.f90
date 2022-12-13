 MODULE hysteresis_empirical
implicit none

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
       if(h_pristine) A=0
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
       if(h_pristine) A=0
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

       increasing =  (h > h_last) 
       if( (increasing .and. .not. increasing_last) .or. &
           (.not. increasing .and. increasing_last)) then
         h_turn = h_last
         h_pristine = .false.
       endif
           
       h_last = h
       increasing_last = increasing
  
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

END MODULE hysteresis_empirical


FUNCTION th_schyster(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  phenomenological approach to sc-hysteresis-curves, with history and B reversal
!  mm
      use theory_description 
      use hysteresis_empirical

      implicit none 
      real    :: th_schyster
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
 
! the reout parameter representation 
 
     double precision :: th
 
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'schyster'
       nparx =        3
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_schyster = 0
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
        parnam ( 2) = 'h_slope '  ! slope of reversal approach (=field expulsion)                                   
        parnam ( 3) = 'h_sense '  ! sense of running through hysteresis loops                                       
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "slope of reversal approach (=field expulsion)" !//cr//parspace//&
        th_param_desc( 3,idesc) = "sense of running through hysteresis loops" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "         > "
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
! 
        th_schyster = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli     =      pa( 1)
      mb_slope  =      pa( 2)
      h_sense   = nint(pa( 3))
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     th  = ampli * mhyster(x)
     th_schyster = th
 
! ---- writing computed parameters to the record >>>  
 
! subroutines and functions entered here are private to this theory and share its variables 
 
 end function th_schyster
