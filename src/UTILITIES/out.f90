 FUNCTION th_fatkim(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  Fatkullin Kimmich expression diffusion on a Gaussian Path
! 
      use theory_description 
      implicit none 
      real    :: th_fatkim
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
     double precision :: wl4_fk     ! Rouse rate                                                                      
     double precision :: a_fk       ! tube step length                                                                
     double precision :: numod      ! exponent modification factor (default=1)                                        
! the recin parameter representation 
     double precision :: q          ! q-value (momentum tranbsfer)                                                    
! the reout parameter representation 
     double precision ::            !                                                                                 
 
     double precision :: th
 
     double precision, parameter :: Navogadro = 6.022140857d23
     double precision :: t
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'fatkim'
       nparx =        4
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_fatkim = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " Fatkullin Kimmich expression diffusion on a Gaussian Path"
       th_citation(idesc)     = ""
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'wl4_fk  '  ! Rouse rate                                                                      
        parnam ( 3) = 'a_fk    '  ! tube step length                                                                
        parnam ( 4) = 'numod   '  ! exponent modification factor (default=1)                                        
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "Rouse rate" !//cr//parspace//&
        th_param_desc( 3,idesc) = "tube step length" !//cr//parspace//&
        th_param_desc( 4,idesc) = "exponent modification factor (default=1)" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q-value (momentum tranbsfer)"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "         > "
! 
        th_fatkim = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      wl4_fk   =      pa( 2)
      a_fk     =      pa( 3)
      numod    =      pa( 4)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: q-value (momentum tranbsfer)
      xh = 
      call parget('q       ',xh,iadda,ier)
      q        = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     t   = x

     th = ampli * AFK(t,q,a_fk,Wl4_fk,numod)

     th_fatkim = th
 
! ---- writing computed parameters to the record >>>  
      call parset('        ',sngl(),iadda,ier)
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 
! Fatkullin-Kimmich-Form von Sinc im local Reptation Regime !
         function AFK(t,Q,a,Wl4,addexp)

         implicit none
         double precision :: AFK
         double precision, parameter :: Pi=4d0*atan(1d0)

         double precision, intent(in) :: Wl4,t,a,AFK,Q,addexp
         double precision             :: arg1, arg2,PhiRouse,SS

         PhiRouse = 2*sqrt(Wl4*t/Pi)
         SS       = PhiRouse/3.0d0
         SS       = SS**addexp
         arg1     = Q**4*a**2*SS/72.0d0
         arg2     = Q**2*a*sqrt(SS)/6.0d0/sqrt(2.0d0)
         ifail    = 1
         AFK      = exp(arg1) * erfc(arg2)

         end function AFK

 end function th_fatkim
