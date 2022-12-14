 FUNCTION th_ntc(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  ntc-kennlinie
!  J.K. Kjems, T. Freloft, D. Richter and S.K. Sinha, Physica 136B (1986) 285-290
      use theory_description 
      implicit none 
      real    :: th_ntc
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
     double precision :: a0         ! a0-coeff (const)                                                                
     double precision :: a1         ! a1 coeff (ln(r))                                                                
     double precision :: a2         ! a2 coeff (ln(r)**2)                                                             
     double precision :: a3         ! a3-coeff (ln(r)**3)                                                             
! the recin parameter representation 
     double precision :: th
 
     double precision   :: r, temp
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'ntc'
       nparx =        4
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_ntc = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " ntc-kennlinie"
       th_citation(idesc)     = " J.K. Kjems, T. Freloft, D. Richter and S.K. Sinha, Physica 136B (1986) 285-290"
!       --------------> set the parameter names --->
        parnam ( 1) = 'a0      '  ! a0-coeff (const)                                                                
        parnam ( 2) = 'a1      '  ! a1 coeff (ln(r))                                                                
        parnam ( 3) = 'a2      '  ! a2 coeff (ln(r)**2)                                                             
        parnam ( 4) = 'a3      '  ! a3-coeff (ln(r)**3)                                                             
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "a0-coeff (const)" !//cr//parspace//&
        th_param_desc( 2,idesc) = "a1 coeff (ln(r))" !//cr//parspace//&
        th_param_desc( 3,idesc) = "a2 coeff (ln(r)**2)" !//cr//parspace//&
        th_param_desc( 4,idesc) = "a3-coeff (ln(r)**3)" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "         > "
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
! 
        th_ntc = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      a0       =      pa( 1)
      a1       =      pa( 2)
      a2       =      pa( 3)
      a3       =      pa( 4)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: 
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     r = x

     th = 1d0/(a0+ a1*log(r) +a2*log(r)**2 +a3*log(r)**3) - 273.15d0


     th_ntc = th
 
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 
 end function th_ntc
