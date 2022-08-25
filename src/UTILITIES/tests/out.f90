 REAL FUNCTION th_thtester(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!====================================================================================
!  a test comment to check the th_template generator and another line for it. With complicated formulae.
!  ref. Plisch und Plum in Grimms Maerchen
      use theory_description 
      implicit none 
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
     double precision :: scope      ! scope scope scope                                                               
     double precision :: knurrrer   ! ! a knurrer                                                                     
! the recin parameter representation 
     double precision :: q          ! Q-VALUE                                                                         
     double precision :: temp       ! Temperatur der Probe                                                            
! the reout parameter representation 
     double precision :: qq         ! square of q                                                                     
 
     double precision :: tau
     integer          :: myint
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'thtester'
       nparx =        3
       IF (npar.lt.nparx) then
          WRITE (6,(' theory: ',a8,' no of parametrs=',i8,' exceeds current max. = ',i8))thnam,nparx,npar
            th_thtester = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " a test comment to check the th_template generator and another line for it. With complicated formulae."
       th_citation(idesc)     = " ref. Plisch und Plum in Grimms Maerchen"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'scope   '  ! scope scope scope                                                               
        parnam ( 3) = 'knurrrer'  ! ! a knurrer                                                                     
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "scope scope scope" !//cr//parspace//&
        th_param_desc( 3,idesc) = "! a knurrer" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > Q-VALUE"
        th_file_param(  2,idesc) = "temp     > Temperatur der Probe"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "qq       > square of q"
! 
        th_thtester = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      scope    =      pa( 2)
      knurrrer =      pa( 3)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: Q-VALUE
      xh =     0.1
      call parget('q       ',xh,iadda,ier)
      q        = xh
! >>> extract: Temperatur der Probe
      xh =     3000
      call parget('temp    ',xh,iadda,ier)
      temp     = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     tau = x
     qq  = q
     th_thtester = ! INSERT RESULTING VALUE OF TH EVALUATION HERE
 
! ---- writing computed parameters to the record >>>  
      call parset('qq      ',sngl(qq),iadda,ier)
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 
  function s2(q)
    s2 = 1
  end function s2
 end program th_thtester
