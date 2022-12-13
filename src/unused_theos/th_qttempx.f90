 FUNCTION th_qttempx(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  general expression allowing the insertion of extracted parameters q**eq * t**et * temp**etem * xpar**ex
! 
      use theory_description 
      implicit none 
      real    :: th_qttempx
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
     double precision :: eq         ! exponent of q                                                                   
     double precision :: et         ! exponent pf t                                                                   
     double precision :: ettemp     ! exponent pf temp                                                                
     double precision :: ex1        ! exponent of xpar1 ( a free parameter at disposition)                            
     double precision :: ex2        ! exponent of xpar2 ( a free parameter at disposition)                            
     double precision :: ex3        ! exponent of xpar3 ( a free parameter at disposition)                            
! the recin parameter representation 
     double precision :: q          ! q-value                                                                         
     double precision :: t          ! t-value                                                                         
     double precision :: temp       ! temp-value                                                                      
     double precision :: xpar1      ! extra parameter free                                                            
     double precision :: xpar2      ! extra parameter free                                                            
     double precision :: xpar3      ! extra parameter free                                                            
! the reout parameter representation 
 
     double precision :: th
 
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'qttempx'
       nparx =        7
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_qttempx = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " general expression allowing the insertion of extracted parameters q**eq * t**et * temp**etem * xpar**ex"
       th_citation(idesc)     = ""
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'eq      '  ! exponent of q                                                                   
        parnam ( 3) = 'et      '  ! exponent pf t                                                                   
        parnam ( 4) = 'ettemp  '  ! exponent pf temp                                                                
        parnam ( 5) = 'ex1     '  ! exponent of xpar1 ( a free parameter at disposition)                            
        parnam ( 6) = 'ex2     '  ! exponent of xpar2 ( a free parameter at disposition)                            
        parnam ( 7) = 'ex3     '  ! exponent of xpar3 ( a free parameter at disposition)                            
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "exponent of q" !//cr//parspace//&
        th_param_desc( 3,idesc) = "exponent pf t" !//cr//parspace//&
        th_param_desc( 4,idesc) = "exponent pf temp" !//cr//parspace//&
        th_param_desc( 5,idesc) = "exponent of xpar1 ( a free parameter at disposition)" !//cr//parspace//&
        th_param_desc( 6,idesc) = "exponent of xpar2 ( a free parameter at disposition)" !//cr//parspace//&
        th_param_desc( 7,idesc) = "exponent of xpar3 ( a free parameter at disposition)" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q-value"
        th_file_param(  2,idesc) = "t        > t-value"
        th_file_param(  3,idesc) = "temp     > temp-value"
        th_file_param(  4,idesc) = "xpar1    > extra parameter free"
        th_file_param(  5,idesc) = "xpar2    > extra parameter free"
        th_file_param(  6,idesc) = "xpar3    > extra parameter free"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
! 
        th_qttempx = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      eq       =      pa( 2)
      et       =      pa( 3)
      ettemp   =      pa( 4)
      ex1      =      pa( 5)
      ex2      =      pa( 6)
      ex3      =      pa( 7)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: q-value
      xh =          1.0
      call parget('q       ',xh,iadda,ier)
      q        = xh
! >>> extract: t-value
      xh =          1.0
      call parget('t       ',xh,iadda,ier)
      t        = xh
! >>> extract: temp-value
      xh =          1.0
      call parget('temp    ',xh,iadda,ier)
      temp     = xh
! >>> extract: extra parameter free
      xh =          1.0
      call parget('xpar1   ',xh,iadda,ier)
      xpar1    = xh
! >>> extract: extra parameter free
      xh =          1.0
      call parget('xpar2   ',xh,iadda,ier)
      xpar2    = xh
! >>> extract: extra parameter free
      xh =          1.0
      call parget('xpar3   ',xh,iadda,ier)
      xpar3    = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     th = q**eq * t**et * temp**ettemp * xpar1**ex1 * xpar2**ex2 * xpar3**ex3
     th_qttempx = th * ampli
 
! ---- writing computed parameters to the record >>>  
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 
 end function th_qttempx
