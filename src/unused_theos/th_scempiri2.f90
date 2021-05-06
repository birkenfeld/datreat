 FUNCTION th_scempiri(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  phenomenological approach to sc-hysteresis-curves
!  mm
      use theory_description 
      implicit none 
      real    :: th_scempiri
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
     double precision :: mmax       ! prefactor       
     double precision :: hoffset    ! offset (peakpos)
     double precision :: hsharp     ! sharpmess1   
     double precision :: hsharp4     ! sharpmess1
     double precision :: hbeta      ! asymptotics                                                                     
! the reout parameter representation 
 
     double precision :: th
 
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'scempiri'
       nparx =        5
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_scempiri = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " phenomenological approach to sc-hysteresis-curves"
       th_citation(idesc)     = " mm"
!       --------------> set the parameter names --->
        parnam ( 1) = 'mmax    '  ! prefactor                                                                       
        parnam ( 2) = 'hoffset '  ! offset (peakpos)                                                                
        parnam ( 3) = 'hsharp  '  ! sharpmess1                                                                      
        parnam ( 4) = 'hbeta   '  ! asymptotics                                                                     
        parnam ( 5) = 'hsharp4 '  ! sharpness4                                                                     
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "offset (peakpos)" !//cr//parspace//&
        th_param_desc( 3,idesc) = "sharpmess1" !//cr//parspace//&
        th_param_desc( 4,idesc) = "asymptotics" !//cr//parspace//&
        th_param_desc( 5,idesc) = "sharpness4" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "         > "
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
! 
        th_scempiri = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      mmax     =      pa( 1)
      hoffset  =      pa( 2)
      hsharp   =      pa( 3)
      hbeta    =      pa( 4)
      hsharp4  =      pa( 5)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: 
!      xh = 
!      call parget('        ',xh,iadda,ier)
!               = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     th  = mmax * (1d0+hsharp*((x-hoffset))**2+hsharp4*(x-hoffset)**4)**(-hbeta)

     th_scempiri = th
 
! ---- writing computed parameters to the record >>>  

 end function th_scempiri
