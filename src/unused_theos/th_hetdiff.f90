 FUNCTION th_hetdif(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  Heterogeneous diffusion modell, average of various diffusion constants form a 
!  log-normal distribution w(D) = 1/N * exp(-[(ln(D/d0)-ln(dc/d0))/sigma]**2 
!  th   = int(w(D) * exp(-q^2 D t)  dD)
! 
      use theory_description 
      implicit none 
      real    :: th_hetdif
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
     double precision :: d0         ! scale of diffusion (reference value)                                            
     double precision :: dc         ! central value of log-norm                                                       
     double precision :: sigma      ! width of log-normal distribution                                                
! the recin parameter representation 
     double precision :: q          ! q-value                                                                         
! the reout parameter representation 
 
     double precision   :: t
!

     double precision   :: adapint, eps=1d-11, erracc
     integer            :: maxit  = 10000
     double precision   :: sigma_multiple = 9d0
     double precision   :: a, b, nrm

! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'hetdif'
       nparx =        3
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_hetdif = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " Heterogeneous diffusion modell, average of various diffusion constants"//cr//&
                                " form a log-normal distribution w(D) = 1/N * exp(-[(ln(D/d0)-ln(dc/d0))/sigma]**2"//cr//&
                                " th   = int(w(D) * exp(-q^2 D t)  dD)"
       th_citation(idesc)     = ""
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'dc      '  ! central value of log-norm                                                       
        parnam ( 3) = 'sigma   '  ! width of log-normal distribution                                                
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "central value of log-norm" !//cr//parspace//&
        th_param_desc( 3,idesc) = "width of log-normal distribution" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q-value"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
! 
        th_hetdif = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      dc       =      pa( 2)
      sigma    =      pa( 3)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: q-value
      xh =  0.1
      call parget('q       ',xh,iadda,ier)
      q        = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     
     a   = - sigma_multiple * sigma
     b   = - a
     t   = 0d0
     nrm = adapint(lognorm_kernel,a,b, eps, maxit, erracc) 
!    write(6,*) t, a,b, nrm
     t   = x
     th_hetdif = ampli * adapint(lognorm_kernel,a,b, eps, maxit, erracc) / nrm
!     write(6,*) th_hetdif
! ---- writing computed parameters to the record >>>  
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 
  function lognorm_kernel(y) result(val)
    implicit none
    double precision, intent(in) :: y
    double precision             :: val
    double precision             :: d


!    val = exp(-((log(y/d0)-log(dc/d0))/sigma)**2) * exp(-q*q*y*t)

    d   = dc * exp(y)

    val = exp(-(y/sigma)**2) * d * exp(-d*q*q*t)

  end function lognorm_kernel
 end function th_hetdif
