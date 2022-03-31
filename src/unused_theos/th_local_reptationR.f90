 FUNCTION th_locrepr(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  generalized local reptation expression along the lines of deGennes, 
!  but with finite summation of integrals and lenthscale, 
!  timescale and fluctuation ratio as parameters
!  Journal de Physique, 1981, 42 (5), pp.735-740. <10.1051/jphys:01981004205073500>
      use theory_description 
      implicit none 
      real    :: th_locrepr
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
     double precision :: b          ! fluctuation intensity (relative)                                                
     double precision :: a          ! length scale                                                                    
     double precision :: tau        ! timescale                                                                       
     double precision :: lz         ! total length                                                                    
! the recin parameter representation 
     double precision :: q          ! q-value                                                                         

     double precision :: local_reptation_lr, sqt, sqt0
     double precision :: teps = 1d-6
 
     double precision :: nu
     integer          :: m, mode
     

     double precision   :: t
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'locrepr'
       nparx =        9
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_locrepr = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " generalized local reptation expression along the lines of deGennes, " //cr//parspace//&
                                " but with finite summation of integrals and lenthscale," //cr//parspace//&
                                " timescale and fluctuation ratio as parameters"
       th_citation(idesc)     = " Journal de Physique, 1981, 42 (5), pp.735-740. <10.1051/jphys:01981004205073500>"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'b       '  ! fluctuation intensity (relative)                                                
        parnam ( 3) = 'a       '  ! length scale                                                                    
        parnam ( 4) = 'tau     '  ! timescale                                                                       
        parnam ( 5) = 'lz      '  ! total length                                                                    
        parnam ( 6) = 'teps    '  ! near zero for t=0 evaluation 
        parnam ( 7) = 'nu      '  ! default 0.5
        parnam ( 8) = 'm       '  ! sharpness of transition 
        parnam ( 9) = 'mode    '  ! 0=std, 1=ring 
      
                                                                         
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "fluctuation intensity (relative)" !//cr//parspace//&
        th_param_desc( 3,idesc) = "length scale" !//cr//parspace//&
        th_param_desc( 4,idesc) = "timescale" !//cr//parspace//&
        th_param_desc( 5,idesc) = "total length" !//cr//parspace//&
        th_param_desc( 6,idesc) = "teps near zero for t=0 evaluation" !//cr//parspace//&
        th_param_desc( 7,idesc) = "chain expansion exponent defaut=0.5" !//cr//parspace//&
        th_param_desc( 8,idesc) = "diffusion mode summation number (only for mode=1)" !//cr//parspace//&
        th_param_desc( 9,idesc) = "mode = 0: std locrep expr, mode=1: ring variant" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q-value"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)   = " "
        th_out_param(1,idesc)   = "a2sqt = a**2/sqrt(tau) "
! 
        th_locrepr = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =         pa( 1)
      b        =   abs ( pa( 2) )
      a        =   abs ( pa( 3) )
      tau      =   abs ( pa( 4) )
      lz       =   abs ( pa( 5) )
      teps     =         pa( 6)
      nu       =   abs ( pa(7)  )
      M        =   nint( pa(8)  )
      mode     =   nint( pa(9)  )
      if (teps <= 0d0) teps = 1d-6
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: q-value
      xh =     0.1d0
      call parget('q       ',xh,iadda,ier)
      q        = xh

      t        = x
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     sqt0 =  local_reptation_lr( q*a, teps,  b, lz, nu, M, mode)
     sqt  =  local_reptation_lr( q*a, t/tau, b, lz, nu, M, mode)
     th_locrepr = ampli * sqt/sqt0


!     write(6,*) t, q, sqt, sqt0
      call parset("a2sqt   ", sngl(a*a/sqrt(tau)),iadda,ier)
                   
 end function th_locrepr




function local_reptation_lr(q, tx, B, L, nu, M, mode) result(val)
 use  ADAPTIVE_INTEGRAL

 implicit none
 double precision, intent(in)  :: q
 double precision, intent(in)  :: tx
 double precision, intent(in)  :: B
 double precision, intent(in)  :: L
 double precision, intent(in)  :: nu
 integer         , intent(in)  :: M
 integer         , intent(in)  :: mode   ! 0=std, 1=ring

 double precision :: val
 double precision :: A, La, Lb, t
 double precision, parameter :: tzero = 1d-6

 double precision :: eps 
 double precision :: erracc
 integer          :: maxit


 t     = max(tzero,tx)
 A     =   1d0
 eps   =   1d-8
 maxit = 10000


 La    = 10d0 * sqrt(t)
 Lb    = L-2*La
 if(Lb <= 0d0) then
  La = L/2
  Lb = 0
 endif


 if(mode == 0) then
   val =       adapint (lr_kernel    , 0d0  , La,    eps, maxit, erracc) 
   val = val + adapint (lr_kernel    , La   , L-La,  eps, maxit, erracc) 
   val = val + adapint (lr_kernel    , L-La , L,     eps, maxit, erracc) 
 else
   val =       adapint (lrring_kernel    , 0d0  , La,    eps, maxit, erracc) 
   val = val + adapint (lrring_kernel    , La   , L-La,  eps, maxit, erracc) 
   val = val + adapint (lrring_kernel    , L-La , L,     eps, maxit, erracc) 
 endif
  


CONTAINS

  function lr_kernel(u) result(val)
    implicit none
    double precision             :: u
    double precision             :: val
    double precision, parameter  :: pi = 4d0*atan(1d0)

!   val = (A + B*exp(-u**2/(4d0*t))/(2d0*sqrt(Pi*t)))*exp(-q**2*abs(u)**(2*nu)/6d0)
    val = (A + B* diff_1d(u)                        )*exp(-q**2*abs(u)**(2*nu)/6d0)
    val = 2*val*(L-u)/L

  end function lr_kernel



  function lrring_kernel(u) result(val)
    implicit none
    double precision             :: u
    double precision             :: val
    double precision, parameter  :: pi = 4d0*atan(1d0)

    val = (A + B*periodic_diff(u))*exp(-q**2*(abs(u)*(1d0-abs(u)/L))**(2*nu)/6d0)
    val = 2*val*(L-u)/L

  end function lrring_kernel


  function periodic_diff(u) result(val)
    implicit none
    double precision             :: u
    double precision             :: val
    double precision, parameter  :: pi = 4d0*atan(1d0)
    integer                      :: k

    val = 0d0
    do k=-M,M
      val = val + cos(2*Pi*u*k/L) * exp(-k*k*t*2d0/5d0)
    enddo 
    val = val / (L)

  end function periodic_diff


  function diff_1d(u) result(val)
    implicit none
    double precision             :: u
    double precision             :: val
    double precision, parameter  :: pi = 4d0*atan(1d0)
    integer                      :: k

    val = exp(-u**2/(4d0*t))/(2d0*sqrt(Pi*t))

  end function diff_1d

end function local_reptation_lr


