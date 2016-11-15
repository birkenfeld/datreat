      FUNCTION th_anisodiff (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!     ===================================================
! anisotropic diffusionincoherent with PFG NMR options
! ref. Hall P-L. and Ross D.K. , Molecular Physics , 1978, p 1549-1554
!
      use theory_description
      implicit none
      CHARACTER(8) thnam, parnam (20)
      real    :: pa (20), qq (3)
      real    :: th_anisodiff
      real    :: x, xh
      integer :: mbuf, nparx, ier, ini, npar
      integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
		real, intent(inout) :: params(mbuf)             ! value des parameters n

      double precision, parameter    :: Pi = 4*atan(1d0)
      double precision :: q, tau, diffusion, lz       ! will be valid in the scope of sqt_ker and used there, keep the names!
      double precision :: adapint, erracc, eps = 1d-12
      integer          :: maxit = 10000
      integer          :: mode

      double precision :: wtextur = 1d0
      double precision :: theta0  = 0d0
      double precision :: wnorm

!
! ----- initialisation -----
      IF (ini.eq.0) then
         thnam = 'anisodi'
         nparx = 8
         IF (npar.lt.nparx) then
            WRITE (6, 1) thnam, nparx, npar
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)
            th_anisodiff = 0
            RETURN
         ENDIF
         npar = nparx
! >>>>> describe theory with    8 parameters >>>>>>>
        idesc = next_th_desc()
        th_identifier(idesc)   = thnam
        th_explanation(idesc)  = "1D anisotropic diffusion incoherent scattering or with PFG NMR options"//cr//&
                                 "           as selected by mode parameter                           "!
 !
        th_citation(idesc) = "ref. Hall P-L. and Ross D.K. , Molecular Physics , 1978, p 1549-1554"
!        --------------> set the number of parameters
         parnam (1) = 'amplit'  ! std prefactor
         parnam (2) = 'diff'    ! diffusion constant
         parnam (3) = 'lz'      ! length of 1D diffusion range
         parnam (4) = 'qval'    ! q (if in data params will be from there)
         parnam (5) = 'tau'     ! tau (if in data params will be taken from there)
         parnam (6) = 'mode'    ! mode 0 x=tau, mode 1 x= q, mode 2 x = q*q*tau ==> PFG NMR data
         parnam (7) = 'theta0'  ! textur center default 0
         parnam (8) = 'wtextur' ! textur weite
! >>>>> describe parameters >>>>>>>
        th_param_desc( 1,idesc) = " prefactor " !//cr//parspace//&
        th_param_desc( 2,idesc) = " diffusion constant in l-units**2/t-units  " !//cr//parspace//&
        th_param_desc( 3,idesc) = " length of 1D diffusion range in l-units   " !//cr//parspace//&
        th_param_desc( 4,idesc) = " q-value if not given in record parameters in 1/l-units " !//cr//parspace//&
        th_param_desc( 5,idesc) = " time-value if not given in record paramters  " !//cr//parspace//&
        th_param_desc( 6,idesc) = " mode 0: x=tau, mode 1: x= q, mode 2: x = q*q*tau ==> PFG NMR " !//cr//parspace//&
        th_param_desc( 7,idesc) = " texture center default 0 " !//cr//parspace//&
        th_param_desc( 8,idesc) = " texture width " !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
       th_file_param(:,idesc) = " " 
       th_file_param(1,idesc) = "q   > q-value   " 
       th_file_param(2,idesc) = "tau > tiem value (only one of those depending on mode) " 
! >>>>> describe record parameters creaqted by this theory >>>>>>>
       th_out_param(:,idesc)  = " "  
       th_out_param(1,idesc)  = "none "  
!


         th_anisodiff = 0.0

         RETURN
      ENDIF
!
! ---- calculate theory here -----
      q   = pa(4)
      CALL getpar ('q       ', xh,nopar ,params,napar,mbuf, ier)
      if(ier.eq.0) q = xh

      tau = pa(5)
      CALL getpar ('tau     ', xh,nopar ,params,napar,mbuf, ier)
      if(ier.eq.0) tau = xh


      diffusion   = abs(pa(2))
      lz          = abs(pa(3))
      mode        = nint(pa(6))

      theta0      = abs(pa(7))
      wtextur     = abs(pa(8))

      select case(mode)
      case (0)       ! x axis is assumed to be time (std sqt), q from params or fit params
         tau = x
         th_anisodiff = pa(1) * 0.5d0 * adapint(sqt_ker,0d0, Pi, eps, maxit, erracc)
      case (1)       ! x axis is assumed to be q
         q = x
         th_anisodiff = pa(1) * 0.5d0 * adapint(sqt_ker,0d0, Pi, eps, maxit, erracc)
      case (2)       ! pfg scale x=q*q*tau
         q = sqrt(x / tau)
         th_anisodiff = pa(1) * 0.5d0 * adapint(sqt_ker,0d0, Pi, eps, maxit, erracc)
      case (3)       ! pfg scale x=q*q*tau mit textur
         q = sqrt(x / tau)
         wnorm        = adapint(w_textur_ker,0d0, Pi, eps, maxit, erracc)
         th_anisodiff = pa(1) * adapint(sqt_textur_ker,0d0, Pi, eps, maxit, erracc) / wnorm
      case default
        write(6,*)"anisodiff mode: ",mode, " not supported"
        th_anisodiff = 0
      end select



contains

double precision function a_n(n,qzl)
! --------------------------------
! Eqns (11)+(12) of ref
  implicit none
  integer         , intent(in)  :: n
  double precision, intent(in)  :: qzl

  double precision, parameter   :: Pi = 4*atan(1d0)
  double precision, parameter   :: eps = 1d-9

  double precision              :: denom

  if(n==0) then
    denom = qzl**2
    if(denom < eps) then
      a_n = 1d0
    else
      a_n = 2/denom * (1d0-cos(qzl))
    endif
    return
  endif

  denom = (qzl**2 - (n*Pi)**2)**2
  if(denom < eps) then
    a_n = 0.5d0
  else
    a_n = (2*qzl)**2/denom * (1d0 - cos(qzl)*(-1d0)**n)
  endif

end function a_n


double precision function lorenzian(dqq, omega)
   implicit none
   double precision,  intent(in)  :: dqq
   double precision,  intent(in)  :: omega
   double precision, parameter    :: Pi = 4*atan(1d0)


   lorenzian   = dqq/Pi/(dqq**2+omega**2)

end function lorenzian

double precision function ftlor(dqq, t)
   implicit none
   double precision,  intent(in)  :: dqq
   double precision,  intent(in)  :: t

   ftlor       = exp(-abs(dqq*t))

end function ftlor


double precision function sqt(diff, qz, l, t)
   implicit none
   double precision, intent(in) :: diff
   double precision, intent(in) :: qz
   double precision, intent(in) :: l
   double precision, intent(in) :: t

   double precision, parameter    :: Pi = 4*atan(1d0)
   double precision :: a, eps = 1d-8
   integer :: i, imax = 200

   sqt  = 0
   imax = nint(sqrt(2d0 * sqrt((qz*l)**2 * eps**3)) / (Pi*eps)) + 1
   do i=0,imax
     a = a_n(i,qz*l)
     sqt = sqt + a * ftlor(diff*(i*Pi/l)**2,t)
   enddo

end function sqt


double precision function sqt_ker(theta)
   implicit none
   double precision, intent(in) ::  theta

   sqt_ker = sqt(diffusion, q * cos(theta), lz, tau) * sin(theta)


end function sqt_ker

double precision function w_textur(theta)
   implicit none
   double precision, intent(in) ::  theta

   w_textur = exp(-(((theta-theta0)/(wtextur))**2)/2d0)

end function w_textur

double precision function w_textur_ker(theta)
   implicit none
   double precision, intent(in) ::  theta

   w_textur_ker = w_textur(theta) * sin(theta)

end function w_textur_ker


double precision function sqt_textur_ker(theta)
   implicit none
   double precision, intent(in) ::  theta

   sqt_textur_ker = w_textur(theta) * sqt(diffusion, q * cos(theta), lz, tau) * sin(theta)

end function sqt_textur_ker



END FUNCTION th_anisodiff

