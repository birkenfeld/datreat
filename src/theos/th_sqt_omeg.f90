 FUNCTION th_sqt_omeg(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  combined s(q,omega) and s(q,t) data from an s(q,t) model
! 
      use theory_description 
      implicit none 
      real    :: th_sqt_omeg
      character(len=8) :: thnam, parnam (*) 
      real    :: pa (*) 
      real    :: x , xh
      integer :: mbuf, nparx, ier, ini, npar, iadda
      integer, intent(inout) :: nopar       
      character(len=80), intent(inout) :: napar(mbuf) 
      real, intent(inout) :: params(mbuf) 
     
      double precision, parameter :: Pi = 4*atan(1d0)
      integer                     :: actual_record_address
      logical                     :: is_in_xaxis

     character(len=8)             :: buf
     integer                      :: np
     
! the internal parameter representation 
     integer, parameter :: n_strexpo = 6   ! maximum number of streched exponentials
     double precision   :: intensit        ! overall prefactor                   
     double precision   :: omega0          ! omega scale zero shift    
     double precision   :: amp(n_strexpo)  ! amplitudes of streched exponentials
     double precision   :: tau(n_strexpo)  ! tau_values of streched exponentials
     double precision   :: nue(n_strexpo)  ! tau q exponent
     double precision   :: beta(n_strexpo) ! streching exponents of streched exponentials
     double precision   :: u_sqr           ! Debeywaller factor generation u**2

     double precision :: epsilon = 1d-8    ! integration accuracy
                                           
! the recin parameter representation 
     double precision :: q          ! q-value    default value    
     double precision :: nse        ! explict switch between sqt and sqomega (if nse=1 ==> sqt)  
     double precision :: bk1level   ! resolution description bgr-levwl    
     double precision :: bk1slope   ! resolution description slope of bgr   
     double precision :: ga1inten   ! resolution description 1. gaussian amplitude                                    
     double precision :: ga1width   ! resolution description 1. gaussian width           
     double precision :: ga1cente   ! resolution description 1. gaussian center          
     double precision :: ga2inten   ! resolution description 2. gaussian amplitude       
     double precision :: ga2width   ! resolution description 2. gaussian width           
     double precision :: ga2cente   ! resolution description 2. gaussian center          
     double precision :: ga3inten   ! resolution description 3. gaussian amplitude       
     double precision :: ga3width   ! resolution description 3. gaussian width           
     double precision :: ga3cente   ! resolution description 3. gaussian center          
     double precision :: ga4inten   ! resolution description 4. gaussian amplitude       
     double precision :: ga4width   ! resolution description 4. gaussian width           
     double precision :: ga4cente   ! resolution description 4. gaussian center          
     double precision :: ga5inten   ! resolution description 5. gaussian amplitude       
     double precision :: ga5width   ! resolution description 5. gaussian width           
     double precision :: ga5cente   ! resolution description 5. gaussian center          
     double precision :: ga6inten   ! resolution description 6. gaussian amplitude       
     double precision :: ga6width   ! resolution description 6. gaussian width           
     double precision :: ga6cente   ! resolution description 6. gaussian center          
     double precision :: ga7inten   ! resolution description 7. gaussian amplitude       
     double precision :: ga7width   ! resolution description 7. gaussian width           
     double precision :: ga7cente   ! resolution description 7. gaussian center          
     double precision :: ga8inten   ! resolution description 8. gaussian amplitude       
     double precision :: ga8width   ! resolution description 8. gaussian width           
     double precision :: ga8cente   ! resolution description 8. gaussian center          

! the reout parameter representation 
     double precision :: xwidth     !



                                  
 
     double precision :: th
 
     integer            :: maxit = 5000
     double precision   :: erraccu
     integer            :: i
     double precision   :: t
     double precision   :: omega
     double precision   :: o0
     double precision   :: rsum
!?     double precision   :: str_tau0, str_beta
!?     double precision   :: str_delta
     double precision   :: a,b
     double precision   :: res
     double precision   :: dwf
     double precision   :: gampli(8)
     double precision   :: gwidth(8)
     double precision   :: gcenter(8)
     double precision   :: str_delta

     double precision   :: adapint

     double precision   :: sqomega_qlimit = 0.00d0   !! for lower q-values automatically s(q,t) is computed
     double precision, parameter :: virtual_zero_t = 1d-6 !! t~=0 (in ns) to avoid potentialk zero divisions

!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'sqt_omeg'
       nparx =       2 + 4*n_strexpo + 1
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_sqt_omeg = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " combined s(q,omega) streched exponentials"
       th_citation(idesc)     = ""
!       --------------> set the parameter names --->

        parnam ( 1) = 'intensit'  
        parnam ( 2) = 'omega0  '  
        np = 2
        do i = 1, n_strexpo 
           np = np + 1
           write(buf,'("amp",i0)') i
           parnam(np) = buf
           np = np + 1
           write(buf,'("tau",i0)') i
           parnam(np) = buf
           np = np + 1
           write(buf,'("nue",i0)') i
           parnam(np) = buf
           np = np + 1
           write(buf,'("beta",i0)') i
           parnam(np) = buf
        enddo
        parnam(2 + 4*n_strexpo + 1) = "u_sqr   "


                                       
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "offset of frequency zero   "
        th_param_desc( 3,idesc) = "amplitude of n-th streched exp   "
        th_param_desc( 4,idesc) = "tau = tau * q**nue  of n-th streched exp   "
        th_param_desc( 5,idesc) = "nue = q-nue  of n-th streched exp (give 0 for no q-dep.)  "

! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q-value    default value"
        th_file_param(  2,idesc) = "nse      > explict switch between sqt and sqomega (if nse=1 ==> sqt)"
        th_file_param(  3,idesc) = "bk1level > resolution description bgr-levwl"
        th_file_param(  4,idesc) = "bk1slope > resolution description slope of bgr"
        th_file_param(  5,idesc) = "ga1inten > resolution description 1. gaussian amplitude"
        th_file_param(  6,idesc) = "ga1width > resolution description 1. gaussian width"
        th_file_param(  7,idesc) = "ga1cente > resolution description 1. gaussian center"
        th_file_param(  8,idesc) = "ga2inten > resolution description 2. gaussian amplitude"
        th_file_param(  9,idesc) = "ga2width > resolution description 2. gaussian width"
        th_file_param( 10,idesc) = "ga2cente > resolution description 2. gaussian center"
        th_file_param( 11,idesc) = "ga3inten > resolution description 3. gaussian amplitude"
        th_file_param( 12,idesc) = "ga3width > resolution description 3. gaussian width"
        th_file_param( 13,idesc) = "ga3cente > resolution description 3. gaussian center"
        th_file_param( 14,idesc) = "ga4inten > resolution description 4. gaussian amplitude"
        th_file_param( 15,idesc) = "ga4width > resolution description 4. gaussian width"
        th_file_param( 16,idesc) = "ga4cente > resolution description 4. gaussian center"
        th_file_param( 17,idesc) = "ga5inten > resolution description 5. gaussian amplitude"
        th_file_param( 18,idesc) = "ga5width > resolution description 5. gaussian width"
        th_file_param( 19,idesc) = "ga5cente > resolution description 5. gaussian center"
        th_file_param( 20,idesc) = "ga6inten > resolution description 6. gaussian amplitude"
        th_file_param( 21,idesc) = "ga6width > resolution description 6. gaussian width"
        th_file_param( 22,idesc) = "ga6cente > resolution description 6. gaussian center"
        th_file_param( 23,idesc) = "ga7inten > resolution description 7. gaussian amplitude"
        th_file_param( 24,idesc) = "ga7width > resolution description 7. gaussian width"
        th_file_param( 25,idesc) = "ga7cente > resolution description 7. gaussian center"
        th_file_param( 26,idesc) = "ga8inten > resolution description 8. gaussian amplitude"
        th_file_param( 27,idesc) = "ga8width > resolution description 8. gaussian width"
        th_file_param( 28,idesc) = "ga8cente > resolution description 8. gaussian center"
        th_file_param( 28,idesc) = "xwidth  the channel width (close to omega=0) "
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
! 
        th_sqt_omeg = 0.0
 
        RETURN
     ENDIF
!

      intensit      = pa( 1)
      omega0        = pa( 2)
      np = 2
      do i = 1, n_strexpo 
         np = np + 1
         amp(i) =     abs(pa(np))
         np = np + 1
         tau(i) =     abs(pa(np))
         np = np + 1
         nue(i) =     pa(np)  
         np = np + 1
         beta(i) =    abs(pa(np))    
      enddo
      
      u_sqr      =  (pa(2 + 4*n_strexpo + 1))



! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: q-value    default value
      xh =      0
      call parget('q       ',xh,iadda,ier)
      q        = xh
! >>> extract: explict switch between sqt and sqomega (if nse=1 ==> sqt)
      xh =      0
      call parget('nse     ',xh,iadda,ier)
      nse      = xh
! >>> extract: explict switch between sqt and sqomega (if nse=1 ==> sqt)
      xh =      0.001
      call parget('_xwidth ',xh,iadda,ier)
      xwidth    = xh
! >>> extract: resolution description bgr-levwl
      xh =      0
      call parget('bk1level',xh,iadda,ier)
      bk1level = xh
! >>> extract: resolution description slope of bgr
      xh =      0
      call parget('bk1slope',xh,iadda,ier)
      bk1slope = xh
! >>> extract: resolution description 1. gaussian amplitude
      xh =      0
      call parget('ga1inten',xh,iadda,ier)
      ga1inten = xh
! >>> extract: resolution description 1. gaussian width
      xh =      0
      call parget('ga1width',xh,iadda,ier)
      ga1width = xh
! >>> extract: resolution description 1. gaussian center
      xh =      0
      call parget('ga1cente',xh,iadda,ier)
      ga1cente = xh
! >>> extract: resolution description 2. gaussian amplitude
      xh =      0
      call parget('ga2inten',xh,iadda,ier)
      ga2inten = xh
! >>> extract: resolution description 2. gaussian width
      xh =      0
      call parget('ga2width',xh,iadda,ier)
      ga2width = xh
! >>> extract: resolution description 2. gaussian center
      xh =      0
      call parget('ga2cente',xh,iadda,ier)
      ga2cente = xh
! >>> extract: resolution description 3. gaussian amplitude
      xh =      0
      call parget('ga3inten',xh,iadda,ier)
      ga3inten = xh
! >>> extract: resolution description 3. gaussian width
      xh =      0
      call parget('ga3width',xh,iadda,ier)
      ga3width = xh
! >>> extract: resolution description 3. gaussian center
      xh =      0
      call parget('ga3cente',xh,iadda,ier)
      ga3cente = xh
! >>> extract: resolution description 4. gaussian amplitude
      xh =      0
      call parget('ga4inten',xh,iadda,ier)
      ga4inten = xh
! >>> extract: resolution description 4. gaussian width
      xh =      0
      call parget('ga4width',xh,iadda,ier)
      ga4width = xh
! >>> extract: resolution description 4. gaussian center
      xh =      0
      call parget('ga4cente',xh,iadda,ier)
      ga4cente = xh
! >>> extract: resolution description 5. gaussian amplitude
      xh =      0
      call parget('ga5inten',xh,iadda,ier)
      ga5inten = xh
! >>> extract: resolution description 5. gaussian width
      xh =      0
      call parget('ga5width',xh,iadda,ier)
      ga5width = xh
! >>> extract: resolution description 5. gaussian center
      xh =      0
      call parget('ga5cente',xh,iadda,ier)
      ga5cente = xh
! >>> extract: resolution description 6. gaussian amplitude
      xh =      0
      call parget('ga6inten',xh,iadda,ier)
      ga6inten = xh
! >>> extract: resolution description 6. gaussian width
      xh =      0
      call parget('ga6width',xh,iadda,ier)
      ga6width = xh
! >>> extract: resolution description 6. gaussian center
      xh =      0
      call parget('ga6cente',xh,iadda,ier)
      ga6cente = xh
! >>> extract: resolution description 7. gaussian amplitude
      xh =      0
      call parget('ga7inten',xh,iadda,ier)
      ga7inten = xh
! >>> extract: resolution description 7. gaussian width
      xh =      0
      call parget('ga7width',xh,iadda,ier)
      ga7width = xh
! >>> extract: resolution description 7. gaussian center
      xh =      0
      call parget('ga7cente',xh,iadda,ier)
      ga7cente = xh
! >>> extract: resolution description 8. gaussian amplitude
      xh =      0
      call parget('ga8inten',xh,iadda,ier)
      ga8inten = xh
! >>> extract: resolution description 8. gaussian width
      xh =      0
      call parget('ga8width',xh,iadda,ier)
      ga8width = xh
! >>> extract: resolution description 8. gaussian center
      xh =      0
      call parget('ga8cente',xh,iadda,ier)
      ga8cente = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!?! write(*,*)"T0: nse=",nse,ga1inten
!?!
!?!    do i=1,n_strexpo
!?!      write(*,*)"T a tau nue beta:",i, amp(i), tau(i), nue(i), beta(i)
!?!    enddo

      dwf        = exp(-u_sqr*q*q/3.0d0)



!     if( nse == 0 .and. is_in_xaxis("omega",iadda)) then
     if( nse == 0 ) then
        gampli(1)   = ga1inten
        gwidth(1)   = ga1width
        gcenter(1)  = ga1cente

        gampli(2)   = ga2inten
        gwidth(2)   = ga2width
        gcenter(2)  = ga2cente

        gampli(3)   = ga3inten
        gwidth(3)   = ga3width
        gcenter(3)  = ga3cente

        gampli(4)   = ga4inten
        gwidth(4)   = ga4width
        gcenter(4)  = ga4cente

        gampli(5)   = ga5inten
        gwidth(5)   = ga5width
        gcenter(5)  = ga5cente

        gampli(6)   = ga6inten
        gwidth(6)   = ga6width
        gcenter(6)  = ga6cente

        gampli(7)   = ga7inten
        gwidth(7)   = ga7width
        gcenter(7)  = ga7cente

        gampli(8)   = ga8inten
        gwidth(8)   = ga8width
        gcenter(8)  = ga8cente

        o0         = x - omega0

        rsum = 0
drs:    do i=1,size(gampli)
!?! write(*,*)"Tdrs=",i,gampli(i)
         if(gampli(i) == 0d0) cycle
         Omega     = o0 - gcenter(i)
         str_delta = abs(gwidth(i))
         a         = 0d0
         b         = 9.0d0/str_delta
         res       = adapint(fqt_kernel,a,b,epsilon,maxit,erraccu)*2
!?! write(*,*)"Tdrs2=",i,rsum

         rsum      = rsum + gampli(i)*res/(2*Pi)*sqrt(Pi)
        enddo drs

        th = intensit * rsum / fqt(virtual_zero_t) * dwf

!?! write(*,*)"Tomx:",x,rsum,th

    else
        t  = x
        th = intensit  * fqt(t)  / fqt(virtual_zero_t) * dwf

!?! write(*,*)"Tnse:",x,th

    endif
     th_sqt_omeg = th
 
! ---- writing computed parameters to the record >>>  
!      call parset('tau_q   ',sngl(tau),iadda,ier)

 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 

       function fqt_kernel(t)
!      -----------------------
       implicit none
       double precision, intent(in) :: t
       double precision             :: fqt_kernel

       fqt_kernel = fqt(t) * exp(-1d0/4d0*(str_delta*t)**2) * cos(t*Omega)*str_delta

       fqt_kernel = fqt(t) * &
                    exp(-1d0/4d0*(str_delta*t)**2) *    &
                    (sin(-t*Omega+0.5d0*t*xwidth)+sin(t*Omega+0.5d0*t*xwidth))/ &   !! this replaces cos(t*Omega)
                    (t*xwidth) * str_delta                                         !! in order to yield the integral

       end function fqt_kernel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! die eigentliche Zeitfunktion                                            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       function fqt(t)
!      ---------------
       implicit none
       double precision, intent(in) :: t
       double precision             :: fqt

       integer          :: j
       double precision :: tau_q

       fqt = 0d0
       do j=1,n_strexpo
         if(amp(j) == 0d0) cycle
         tau_q = tau(j)*(q**nue(j))
         fqt = fqt + amp(j) * exp( -(abs(t)/tau_q)**beta(j))
       enddo

       end function fqt

 end function th_sqt_omeg
