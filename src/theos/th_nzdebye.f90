 FUNCTION th_nzdebye(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  two-setp approximation with Gaussian walks bolb(N,l) * coil of blobs(Z,sqrt(N)*l)
! 
      use theory_description 
      implicit none 
      real    :: th_nzdebye
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
     double precision :: n          ! number of segments in blob                                                      
     double precision :: z          ! number of blobs in chain                                                        
     double precision :: l          ! effective segment length                                                        
     double precision :: nu         ! exponent                                                                        
! the recin parameter representation 
                                                                           
! the reout parameter representation 
     double precision :: rg         ! radius of gyration                                                              
 
     double precision :: th
 
   double precision :: pq, q
   integer          :: ni, nz
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'nzdebye'
       nparx =        5
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_nzdebye = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " two-setp approximation with Gaussian walks bolb(N,l) * coil of blobs(Z,sqrt(N)*l)"
       th_citation(idesc)     = ""
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'n       '  ! number of segments in blob                                                      
        parnam ( 3) = 'z       '  ! number of blobs in chain                                                        
        parnam ( 4) = 'l       '  ! effective segment length                                                        
        parnam ( 5) = 'nu      '  ! exponent                                                                        
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "number of segments in blob" !//cr//parspace//&
        th_param_desc( 3,idesc) = "number of blobs in chain" !//cr//parspace//&
        th_param_desc( 4,idesc) = "effective segment length" !//cr//parspace//&
        th_param_desc( 5,idesc) = "exponent" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "         > "
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "rg       > radius of gyration"
! 
        th_nzdebye = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      n        =      pa( 2)
      z        =      pa( 3)
      l        =      pa( 4)
      nu       =      pa( 5)
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
     q   = x
     ni  = nint(n)
     nz  = nint(z)
     rg  =  l * dble(ni*nz)**nu / sqrt(6d0)
     pq  = nndebye(q, l, nu, ni ) * nndebye(q, l*sqrt(n), nu, nz )

     th_nzdebye = pq
 
! ---- writing computed parameters to the record >>>  
      call parset('rg      ',sngl(rg),iadda)
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 


function nndebye(q, l, nue,n ) result(val)
   implicit none
   double precision ,intent(in) :: q           ! Q-value
   double precision ,intent(in) :: l           ! effective segment length
   double precision ,intent(in) :: nue         ! exponent
   integer          ,intent(in) :: n           ! number of segments
   double precision             :: val

   integer          :: i, j
   double precision :: eterms(0:n-1)
   double precision :: ql2, mu, Rg

   ql2 = (1d0/6d0) * (q*l)**2
   mu  = 2d0 * nue
!  Rg  = sqrt((l*l*dble(n)**mu) / 6d0) )
!$OMP PARALLEL DO
   do i=0,n-1
     eterms(i) = exp(-ql2*dble(i)**mu)
   enddo
!$OMP END PARALLEL DO
   val = 0
!$OMP PARALLEL DO REDUCTION(+:val)
   do i=1,n
    do j=1,n
      val = val + eterms(abs(i-j))
    enddo
   enddo
!$OMP END PARALLEL DO

   val = val / (n*n)

end function nndebye



 end function th_nzdebye
