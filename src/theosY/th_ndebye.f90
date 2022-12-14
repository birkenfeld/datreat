 FUNCTION th_ndebye(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  random chain structure factor P(Q) (P(Q==)=1) by direct summation over N segments chain statistics is expressed ny nu (nu=0.5 ==> Gaussian chain) in addition a chi parameter my be set: S(Q) = 1/(1/(N*P(Q) + 2*chi)
!  J.K. Kjems, T. Freloft, D. Richter and S.K. Sinha, Physica 136B (1986) 285-290
      use theory_description 
      implicit none 
      real    :: th_ndebye
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
     double precision :: n          ! number of segments (do not fit, its integer)                                    
     double precision :: l          ! effective segment length Rg = l * N**nu /sqrt(6)                                
     double precision :: nu         ! expansion parameter (Gaussian random walk =0.5)                                 
     double precision :: chi        ! chi parameter                                                                   
 
! the reout parameter representation 
     double precision :: rg         ! radius of gyration                                                              
     double precision :: re         ! end-end radius of gyration                                                              
 
     double precision :: th
 
   double precision :: q, pq, sq
   integer          :: ni
!
   double precision :: sq_molecule_dendrimer
   integer          :: withsq


! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'ndebye'
       nparx =        6
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_ndebye = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " random chain structure factor P(Q) (P(Q==)=1) by direct summation over N segments chain statistics is expressed ny nu (nu=0.5 ==> Gaussian chain) in addition a chi parameter my be set: S(Q) = 1/(1/(N*P(Q) + 2*chi)"
       th_citation(idesc)     = " J.K. Kjems, T. Freloft, D. Richter and S.K. Sinha, Physica 136B (1986) 285-290"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'n       '  ! number of segments (do not fit, its integer)                                    
        parnam ( 3) = 'l       '  ! effective segment length Rg = l * N**nu /sqrt(6)                                
        parnam ( 4) = 'nu      '  ! expansion parameter (Gaussian random walk =0.5)                                 
        parnam ( 5) = 'chi     '  ! chi parameter                                                                   
        parnam ( 6) = 'withsq  '  ! switch sq of molecule
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "number of segments (do not fit, its integer)" !//cr//parspace//&
        th_param_desc( 3,idesc) = "effective segment length Rg = l * N**nu /sqrt(6)" !//cr//parspace//&
        th_param_desc( 4,idesc) = "expansion parameter (Gaussian random walk =0.5)" !//cr//parspace//&
        th_param_desc( 5,idesc) = "chi parameter" !//cr//parspace//&
        th_param_desc( 6,idesc) = "switch on/off (1/0) the molecular formfactor of sq_dendrimer" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "         > "
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "rg       > radius of gyration for nu .ne. 1/2 this differs from Re/sqrt(6)"
        th_out_param(  2,idesc) = "re       > end-to-end radius"
! 
        th_ndebye = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      n        =      pa( 2)
      l        =  abs(pa( 3))
      nu       =  abs(pa( 4))
      chi      =      pa( 5)
      withsq   = nint(pa( 6))
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: 
      xh = n 
      call parget('nseg     ',xh,iadda,ier)
      n         = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     q   = x
     ni  = nint(n)
     re  =  l * dble(ni)**nu 
     pq  = nndebye(q, l, nu, ni, rg )
     sq  = 1d0 / ( 1d0/(ni*pq) + 2d0 * chi )
     th  = ampli * sq
     if(withsq == 1) then
        th = th * sq_molecule_dendrimer(q)
     endif

     th_ndebye = th
 
! ---- writing computed parameters to the record >>>  
      call parset('rg      ',sngl(rg),iadda,ier)
      call parset('re      ',sngl(re),iadda,ier)
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 


function nndebye(q, l, nue,n ,Rg ) result(val)
   implicit none
   double precision ,intent(in) :: q           ! Q-value
   double precision ,intent(in) :: l           ! effective segment length
   double precision ,intent(in) :: nue         ! exponent
   integer          ,intent(in) :: n           ! number of segments
   double precision ,intent(out):: Rg          ! radius of gyration
   double precision             :: val

   integer          :: i, j
   double precision :: eterms(0:n-1)
   double precision :: ql2, mu

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

!!! ---- consistent Rg determination ----
!$OMP PARALLEL DO
         do i=0,n-1
           eterms(i) = dble(i)**mu
         enddo
!$OMP END PARALLEL DO

         Rg = 0
!$OMP PARALLEL DO REDUCTION(+:Rg)
         do i=1,n
          do j=1,n
            Rg = Rg + eterms(abs(i-j))
          enddo
         enddo
!$OMP END PARALLEL DO
        Rg    = sqrt( (0.5d0 * l**2) * Rg/(n**2) )
 

end function nndebye



 end function th_ndebye
