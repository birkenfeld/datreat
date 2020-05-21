 FUNCTION th_ringform(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  ring chain structure factor P(Q) (P(Q==)=1) by direct summation over N segments with change statistics as a function of distance along the chain 	chain statistics is expressed by nu (nu=0.5 ==> Gaussian chain) 	the closure is taken Bensafi et al. in addition a chi parameter may be set: S(Q) = 1/(1/(N*P(Q) + 2*chi) or left out
!  J.K. Kjems, T. Freloft, D. Richter and S.K. Sinha, Physica 136B (1986) 285-290
      use theory_description 
      implicit none 
      real    :: th_ringform
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
     double precision :: n1         ! Distance along chain where Gaussian statistics is valid                         
     double precision :: l          ! effective segment length Rg = l * N**nu /sqrt(6)                                
     double precision :: nu         ! expansion parameter (Gaussian random walk =0.5)                                 
     double precision :: nuwidth    ! width of the exponent transition                                                
     double precision :: chi        ! chi parameter                                                                   
! the recin parameter representation 
     double precision :: conc       ! monomer conc.                                                                   
! the reout parameter representation 
     double precision :: rg         ! radius of gyration                                                              
 
     double precision :: th
 
   double precision :: pq, sq, q
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'ringform'
       nparx =        7
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_ringform = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " ring chain structure factor P(Q) (P(Q==)=1) by direct summation over N segments with change statistics as a function of distance along the chain 	chain statistics is expressed by nu (nu=0.5 ==> Gaussian chain) 	the closure is taken Bensafi et al. in addition a chi parameter may be set: S(Q) = 1/(1/(N*P(Q) + 2*chi) or left out"
       th_citation(idesc)     = " J.K. Kjems, T. Freloft, D. Richter and S.K. Sinha, Physica 136B (1986) 285-290"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'n       '  ! number of segments (do not fit, its integer)                                    
        parnam ( 3) = 'n1      '  ! Distance along chain where Gaussian statistics is valid                         
        parnam ( 4) = 'l       '  ! effective segment length Rg = l * N**nu /sqrt(6)                                
        parnam ( 5) = 'nu      '  ! expansion parameter (Gaussian random walk =0.5)                                 
        parnam ( 6) = 'nuwidth '  ! width of the exponent transition                                                
        parnam ( 7) = 'chi     '  ! chi parameter                                                                   
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "number of segments (do not fit, its integer)" !//cr//parspace//&
        th_param_desc( 3,idesc) = "Distance along chain where Gaussian statistics is valid" !//cr//parspace//&
        th_param_desc( 4,idesc) = "effective segment length Rg = l * N**nu /sqrt(6)" !//cr//parspace//&
        th_param_desc( 5,idesc) = "expansion parameter (Gaussian random walk =0.5)" !//cr//parspace//&
        th_param_desc( 6,idesc) = "width of the exponent transition" !//cr//parspace//&
        th_param_desc( 7,idesc) = "chi parameter" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "conc     > monomer conc."
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "rg       > radius of gyration"
! 
        th_ringform = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      n        =      pa( 2)
      n1       =      pa( 3)
      l        =      pa( 4)
      nu       =      pa( 5)
      nuwidth  =      pa( 6)
      chi      =      pa( 7)
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: monomer conc.
      xh =     0.1
      call parget('conc    ',xh,iadda,ier)
      conc     = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     q   = x
     pq  = nring(q, l, nu, nuwidth, n, n1)
     sq  = 1d0 / ( 1d0/(n*pq) + 2d0 * chi )
     th  = ampli * sq


     th_ringform = th
 
! ---- writing computed parameters to the record >>>  
      call parset('Rg      ',sngl(rg),iadda,ier)
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 


function nring(q, l, nu, nuwidth, n, n1 ) result(val)
   implicit none
   double precision ,intent(in) :: q           ! Q-value
   double precision ,intent(in) :: l           ! effective segment length
   double precision ,intent(in) :: nu          ! exponent
   double precision ,intent(in) :: nuwidth
   double precision ,intent(in) :: n           ! number of segments
   double precision ,intent(in) :: n1
   double precision             :: val

   integer          :: i, j
   double precision :: eterms(0:nint(n)-1)
   double precision :: rterms(0:nint(n)-1)
   double precision :: Rg0, rr, mu=0.5d0


!!$OMP PARALLEL DO PRIVATE(rr)   
 do i=0,nint(n)-1
  rr = l**2 * (i*(1d0-i/N))**(nucross(dble(i),n1,nuwidth,nu,mu)*2)
  rterms(i) = 0.5d0 * rr 
  eterms(i) = exp(-q*q*rr/6d0)
 enddo
 val = nint(n) * eterms(0)
 rg0 = nint(n) * rterms(0)
!!$OMP PARALLEL DO REDUCTION(+:val,rg0)
   do i=1,nint(n)-1
      val = val + eterms(i)*(nint(n)-i) * 2
      rg0 = rg0 + rterms(i)*(nint(n)-1) * 2
   enddo
!!$OMP END PARALLEL DO

   val = val / (nint(n)**2)
   rg  = sqrt(rg0 / (nint(n)**2))

end function nring

function nucross(x,x0,width,nu,mu) result(nuc)
  double precision, intent(in) :: x, x0, width, nu, mu
  double precision :: nuc
  nuc = fermi(x,x0,width)*mu + (1d0-  fermi(x,x0,width))*nu
end function nucross


function fermi(x,x0,width) result(ff)
  double precision, intent(in) :: x, x0, width
  double precision :: ff
  ff = 1d0/(1d0+exp((x-x0)/width))
  
end function fermi



 end function th_ringform
