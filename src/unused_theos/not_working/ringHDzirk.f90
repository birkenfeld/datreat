 FUNCTION th_ringHDzirk(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  ring chain structure factor P(Q) (P(Q==)=1) by direct summation over N segments with change statistics as a function of distance along the chain 	chain statistics is expressed by nu (nu=0.5 ==> Gaussian chain) 	the closure is taken Bensafi et al. in addition a chi parameter may be set: S(Q) = 1/(1/(N*P(Q) + 2*chi) or left out
!  J.K. Kjems, T. Freloft, D. Richter and S.K. Sinha, Physica 136B (1986) 285-290
      use theory_description 
      implicit none 
      real    :: th_ringHDzirk
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

     double precision :: llin
     double precision :: nlin
     double precision :: nulin

                                         
     double precision :: chi, chi12, chi13, chi23
     double precision :: v11, v22, v12    

     double precision :: chi_aa   
     double precision :: chi_ab   
     double precision :: chi_aab   
     double precision :: chi_gamh 
     double precision :: chi_siga 
     double precision :: chi_sigb 
     double precision :: chi_siab 
     double precision :: chi_xi   



! the recin parameter representation 
     double precision :: conc       ! monomer conc.                                                                   
! the reout parameter representation 
     double precision :: rgring         ! radius of gyration   
     double precision :: rglin          ! radius of gyration   
    
     double precision :: phiH, phiD
                                                               
 
     double precision :: th
 
     double precision :: pq, pqlin, q
     double precision :: s11, s011, s022, s033
 
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'rihdzrpa'
       nparx =        17
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_ringHDzirk = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " ring chain structure factor P(Q) (P(Q==)=1) by direct summation over N segments with change statistics as a function of distance along the chain 	chain statistics is expressed by nu (nu=0.5 ==> Gaussian chain) 	the closure is taken Bensafi et al. in addition a chi parameter may be set: S(Q) = 1/(1/(N*P(Q) + 2*chi) or left out" //cr//parspace//&
       "h,d-ring mix melt in linear polymer with chi parameter, volfrac of h/d-ring must be a record param phiringh/d: " //cr//parspace//&
       " interaction:        chi   (A9) from Zirkel paper:A. Zirkel et al. Macromolecules 2002, 35, 7375-7386" 


       th_citation(idesc)     = " Hammouda, Multicomponent RPA, chapter 34 section 5 Eq. (16)"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'n       '  ! number of segments (do not fit, its integer)                                    
        parnam ( 3) = 'n1      '  ! Distance along chain where Gaussian statistics is valid                         
        parnam ( 4) = 'l       '  ! effective segment length Rg = l * N**nu /sqrt(6)                                
        parnam ( 5) = 'nu      '  ! expansion parameter (Gaussian random walk =0.5)                                 
        parnam ( 6) = 'nuwidth '  ! width of the exponent transition 
        parnam ( 7) = 'nlin    '  ! number of segments (do the linear component   
        parnam ( 8) = 'llin    '  ! effective segment length Rg = l * N**nu /sqrt(6)                            
        parnam ( 9) = 'nulin   '  ! expansion parameter (Gaussian random walk =0.5)   
        parnam (10) = 'chi_aa  '  ! chi parameter   a_A         
        parnam (11) = 'chi_ab  '  ! chi parameter   a_B 
        parnam (12) = 'chi_aab '  ! chi parameter   a_AB 
        parnam (13) = 'chi_gamh'  ! chi parameter   Gamma_h 
        parnam (14) = 'chi_siga'  ! chi parameter   sigma_a
        parnam (15) = 'chi_sigb'  ! chi parameter   sigma_b
        parnam (16) = 'chi_siab'  ! chi parameter   sigma_ab
        parnam (17) = 'chi_xi  '  ! chi parameter   xi          

! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "number of segments (do not fit, its integer)" !//cr//parspace//&
        th_param_desc( 3,idesc) = "Distance along chain where Gaussian statistics is valid" !//cr//parspace//&
        th_param_desc( 4,idesc) = "effective segment length Rg = l * N**nu /sqrt(6)" !//cr//parspace//&
        th_param_desc( 5,idesc) = "expansion parameter (Gaussian random walk =0.5)" !//cr//parspace//&
        th_param_desc( 6,idesc) = "width of the exponent transition" !//cr//parspace//&
        th_param_desc( 7,idesc) = "number of segments linear " !//cr//parspace//&
        th_param_desc( 8,idesc) = "segemnt length of linear" !//cr//parspace//&
        th_param_desc( 9,idesc) = "expansion factor of linear" !//cr//parspace//&
        th_param_desc( 10,idesc) = "chi parameter   a_A    "
        th_param_desc( 11,idesc) = "chi parameter   a_B "
        th_param_desc( 12,idesc) = "chi parameter   a_AB "
        th_param_desc( 13,idesc) = "chi parameter   Gamma_h "
        th_param_desc( 14,idesc) = "chi parameter   sigma_a"
        th_param_desc( 15,idesc) = "chi parameter   sigma_b"
        th_param_desc( 16,idesc) = "chi parameter   sigma_ab"
        th_param_desc( 17,idesc) = "chi parameter   xi"
! >>>>> describe record parameters       used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "phiringh     > h-ring volume fraction."
        th_file_param(  2,idesc) = "phiringd     > d-ring volume fraction."
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "rgring       > radius of gyration"
        th_out_param(  2,idesc) = "rglin        > radius of gyration"
! 
        th_ringHDzirk = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      n        = abs( pa( 2)) 
      n1       = abs( pa( 3))
      l        = abs( pa( 4))
      nu       = abs( pa( 5))
      nuwidth  = abs( pa( 6))
      nlin     = abs( pa( 7))
      llin     = abs( pa( 8))
      nulin    = abs( pa( 9))
      chi_aa   = pa(10)! chi parameter   a_A         
      chi_ab   = pa(11)! chi parameter   a_B 
      chi_aab  = pa(12)! chi parameter   a_AB 
      chi_gamh = pa(13)! chi parameter   Gamma_h 
      chi_siga = abs( pa(14))  ! chi parameter   sigma_a
      chi_sigb = abs( pa(15))  ! chi parameter   sigma_b
      chi_siab = abs( pa(16))  ! chi parameter   sigma_ab
      chi_xi   = abs( pa(17))  ! chi parameter   xi          


! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: monomer conc.
      call parget('phiringh',xh,iadda,ier)
      phiH = xh

      if(ier .ne. 0) then
       write(*,*)"WARNING: parameter phiringh missing in record parameters"
       write(*,*)"WARNING: using default phiringh = 0.01"
       phiH = 0.01d0
      endif

      call parget('phiringd',xh,iadda,ier)
      phiD     = xh

      if(phid > 1d-4) then
       write(*,*)"WARNING: parameter phiringd .ne. is inot supported !"
       phiD = 0.00d0
      endif
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     q      = x
     pq     = nring(q, l, nu, nuwidth, n, n1, rgring )
     pqlin  = nndebye(q, llin, nulin,nint(nlin) ,Rglin )

     s011   = phiH            * n    * pq                !! H-ring component=1
     s022   = phiD            * n    * pq                !! D-ring component=2
     s033   = (1d0-phiH-phiD) * nlin * pqlin             !! linear matrix background coponent=3

     chi    = 0.5d0 * ( chi_aa*(4*Pi)/3*chi_siga**3*(1-0.1d0*(q*chi_siga)**2) + &
                        chi_ab*(4*Pi)/3*chi_sigb**3*(1-0.1d0*(q*chi_sigb)**2) - &
                        chi_aab*(4*Pi)/3*chi_siab**3*(1-0.1d0*(q*chi_siab)**2)& 
              ) + &
             chi_gamh*exp(-chi_siab/chi_xi)/(1+(chi_xi*q)**2) * &
              (sin(q*chi_siab)/(q*chi_xi)+cos(q*chi_siab)) 
               

!!prev??     chi12  = 0
!!prev??     chi13  = chi
!!prev??     chi23  = chi
!!prev??
!!prev??     v11    = 1d0/s033         - 2d0 * chi13
!!prev??     v22    = 1d0/s033         - 2d0 * chi23
!!prev??     v12    = 1d0/s033 + chi12 - chi13 - chi23
!!prev??
!!prev??     s11    = (s011*(1d0+v22*s022)) / &
!!prev??              ( (1d0+v11*s011)*(1d0+v22*s022) - v12**2*s011*s022 )

!! Zirkel Paper (3.1)
     s11 =1d0 / ( 1d0/s011 + 1d0/s033 - 2d0 * chi )  !! eigentlich chi/v0 aber vo in den anderen faktore..?!

!write(*,*)"T",q,s011,s033,chi,s11

!     sq  = phiring / ( 1d0/(phiring*n*pq) + 1d0/((1-phiring)*nlin*pqlin)- 2d0 * chi )
     th  = ampli * s11

!!TEST!!
!T! write(*,'("T:",f8.4,2x,2f7.4,2x,2f6.3,2x,2f6.0,2x,9f11.5,2x,f11.5)') &
!T!         q,pq,pqlin,phiH,phiD,n,nlin,s011,s022,s033,v11,v22,v12,&
!T!         (1d0+v22*s022),(1d0+v11*s011)*(1d0+v22*s022),v12**2*s011*s022 , s11
!T!


     th_ringHDzirk = th
 
! ---- writing computed parameters to the record >>>  
      call parset('rgring  ',sngl(rgring),iadda,ier)
      call parset('rgrlin  ',sngl(rglin),iadda,ier)
 
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




function nring(q, l, nu, nuwidth, n, n1, Rg ) result(val)
   implicit none
   double precision ,intent(in) :: q           ! Q-value
   double precision ,intent(in) :: l           ! effective segment length
   double precision ,intent(in) :: nu          ! exponent
   double precision ,intent(in) :: nuwidth
   double precision ,intent(in) :: n           ! number of segments
   double precision ,intent(in) :: n1
   double precision ,intent(out):: Rg          ! radius of gyration
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
      rg0 = rg0 + rterms(i)*(nint(n)-i) * 2
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



 end function th_ringHDzirk
