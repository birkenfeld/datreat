 FUNCTION th_sqt_om1(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  combined s(q,omega) and s(q,t) data from an s(q,t) model
! 
      use theory_description 
      implicit none 
      real    :: th_sqt_om1
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
     
! the internal parameter representation 
     double precision :: intensit   ! prefactor                   
     double precision :: omega0     ! omega scale zero shift    
     double precision :: u_sqr      ! <u**2> value for Debye-Waller-Factor
     double precision :: r_jump     ! jump length (if applicable)   
     double precision :: tau0_inc   ! tau_inc = tau_qinc * q**nue_inc + tau0_inc       
     double precision :: tauq_inc   ! tau_inc = tau_qinc * q**nue_inc + tau0_inc       
     double precision :: nue_inc    ! tau_inc = tau_qinc * q**nue_inc + tau0_inc   
     double precision :: bet0_inc   ! beta_inc   = bet_qinc * q**mue_inc + bet0_inc 
     double precision :: bet_qinc   ! beta_inc   = bet_qinc * q**mue_inc + bet0_inc 
     double precision :: mue_inc    ! beta_inc   = bet_qinc * q**mue_inc + bet0_inc 
     double precision :: tau0_coh   ! tau_coh = tau_qcoh * q**nue_coh + tau0_coh       
     double precision :: tauq_coh   ! tau_coh = tau_qcoh * q**nue_coh + tau0_coh       
     double precision :: nue_coh    ! tau_coh = tau_qcoh * q**nue_coh + tau0_coh   
     double precision :: bet0_coh   ! beta_coh   = bet_qcoh * q**mue_coh + bet0_coh 
     double precision :: bet_qcoh   ! beta_coh   = bet_qcoh * q**mue_coh + bet0_coh 
     double precision :: mue_coh    ! beta_coh   = bet_qcoh * q**mue_coh + bet0_coh 
     double precision :: a1_0       ! A1 = a1_0  +  a1_qnue * q**nue_a1 + a2_qnue * q**nue_a2
     double precision :: a1_qnue    ! A1 = a1_0  +  a1_qnue * q**nue_a1 + a2_qnue * q**nue_a2
     double precision :: nue_a1     ! A1 = a1_0  +  a1_qnue * q**nue_a1 + a2_qnue * q**nue_a2
     double precision :: a2_qnue    ! A1 = a1_0  +  a1_qnue * q**nue_a1 + a2_qnue * q**nue_a2
     double precision :: nue_a2     ! A1 = a1_0  +  a1_qnue * q**nue_a1 + a2_qnue * q**nue_a2
     double precision :: Bcoh0      ! Bcoh = Bcoh0 +  b1_qnue * q**nue_b1 + b2_qnue * q**nue_b2
     double precision :: b1_qnue    ! Bcoh = Bcoh0 +  b1_qnue * q**nue_b1 + b2_qnue * q**nue_b2
     double precision :: nue_b1     ! Bcoh = Bcoh0 +  b1_qnue * q**nue_b1 + b2_qnue * q**nue_b2
     double precision :: b2_qnue    ! Bcoh = Bcoh0 +  b1_qnue * q**nue_b1 + b2_qnue * q**nue_b2
     double precision :: nue_b2     ! Bcoh = Bcoh0 +  b1_qnue * q**nue_b1 + b2_qnue * q**nue_b2


     double precision :: epsilon = 1d-9
                                           
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
     double precision :: tau_inc    !                              
     double precision :: tau_coh    !    

     double precision :: beta_inc   !                              
     double precision :: beta_coh   !

     double precision :: EISF       !  
     double precision :: A1         !  
     double precision :: Bcoh       !

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

!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'sqt_om1'
       nparx =       26
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_sqt_om1 = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " combined s(q,omega) and s(q,t) data from an s(q,t) model"
       th_citation(idesc)     = ""
!       --------------> set the parameter names --->

        parnam ( 1) = 'intensit'  
        parnam ( 2) = 'omega0  '  
        parnam ( 3) = 'u_sqr   '  
        parnam ( 4) = 'r_jump  '  
        parnam ( 5) = 'tau0_inc'  
        parnam ( 6) = 'tauq_inc'  
        parnam ( 7) = 'nue_inc '  
        parnam ( 8) = 'bet0_inc'  
        parnam ( 9) = 'bet_qinc'  
        parnam (10) = 'mue_inc '  
        parnam (11) = 'tau0_coh'  
        parnam (12) = 'tauq_coh'  
        parnam (13) = 'nue_coh '  
        parnam (14) = 'bet0_coh'  
        parnam (15) = 'bet_qcoh'  
        parnam (16) = 'mue_coh '  
        parnam (17) = 'a1_0    '  
        parnam (18) = 'a1_qnue '  
        parnam (19) = 'nue_a1  '  
        parnam (20) = 'a2_qnue '  
        parnam (21) = 'nue_a2  '  
        parnam (22) = 'Bcoh0   '  
        parnam (23) = 'b1_qnue '  
        parnam (24) = 'nue_b1  '  
        parnam (25) = 'b2_qnue '  
        parnam (26) = 'nue_b2  ' 


                                       
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "offset of frequency zero   "
        th_param_desc( 3,idesc) = "<u**2> for Debye-Waller factor    "
        th_param_desc( 4,idesc) = "jumpdistance for EISF compuation    "
        th_param_desc( 5,idesc) = "q-independnt part of the incoherent streched exp-tau  "
        th_param_desc( 6,idesc) = "q-dpendent part (factor to q**nue_inc) of tau incoheren1 "
        th_param_desc( 7,idesc) = "q-exponent of the tau_inc q-dependencen "
        th_param_desc( 8,idesc) = "constant part of streching exponent beta of the incoherent sterched exp. sqt  "
        th_param_desc( 9,idesc) = "q-dependent part of beta_inc betainc =bet0_inc+bet_qinc*q*nue_inc "
        th_param_desc(10,idesc) = "q_exponent of beta_inc  "                                                       
        th_param_desc(11,idesc) = "q-independnt part of the coherent streched exp-tau  "
        th_param_desc(12,idesc) = "q-dpendent part (factor to q**nue_inc) of tau coherent "
        th_param_desc(13,idesc) = "q-exponent of the tau_coh q-dependencen "
        th_param_desc(14,idesc) = "constant part of streching exponent beta of the cohoherent sterched exp. sqt  "
        th_param_desc(15,idesc) = "q-dependent part of beta_coh betacoh =bet0_coh+bet_qcoh*q*nue_coh "  
        th_param_desc(16,idesc) = "q_exponent of beta_coh  "                                          
        th_param_desc(17,idesc) = "constant part of the A1 factor     "
        th_param_desc(18,idesc) = "coefficient21 of q-dependence of A1  "
        th_param_desc(19,idesc) = "q-exponent1 of A1 "
        th_param_desc(20,idesc) = "coefficient 2 of A1 =a1_0+a1_qnue*q**nue_a1+a2_nue*q**nue_a2 "
        th_param_desc(21,idesc) = "q_exponent 2 for A1   "
        th_param_desc(22,idesc) = "constant part of coherent inetnsity B    "
        th_param_desc(23,idesc) = "q-coefficient 1 of Bcoh   "
        th_param_desc(24,idesc) = "q_exponent 1 of Bcoh   "
        th_param_desc(25,idesc) = "q_coefficient 2 of Bcoh   "
        th_param_desc(26,idesc) = "q-exponent 2 of Bcoh    "







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
        th_out_param(  1,idesc) = "tau0_q   > scaled: str_tau0 * (1+jlen*qz**qexpt0)/(qz**qexpt0)"
        th_out_param(  2,idesc) = "beta_q   > scaled: str_beta * (qz**qexpbeta + beta0)"
! 
        th_sqt_om1 = 0.0
 
        RETURN
     ENDIF
!

      intensit      = pa( 1)
      omega0        = pa( 2)
      u_sqr         = abs(pa( 3))
      r_jump        = abs(pa( 4))
      tau0_inc      = pa( 5)
      tauq_inc      = pa( 6)
      nue_inc       = pa( 7)
      bet0_inc      = pa( 8)
      bet_qinc      = pa( 9)
      mue_inc       = pa(10)
      tau0_coh      = pa(11)  
      tauq_coh      = pa(12)  
      nue_coh       = pa(13)  
      bet0_coh      = pa(14)
      bet_qcoh      = pa(15)
      mue_coh       = pa(16)
      a1_0          = pa(17)
      a1_qnue       = pa(18)
      nue_a1        = pa(19)
      a2_qnue       = pa(20)
      nue_a2        = pa(21)
      Bcoh0         = pa(22)
      b1_qnue       = pa(23)
      nue_b1        = pa(24)
      b2_qnue       = pa(25)
      nue_b2        = pa(26)

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
!! some generic q-dependencies, to be improved according to model
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      tau_inc    = tauq_inc * q**nue_inc + tau0_inc
      tau_coh    = tauq_coh * q**nue_coh + tau0_coh
      
      beta_inc   = bet_qinc * q**mue_inc + bet0_inc 
      beta_coh   = bet_qcoh * q**mue_coh + bet0_coh

     
      EISF       = (3d0  * besj1(q*r_jump) / (q * r_jump) )**2
      A1         = a1_0  +  a1_qnue * q**nue_a1 + a2_qnue * q**nue_a2
      Bcoh       = Bcoh0 +  b1_qnue * q**nue_b1 + b2_qnue * q**nue_b2

!! add hier S(Q,inf) and(?) S(Q,0) of crosslinker blobs
!!    ScrossInf =
!!    Scross0   =

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     if( nse == 0 .or. q > sqomega_qlimit .or. .not. is_in_xaxis("ns",iadda)) then
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
         if(gampli(i) == 0d0) cycle
         Omega     = o0 - gcenter(i)
         str_delta = abs(gwidth(i))
         a         = 0d0
         b         = 9.0d0/str_delta
         res       = adapint(fqt_kernel,a,b,epsilon,maxit,erraccu)*2
         rsum      = rsum + gampli(i)*res/(2*Pi)*sqrt(Pi)
        enddo drs

        dwf        = exp(-u_sqr*q*q/3.0d0)

        th = intensit * dwf * rsum    

! write(*,*)"Tomx:",x,dwf,rsum,th

    else
        t  = x
        dwf        = exp(-u_sqr*q*q/3.0d0)
        th = intensit * dwf * fqt(t) 
! write(*,*)"Tnse:",x,dwf,th

    endif
     th_sqt_om1 = th
 
! ---- writing computed parameters to the record >>>  
      call parset('tau_inc   ',sngl(tau_inc),iadda,ier)
      call parset('tau_coh   ',sngl(tau_coh),iadda,ier)
      call parset('beta_inc  ',sngl(beta_inc),iadda,ier)
      call parset('beta_coh  ',sngl(beta_coh),iadda,ier)
      call parset('EISF      ',sngl(EISF),iadda,ier)
      call parset('A1        ',sngl(A1),iadda,ier)
      call parset('Bcoh      ',sngl(Bcoh),iadda,ier)
 
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

!       fqt =  exp(-(t/str_tau0)**str_beta)


       fqt= (EISF+(1-EISF)*exp(-(t/tau_inc)**beta_inc)*A1+(1-A1)*exp(-(t/tau_coh)**beta_coh))+Bcoh

       end function fqt

! 
! 
! schnittkugel1[q_, 
!   d_] := (4*Pi/3*(2*Pi/q)^3 - Pi*(2*Pi/q)^2*d + Pi/12*d^3)/(4*
!     Pi/3*(2*Pi/q)^3)
!        end function Bcoh
! 
! 

 end function th_sqt_om1
