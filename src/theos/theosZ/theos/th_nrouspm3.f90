 FUNCTION th_nrouspm3(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  Rouse by discrete summation with mode restriction
!  Rouse, Doi_Edwards
      use theory_description 
      implicit none 
      real    :: th_nrouspm3
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
     double precision :: amplitu    ! prefactor                                                                       
     double precision :: wl4        ! rouse rate                                                                      
     integer          :: n_segmen   ! effective N for summation                                                       
     double precision :: re         ! end-to-end radius of Gaussian coil                                              
     double precision :: temp       ! temperature if not in parameters                                                
     double precision :: com_diff   ! center of mass diffusion (if 0, use Rouse default)                              
 
! the recin parameter representation 
     double precision :: q          ! momentum transfer                                                               
                                                              
! the reout parameter representation 
     double precision :: Dr         ! diffusion   
     double precision :: W          ! rouse     W
     double precision :: l          ! sqgment length    
                                                                    
     double precision :: sqt, sq, t
     double precision :: labelarray(1000)
     double precision :: modeamp   (size(labelarray))

     double precision :: label1_s, label2_s, label3_s
     double precision :: label1_a, label2_a, label3_a, label4_a
     integer          :: nl1, nl2, nl3

     double precision, save :: labelstat(7)
     character(len=8) :: buf
     integer          :: i

     double precision :: rr, a_cross, nu_subdiff, r02

     double precision :: pfix      ! selects fixed chain spectrum if =-0.5 else 0=normal rouse
     integer          :: fixend
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'nrouspm3'
       nparx =        20
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_nrouspm3 = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " Rouse by discrete summation with 4 zone labelled chain and mode modification "//cr//parspace//&
                                " modeamplitudes amodx = arctan(amtanx)*2/Pi "

       th_citation(idesc)     = " Rouse, Doi_Edwards"
!       --------------> set the parameter names --->
        parnam ( 1) = 'amplitu '  ! prefactor                                                                       
        parnam ( 2) = 'wl4     '  ! rouse rate                                                                      
        parnam ( 3) = 'n_segmen'  ! effective N for summation                                                       
        parnam ( 4) = 're      '  ! end-to-end radius of Gaussian coil                                              
        parnam ( 5) = 'temp    '  ! temperature if not in parameters                                                
        parnam ( 6) = 'com_diff'  ! center of mass diffusion (if 0, use Rouse default)                              
        parnam ( 7) = 'nu_subd '  ! subdiff nu                             
        parnam ( 8) = 'r02     '  ! transition r**2 to com diff                            
        parnam ( 9) = 'a_cross '  ! cross oover a                              
        parnam (10) = 'amatan1   '  ! minimum mode to be included                     
        parnam (11) = 'amatan2   '  ! minimum mode to be included                     
        parnam (12) = 'amatan3   '  ! minimum mode to be included                     
        parnam (13) = 'amatan4   '  ! minimum mode to be included                     
        parnam (14) = 'amatan5   '  ! minimum mode to be included                     
        parnam (15) = 'amatan6   '  ! minimum mode to be included                     
        parnam (16) = 'amatan7   '  ! minimum mode to be included                     
        parnam (17) = 'amatan8   '  ! minimum mode to be included                     
        parnam (18) = 'amatan9   '  ! minimum mode to be included                     
        parnam (19) = 'amatan10  '  ! minimum mode to be included                     
        parnam (20) = 'amatan11  '  ! minimum mode to be included                     
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "rouse rate" !//cr//parspace//&
        th_param_desc( 3,idesc) = "effective N for summation" !//cr//parspace//&
        th_param_desc( 4,idesc) = "end-to-end radius of Gaussian coil" !//cr//parspace//&
        th_param_desc( 5,idesc) = "temperature if not in parameters" !//cr//parspace//&
        th_param_desc( 6,idesc) = "center of mass diffusion " !//cr//parspace//&
        th_param_desc( 7,idesc) = "sub diffusion strech exponent" !//cr//parspace//&
        th_param_desc( 8,idesc) = "sub diffusion to com r**2" !//cr//parspace//&
        th_param_desc( 9,idesc) = "transition sharpness exponent" !//cr//parspace//&
        th_param_desc(10,idesc) = "tan(mode amplitude) " !//cr//parspace//&
        th_param_desc(11,idesc) = "tan(mode amplitude) " !//cr//parspace//&
        th_param_desc(12,idesc) = "tan(mode amplitude) " !//cr//parspace//&
        th_param_desc(13,idesc) = "tan(mode amplitude) " !//cr//parspace//&
        th_param_desc(14,idesc) = "tan(mode amplitude) " !//cr//parspace//&
        th_param_desc(15,idesc) = "tan(mode amplitude) " !//cr//parspace//&
        th_param_desc(16,idesc) = "tan(mode amplitude) " !//cr//parspace//&
        th_param_desc(17,idesc) = "tan(mode amplitude) " !//cr//parspace//&
        th_param_desc(18,idesc) = "tan(mode amplitude) " !//cr//parspace//&
        th_param_desc(19,idesc) = "tan(mode amplitude) " !//cr//parspace//&
        th_param_desc(20,idesc) = "tan(mode amplitude) " !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > momentum transfer"
        th_file_param(  2,idesc) = "temp     > temperature"
        th_file_param(  2,idesc) = "fixend   > 1=fixed end chain, 0=normal rouse"
        th_file_param(  3,idesc) = "label1_s > relative (0..1) position of label l1=1..N*l1s, l2=..N*l2s, ..."
        th_file_param(  4,idesc) = "label2_s > relative (0..1) position of label l1=1..N*l1s, l2=..N*l2s, ..."
        th_file_param(  5,idesc) = "label3_s > relative (0..1) position of label l1=1..N*l1s, l2=..N*l2s, ..."
        th_file_param(  6,idesc) = "label1_a > contrast of first label section I (1..N+ls1)"
        th_file_param(  7,idesc) = "label2_a > contrast of first label section II "
        th_file_param(  8,idesc) = "label3_a > contrast of first label section III"
        th_file_param(  9,idesc) = "label4_a > contrast of first label section IV"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "diff     > diffusion in cm**2/s"
        th_out_param(  2,idesc) = "l        > effective segment length"
        th_out_param(  3,idesc) = "w        > rouse rate: W"
        th_out_param(  4,idesc) = "wl4      > inferred value of Wl4"
        th_out_param(  4,idesc) = "pfix     > used pfix value -0.5=fixed end, 0=normal"
! 
        th_nrouspm3 = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      amplitu  =      pa( 1)
      wl4      =  abs(pa( 2))
      n_segmen = nint(pa( 3))
      re       =  abs(pa( 4))
      temp     =      pa( 5)
      com_diff =  abs(pa( 6))
      nu_subdiff= abs(pa( 7))
      r02      =  abs(pa( 8))
      a_cross  =  abs(pa( 9))
      modeamp        =  1d0            !! default: ALL = 1 
      modeamp(1:11)  =  atan(abs(pa(10:20)))*2/Pi !! possibly modify the first 11

! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: momentum transfer
      xh =     0.1d0 
      call parget('q       ',xh,iadda,ier)
      q        = xh
! >>> extract: temperature
      xh =    300d0 
      call parget('temp    ',xh,iadda,ier)
      if( xh > 0d0) temp     = xh
!! >>> extract: labelzones fractional limits
      xh =    1d0 
      call parget('label1_s ',xh,iadda,ier)
      label1_s     = xh
      xh =    1d0 
      call parget('label2_s ',xh,iadda,ier)
      label2_s     = xh
      xh =    1d0 
      call parget('label3_s ',xh,iadda,ier)
      label3_s     = xh
!! >> label contrasts
      xh =    1d0 
      call parget('label1_a ',xh,iadda,ier)
      label1_a     = xh
      xh =    1d0 
      call parget('label2_a ',xh,iadda,ier)
      label2_a     = xh
      xh =    1d0 
      call parget('label3_a ',xh,iadda,ier)
      label3_a     = xh
      xh =    1d0 
      call parget('label4_a ',xh,iadda,ier)
      label4_a     = xh

      xh =    0d0
      call parget('fixend   ',xh,iadda,ier)
      fixend    = nint(xh)
      if(fixend == 1) then 
         pfix = -0.5d0
      else
         pfix = 0d0
      endif

! distribute contrasts to label-zone array
      if(n_segmen > size(labelarray)) n_segmen = size(labelarray)
      labelarray = 1d0
      nl1 = max(nint(n_segmen * label1_s),1)
      nl2 =     nint(n_segmen * label2_s)
      nl3 =     nint(n_segmen * label3_s)
      labelarray(1:nl1-1)        = label1_a 
      labelarray(nl1:nl2-1)      = label2_a 
      labelarray(nl2:nl3-1)      = label3_a 
      labelarray(nl3:n_segmen)   = label4_a 
      
     
!! TEST
!! write(*,*)"Tls:", n_segmen, label1_s, nl1,  label2_s, nl2,  label3_s, nl3 
!! write(*,*)"Tla:",  label1_a, label2_a,   label3_a, label4_a 
      if(maxval(abs(labelstat-[dble(nl1),dble(nl2),dble(nl3),label1_a,label2_a,label3_a,label4_a])) > 0.001d0) then
        labelstat =           [dble(nl1),dble(nl2),dble(nl3),label1_a,label2_a,label3_a,label4_a]
        write(*,'(a)')"Label status:"
        write(*,'(200i1)') nint(labelarray(1:n_segmen))
      endif
   
!! write(*,'(a,12f8.3,a)')"Tmod:",modeamp(1:12)," ...."

 
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     t   = x

!     Dr  = com_diff * 1d-9 / 1d-16  ! in A**2/ns
     Dr  = 1d-16  ! eff=0 so nrouse is without diff, which is added via rr

     call nrouspm3(q,t,temp,Dr,wl4,n_segmen,Re, W, l,labelarray, modeamp, Sq,Sqt)


     rr = ((exp(-log(r02/com_diff/6)*nu_subdiff)*r02*t** nu_subdiff)** a_cross +&
                (6*com_diff*t)**a_cross)**(1d0/a_cross)

     if(fixend == 1) rr = 0d0

     th_nrouspm3 = amplitu * sqt/sq  * exp( -(q*q*rr/6d0) )


 
! ---- writing computed parameters to the record >>>
      call        parset('l       ',sngl(l),iadda,ier)      ! in ns A units
      call        parset('w       ',sngl(W),iadda,ier)      !     "
      call        parset('wl4     ',sngl(W*l**4),iadda,ier) !     "
      call        parset('pfix    ',sngl(pfix),iadda,ier)   !     "
      do i=1,12
         write(buf,'("amod",i0)') i
         call        parset(buf,sngl(modeamp(i)),iadda,ier) 
      enddo
  
!      Dr        = Dr /( 1d-9 / 1d-16 ) ! in cm**2/s
!      call parset('diff    ',sngl(Dr),iadda,ier)
 
 CONTAINS 
 

       subroutine nrouspm3(q,t,temp,Dr,wl4,N,R, W, l,labelarray, modeamplitudes, Sq,Sqt)
!      ========================================================
!
! Rouse expression for a chain of finite length:
! Input parameters:
!    q     ----> momentum transfer in A**-1
!    t     ----> time in nano-sec
!    temp  ----> temperature in K
!    Dr    ----> center of mass diffusion constant in A**2/ns, if 0 <-- Rouse-expectation
!    wl4   ----> friction coefficient in A**4/ns
!    N     ----> number of chain segments
!    R     ----> end-to-end distance of the polymer molecule
!    labelarray      ----> labels along chain = segemntwise contrast array
!    modeamplitudes  ----> mode amplitude modification factors  
! Output parameters:
!    W     <--- "Rouse factor" 3kT/(xi*l**2); R**2=N*l**2
!    l     <--- "Segment length l"
!    Sq    <--- S(Q)
!    Sqt   <--- S(Q,t)
! ------------------------------------------------------------
!
       implicit none

       double precision kb, pi
       parameter(kb=1.380662d-23)
       parameter(pi=3.141592654d0)

       double precision :: q,t,temp,Dr,xi,R, W,Sq,Sqt, wl4
       integer          :: N, nn,mm,ifix,ip
       double precision :: labelarray(N), modeamplitudes(N)

       double precision l, tau_p, kbt, Sq0, arg1, arg2
       double precision a0,e0, ff2, ffc,    arg10,arg20
       double precision aa1 , aa2
       double precision p, p0fix, pfac

       double precision :: cosarray(N,N), ewfac(N)

       integer :: ipmin, ipmax, i

!       integer iout
       
       if(N.le.0) then
         W  = 999
         Sq = 999
         Sqt= 999
         write(6,*)'Error Number of chain segments is <= 0!',N
         return
       endif

! ---- determine the segment length l ----
       l = sqrt(R**2/N)       
       
! ---- and the Rousefactor ----
       kbt = temp*kb            ! in Joule = kg*m**2/s**2
       kbt = kbt * 100          ! in         kg*A**2/ns**2
       xi  = 3*kbt*l**2 / wl4
       W   = 3*kbt/(xi*(l**2))  ! in 1/ns


! ---- set the diffusion constant if input is zero --- !
       if(Dr.eq.0.0d0) then
         Dr = kbt/(N*xi)
       endif

!$OMP PARALLEL DO     
       do nn=1,N
        do ip=1,N
         cosarray(nn,ip) = cos((pi*(ip+pfix)*nn)/dfloat(N)) / (ip+pfix)  
        enddo
       enddo
!$OMP END PARALLEL DO   

!$OMP PARALLEL DO    
       do i=1,N
         ewfac(i) = (1d0-exp(-2*W*(1-cos((pi*(i+pfix))/dfloat(N)))*t)) * modeamplitudes(i)
       enddo
!$OMP END PARALLEL DO    

       ipmin = 1
       ipmax = N

! ---- init sums ----
       Sq0 = 0
       Sq  = 0
       Sqt = 0
       ff2  = -2*N*(l*q)**2/(3*pi**2)

! ---- Do the sums -----

!$OMP PARALLEL DO REDUCTION(+:Sq,Sqt)
       do nn = 1,N
        do mm = 1,N

          Sq  = Sq  + (exp(-(q**2)*(       abs(nn-mm)*(l**2)/6.0d0)) &
                      ) * labelarray(nn)*labelarray(mm)
          Sqt = Sqt + (exp(-(q**2)*(Dr*t + abs(nn-mm)*(l**2)/6.0d0) + &
                ff2 * sum(cosarray(nn,ipmin:ipmax) * cosarray(mm,ipmin:ipmax) *  ewfac(ipmin:ipmax) )) &
                      ) * labelarray(nn)*labelarray(mm) 

        enddo
       enddo
!$OMP END PARALLEL DO

       Sq  = Sq /N
       Sqt = Sqt/N

!!       write(6,'(1x,a,6E14.6)')'q,t,Sq,Sqt, Sqt/Sq, w=', q,t,Sq,Sqt, Sqt/Sq, w 

       return
       end
 

 end function th_nrouspm3
