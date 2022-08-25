 FUNCTION th_nvilgis(x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!================================================================================
!  Scattering factor of a Rouse  chain in a harmonic potential. Using direct summation over an effective number of beads, however, using the infinite chain rms distance form Eq(24) of the publication.
!  T.A. Vilgis and F. Boue, Journal of Polymer Science Part B Polymer Physics 26, 2291-2301 (1988)
      use theory_description 
      implicit none 
      real    :: th_nvilgis
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
     double precision :: re         ! effective segment length via Gaussian end-to-end re
     double precision :: wl4        ! Rouse rate                                                                      
     double precision :: rmesh      ! effective potential parameter in terms mesh size, see paper                     
     double precision :: tmesh      ! lifetime of the rmesh constraint (experimental)                                 
! the recin parameter representation 
     double precision :: q          ! q-value    default value                                                        
! the reout parameter representation 
     double precision :: rg         ! radius of gyration                                                              
 
     double precision :: th
 
     double precision :: q0, l, Sq, Sqt, t
     integer          :: ni
!
! ----- initialisation ----- 
    IF (ini.eq.0) then     
       thnam = 'nvilgis'
       nparx =        6
       IF (npar.lt.nparx) then
           WRITE (6,*)' theory: ',thnam,' no of parametrs=',nparx,' exceeds current max. = ',npar
          th_nvilgis = 0
          RETURN
       ENDIF
       npar = nparx
! >>>>> describe theory with >>>>>>> 
       idesc = next_th_desc()
       th_identifier(idesc)   = thnam
       th_explanation(idesc)  = " Scattering factor of a Rouse  chain in a harmonic potential. Using direct summation over an effective number of beads, however, using the infinite chain rms distance form Eq(24) of the publication."
       th_citation(idesc)     = " T.A. Vilgis and F. Boue, Journal of Polymer Science Part B Polymer Physics 26, 2291-2301 (1988)"
!       --------------> set the parameter names --->
        parnam ( 1) = 'ampli   '  ! prefactor                                                                       
        parnam ( 2) = 'n       '  ! number of segments (do not fit, its integer)                                    
        parnam ( 3) = 're      '  ! effective segment length Re = l * N**nu /sqrt(6)                                
        parnam ( 4) = 'wl4     '  ! Rouse rate                                                                      
        parnam ( 5) = 'rmesh   '  ! effective potential parameter in terms mesh size, see paper                     
        parnam ( 6) = 'tmesh   '  ! lifetime of the rmesh constraint (experimental)                                 
! >>>>> describe parameters >>>>>>> 
        th_param_desc( 1,idesc) = "prefactor" !//cr//parspace//&
        th_param_desc( 2,idesc) = "number of segments (do not fit, its integer)" !//cr//parspace//&
        th_param_desc( 3,idesc) = "effective segment length Rg = l * N**nu /sqrt(6)" !//cr//parspace//&
        th_param_desc( 4,idesc) = "Rouse rate" !//cr//parspace//&
        th_param_desc( 5,idesc) = "effective potential parameter in terms mesh size, see paper" !//cr//parspace//&
        th_param_desc( 6,idesc) = "lifetime of the rmesh constraint (experimental)" !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
        th_file_param(:,idesc) = " " 
        th_file_param(  1,idesc) = "q        > q-value    default value"
! >>>>> describe record parameters creaqted by this theory >>>>>>> 
        th_out_param(:,idesc)  = " "
        th_out_param(  1,idesc) = "rg       > radius of gyration"
! 
        th_nvilgis = 0.0
 
        RETURN
     ENDIF
!
! ---- transfer parameters -----
      ampli    =      pa( 1)
      n        = ABS( pa( 2) )
      re       = ABS( pa( 3) )
      wl4      = ABS( pa( 4) )
      rmesh    = ABS( pa( 5) )
      tmesh    = ABS( pa( 6) )
! ---- extract parameters that are contained in the present record under consideration by fit or thc ---
      iadda = actual_record_address()
! >>> extract: q-value    default value
      xh =      0
      call parget('q       ',xh,iadda,ier)
      q        = xh
! 
! ------------------------------------------------------------------
! ----------------------- implementation ---------------------------
! ------------------------------------------------------------------
! 
     t   = x
     ni  = nint(n)
     l   = re/sqrt(n)
     rg  =  l * dble(ni) / sqrt(6d0)
     q0  = (1d0/rmesh**2) * (l**2) / 3d0
     call  NrouseV(q,t,wl4,l,q0,ni, Sq,Sqt)

     th  = ampli * Sqt/Sq


     th_nvilgis = th
 
! ---- writing computed parameters to the record >>>  
      call parset('rg      ',sngl(rg),iadda,ier)
 
 CONTAINS 
 
! subroutines and functions entered here are private to this theory and share its variables 
 


       subroutine NrouseV(q,t,wl4,l,q0,N, Sq,Sqt)
!      ==========================================
!
! Vilgis Boue: Rousechain in a potential
!
       implicit none

       double precision, intent(in)   :: q             ! momentum transfer
       double precision, intent(in)   :: t             ! time
       double precision, intent(in)   :: wl4           ! Rouse rate
       double precision, intent(in)   :: l             ! (effective) segment length
       double precision, intent(in)   :: q0            ! potential parameter as defined in Vilgis Boue
       integer         , intent(in)   :: N             ! number of effective segments

       double precision, intent(out)  :: Sq            ! Sqt=0
       double precision, intent(out)  :: Sqt           ! Sqt

!       double precision, parameter  :: kb=1.380662d-23
       double precision, parameter :: pi=3.141592654d0

       double precision :: wlt05, w
 
       double precision :: rrnnt(0:N-1)
       double precision :: rrnn0(0:N-1)

       integer :: nn, i

!       integer iout

       if(N.le.0) then
         Sq = 999
         Sqt= 999
         write(6,*)'Error Number of chain segments is <= 0!',N
         return
       endif

! ---- and the Rousefactor ----
!       xi  = 3*kbt*l**2 / wl4
!       W   = 3*kbt/(xi*(l**2))  ! in 1/ns

       w     = wl4 / l**4
       wlt05 = sqrt(w*t)

iq0:   if(q0 > 1d-2) then

!$OMP PARALLEL DO
          do nn=0,N-1
             rrnnt(nn) = l**2/q0 &
                        * (1d0- 0.5d0 &
                        *  ( 2d0 * cosh(q0*nn)  &
                        -    ( &
                                exp(-q0*nn) * erf(q0*wlt05 - nn/(2d0*wlt05)) &
                              + exp( q0*nn) * erf(q0*wlt05 + nn/(2d0*wlt05)) &
                             ) &
                           ) &
                          )
   ! limit for  t==> 0
             rrnn0(nn)= (l**2*exp(q0*nn)**2-0.2D1*l**2*cosh(q0*nn)*exp(q0*nn)+0.2D1*l**2*exp(q0*nn)-l**2)/exp(q0*nn)/q0/0.2D1
          enddo
!$OMP END PARALLEL DO
       else  ! series expansion with respect to q0 up to q0**3
!$OMP PARALLEL DO
          do nn=0,N-1
             rrnnt(nn) = &
                      l**2*(nn*erf(nn/wlt05/2d0)*sqrt(Pi)+2d0*exp(-nn**2/wlt05**2/4d0)*wlt05)*Pi**(-1d0/2d0)- &
                      l**2*nn**2*q0/2d0+l**2*(nn**3*erf(nn/wlt05/2d0)*sqrt(Pi)+0.2D1*exp(-nn**2/wlt05**2/4d0)*nn**2*wlt05- &
                      4d0*exp(-nn**2/wlt05**2/4d0)*wlt05**3)*Pi**(-1d0/2d0)*q0**2/6d0-l**2*nn**4*q0**3/24d0
  ! limit for  t==> 0
             rrnn0(nn) = l**2*nn-l**2*nn**2*q0/2d0+l**2*nn**3*q0**2/6d0-l**2*nn**4*q0**3/24d0
          enddo
!$OMP END PARALLEL DO
       endif iq0

! ---- init sums ----
       Sq  = 0
       Sqt = 0

! ---- Do the sums -----

! !$OMP PARALLEL DO REDUCTION(+:Sq,Sqt)
!        do nn = 1,N
!         do mm = 1,N
! 
!           Sq  = Sq  + exp(-(q**2)*rrnn0(abs(nn-mm))/6.0d0)
!           Sqt = Sqt + exp(-(q**2)*rrnn1(abs(nn-mm))/6.0d0)
! 
!         enddo
!        enddo
! !$OMP END PARALLEL DO


!$OMP PARALLEL DO REDUCTION(+:Sq,Sqt)
       do nn = 1,N-1
          Sq  = Sq  + (N-nn)*exp(-(q**2)*rrnn0(nn)/6.0d0)
          Sqt = Sqt + (N-nn)*exp(-(q**2)*rrnnt(nn)/6.0d0)
       enddo
!$OMP END PARALLEL DO
       Sq  = 2*Sq  + N*exp(-(q**2)*rrnn0(0)/6.0d0)
       Sqt = 2*Sqt + N*exp(-(q**2)*rrnnt(0)/6.0d0)


       Sq  = Sq /N
       Sqt = Sqt/N

!       write(6,'(1x,a,6E14.6)')'q,t,Sq,Sqt, Sqt/Sq ', q,t,Sq,Sqt, Sqt/Sq

       return
       end





 end function th_nvilgis
