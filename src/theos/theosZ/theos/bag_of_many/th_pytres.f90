      FUNCTION th_pythres (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
! -------> Resistance of Zones of Pythargoras coil with finite cuts and deltah <----                                               
!                                                                       
!                                                                       
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20)
		integer     , intent(inout) :: mbuf
		integer     , intent(inout) :: nopar                    ! Anzahl der Parameter data
		character*80, intent(inout) :: napar(mbuf)              ! name des parameters n
		real        , intent(inout) :: params(mbuf)             ! value des parameters n
		DATA zpi / 6.283185 / 
!                 

      double precision n, alpha, beta, rho, L, dh, delta, scale
                                                      
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'pythres' 
         nparx = 7
		if(npar.lt.nparx) then 
           write(6,"(' theory: ',a8,' no of parametrs=',i8,' exceeds current max. = ',i8)")thnam,nparx,npar
           th_pythres = 0
           return
        endif
         npar = nparx 
!        --------------> set the number of parameters                   

         parnam (1) = 'scale'    ! Skalierungsfaktor (Strom fuer Spannungsabfal o.ae)
         parnam (2) = 'rho'      ! Leitfaehigkeit
         parnam (3) = 'alpha'    ! Nominalbreite der 1. Zone
         parnam (4) = 'beta'     ! Steigung des Keils (8mm/50mm)
         parnam (5) = 'L'        ! Laenge der Messstrecke
         parnam (6) = 'dcuts'    ! Breite der Schnitte
         parnam (7) = 'hoffset'  ! Offsetdicke
         
!                                                                       
         th_pythres = 0
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      
      n = x
      scale = pa(1) 
      rho   = pa(2)
      alpha = pa(3) 
      beta  = pa(4)
      L     = pa(5) 
      delta = pa(6) 
      dh    = pa(7)

      

      th_pythres = pa(1)*L*rho/(alpha*(sqrt(n+1)-sqrt(n))-Delta)/  &
                 & ((alpha*sqrt(n)+alpha*(sqrt(n+1)-sqrt(n))/2)*beta+dh)
!                                                                       
      RETURN 
      END FUNCTION th_pythres
