                                                       
                                                    
                                                                    
      FUNCTION echo (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
!    nse echo function 
!                                                                       
!                                                                       
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
			integer :: mbuf
			integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
			real, intent(inout) :: params(mbuf)             ! value des parameters n
			      DATA zpi / 6.283185 / 
      DATA gam_larmor / 183.0333e6 / 
      DATA amass_n / 0.16749286e-26 / 
      DATA hplanck / 0.66260755e-33 / 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'echo ' 
         nparx = 12 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            echo = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'amplitu' 
         parnam (2) = 'average' 
         parnam (3) = 'derotat' 
         parnam (4) = 'lambda0' 
         parnam (5) = 'lamfwhm' 
         parnam (6) = 'ishift ' 
         parnam (7) = 'nturns ' 
      parnam (8)  = 'rcoil  ' 
      parnam (9)  = 'lcoil  ' 
      parnam (10)  = 'acoil  ' 
      parnam (11)  = 'bcoil  ' 
      parnam (12)  = 'a2coil  ' 
                                                                        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Zugrunde liegende Geometrie                                          
!!                                                                      
!!                                                                      
!!  a2coil                 acoil                   bcoil                
!!                                 lcoil                                
!!  |----------------------|-------====------------|                    
!!                                 ^                                    
!!  pi/2                   pi      0-pkt           pi/2                 
!!                                                                      
!!                          <acoil><====bcoil=====>                     
!!  <=====acoil2==================>                                     
!!                                 <lcoil>                              
!!                                                                      
!!                                                                      
!!  rcoil: effektiver Spulenradius                                      
!!                                                                      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                        
!                                                                       
         echo = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      xs = x - pa (6) 
      alam = pa (4) * 1e-10 
      dlam = pa (5) * 1e-10 * 0.5 / sqrt (alog (2.0) ) 
      xtur = pa (7) 
      rcoil = pa (8) 
      coill = pa (9) 
      acoil = pa (10) 
      bcoil = pa (11) 
      acoil2 = pa (12) 
      omeg = gam_larmor * amass_n / hplanck * (2 * zpi * 1e-7) 
      omeg = omeg * xtur * (blinteg (acoil, bcoil, coill, rcoil)        &
      - blinteg (acoil2, acoil, coill, rcoil) ) / 2                     
                                                                        
                                                                        
!      write(6,*)'omeg: ',omeg,omeg*alam,omeg*alam*57.29578             
      CALL setpar ('pha_sens ', sngl (omeg * alam * 57.29578d0) ,nopar ,params,napar,mbuf, ier)  
                                                                        
      echo = pa (2) + pa (1) * exp ( - (xs * omeg * dlam / 2) **2)       &
      * cos (omeg * alam * (xs + pa (3) ) )                             
!                                                                       
      RETURN 
      END FUNCTION echo                              
                                                                        
                                                                        
      FUNCTION blinteg (a, b, xl, r) 
                                                                        
      blinteg = sqrt ( (r / xl) **2 + ( (xl - a) / xl) **2) - sqrt ( (r &
      / xl) **2 + ( (xl - b) / xl) **2) + sqrt ( (r / xl) **2 + (b / xl)&
      **2) - sqrt ( (r / xl) **2 + (a / xl) **2)                        
      RETURN 
      END FUNCTION blinteg                          
!                                                                       
!                                                                       
!*ds          
