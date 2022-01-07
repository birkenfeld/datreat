                                                       
                                                    
                                                                    
      FUNCTION echotri3 (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
!    nse echotri3 function 
!                                                                       
!                                                                       
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
            integer :: mbuf
            integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
            character(len=80), intent(inout) :: napar(mbuf)      ! name des parameters n
            real, intent(inout) :: params(mbuf)             ! value des parameters n
                  DATA zpi / 6.283185 / 
      DATA gam_larmor / 183.0333e6 / 
      DATA amass_n / 0.16749286e-26 / 
      DATA hplanck / 0.66260755e-33 / 

      double precision :: fwhm, z, wx, dlam2, beta, p
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'echotri3 ' 
         nparx = 7 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            echotri3 = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'amplitu' 
         parnam (2) = 'average' 
         parnam (3) = 'lslope'    ! spectral slope 
         parnam (4) = 'lambda0' 
         parnam (5) = 'lamfwhm' 
         parnam (6) = 'lamwid2' 
         parnam (7) = 'jshift ' 
                                                                        
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
         echotri3 = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      xs = x - pa (7)   !! delta j in Tm
      beta = pa(3)  * 1e10
      alam = pa (4) * 1e-10 
      fwhm = pa (5) * 1e-10
      dlam2= pa (6) * 1d-10
      omeg = gam_larmor * amass_n / (hplanck) 
!      omeg = gam_larmor * amass_n / hplanck * (2 * zpi * 1e-7) 
!     omeg = omeg * xtur * (blinteg (acoil, bcoil, coill, rcoil)        &
!      - blintegT (acoil2, acoil, coill, rcoil) ) / 2                     
                                                                        
                                                                        
!      write(6,*)'omeg: ',omeg,omeg*alam,omeg*alam*57.29578             
      CALL setpar ('pha_sens ', sngl (omeg * alam * 57.29578d0) ,nopar ,params,napar,mbuf, ier)  
                                                                        
!      echotri3 = pa (2) + pa (1) * exp ( - (xs * omeg * dlam / 2) **2)       &
!      * cos (omeg * alam * (xs + pa (3) ) )                             
!                     
       z    = xs * omeg * fwhm
       p    = xs * omeg

! write(*,*)"z=",z,"  p=",p," omeg=",omeg

       if( abs(z) > 1d-5 ) then 
         wx   = cos(p*alam) * 2d0*(1d0-cos(z))/z**2  &
              + beta * fwhm* ((2d0*(z*sin(z)+2d0*cos(z)-2d0))*sin(alam*p)/(z**3))
       else
         wx = 1d0
       endif

       wx   = wx *  exp ( - (xs * omeg * dlam2 / 2) **2)  

!       write(6,*) x, wx
!       wx   = 2d0*(1d0-cos(z))/z**2 *  exp ( - (xs * omeg * dlam2 / 2) **2)  
!       echotri3 = pa(2)+pa(1)*cos(omeg*alam*(xs))*wx 

       echotri3 = pa(2)+pa(1)*wx 

                                                  
      RETURN 
      END FUNCTION echotri3                              
       
