                                                                 
      FUNCTION th_moderators(x, pa, thnam, parnam, npar,  ini, nopar ,params,napar,mbuf)
!     --------------------------------------------------------------------------------
! Moderator Brilliance models with max three maxwellians + epithermal
! assumed units: Brilliance in  n/cm**2/s/sterad/A
!                       
      implicit none
                                           
      CHARACTER(8) thnam, parnam (20) 
      real*4            :: pa (20)
			integer :: mbuf
			integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
			real, intent(inout) :: params(mbuf)             ! value des parameters n

      integer           :: ini, npar, nparx

      real*4            :: th_moderators, x

      double precision  :: Temp1, Temp2, Temp3
      double precision  :: Flux1, Flux2, Flux3, Flux_epithermal
      double precision  :: a1_epithermal, a2_epithermal
      double precision  :: lambda                        ! in A
  
      double precision  :: Maxwell

 

!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'moderato' 
         nparx = 9
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th_moderators = 0
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'flux1' 
         parnam (2) = 'temp1' 
         parnam (3) = 'flux2' 
         parnam (4) = 'temp2' 
         parnam (5) = 'flux3' 
         parnam (6) = 'temp3' 
         parnam (7) = 'flux_epi' 
         parnam (8) = 'a1_epi' 
         parnam (9) = 'a2_epi' 
       
!                                                                       
         th_moderators = 0
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----    
      flux1      = abs(pa(1))
      temp1      = abs(pa(2))
      flux2      = abs(pa(3))
      temp2      = abs(pa(4))
      flux3      = abs(pa(5))
      temp3      = abs(pa(6))
      Flux_epithermal = abs(pa(7))
      a1_epithermal   = pa(8)
      a2_epithermal   = pa(9)

                                  
      lambda = x 
      th_moderators = flux1 * Maxwell(lambda, temp1) &
                 + flux2 * Maxwell(lambda, temp2) &        
                 + flux3 * Maxwell(lambda, temp3) &     
          +  Flux_epithermal/(1d0+exp(a1_epithermal*lambda-a2_epithermal))/lambda
  

      ! write(6,*)'flux1           ',flux1          
      ! write(6,*)'temp1           ',temp1          
      ! write(6,*)'flux2           ',flux2          
      ! write(6,*)'temp2           ',temp2          
      ! write(6,*)'flux3           ',flux3          
      ! write(6,*)'temp3           ',temp3          
      ! write(6,*)'Flux_epithermal ',Flux_epithermal
      ! write(6,*)'a1_epithermal   ',a1_epithermal  
      ! write(6,*)'a2_epithermal   ',a2_epithermal  
      ! write(6,*)'max1            ', Maxwell(lambda, temp1) 
      ! write(6,*)'max2            ', Maxwell(lambda, temp2) 
      ! write(6,*)'max2            ', Maxwell(lambda, temp3)
      ! write(6,*)'epi             ', Flux_epithermal/(1d0+exp(a1_epithermal*lambda-a2_epithermal))/lambda
      ! write(6,*)'moderators      ',th_moderators
     


!                                                                       
      RETURN 
      END FUNCTION th_moderators




  double precision function Maxwell(lambda, Temp)
! -----------------------------------------------
!
! normalized Maxwell Spectrum 
!
  implicit none
  double precision, intent(in)    ::  lambda      ! wavelength in m
  double precision, intent(in)    ::  Temp        ! Moderator Temp in K
 
  double precision                :: lambda_A
  double precision                :: A
  double precision, parameter     :: A0=949.0d0
  double precision                :: x, y

  lambda_A = lambda   !/ unit_Angstroem

  if(lambda_A.lt.1d-6) lambda_A = 1d-6
  
  if( Temp > 0 ) then 
    A = A0/Temp
  else
    Maxwell = 0
    return
  endif

  y = -A/lambda_A**2

  if( y.lt.-100.0d0) then
    y = -100.0d0
  endif

  Maxwell = 2 * A**2/(lambda_A**5) * exp(y)
 
  return
  end function Maxwell
