      FUNCTION th_cloizeaux (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf) 
!     ================================================================================               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Formel (23) aus Cloizeaux, J. de Physique I, 3 (1993) p.1533 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none
      
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20)
		integer     , intent(inout) :: mbuf
		integer     , intent(inout) :: nopar                    ! Anzahl der Parameter data
		character*80, intent(inout) :: napar(mbuf)              ! name des parameters n
		real        , intent(inout) :: params(mbuf)             ! value des parameters n
		DATA zpi / 6.283185 / 


         DOUBLE PRECISION Pi
         Parameter       (Pi=3.141592654d0)


         REAL            th_cloizeaux, x, pa, qq, zpi, xh, vol_frac
         INTEGER         ini, npar, nparx
         DIMENSION       pa(20),qq(3)
         DATA            zpi/6.283185/

         double precision w,d,q,n,ne,cloizeaux,t,a,wl4


 
         INTEGER iadda, ier, iout, iot                
         COMMON/thiadd/iadda         !! transfer of address for parget

!
! ----- initialisation -----
       if(ini.eq.0) then
         thnam = 'cloizeaux'
         nparx =  3
         if(npar.lt.nparx) then
           write(6,1)thnam,nparx,npar
1          format(' theory: ',a8,' no of parametrs=',i8,' exceeds current max. = ',i8)
           th_cloizeaux   = 0
           return
         endif
         npar = nparx
!        --------------> set the number of parameters
         parnam(1) = 'amplitud'          
         parnam(2) = 'wl4     '
         parnam(3) = 'dtube   '
!
         th_cloizeaux  = 0
         return
       endif
!
! ---- calculate theory here -----
        a       = pa(1)   ! amplitude (should be 1)
        wl4     = pa(2)   ! Rouse rate W*l**4
        d       = pa(3)   ! Tube diameter

        call parget('q       ',xh ,iadda,ier)
        if(ier.eq.0) then
          q       = xh
        else
          q       = 1.0  
        endif
 
        t = x 
 
        th_cloizeaux = a * cloizeaux(t,q,d,Wl4) 

         return
        end 

