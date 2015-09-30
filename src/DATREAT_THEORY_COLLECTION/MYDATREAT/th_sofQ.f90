      FUNCTION th20 (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
! -------> sofq <--------                                               
!                                                                       
!                                                                       
      PARAMETER (pi = 3.141592654) 
      CHARACTER(8) thnam, parnam (20) 
      DIMENSION pa (20), qq (3) 
      integer :: mbuf
      integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
      real, intent(inout) :: params(mbuf)             ! value des parameters n
      REAL(8) qqr (1), sq, eta, scl, gamma, r 
      REAL(8) d, den 
                                                        !! aix.sp extchk
      REAL(8) pschulz, pschj1, betaj1, adapint, peryev 
                                                        !! aix          
      REAL(8) q 
       
      double precision sofq_table

      integer new
      common/stabnew/new 

      DATA new /0/

      DATA zpi / 6.283185 / 
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'sofq' 
         nparx = 5 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th20 = 0 
            RETURN 
         ENDIF 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'volfrac' 
         parnam (2) = 'scalelen' 
         parnam (3) = 'gamma' 
         parnam (4) = 'r' 
         parnam (5) = 'mode' 
!                                                                       
         th20 = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
      eta = pa (1) 
      scl = pa (2) 
      gamma = pa (3) 
      r = pa (4) 
      mode = pa (5) + 0.1 
                                                                        
      IF (mode.ge.1.and.mode.le.5) then 
         qqr (1) = x * r * 2.0 
         CALL sofq (qqr, sq, 1, eta, scl, gamma, r, mode, ierr) 
         IF (ierr.lt.0) write (6, * ) 'sofq: ierr=', ierr 
      ELSE 
         IF (mode.eq.6) then 
            q = x 
            d = 2 * r 
            den = 6.0d0 * eta / (pi * d * d * d) 
            sq = peryev (q, r, den, - 1d-7) 
         ELSE 
            WRITE (6, * ) 'sofq: mode=', mode, ' is out of range' 
            sq = 1.0 
         ENDIF 
      ENDIF 
  
      if(mode.le.0) then
        if(mode.lt.0) new=0
        qqr (1) = x * r * 2
        th20 = sofq_table(qqr(1))
      endif 

      th20 = sq 
!                                                                       
      RETURN 
      END FUNCTION th20          







       real*8 function sofq_table(y)
!!     ==============================
!!                                y = Q*R
        implicit none

        integer ntab, i
        integer new
        integer mtab
        Parameter(mtab=10000) 
        real*8 dgy, qr(mtab), sofq_tab(mtab)
        common/csofqtab/dgy, sofq_tab
        common/stabnew/new 


        real*8 ya,y,h, yout,arg
        integer i1,i2,i3


        if(new.eq.0) then
          write(6,*)'Reading sofq.dat .....'
          open(19,file='sofq.dat',status='old')
          do i=1,ntab
            read(19,*,end=999)qr(i),sofq_tab(i) 
          enddo
999       continue
          close(19)
          ntab = i-1
          new  = 1
          write(6,*)'..done no pts= ', ntab
        endif

        dgy = qr(2)-qr(1)
     
        ya = dabs(y)
        i1 = ya/dgy + 1
        i2 = i1+1
        i3 = i1+2

        if(i2.ge.ntab) then
          sofq_table = 1.0d0
          return
        endif

        
     
        yout =  sofq_tab(i1)*(ya-i3*dgy)*(ya-i2*dgy)/((i1-i3)*(i1-i2)) &
               +sofq_tab(i2)*(ya-i1*dgy)*(ya-i3*dgy)/((i2-i1)*(i2-i3)) &
               +sofq_tab(i3)*(ya-i1*dgy)*(ya-i2*dgy)/((i3-i1)*(i3-i2))
        yout = yout/(dgy*dgy)

        sofq_table = yout

       return
       end

                
