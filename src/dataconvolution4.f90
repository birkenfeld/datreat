      SUBROUTINE datconv (idat, xio, yio, nin, xmaster, nout) 
!      -------------------------------------------------                
!                                                                       
! Variante, die mit den Scheiderschen Parameteren arbeitet              
! und den vollen Radialmittelungskertn enthaelt                         
!      -------------------------------------------------                
!      | convolution routine to be used with datreatz  |                
!      | userdefined convolution procedure(s)          |                
!      | to be supplied by the user on his own disk    |                
!      |                                               |                
!      | idat -- datareference no. of datreat item     |                
!      |         to be used as adress for parget       |                
!      |         to extract any resolution parameters  |                
!      |         needed by this subroutine.            |                
!      | xio  -- input x-vector  ---> result on return |                
!      | yio  -- input y-vector  ---> result on return |                
!      |         these vectors contain data to be conv.|                
!      | nin  -- no. of elements in xio and yio        |                
!      |         will be set to nout on output         |                
!      | xmaster x-vector containing the x-values to be|                
!      |         supplied on output                    |                
!      | nout -- no. of xmaster elememnts              |                
!      |-----------------------------------------------|                
!                                                                       
      IMPLICIT none 
      INTEGER maxbuf 
      PARAMETER (maxbuf = 4000) 
                                                                        
      INTEGER idat, nin, nout 
      REAL xio, yio, xmaster 
      DIMENSION xio ( * ), yio ( * ), xmaster ( * ) 
                                                                        
      REAL ybuff 
      DIMENSION ybuff (maxbuf) 
                                                                        
      REAL sa, sbq 
      REAL resker 
      REAL q, q0, sigmasqr 
      REAL xn, sum, rk 
      INTEGER iout, iot, nbuf, ier 
      INTEGER i, j, j1, janf 
                                                                        
      REAL fwhm, r_c, r_s, l_c, l_d, d_qav, d_qd, xh 
      REAL lambda, k0, dtheta 
                                                                        
                                                                        
      LOGICAL mark 
      SAVE sa, sbq 
                                                                        
      REAL pi, sqpi, ln2 
      PARAMETER (pi = 3.141592654, sqpi = 1.772453851, ln2 = 0.69314718) 
                                                                        
                                                                        
      DATA sa / 0.15 /, sbq / 0.1e-5 / 
!                                                                       
!                                                                       
!      ----- standard-preliminaries: do not change ! -------            
!            ---------------------------------------                    
                                                                        
      nbuf = nout 
                                                                        
      IF (nout.gt.maxbuf) then 
         CALL errsig (99, 'datconv: dimension nout exceeds buffer') 
         nbuf = maxbuf 
      ENDIF 
      iot = iout () 
                                                                        
!     ------- body 1: extract the parameters ------------------------   
!             include here parget calls for all resolution parameters   
!             that are needed and expected                              
!     ---------------------------------------------------------------   
!                                                                       
! here we try to compute sa and sbq from other information              
! Vgl. Doktorarbeit D. Schneiders, p.44                                 
                                                                        
      ier = - 1 
      lambda = 7.0d0 
      CALL parget ('lambda  ', lambda, idat, ier) 
      ier = - 1 
      fwhm = 0.2d0 
      CALL parget ('fwhm_sel', fwhm, idat, ier) 
      ier = - 1 
                        ! Radius collimation hole                       
      r_c = 0.015 
      CALL parget ('diam_c  ', xh, idat, ier) 
      IF (ier.eq.0) r_c = xh / 2 
      ier = - 1 
                        ! Radius sample area                            
      r_s = 0.005 
      CALL parget ('diam_s  ', xh, idat, ier) 
      IF (ier.eq.0) r_s = xh / 2 
      ier = - 1 
                        ! collimation length                            
      l_c = 8.0 
      CALL parget ('l_c     ', l_c, idat, ier) 
      ier = - 1 
                        ! detektor distance                             
      l_d = 8.0 
      CALL parget ('l_d     ', l_d, idat, ier) 
      ier = - 1 
                        ! detektor element dq(fwhm)                     
      d_qd = 0.0 
      CALL parget ('d_qd    ', d_qd, idat, ier) 
      ier = - 1 
                        ! detektor averaging ring width dq(fwhm)        
      d_qav = 0.0 
      CALL parget ('d_qav   ', d_qav, idat, ier) 
                                                                        
!! > Bug fix 22.5.02 !! l_c <-> l_d !!                                  
      k0 = 2 * Pi / lambda 
      IF (r_c / (l_c + l_d) .gt.r_s / l_d) then 
         dtheta = 2 * r_c / l_c 
      ELSE 
         dtheta = 2 * r_s * (1 / l_c + 1 / l_d) 
      ENDIF 
                                                                        
      sbq = ( (k0 * dtheta) **2 + d_qd**2 + d_qav**2) / (8 * ln2) 
      sa = fwhm / sqrt (8 * ln2) 
                                                                        
                                                                        
      CALL parget ('sa      ', sa, idat, ier) 
      IF (ier.ne.0) then 
!             call errsig(98,'datconv replace sa     by default')       
      CALL parset ('sa      ', sa, idat) 
      ENDIF 
                                                                        
      CALL parget ('sbq     ', sbq, idat, ier) 
      IF (ier.ne.0) then 
!             call errsig(98,'datconv replace sbq    by default')       
      CALL parset ('sbq     ', sbq, idat) 
      ENDIF 
                                                                        
                                                                        
! ------- now do the convolution ! --------------------------------     
!         store result on ybuff  !                                      
!         compute for all x-values from xmaster !                       
!         ---------------------------------------                       
                                                                        
      IF (iot.gt.0) then 
         WRITE (6, * ) 'convoluting:' 
         WRITE (6, * ) 'theory-field: x1=', xio (1) , ' x2=', xio (nin) &
         , ' nin =', nin                                                
         WRITE (6, * ) 'result-field: x1=', xmaster (1) , ' x2=',       &
         xmaster (nout) , ' nout=', nbuf                                
      WRITE (6,  * ) 'sa       = ', sa 
      WRITE (6,  * ) 'sbq      = ', sbq 
      ENDIF 
!                                                                       
      janf = 1 
      j1 = 1 
      DO i = 1, nbuf 
      Q0 = xmaster (i) 
      sigmasqr = (Q0 * sa) **2 + abs (sbq) 
                                                                        
                                                                        
      sum = 0.0 
      xn = 0.0 
      mark = .false. 
      DO j = janf, nin 
      q = xio (j) 
      rk = resker (q, q0, sigmasqr) 
!           ---- mark if nonzero part of resolution is there --         
      IF (.not.mark) then 
         IF (rk.ne.0.0) then 
            mark = .true. 
            j1 = j 
         ENDIF 
      ENDIF 
!           ---- finish if resolution is done ---                       
      IF (mark) then 
         IF (rk.eq.0.0) goto 200 
      ENDIF 
!           ----                                                        
      sum = sum + rk * yio (j) 
      xn = xn + rk 
      IF (iot.gt.3) write (6, 111) i, j, q0, q, rk, yio (j) 
  111 FORMAT        (2i4,' q0=',f12.6,' q=',f12.6,' rk=',e13.6,         &
     &                   ' y=',e13.6)                                   
      enddo 
  200 CONTINUE 
      janf = j1 
      IF (xn.gt.0.0) then 
         ybuff (i) = sum / xn 
      ELSE 
         ybuff (i) = 0.0 
      ENDIF 
      IF (iot.gt.2) write (6, * ) 'y(', i, ')=', ybuff (i) 
      enddo 
                                                                        
! ------- put the buffer and xmaster onto yio and xio ! -----------     
!         set nin to nbuf !                                             
!         and return !                                                  
!         ----------------------------------------------                
      nin = nbuf 
      DO i = 1, nbuf 
      xio (i) = xmaster (i) 
      yio (i) = ybuff (i) 
      enddo 
                                                                        
      RETURN 
      END SUBROUTINE datconv                        
                                                                        
      FUNCTION resker (q, q0, sigmasqr) 
!      ------------------------------------------------                 
!      | customized resolution kernel for small-angle |                 
!      | scattering in gaussian approximation         |                 
!      | q0      -- momentum transfer(data)           |                 
!      | q       -- scattering angle / 2(theo)        |                 
!      | sigmasqr-- sigma**2                          |                 
!      ------------------------------------------------                 
      IMPLICIT none 
                                                                        
      REAL resker, q, q0, sigmasqr 
      REAL qs, arg_exp, arg2 
      REAL arg_I0 
                                                    !<- IMSL            
      REAL BSI0, BSI0E 
      REAL(8) DBSI0, DBSI0E 
                                                                        
      qs = q / sigmasqr 
      arg_exp = - 0.5 * (q**2 + q0**2) / sigmasqr 
      arg_I0 = abs (q * q0 / sigmasqr) 
      arg2 = - ( (q - q0) **2 / (2 * sigmasqr) ) 
                                                                        
      IF (arg2.lt. - 15.0) then 
         resker = 0.0 
      ELSE 
         resker = qs * exp (arg_exp + arg_I0) * BSI0E (arg_I0) 
!         resker = qs*   exp(arg2)                                      
      ENDIF 
!                                                                       
      RETURN 
      END FUNCTION resker                           
