      FUNCTION qspot (t, p, x1, x2) 
!      =========================                                        
!                                                                       
      DIMENSION p (10), xlim (2) 
! ** common-block to hold information for theories **                   
!                  ------> range of data to be used                     
!                                                                       
!                                                                       
      COMMON / gai / a (16), x (16), sxphi (320), cxphi (320), aphi (   &
      320), cxtht (320), sxtht (320), atht (320), ngs                   
! ** holds arguments and weightfactors and sin/cos of arguments         
! ** for multiple 32-pts gaussian integration **                        
      COMMON / qspt / pold (10), y (300), ntbl 
! ** holds last version of parameters,                                  
! ** table of function values from numerical integration                
! ** length of this table                                               
                                                                     !! 
!                                                                       
! ** preliminary test with 6-pts. gauss-integration**                   
!      data a/0.1713245,0.3607616,0.4679139,13*0./                      
!      data x/0.93247,0.6612094,0.23866192,13*0./                       
!      data ngs/3/                                                      
! ** contains data for standard  6-pts gaussian integration **          
! ** in b l o c k  d a t a  **                                          
!                                                                       
! ** meaning of parameters **                                           
!    p(1) = hight (multiplicative factor)                               
!    p(2) = 2theta of ideal peak position                               
!    p(3) = dqx  width of gaussian ellipsoid in x-direction             
!    p(4) = dqy   "         "        "       "  y   "                   
!    p(5) = dqz   "         "        "       "  z   "                   
!    p(6) = wavelength                                                  
!    p(7) = mphi = number of 32-pts integrations for phi-integration    
!    p(8) = ntht =   "         "          "       "  tht     "          
!    p(9) = thtstep = increment of tabel of integrated values           
!           for interpolation                                           
!    p(10) = print-interpolationtable                                   
!                                                                       
! ** xlim in common theo gives the limits of theta                      
!    between which the interpolation tabel is established **            
!                                                                       
      DATA rad / 0.01745329 /, pi / 3.1415926 / 
!                                                                       
 
! ** are parametrs (except height) changed **                           
      DO 10 i = 2, 9 
         IF (p (i) .ne.pold (i) ) goto 2000 
   10 END DO 
! ** ok no parametr change since last call ==>                          
!    take old interpolation table y**                                   
!                                                                       
 1000 CONTINUE 
! ** interpolation from table **                                        
      t0 = xlim (1) 
      ts = p (9) 
      n = (t - t0) / ts + 1 
      IF (n.gt.ntbl - 2) then 
         WRITE (6, 999) t, n, ntbl 
  999 FORMAT   (' qspot >>> call is out of limits ||'/                  &
     &          ' t =',f12.6/                                           &
     &          ' n =',i10/                                             &
     &          ' nt=',i10/)                                            
         qspot = 0.0 
         RETURN 
      ENDIF 
!                                                                       
      t1 = t0 + (n - 1) * ts 
      t2 = t1 + ts 
      t3 = t2 + ts 
      y1 = y (n) 
      y2 = y (n + 1) 
      y3 = y (n + 2) 
!                                                                       
      qspot = y1 * (t3 - t) * (t2 - t) / ( (t3 - t1) * (t2 - t1) )      &
      + y2 * (t1 - t) * (t3 - t) / ( (t1 - t2) * (t3 - t2) ) + y3 *     &
      (t1 - t) * (t2 - t) / ( (t1 - t3) * (t2 - t3) )                   
      qspot = qspot * p (1) * cos (t * rad / 2) 
!                                                                       
      RETURN 
!                                                                       
! ** calculate tabel of function by numerical integration **            
!                                                                       
 2000 CONTINUE 
      WRITE (6, * ) 'qspot generating table ...' 
      xlim (1) = x1 
      xlim (2) = x2 
!                                                                       
      zq0 = 4 * pi / p (6) 
      qhkl = zq0 * sin (p (2) * rad / 2) 
! ** width of ellipsoid **                                              
      a1 = 1 / p (3) **2 
      a2 = 1 / p (4) **2 
      a3 = 1 / p (5) **2 
! ** length of table **                                                 
      ntbl = (xlim (2) - xlim (1) ) / p (9) + 3 
      WRITE (6, * ) 'ntbl=', ntbl, xlim, p (9) 
      IF (ntbl.gt.300) then 
         WRITE (6, 9999) ntbl 
 9999 FORMAT (' qspot >>> length of table :',i10/                       &
     & ' exceeds max length (300) ')                                    
         qspot = 0 
         RETURN 
      ENDIF 
! ** store old parametervalues **                                       
      DO 2005 i = 1, 9 
 2005 pold (i) = p (i) 
! ** prepare integration-coefficients **                                
!                                                                       
! ** no. of 32 pts-integrations **                                      
      m = p (7) + 0.5 
      n = p (8) + 0.5 
      maxm = 320 / (2 * ngs) 
      IF (m.gt.maxm) m = maxm 
      IF (n.gt.maxm) n = maxm 
!                                                                       
! ** phi integration (0 - 2*pi) **                                      
      w = pi / m 
      iadd = 0 
      DO 2001 i = 1, m 
         DO 2001 j = 1, ngs 
            iadd = iadd+1 
            aphi (iadd) = w * a (j) 
            xh = w * x (j) + w + (i - 1) * 2 * w 
            sxphi (iadd) = sin (xh) 
            cxphi (iadd) = cos (xh) 
            iadd = iadd+1 
            aphi (iadd) = w * a (j) 
            xh = - w * x (j) + w + (i - 1) * 2 * w 
            sxphi (iadd) = sin (xh) 
            cxphi (iadd) = cos (xh) 
 2001 CONTINUE 
      mphi = iadd 
! ** tht integration (0 - pi) **                                        
      w = pi / (2 * n) 
      iadd = 0 
      DO 2002 i = 1, n 
         DO 2002 j = 1, ngs 
            iadd = iadd+1 
            atht (iadd) = w * a (j) 
            xh = w * x (j) + w + (i - 1) * 2 * w 
            sxtht (iadd) = sin (xh) 
            cxtht (iadd) = cos (xh) 
            iadd = iadd+1 
            atht (iadd) = w * a (j) 
            xh = - w * x (j) + w + (i - 1) * 2 * w 
            sxtht (iadd) = sin (xh) 
            cxtht (iadd) = cos (xh) 
 2002 CONTINUE 
      ntht = iadd 
!                                                                       
! ** loop over elements of table y **                                   
      DO 2010 i = 1, ntbl 
         th = (i - 1) * p (9) + xlim (1) 
         q = zq0 * sin (th * rad / 2) 
         sum = 0 
! ** loop for phi-integration **                                        
         DO 2020 j = 1, mphi 
            ap = aphi (j) 
            sp = sxphi (j) 
            cp = cxphi (j) 
! ** loop for theta-integration **                                      
            DO 2030 l = 1, ntht 
               st = sxtht (l) 
               ct = cxtht (l) 
               at = atht (l) 
               qx = q * st * sp 
               qy = q * st * cp 
               qz = q * ct - qhkl 
               arg = - (qx * a1 * qx + qy * a2 * qy + qz * a3 * qz) 
               IF (arg.gt. - 20) then 
                  sum = sum + exp (arg) * st * at * ap 
               ENDIF 
 2030       END DO 
 2020    END DO 
         y (i) = sum 
 2010 END DO 
!                                                                       
! ** printout of interpolationtable **                                  
! ** if p(10) = 1 **                                                    
!                                                                       
      ipr = p (10) + 0.5 
      IF (ipr.eq.1) then 
         WRITE (6, 3000) xlim, p, ntbl, (y (i), i = 1, ntbl) 
 3000 FORMAT (' qspot >> new interploation table <<'/    &
     & '    2theta-limits =',f12.4,' - ',f12.4,/         &
     & ' parametrs are :', 1x,10f12.6/                   &
     & ' length of table =',i5/                          &
     & ' table :'(1x,10f12.6) )                                                    
      ENDIF 
!                                                                       
! ** now interpolate value **                                           
      GOTO 1000 
!                                                                       
      END FUNCTION qspot                            
!                                                                       
!                                                                       
      BLOCKDATA 
      COMMON / gai / a (16), x (16), sxphi (320), cxphi (320), aphi (   &
      320), cxtht (320), sxtht (320), atht (320), ngs                   
! ** preliminary test with 6-pts. gauss-integration**                   
      DATA a / 0.1713245, 0.3607616, 0.4679139, 13 * 0. / 
      DATA x / 0.93247, 0.6612094, 0.23866192, 13 * 0. / 
      DATA ngs / 3 / 
! ** contains data for standard  6-pts gaussian integration **          
      END BLOCKDATA                                 
