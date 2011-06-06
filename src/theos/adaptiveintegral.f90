RECURSIVE FUNCTION adapint (f, a, b, epsilon, maxiter, erroraccu) 
!      =================================================                
!                                                                       
!      lokal adaptives integrationsverfahren 2te ver. wg. rekur.        
!                                                                       
      use outlev
      IMPLICIT real (8)(a - h, o - z) 
!                                                                       
      PARAMETER (maxstack = 200) 
!                 --------------->  stacktiefe                          
      DIMENSION s (0:maxstack), xa (0:maxstack), xb (0:maxstack) 
!                                                                       
!      LOGICAL cray 
!!     common/outlev/iot,ibild,ierrs,inka, cray                         
                                                                     !! 
!      REAL xxxx, yyyy, ptxf 
                                                                     !! 
!      COMMON / outlev / iot, ibild, ierrs, inka, iibuf, xxxx, yyyy,     &
!      ptxf (20)                                                         
!                                                                       
      EXTERNAL f 
!                                                                       
!      -----------------------------------------------------------------
                                                                        
                                                                        
      iterationcounter = 0 
      erg = 0.0d0 
      itop = 0 
      xa (itop) = a 
      xb (itop) = b 
      s (itop) = gaus8a (f, a, b) 
                                                                        
      DO i = 1, maxiter 
      IF (itop.ge. (maxstack - 2) ) then 
         erg = erg + s (itop) 
         xbb = xb (itop) 
         itm = itop - 1 
         DO itop = itm, 0, - 1 
         IF (xa (itop) .eq.xbb) goto 1 
         enddo 
    1    CONTINUE 
         WRITE (6, * ) 'warning! adaptint stack overflow!' 
      ELSE 
         iterationcounter = iterationcounter + 1 
         itop = itop + 1 
         xa (itop) = (xa (itop - 1) + xb (itop - 1) ) * 0.5d0 
         xb (itop) = xb (itop - 1) 
         s (itop) = gaus8a (f, xa (itop), xb (itop) ) 
         itop = itop + 1 
         xa (itop) = xa (itop - 2) 
         xb (itop) = xa (itop - 1) 
         s (itop) = gaus8a (f, xa (itop), xb (itop) ) 
         error = dabs ( (s (itop) + s (itop - 1) ) - s (itop - 2) ) 
                                                                        
!         IF (iot.gt.2) then 
!            WRITE (6, '(1x,i3,i5,4e13.6)') itop, iterationcounter, xa ( &
!            itop) , xb (itop) , (s (itop) + s (itop - 1) ) , error      
!         ENDIF 
                                                                        
         IF (error.lt.epsilon) then 
            erg = erg + s (itop) + s (itop - 1) 
            erroraccu = erroraccu + error 
            itop = itop - 1 
            xbb = xb (itop) 
            itop = itop - 1 
            itm = itop - 1 
            DO itop = itm, 0, - 1 
            IF (xa (itop) .eq.xbb) goto 2 
            enddo 
    2       CONTINUE 
         ENDIF 
      ENDIF 
      IF (itop.le.0) goto 3 
      enddo 
      WRITE (6, * ) 'adapint fatal error iterationnumber exceeded!' 
    3 CONTINUE 
                                                                        
                                                                        
      adapint = erg 
                                                                        
      RETURN 
                                                                        
      END FUNCTION adapint                          

RECURSIVE  FUNCTION a2dapint (f, a, b, epsilon, maxiter, erroraccu) 
!      =================================================                
!                                                                       
!      lokal adaptives integrationsverfahren 2te ver. wg. rekur.        
!                                                                       
      use outlev
      IMPLICIT real (8)(a - h, o - z) 
!                                                                       
      PARAMETER (maxstack = 200) 
!                 --------------->  stacktiefe                          
      DIMENSION s (0:maxstack), xa (0:maxstack), xb (0:maxstack) 
!                                                                       
!      LOGICAL cray 
                                                      !! aix            
!      REAL xxxx, yyyy, ptxf 
                                                                    !! a
!      COMMON / outlev / iot, ibild, ierrs, inka, iibuf, xxxx, yyyy,     &
!      ptxf (20)                                                         
!                                                                       
      EXTERNAL f 
!                                                                       
!      -----------------------------------------------------------------
                                                                        
                                                                        
      iterationcounter = 0 
      erg = 0.0d0 
      itop = 0 
      xa (itop) = a 
      xb (itop) = b 
      s (itop) = g2aus8a (f, a, b) 
                                                                        
      DO i = 1, maxiter 
      IF (itop.ge. (maxstack - 2) ) then 
         erg = erg + s (itop) 
         xbb = xb (itop) 
         itm = itop - 1 
         DO itop = itm, 0, - 1 
         IF (xa (itop) .eq.xbb) goto 1 
         enddo 
    1    CONTINUE 
         WRITE (6, * ) 'warning! a2daptint stack overflow!' 
      ELSE 
         iterationcounter = iterationcounter + 1 
         itop = itop + 1 
         xa (itop) = (xa (itop - 1) + xb (itop - 1) ) * 0.5d0 
         xb (itop) = xb (itop - 1) 
         s (itop) = g2aus8a (f, xa (itop), xb (itop) ) 
         itop = itop + 1 
         xa (itop) = xa (itop - 2) 
         xb (itop) = xa (itop - 1) 
         s (itop) = g2aus8a (f, xa (itop), xb (itop) ) 
         error = dabs ( (s (itop) + s (itop - 1) ) - s (itop - 2) ) 
                                                                        
!         IF (iot.gt.2) then 
!            WRITE (6, '(1x,i3,i5,4e13.6)') itop, iterationcounter, xa ( &
!            itop) , xb (itop) , (s (itop) + s (itop - 1) ) , error      
!         ENDIF 
                                                                        
         IF (error.lt.epsilon) then 
            erg = erg + s (itop) + s (itop - 1) 
            erroraccu = erroraccu + error 
            itop = itop - 1 
            xbb = xb (itop) 
            itop = itop - 1 
            itm = itop - 1 
            DO itop = itm, 0, - 1 
            IF (xa (itop) .eq.xbb) goto 2 
            enddo 
    2       CONTINUE 
         ENDIF 
      ENDIF 
      IF (itop.le.0) goto 3 
      enddo 
      WRITE (6, * ) 'a2dapint fatal error iterationnumber exceeded!' 
    3 CONTINUE 
                                                                        
                                                                        
      a2dapint = erg 
                                                                        
      RETURN 
                                                                        
      END FUNCTION a2dapint         
