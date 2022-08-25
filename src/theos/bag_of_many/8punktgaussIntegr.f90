      FUNCTION g2aus8a (f, xu, xo) 
!-----------------------------------------------------------------------
!      8-punkte gauss-integration  : int(f,xu..xo)  
!-----------------------------------------------------------------------
      PARAMETER (ndim = 8) 
      PARAMETER (ndim2 = ndim / 2) 
      IMPLICIT real (8)(a - h, o - z) 
      DIMENSION a (ndim2), x (ndim2) 
      DATA a / 0.362683783378362d0, 0.313706645877887d0,                &
      0.222381034453374d0, 0.101228536290376d0 /                        
      DATA x / 0.183434642495650d0, 0.525532409916329d0,                &
      0.796666477413627d0, 0.960289856497536d0 /                        
                                                                        
      xave = (xo + xu) * 0.5d0 
      range = (xo - xu) * 0.5d0 
      sum = 0.d0 
      DO i = 1, ndim2 
      sum = sum + a (i) * (f (xave+range * x (i) ) + f (xave-range * x (&
      i) ) )                                                            
      enddo 
      g2aus8a = sum * range 
                                                                        
      RETURN 
      END FUNCTION g2aus8a  
		
		FUNCTION gaus8a (f, xu, xo) 
!-----------------------------------------------------------------------
!      8-punkte gauss-integration  : int(f,xu..xo)  2te kopie wg. rekur.
!-----------------------------------------------------------------------
      PARAMETER (ndim = 8) 
      PARAMETER (ndim2 = ndim / 2) 
      IMPLICIT real (8)(a - h, o - z) 
      DIMENSION a (ndim2), x (ndim2) 
      DATA a / 0.362683783378362d0, 0.313706645877887d0,                &
      0.222381034453374d0, 0.101228536290376d0 /                        
      DATA x / 0.183434642495650d0, 0.525532409916329d0,                &
      0.796666477413627d0, 0.960289856497536d0 /                        
                                                                        
      xave = (xo + xu) * 0.5d0 
      range = (xo - xu) * 0.5d0 
      sum = 0.d0 
      DO i = 1, ndim2 
      sum = sum + a (i) * (f (xave+range * x (i) ) + f (xave-range * x (&
      i) ) )                                                            
      enddo 
      gaus8a = sum * range 
                                                                        
      RETURN 
      END FUNCTION gaus8a   