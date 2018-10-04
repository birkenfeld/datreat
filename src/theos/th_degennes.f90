      FUNCTION th38 (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf) 
!  degennes 
!                                                                       
      use theory_description                                                                   
      IMPLICIT none 
                                                                        
      DOUBLEPRECISION Pi 
      Parameter (Pi = 3.141592654d0) 
                                                                        
                                                                        
      CHARACTER(8) thnam, parnam (20) 
      REAL th38, x, pa, qq, zpi, xh, vol_frac 
      INTEGER ini, npar, nparx 
      DIMENSION pa (20), qq (3) 
      integer :: mbuf
      integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
      real, intent(inout) :: params(mbuf)             ! value des parameters n
      DATA zpi / 6.283185 / 
                                                                        
      DOUBLE PRECISION :: w, d, q, n, ne, degennes, t, a, l 
                                                                        
                                                                        
                                                                        
      INTEGER  ier, iot 
                                                                        
!                                                                       
! ----- initialisation -----                                            
      IF (ini.eq.0) then 
         thnam = 'degennes ' 
         nparx = 5 
         IF (npar.lt.nparx) then 
            WRITE (6, 1) thnam, nparx, npar 
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)                              
            th38 = 0 
            RETURN 
         ENDIF

        idesc = next_th_desc()
        th_identifier(idesc)   = thnam
        th_explanation(idesc)  =  &
         "deGennes expression for reptation S(q,t):                   "//cr//&
         " with tau0  = 36 / (W * ( (q * l) **4) )                    "//cr//&
         "      taud  = 3 * (n**3) * (l**2) / ( (pi**2) * W * (d**2)) "//cr//&
         "      td    = t / taud ;  t0 = t / tau0                     "//cr//&
         "      eqd   = exp ( - ( (q * d / 6) **2) )                  "//cr//&
         "      sum   = sum(i) of { exp ( - i * i * td) / (i * i) }   "//cr//&
         "      creep = sum * (8 / (pi**2) ) * eqd                    "//cr//&
         " degennes = (1 - eqd) * exp (t0) * erfc (sqrt (t0) ) + creep"
  
       th_citation(idesc)     = " e.g. P.Schleger et. al. PRL 81, p.124 (1998) and references therein" 


 
         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'amplitud' 
         parnam (2)  = 'w       ' 
         parnam (3)  = 'dtube   ' 
         parnam (4)  = 'n       ' 
         parnam (5)  = 'l       ' 
! 
         th_param_desc(1,idesc) = "prefactor "  
         th_param_desc(2,idesc) = "Rouse rate (the W part of wl4) typically in units of 1/ns"  
         th_param_desc(3,idesc) = "repation tube diameter (in THE LENGTH UNITS for q, l..., typ. Angstroem) "  
         th_param_desc(4,idesc) = "number of segments in the cains (only infulences the creep term, taud) "  
         th_param_desc(5,idesc) = "segment length in the length units of q, d"  

         th_file_param(:,idesc) = " "
         th_file_param(1,idesc) = " q  =  q value of the S(q,t), units THE LENGTH UNITS"

         th_out_param(:,idesc)  = " "
         th_out_param(1,idesc)  = " ne = (dtube/l)**2 = number of segments l in entanglement"
                                                                    
         th38 = 0 
         RETURN 
      ENDIF 
!                                                                       
! ---- calculate theory here -----                                      
                          ! amplitude (should be 1)                     
      a = pa (1) 
                          ! Rouse rate                                  
      w = pa (2) 
                          ! Tube diameter                               
      d = pa (3) 
                          ! no of segments of polymer                   
      n = pa (4) 
                          ! segment length                              
      l = pa (5) 
                                                                        
      CALL getpar ('q       ', xh,nopar ,params,napar,mbuf, ier)  
      IF (ier.eq.0) then 
         q = xh 
      ELSE 
         q = 1.0 
      ENDIF 
                                                                        
      t = x 
                                                                        
      th38 = degennes (t, q, d, W, n, l, ne) 
                                                                        
      th38 = th38 * a 
                                                                        
      CALL setpar ('ne      ', sngl (ne) ,nopar ,params,napar,mbuf, ier)  
                                                                        
      RETURN 
      END FUNCTION th38                             
                                                                        
                                                                        
      DOUBLEPRECISION function degennes (t, q, d, W, n, l, ne) 
!      ----------------------------------------------------             
!                                                                       
! siehe P.Schleger et. al. PRL 81, p.124 (1998)   
!
! NOTE:
! NOTE:
! NOTE: 
!  > update the description provided in th_secription in th_degennes 
!  > init section!                      
!                                                                       
!                                                                       
      IMPLICIT none 
      DOUBLEPRECISION t, q, d, W 
      DOUBLEPRECISION derfc, pi 
      DOUBLEPRECISION tau0, taud, t0, td 
      DOUBLEPRECISION n, ne, l, y, sum, eqd 
      INTEGER m, i, j 
                                                                        
      PARAMETER (pi = 3.141592654d0) 
                                                                        
!      l    =   d/sqrt(ne)                                              
      ne = (d / l) **2 
      tau0 = 36 / (W * ( (q * l) **4) ) 
      taud = 3 * (n**3) * (l**2) / ( (pi**2) * W * (d**2) ) 
                                                                        
      td = t / taud 
      t0 = t / tau0 
      eqd = exp ( - ( (q * d / 6) **2) )

!      write(6,*)'td,t0,eqd ',td, t0, eqd
                                                                        
      m = 10 / sqrt (td+0.001D0) + 2 
      sum = 0 
      DO i = 1, m, 2 
      sum = sum + exp ( - i * i * td) / (i * i) 
      enddo 
      sum = sum * (8 / (pi**2) ) * eqd
!     write(6,*)'m, sum ',m,sum
                                                                        
                                                                        
      degennes = (1 - eqd) * exp (t0) * erfc (sqrt (t0) ) + sum 

!      write(6,*)'degennes ',degennes, erfc(sqrt(t0))
      
      RETURN 
      END FUNCTION degennes                         
