      FUNCTION sofqq (x, pa, thnam, parnam, npar, ini, nopar ,params,napar,mbuf) 
!     ===================================================               
!                                                                       
!      sofq
!      calculates the structure factor according to
!      mode=   1       H.-P. form
!              2       Shieu-Chen form
!              3       Sharma-Sharma form
!              4       Critical diverg.
!              5       none (sq=1.)
!              6       PercusYevik                             
!   
      use theory_description
                                                                    
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
            sofqq = 0 
            RETURN 
         ENDIF 

        idesc = next_th_desc()
        th_identifier(idesc)   = thnam
        th_explanation(idesc)  =  &
         "      sofq  (typically used with attribute: multiply)         "//cr//&
         "      calculates the structure factor according to            "//cr//&
         "      mode=   1       Hayter-Penfold form                     "//cr//&
         "              2       Shieu-Chen form                         "//cr//&
         "              3       Sharma-Sharma form                      "//cr//&
         "              4       Critical diverg.                        "//cr//&
         "              5       none (sq=1.)                            "//cr//&
         "              6       PercusYevik                             "
  
       th_citation(idesc)     = "Modified version of the Hayter-Penfold program  bye E.Sheu      (Exxon Res. & Eng.)" 



         npar = nparx 
!        --------------> set the number of parameters                   
         parnam (1) = 'volfrac' 
         parnam (2) = 'scalelen' 
         parnam (3) = 'gamma' 
         parnam (4) = 'r' 
         parnam (5) = 'mode' 
!
         th_param_desc(1,idesc) = "volume fraction of (spherical particles) "  
         th_param_desc(2,idesc) = "characteristic length of the interaction screening length for mode=1,2"//cr//&
                        parspace//"the width of the potential well for mode=3"  
         th_param_desc(3,idesc) = "mode=1,2: gamma*exp(-r/scl) contact pot. "//cr//&    
                        parspace//"mode=3:   depth of the potential well in units of kT, + = repulsive, - =attractive"                
         th_param_desc(4,idesc) = "radius of the spheres"  
         th_param_desc(5,idesc) = "selector (integer) of S(Q) type.   DO NOT FIT!"  


         th_file_param(:,idesc) = " "

         th_out_param(:,idesc)  = " "
 
!                                                                       
         sofqq = 0 
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
      sofqq = sq 
!                                                                       
      RETURN 
      END FUNCTION sofqq          
