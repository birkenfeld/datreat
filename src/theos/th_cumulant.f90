      FUNCTION th_cumulant (x, pa, thnam, parnam, npar,ini, nopar ,params,napar,mbuf)
!     ===================================================
!  cumulant analysis
!  cumulant= ampli * exp ( - sum_i mui * x^i/i! )
!
      use theory_description
      CHARACTER(8) thnam, parnam (20)
      DIMENSION pa (20), qq (3)
      			integer :: mbuf
			integer, intent(inout) :: nopar                 ! Anzahl der Parameter data
      character*80, intent(inout) :: napar(mbuf)      ! name des parameters n
			real, intent(inout) :: params(mbuf)             ! value des parameters n
			REAL mu, ampli, sum, ff
      real   :: th_cumulant
      DIMENSION mu (20)
      DATA zpi / 6.283185 /
!
! ----- initialisation -----
      IF (ini.eq.0) then
         thnam = 'cumulant'
         nparx = 7
         IF (npar.lt.nparx) then
            WRITE (6, 1) thnam, nparx, npar
    1 FORMAT     (' theory: ',a8,' no of parametrs=',i8,                &
     &      ' exceeds current max. = ',i8)
            th_cumulant= 0
            RETURN
         ENDIF
         npar = nparx
! >>>>> describe theory with    7 parameters >>>>>>>
        idesc = next_th_desc()
        th_identifier(idesc)   = thnam
        th_explanation(idesc)  = "cumulant expression                                "//cr//&
                                 "th = A * exp( -sum[i=1..6](mu_i/i! * t**i))        "!
 !
        th_citation(idesc) = "  "
!        --------------> set the number of parameters
         parnam (1) = 'amplitu'
         parnam (2) = 'mu1'
         parnam (3) = 'mu2'
         parnam (4) = 'mu3'
         parnam (5) = 'mu4'
         parnam (6) = 'mu5'
         parnam (7) = 'mu6'
! >>>>> describe parameters >>>>>>>
        th_param_desc( 1,idesc) = " amplitude prefactor A " !//cr//parspace//&
        th_param_desc( 2,idesc) = " mu_1                  " !//cr//parspace//&
        th_param_desc( 3,idesc) = " mu_2                  " !//cr//parspace//&
        th_param_desc( 4,idesc) = " mu_3                  " !//cr//parspace//&
        th_param_desc( 5,idesc) = " mu_4                  " !//cr//parspace//&
        th_param_desc( 6,idesc) = " mu_5                  " !//cr//parspace//&
        th_param_desc( 7,idesc) = " mu_6                  " !//cr//parspace//&
! >>>>> describe record parameters used >>>>>>>
       th_file_param(:,idesc) = " " 
 ! >>>>> describe record parameters creaqted by this theory >>>>>>>
       th_out_param(:,idesc)  = " "  
 !
         th_cumulant= 0
         RETURN
      ENDIF
!
! ---- calculate theory here -----
      ampli = pa (1)
      DO i = 2, nparx
         mu (i + 1 - 2) = pa (i)
      enddo

      ff = 1.0
      sum = 0.0
      DO i = 1, 6
        ff = ff * i
        sum = sum + mu (i) * (x**i) / ff
      enddo

      th_cumulant= ampli * exp ( - sum)
!
      RETURN
      END FUNCTION th_cumulant
