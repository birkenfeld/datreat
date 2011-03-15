CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   FFTPACK 5.0 
C   Copyright (C) 1995-2004, Scientific Computing Division,
C   University Corporation for Atmospheric Research
C   Licensed under the GNU General Public License (GPL)
C
C   Authors:  Paul N. Swarztrauber and Richard A. Valent
C
C   $Id: sinq1i.f,v 1.2 2004/06/15 21:29:20 rodney Exp $
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE SINQ1I (N, WSAVE, LENSAV, IER)
      INTEGER    N, LENSAV, IER
      REAL       WSAVE(LENSAV)
C
      IER = 0
C
      IF (LENSAV .LT. 2*N + INT(LOG(REAL(N))) +4) THEN
        IER = 2
        CALL XERFFT ('SINQ1I', 3)
        GO TO 300
      ENDIF
C
      CALL COSQ1I (N, WSAVE, LENSAV, IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('SINQ1I',-5)
      ENDIF
  300 RETURN
      END
