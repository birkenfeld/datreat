CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   FFTPACK 5.0 
C   Copyright (C) 1995-2004, Scientific Computing Division,
C   University Corporation for Atmospheric Research
C   Licensed under the GNU General Public License (GPL)
C
C   Authors:  Paul N. Swarztrauber and Richard A. Valent
C
C   $Id: cost1b.f,v 1.2 2004/06/15 21:14:57 rodney Exp $
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE COST1B ( N, INC, X, LENX, WSAVE, LENSAV, 
     1                   WORK, LENWRK, IER)
      INTEGER    N, INC, LENX, LENSAV, LENWRK, IER
      REAL       X(INC,*), WSAVE(LENSAV), WORK(LENWRK)
C
      IER = 0
      IF (LENX .LT. INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('COST1B', 6)
        GO TO 100
      ELSEIF (LENSAV .LT. 2*N + INT(LOG(REAL(N))) +4) THEN
        IER = 2
        CALL XERFFT ('COST1B', 8)
        GO TO 100
      ELSEIF (LENWRK .LT. N-1) THEN
        IER = 3
        CALL XERFFT ('COST1B', 10)
        GO TO 100
      ENDIF
C
      CALL COSTB1 (N,INC,X,WSAVE,WORK,IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('COST1B',-5)
      ENDIF
C
  100 CONTINUE
      RETURN
      END
