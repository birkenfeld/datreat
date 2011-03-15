CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   FFTPACK 5.0 
C   Copyright (C) 1995-2004, Scientific Computing Division,
C   University Corporation for Atmospheric Research
C   Licensed under the GNU General Public License (GPL)
C
C   Authors:  Paul N. Swarztrauber and Richard A. Valent
C
C   $Id: cosqmb.f,v 1.2 2004/06/15 21:14:57 rodney Exp $
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE COSQMB (LOT, JUMP, N, INC, X, LENX, WSAVE, LENSAV, 
     1                   WORK, LENWRK, IER)
      INTEGER    LOT, JUMP, N, INC, LENX, LENSAV, LENWRK, IER
      REAL       X(INC,*), WSAVE(LENSAV), WORK(LENWRK)
      LOGICAL    XERCON
C
      IER = 0
C
      IF (LENX .LT. (LOT-1)*JUMP + INC*(N-1) + 1) THEN
        IER = 1
        CALL XERFFT ('COSQMB', 6)
        GO TO 300
      ELSEIF (LENSAV .LT. 2*N + INT(LOG(REAL(N))) +4) THEN
        IER = 2
        CALL XERFFT ('COSQMB', 8)
        GO TO 300
      ELSEIF (LENWRK .LT. LOT*N) THEN
        IER = 3
        CALL XERFFT ('COSQMB', 10)
        GO TO 300
      ELSEIF (.NOT. XERCON(INC,JUMP,N,LOT)) THEN
        IER = 4
        CALL XERFFT ('COSQMB', -1)
        GO TO 300
      ENDIF
C
      LJ = (LOT-1)*JUMP+1
      IF (N-2) 101,102,103
 101  DO 201 M=1,LJ,JUMP
      X(M,1) = X(M,1)
 201  CONTINUE
      RETURN
 102  SSQRT2 = 1./SQRT(2.)
      DO 202 M=1,LJ,JUMP
      X1 = X(M,1)+X(M,2)
      X(M,2) = SSQRT2*(X(M,1)-X(M,2))
      X(M,1) = X1
 202  CONTINUE
      RETURN
  103 CALL MCSQB1 (LOT,JUMP,N,INC,X,WSAVE,WORK,IER1)
      IF (IER1 .NE. 0) THEN
        IER = 20
        CALL XERFFT ('COSQMB',-5)
      ENDIF
C
  300 CONTINUE
      RETURN
      END
