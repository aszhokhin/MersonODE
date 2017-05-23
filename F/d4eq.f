*  THE GENERAL SOLUTION
      PROGRAM INT4G
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(KD=8000,LD=4)
      REAL*8 U(KD,LD),X(LD),X0(LD),Y(LD,LD),V(LD),W(LD)
      CHARACTER UK*7
      EXTERNAL D4EQ
      X0(1)=.8
      X0(2)=.9
      X0(3)=.8
      X0(4)=.9
      CALL OPN6
      T1=130.
      T=160.
      DO1 I=1,LD
    1 X(I)=X0(I)
      CALL SOLEQ(LD,KD,X,U,T1,T,D4EQ)
      CALL RESULTATES(LD,KD,T1,T,U)
      STOP
      END

      SUBROUTINE OPN6
      IMPLICIT REAL*8 (A-H,O-Z)
      OPEN(6,FILE='D4EQ.RES',STATUS='UNKNOWN')
      RETURN
      END

      SUBROUTINE RESULTATES(N,M,T1,T,U)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 U(M,N)
      OPEN(3,FILE='D4EQ.DAT',STATUS='UNKNOWN')
      DO1 I=1,M
      TI=(I-1)*(T-T1)/(M-1)+T1
    1 WRITE(3,2)TI,U(I,1),U(I,2),U(I,3),U(I,4)
    2 FORMAT(5E15.5)
      RETURN
      END

      SUBROUTINE D4EQ(T,X,F)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 F(4),X(4)
      F(2)=X(4)/(1+X(4))*X(1)/(1+X(1)+X(2))
      F(3)=X(2)*X(3)/(1+X(2)+X(3))
      F(1)=1-200*F(2)
      F(2)=400*F(2)-400*F(3)
      F(3)=230*F(3)-43.565*X(3)
      F(4)=1100*X(3)*X(1)/(1+X(1))/(1+X(2)*20)-20*X(4)
      RETURN
      END

      SUBROUTINE SOLEQ(N,M,X,U,Z,V,DTU)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(N),U(M,N)
      LOGICAL OK
      EXTERNAL DTU
      K=0
      A=.00001
      H=.01
      W=.00001
      B=5000.
      K=0
      L=1
      J=0
      T=0.
      CALL MERSON(T,Z,X,N,A,H,W,J,OK,DTU)
      DO2 I=1,N
    2 U(1,I)=X(I)
      S=(V-Z)/(M-1)
      T=Z
      DO4 L=2,M  
      T1=T+S    
      CALL MERSON(T,T1,X,N,A,H,W,J,OK,DTU)
      T=T1
      DO3 I=1,N
    3 U(L,I)=X(I)
    4 CONTINUE
      RETURN
      END

      SUBROUTINE MERSON(T,Q,Y,N,A,H,O,J,L,DTU)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 D(500)
      LOGICAL L
      REAL*8 Y(N)
      DATA R/1.E-13/
      L=.TRUE.
      N4=4*N
      N3=3*N
      N2=2*N
      N41=N4+1
      N31=N3+1
      DO1 K=1,N
    1 D(K+N4)=Y(K)
      Z=T
      S=H
      IS=0
    2 P=S
      C=Q-Z
      IF(ABS(S).LT.ABS(C)) GOTO3
      S=C
      IF(ABS(C/P).LT.R) GOTO11
      IS=1
    3 DO4 K=1,N
    4 D(K)=D(K+N4)
      F=S/3
      CALL DTU(Z,D(N41),D(N31))
      Z=Z+F
      DO5 K=1,N
      KN=K+N
      KN3=K+N3
      KN4=K+N4
      D(KN)=F*D(KN3)
    5 D(KN4)=D(KN)+D(K)
      CALL DTU(Z,D(N41),D(N31))

      DO6 K=1,N
      KN=K+N
      KN3=K+N3
      KN4=K+N4
      D(KN)=.5*D(KN)
    6 D(KN4)=.5*F*D(KN3)+D(KN)+D(K)
      CALL DTU(Z,D(N41),D(N31))

      Z=Z+.5*F
      DO7 K=1,N
      KN=K+N
      KN2=K+N2
      KN3=K+N3
      KN4=K+N4
      D(KN2)=4.5*F*D(KN3)
    7 D(KN4)=.25*D(KN2)+.75*D(KN)+D(K)
      CALL DTU(Z,D(N41),D(N31))
      Z=Z+.5*S
      DO8 K=1,N
      KN=K+N
      KN2=K+N2
      KN3=K+N3
      KN4=K+N4
      D(KN)=2*F*D(KN3)+D(KN)
    8 D(KN4)=3*D(KN)-D(KN2)+D(K)
      CALL DTU(Z,D(N41),D(N31))
      DO9 K=1,N
      KN=K+N
      KN2=K+N2
      KN3=K+N3
      KN4=K+N4
      D(KN2)=-.5*F*D(KN3)-D(KN2)+2*D(KN)
      D(KN4)=D(K+N4)-D(KN2)
      D(KN)=ABS(.5*A*D(KN4))
      D(KN2)=ABS(D(KN2))
      IF(ABS(D(KN4)).LE.R) GOTO9
      IF(D(KN2).GT.D(KN)) GOTO13
    9 CONTINUE
      IF(IS.EQ.1) GOTO11
      DO10 K=1,N
      KN=K+N
      KN2=K+N2
      IF(D(KN2).GT.D(KN)/32) GOTO2
   10 CONTINUE
      S=S*2
      GOTO2
   11 H=P
      T=Z
      DO12 K=1,N
      KN4=K+N4
   12 Y(K)=D(KN4)
      RETURN
   13 C=.5*S
      IF(ABS(C).GE.O) GOTO14
      IF(J.EQ.0) GOTO16
      S=O
      IF(P.LT.0.) S=-S
      IF(IS.EQ.1) GOTO11
      GOTO2
   14 DO15 K=1,N
      KN4=K+N4   
   15 D(KN4)=D(K)
      Z=Z-S
      S=C
      IS=0
      GOTO2
   16 L=.FALSE.
      GOTO11
      END

