      PROGRAM FITBRANDT

C General case of fit for correlated variables 
C Both variables with errors
C Uses procedure LSQGEN (Siegmund Brandt, Data Analysis)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(MAXN=2000,MAXNR=4,MAXM=1000)
      DIMENSION X(MAXNR),Y(MAXN),CY(MAXN,MAXN)
      DIMENSION GY(MAXN,MAXN),CX(MAXNR,MAXNR)
      DIMENSION F(MAXN+MAXNR,MAXN+MAXNR)
      DIMENSION E(MAXM,MAXN+MAXNR)
      DIMENSION A2(MAXN+MAXNR,MAXN+MAXNR),LIST(4)
      DIMENSION T(1000),S(1000),DT(1000),DS(1000),RHO(1000),OT(1000)
      DIMENSION u(20000),os(20000)

      INTEGER it(20000)

      CHARACTER*80 c_line, infile

      COMMON /oshape/ u,os,it

      OPEN(21,FILE='fitbrandt.res')
      OPEN(22,FILE='fitbrandt.dat')
      OPEN(23,FILE='fitbrandt.test')

      NP = 0
      NRD = 0

      ersk1 = 2.0
      ersk2 = 0.01
      roh1 = 0.0
      roh2 = 0.0

      CALL GETARG(1,c_line)
      WRITE(*,*) c_line
      READ(c_line,*) infile
      CALL GETARG(2,c_line)
      READ(c_line,*) NP1
      CALL GETARG(3,c_line)
      READ(c_line,*) NP2

C     Initial signal - hybrid gaussian
      x0 = 98.0
      ao = 280.0
      xpa = x0 + 4.82
      patif = 1.4
      tau = 0.9
      dx = 0.01

      DO 1 i = 1, 10001
        u(i) = x0 + dx*REAL(i-1)
        fas = (u(i)-xpa)*tau
        ex2 = 0.0
        IF (2*patif*patif + fas > 0.0) THEN
          ex2 = dexp(-(u(i)-xpa)*(u(i)-xpa)/(2.0*patif*patif + fas))
        END IF
        os(i) = ao*ex2
        WRITE(23,*) i,u(i),os(i)
  1   CONTINUE     


      open(11,file=infile)

      do 7 i = 1, 1000

      read(11,*,end=8)  TTEMP, STEMP

C      if(ISK.eq.0) go to 4

      NRD = NRD + 1

      IF (NRD .LT. NP1 .OR. NRD .GT. NP2) GO TO 9

      NP = NP + 1

      T(NP) = TTEMP

      S(NP) = STEMP

      DT(NP) = 0.5

      DS(NP) = 0.05*STEMP

      RHO(NP) = 0.0

  9   continue
  7   continue
  8   continue

C     Find integration limit for each input time
      DO 2 i = 1, NP
        DO 3 k = 1, 10001
          IF ( T(i) .GE. u(k) .AND. T(i) .LT. u(k+1) ) THEN
            it(i) = k
            IF ( DABS(T(i)-u(k)) .GT. DABS(u(k+1)-T(i)) ) it(i) = k+1
          END IF
  3     CONTINUE
      WRITE(*,'(I3,1X,F10.3,2X,I4,1X,F10.3)') i,T(i),it(i),u(it(i))
  2   CONTINUE

      WRITE(*,*) it(315),u(315),os(315)

C identify program to user
      WRITE(21,*)' Program E7LSQ demonstrates use of LSQGEN'
      WRITE(21,*)' '
C write table of data
      WRITE(21,*)'    T         S         DT        DS        RHO'
      WRITE(*,*)'    T         S         DT        DS        RHO'
      DO 5 I=1,NP
        WRITE(21,'(5F10.3)')T(I),S(I),DT(I),DS(I),RHO(I)
        WRITE(*,'(5F10.3)')T(I),S(I),DT(I),DS(I),RHO(I)
5     CONTINUE
      WRITE(2,'(/)')
C set up data for input to LSQGEN
      N=2*NP
      M=NP
      NR=4
      DO 20 NRED=4,1,-1
        CALL MTXUNT(CY,N)
        DO 10 I=1,NP
          Y((I-1)*2+1)=T(I)
          Y((I-1)*2+2)=S(I)
          CY((I-1)*2+1,(I-1)*2+1)=DT(I)**2
          CY((I-1)*2+2,(I-1)*2+2)=DS(I)**2
          CY((I-1)*2+1,(I-1)*2+2)=RHO(I)*DS(I)*DT(I)
          CY((I-1)*2+2,(I-1)*2+1)=RHO(I)*DS(I)*DT(I)
10      CONTINUE
C determine first approximation
        IF(NRED.EQ.4) THEN
          LIST(1)=1
          LIST(2)=1
          LIST(3)=1
          LIST(4)=1
          X(1) = 2.5
          X(2) = 2.0
          X(3) = 0.034
          X(4) = 3.3
        ELSE IF(NRED.EQ.3) THEN
          LIST(1)=1
          LIST(2)=1
          LIST(3)=1
          LIST(4)=0
          X(1) = 2.5
          X(2) = 2.0
          X(3) = 1.0
          X(4) = 2.7
        ELSE IF(NRED.EQ.2) THEN
          LIST(1)=1
          LIST(2)=1
          LIST(3)=0
          LIST(4)=0
          X(1) = 2.5
          X(2) = 2.0
          X(3) = 1.0
          X(4) = 2.7
        ELSE IF(NRED.EQ.1) THEN
          LIST(1)=1
          LIST(2)=0
          LIST(3)=0
          LIST(4)=0
          X(1) = 2.5
          X(2) = 2.0
          X(3) = 1.0
          X(4) = 2.7
        END IF
C header for output of results
        WRITE(21,*)' performing fit with LSQGEN'
        WRITE(21,'(3(A,I2),A,4I2)')' N = ',N,', NR = ',NR,
     +  ', NRED = ',NRED,', LIST = ',LIST
        WRITE(21,'(A,4(G12.5,2x))')' first approx.: X = ',X

        WRITE(*,*)' performing fit with LSQGEN'
        WRITE(*,'(3(A,I2),A,4I2)')' N = ',N,', NR = ',NR,
     +  ', NRED = ',NRED,', LIST = ',LIST
        WRITE(*,'(A,4(G12.5,2x))')' first approx.: X = ',X

        NSTEP=100
        CALL LSQGEN(Y,CY,GY,F,E,M,N,NR,NRED,LIST,X,CX,R,A2,NSTEP)
C output of results
        WRITE(21,'(A,G12.4,A,I3,/,A,(4(G12.5,2x)))')
     +  ' result of fit: R = ',R,', NSTEP =',NSTEP,' X =',X

        WRITE(*,'(A,G12.4,A,I3,/,A,(4(G12.5,2x)))')
     +  ' result of fit: R = ',R,', NSTEP =',NSTEP,' X =',X

        IF(NRED.GT.0) THEN
          WRITE(21,*)' covariance matrix CX = '
          WRITE(*,*)' covariance matrix CX = '
          CALL MTXWRT(CX,NRED,NRED)

          sig1 = sqrt(CX(1,1))
          sig2 = sqrt(CX(2,2))
          sig3 = sqrt(CX(3,3))
          sig4 = sqrt(CX(4,4))

          coref12 = CX(2,1)/(sig1*sig2)
          coref13 = CX(3,1)/(sig1*sig3)
          coref23 = CX(3,2)/(sig2*sig3)

          write(21,121) sig1, sig2, sig3, sig4
          write(21,122) coref12, coref13, coref23

        END IF
        WRITE(21,'(/)')
        WRITE(*,'(/)')

      IF (NRED .EQ. 4) THEN
      X(1) = 2.5
      X(2) = 2.0
      X(3) = 0.034
      X(4) = 3.3
      DO 30 j=1,10001
        oco = 0.0
        DO 29 i = 51, j
          xu = u(i)-98.5
C          etau = (X(1)+X(2)*xu+X(3)*xu**2.0)*DEXP(-xu/X(4))
          etau = (X(1)+X(2)*xu)*DEXP(-xu/X(4))
          oco = oco + dx*etau*os(j-i+51)
C      IF (j .EQ. 222) WRITE(*,*) i,u(i),xu,etau,oco
  29    CONTINUE
C        funco = os(j) - oco
        funco = X(3)*oco
        WRITE(22,*) u(j),funco
30    CONTINUE       

      END IF

20    CONTINUE


  121 format(1x,'SIGMA1 = ',g12.5,2x,'SIGMA2 = ',g12.5,2x,'SIGMA3 = ',g1
     &2.5,2x,'SIGMA4 = ',g12.5)
  122 format(1x,'CORR.COEF.12 = ',g12.5,2x,'CORR.COEF13 = ',g12.5,2x,'CO
     &RR.COEF.23 = ',g12.5)

      END


      SUBROUTINE LSQGEN(Y,CY,FY,F,E,M,N,NR,NRED,LIST,X,CX,R,
     +                  A2,NSTEP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION LSQGFN
      PARAMETER(LMAX=1050)
      DIMENSION Y(N),CY(N,N),X(NR),CX(NRED,NRED),LIST(NR)
      DIMENSION E(M,N+NRED),FY(N,N),F(N+NRED,N+NRED),A2(N+NRED,*)
      DIMENSION D(LMAX),T(LMAX),B(LMAX),U(LMAX)
      COMMON /DASV03/ D,T,B,U
      PARAMETER(MAXSTP=100,EPSILN=1.D-8,TT=1.D-15,ZERO=0.D0)
      LOGICAL OK,COVMAT
C general case of least squares fitting
      OK=.TRUE.
      COVMAT=.TRUE.
      IF(NSTEP.LT.0) THEN
        COVMAT=.FALSE.
        NSTEP=ABS(NSTEP)
      END IF
      IF(NSTEP.LT.1) NSTEP=MAXSTP
      IF(NR.GT.0) THEN
C For NR=NRED :  set LIST
        IF(NR.EQ.NRED) THEN
          DO 10 I=1,NR
            LIST(I)=1
   10     CONTINUE
        END IF
      END IF
      L=N+NRED
      CALL MTXZRV(T,L)
      CALL MTXTRA(CY,F,N,N)
      CALL MTXCHI(CY,FY,N)
      CALL MTXCHL(CY,FY,N)
      CALL MTXTRA(F,CY,N,N)
C start iteration
      R=ZERO
      DO 50 ISTEP=1,NSTEP
        WRITE(*,*) ISTEP,X
        RLST=R
        CALL MTXZER(F,L,L)
        CALL MTXPSM(F,FY,L,L,N,N,NRED+1,NRED+1)
        DO 20 K=1,M
          D(K)=-LSQGFN(Y,X,N,NR,K)
   20   CONTINUE
C Numerical Derivatives
        CALL AUXDRG(X,Y,M,N,NR,NRED,LIST,E,OK)
        IF(.NOT.OK) THEN
          NSTEP=-3
          GO TO 60
        END IF
        CALL MTXCHM(F,T,B,L,1)
        CALL MTXMSV(B,B,DBLE(-1.),L)
        CALL MTXLSC(F,B,E,D,U,R,A2,L,L,M,DBLE(0.),OK)
        WRITE(*,*) 'After MTXLSC', U(1),U(2),U(3),U(4),OK
        IF(.NOT.OK) THEN
          NSTEP=-1
          GO TO 60
        END IF
        IF(NRED.GT.0) THEN
          IRED=0
          DO 30 I=1,NR
            IF(LIST(I).NE.0) THEN
              IRED=IRED+1
        WRITE(*,*) X(I),U(IRED)
              X(I)=X(I)+U(IRED)
            END IF
   30     CONTINUE
        END IF
        DO 40 I=1,N
          Y(I)=Y(I)+U(I+NRED)
          T(I+NRED)=T(I+NRED)+U(I+NRED)
   40   CONTINUE
C test for convergence
        IF(ISTEP.GT. 1.AND. ABS(R-RLST).LT.EPSILN*R+TT) THEN
          NSTEP=ISTEP
          IF(COVMAT) THEN
C compute matrix GB
            CALL AUXDRG(X,Y,M,N,NR,NRED,LIST,E,OK)
            IF(.NOT.OK) THEN
              NSTEP=-3
              GO TO 60
            END IF
            CALL MTXGSM(E,A2,M,L,M,N,1,NRED+1)
            CALL MTXMBT(CY,A2,F,N,N,M)
            CALL MTXMLT(A2,F,FY,M,N,M)
            CALL MTXCHI(FY,F,M)
C array FY now contains matrix GB
            IF(NRED.GT.0) THEN
              CALL MTXGSM(E,A2,M,L,M,NRED,1,1)
              CALL MTXMAT(A2,FY,F,NRED,M,M)
              CALL MTXMLT(F,A2,CX,NRED,M,NR)
              CALL MTXCHI(CX,F,NRED)
C array CX now contains covariance matrix of unknowns
            ELSE
              CALL MTXMLT(A2,CY,F,M,N,N)
              CALL MTXMLT(FY,F,A2,M,M,N)
              CALL MTXMAT(F,A2,FY,N,M,N)
              CALL MTXSUB(CY,FY,CY,N,N)
C array CY now contains covariance matrix of 'improved'
C measurements
            END IF
          END IF
          GO TO 60
        END IF
   50 CONTINUE
      NSTEP=-2
   60 RETURN
      END


      DOUBLE PRECISION FUNCTION LSQGFN(ETA,X,N,NR,K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(NR),ETA(N)
      DIMENSION u(20000),os(20000)
      INTEGER it (20000)
C     Signal shape model after filter
      COMMON /oshape/ u,os,it
      dx = 0.01
      oco = 0.0
C      WRITE(*,*) it(315),u(315),os(315)
      DO 10 i = 51, it(K)
        xu = u(i)-98.5
C        etau = (X(1)+X(2)*xu+X(3)*xu**2.0)*DEXP(-xu/X(4))
        etau = (X(1)+X(2)*xu)*DEXP(-xu/X(4))
        oco = oco + dx*etau*os(it(K)-i+51)
C      IF (K .EQ. 3) WRITE(*,*) i,u(i),xu,os(it(K)-i+51),etau,oco
  10  CONTINUE
      funco = X(3)*oco
C      WRITE(*,*) it(315),u(315),os(315)

      LSQGFN = ETA(2*K) - funco
C      WRITE(*,*) K,ETA(2*K-1),ETA(2*K),it(K),os(it(K)),oco,funco,LSQGFN
C      WRITE(*,*) X(1),X(2),X(3),X(4)
C      LSQGFN=ETA(2*K)-X(1)*DEXP(-(ETA(2*K-1)-X(2))*(ETA(2*K-1)-X(2))/(2.
C     *D0*X(3)*X(3)))
      END


      SUBROUTINE MTXZRV(U,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(N)
      DO 10 I=1,N
        U(I)=0.
   10 CONTINUE
      END


      SUBROUTINE MTXTRA(A,R,M,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,N),R(M,N)
      DO 20 J=1,N
        DO 10 I=1,M
          R(I,J)=A(I,J)
   10   CONTINUE
   20 CONTINUE
      END


      SUBROUTINE MTXCHI(A,U,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N,N),U(N,N)
C Step 1: Cholesky decomposition
      CALL MTXCHL(A,U,N)
      DO 50 I=1,N
C Step 2: Forward Substitution
        DO 20 L=I,N
          IF(L.EQ.I) THEN
            A(N,L)=DBLE(1.)/U(L,L)
          ELSE
            A(N,L)=DBLE(0.)
            DO 10 K=I,L-1
              A(N,L)=A(N,L)-U(K,L)*A(N,K)
   10       CONTINUE
            A(N,L)=A(N,L)/U(L,L)
          END IF
   20   CONTINUE
C Step 3: Back Substitution
        DO 40 L=N,I,-1
          IF(L.EQ.N) THEN
            A(I,L)=A(N,L)/U(L,L)
          ELSE
            A(I,L)=A(N,L)
            DO 30 K=N,L+1,-1
              A(I,L)=A(I,L)-U(L,K)*A(I,K)
   30       CONTINUE
            A(I,L)=A(I,L)/U(L,L)
          END IF
   40   CONTINUE
   50 CONTINUE
C Fill lower triangle symmetrically
      IF(N.GT.1) THEN
        DO 70 I=1,N
          DO 60 L=1,I-1
            A(I,L)=A(L,I)
   60     CONTINUE
   70   CONTINUE
      END IF
      END


      SUBROUTINE MTXCHL(A,U,N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(N,N),U(N,N)
      CALL MTXZER(U,N,N)
      DO 30 K=1,N
        S=0.
        DO 20 J=K,N
          IF(K.GT.1) THEN
            S=0.
            DO 10 L=1,K-1
              S=S+U(L,K)*U(L,J)
   10       CONTINUE
          END IF
          U(K,J)=A(K,J)-S
          IF(K.EQ.J) THEN
            U(K,J)=SQRT(ABS(U(K,J)))
          ELSE
            U(K,J)=U(K,J)/U(K,K)
          END IF
   20   CONTINUE
   30 CONTINUE
      END


      SUBROUTINE MTXZER(R,N,M)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(N,M)
      DO 20 J=1,M
        DO 10 I=1,N
          R(I,J)=0.
   10   CONTINUE
   20 CONTINUE
      END


      SUBROUTINE MTXPSM(A,S,M,N,K,L,M1,N1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,N), S(K,L)
      DO 20 I=1,K
        DO 10 J=1,L
          A(M1-1+I,N1-1+J)=S(I,J)
   10   CONTINUE
   20 CONTINUE
      END


      SUBROUTINE AUXDRG(X,ETA,MM,N,NR,NRED,LIST,E,OK)
C Computes the derivative f'(x) of f(x) at x = X. Based on
C H. Rutishauser, Ausdehnung des Rombergschen Prinzips
C (Extension of Romberg's Principle), Numer. Math. 5 (1963) 48-54
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (EPS = 5.D-8, EPSI=1.D-10, DELTA=10.D0, S=1.D-1)
      DIMENSION X(NR),ETA(N),E(MM,NRED+N),LIST(NR)
      DIMENSION DX(0:9),W(0:9,3),T(0:9,0:9),A(0:9)
      LOGICAL LEV(0:9),LMT
      LOGICAL OK,LX
      DOUBLE PRECISION LSQGFN
      DIMENSION u(20000),os(20000)
      INTEGER it (20000)
      COMMON /oshape/ u,os,it
      DATA DX /0.0256D0, 0.0192D0, 0.0128D0, 0.0096D0, 0.0064D0,
     +0.0048D0, 0.0032D0, 0.0024D0, 0.0016D0, 0.0012D0/
      DATA (LEV(K),K=0,8,2) /5*.TRUE./
      DATA (LEV(K),K=1,9,2) /5*.FALSE./
      DATA W(1,1) /1.33333 33333 333333D+00/
      DATA W(3,1) /1.06666 66666 666667D+00/
      DATA W(5,1) /1.01587 30158 730159D+00/
      DATA W(7,1) /1.00392 15686 274510D+00/
      DATA W(2,1) /3.33333 33333 333333D-01/
      DATA W(4,1) /6.66666 66666 666667D-02/
      DATA W(6,1) /1.58730 15873 015873D-02/
      DATA W(8,1) /3.92156 86274 509804D-03/
      DATA W(0,2) /2.28571 42857 142857D+00/
      DATA W(2,2) /1.16363 63636 363636D+00/
      DATA W(4,2) /1.03643 72469 635628D+00/
      DATA W(6,2) /1.00886 69950 738916D+00/
      DATA W(8,2) /1.00220 21042 329337D+00/
      DATA W(1,2) /1.28571 42857 142857D+00/
      DATA W(3,2) /1.63636 36363 636364D-01/
      DATA W(5,2) /3.64372 46963 562753D-02/
      DATA W(7,2) /8.86699 50738 916256D-03/
      DATA W(9,2) /2.20210 42329 336922D-03/
      DATA W(0,3) /1.80000 00000 000000D+00/
      DATA W(2,3) /1.12500 00000 000000D+00/
      DATA W(4,3) /1.02857 14285 714286D+00/
      DATA W(6,3) /1.00699 30069 930070D+00/
      DATA W(8,3) /1.00173 91304 347826D+00/
      DATA W(1,3) /8.00000 00000 000000D-01/
      DATA W(3,3) /1.25000 00000 000000D-01/
      DATA W(5,3) /2.85714 28571 428571D-02/
      DATA W(7,3) /6.99300 69930 069930D-03/
      DATA W(9,3) /1.73913 04347 826087D-03/
      OK=.TRUE.
C      WRITE(*,*) it(3),u(315),os(315)
      L=NR+N
      DO 90 IM=1,MM
        I2=0
        DO 80 IL=1,L
          DEL=DELTA
          IF(IL.LE.NR) THEN
            IF(LIST(IL).NE.0) THEN
              I2=I2+1
              LX=.TRUE.
              XSAV=X(IL)
            ELSE
              GO TO 80
            END IF
          ELSE
            I2=I2+1
            LX=.FALSE.
            IY=IL-NR
            XSAV=ETA(IY)
          END IF
          DO 40 I=1,10
            DEL=S*DEL
            WRITE(*,*) X
            IF(I.EQ. 10.OR. ABS(XSAV+DEL*DX(9)-XSAV) .LT. EPS) THEN
            WRITE(*,*) I, DEL, DX(9),XSAV,EPS,OK
             OK=.FALSE.
              RETURN
            END IF
            DO 10 K = 0,9
              H=DEL*DX(K)
              IF(LX) THEN
                X(IL)=XSAV+H
              ELSE
                ETA(IY)=XSAV+H
              END IF
              FPLUS=LSQGFN(ETA,X,N,NR,IM)
              IF(LX) THEN
                X(IL)=XSAV-H
              ELSE
                ETA(IY)=XSAV-H
              END IF
              FMINUS=LSQGFN(ETA,X,N,NR,IM)
              IF(LX) THEN
                X(IL)=XSAV
              ELSE
                ETA(IY)=XSAV
              END IF
              T(K,0)=(FPLUS-FMINUS)/(H+H)
              A(K)=T(K,0)
   10       CONTINUE
            IF(A(0) .GE. A(9)) THEN
              DO 20 K = 0,9
                A(K)=-A(K)
   20         CONTINUE
            END IF
            LMT=.TRUE.
            DO 30 K = 1,9
              H=A(K-1)-A(K)
              LMT=LMT .AND. (H.LE.EPSI .OR. ABS(H).LE.EPS*ABS(A(K))
     +        +EPSI)
   30       CONTINUE
            IF(LMT) GO TO 50
   40     CONTINUE
   50     CONTINUE
          DO 70 M = 1,9
            DO 60 K = 0,9-M
              IF(LEV(M)) THEN
                T(K,M)=W(M-1,1)*T(K+1,M-1)-W(M,1)*T(K,M-1)
              ELSE IF(LEV(K)) THEN
                T(K,M)=W(M-1,2)*T(K+1,M-1)-W(M,2)*T(K,M-1)
              ELSE
                T(K,M)=W(M-1,3)*T(K+1,M-1)-W(M,3)*T(K,M-1)
              END IF
   60       CONTINUE
   70     CONTINUE
          E(IM,I2)=T(0,9)
   80   CONTINUE
   90 CONTINUE
      END


      SUBROUTINE MTXCHM(U,A,R,M,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(M,M),A(M,N),R(M,N)
      PARAMETER(Z=0.D0)
      DO 30 I=1,M
        DO 20 K=1,N
          R(I,K)=Z
          DO 10 L=I,M
            R(I,K)=R(I,K)+U(I,L)*A(L,K)
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
      END


      SUBROUTINE MTXMSV(U,V,S,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION U(N),V(N)
      DO 10 I=1,N
        V(I)=S*U(I)
   10 CONTINUE
      END


      SUBROUTINE MTXLSC(A,B,E,D,X,R,A2,M,N,L,FRAC,OK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,N),B(M),E(L,N),D(L),X(N)
      LOGICAL OK
      PARAMETER(LMAX=1050,NMAX=1000)
      DIMENSION UP(LMAX),BB(LMAX),A2(M,*),P2(NMAX),S(NMAX),V(NMAX)
      COMMON /DASV02/ UP,BB,P2,S,V
C step 1
      NMINL=N-L
      DO 40 I=1,L
        CALL MTXGRW(E,V,L,N,I)
        CALL MTXHSD(V,UP(I),BB(I),N,I,I+1)
        DO 20 J=I,L
          CALL MTXGRW(E,S,L,N,J)
          CALL MTXHST(V,UP(I),BB(I),S,N,I,I+1)
          IF(J.EQ.I .AND. N.GT.I) THEN
            DO 10 K=I+1,N
              S(K)=V(K)
   10       CONTINUE
          END IF
          CALL MTXPRW(E,S,L,N,J)
   20   CONTINUE
        DO 30 J=1,M
          CALL MTXGRW(A,S,M,N,J)
          CALL MTXHST(V,UP(I),BB(I),S,N,I,I+1)
          CALL MTXPRW(A,S,M,N,J)
   30   CONTINUE
   40 CONTINUE
C step 2
      X(1)=D(1)/E(1,1)
      IF(L.GT.1) THEN
        DO 60 J=2,L
          X(J)=D(J)
          DO 50 K=1,J-1
            X(J)=X(J)-E(J,K)*X(K)
   50     CONTINUE
          X(J)=X(J)/E(J,J)
C          WRITE(*,*) J,X(J)
   60   CONTINUE
      END IF
C step 3
      DO 80 J=1,M
        DO 70 K=1,L
          B(J)=B(J)-A(J,K)*X(K)
   70   CONTINUE
   80 CONTINUE
C step 4
      L2=1
      CALL MTXGSM(A,A2,M,N,M,NMINL,1,L+1)
      CALL MTXSVD(A2,B,P2,R,M,NMINL,L2,FRAC,OK)
      IF(OK) THEN
        CALL MTXPSM(X,P2,N,1,NMINL,1,L+1,1)
        WRITE(*,*) 'After MTXPSM', X(1),X(2),X(3),X(4)
        DO 90 I=L,1,-1
          CALL MTXGRW(E,V,L,N,I)
          CALL MTXHST(V,UP(I),BB(I),X,N,I,I+1)
          WRITE(*,*) 'After MTXHST', X(1),X(2),X(3),X(4)
   90   CONTINUE
      END IF
      END


      SUBROUTINE MTXGSM(A,S,M,N,K,L,M1,N1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,N), S(K,L)
      DO 20 I=1,K
        DO 10 J=1,L
          S(I,J)=A(M1-1+I,N1-1+J)
   10   CONTINUE
   20 CONTINUE
      END


      SUBROUTINE MTXMBT(A,B,R,M,L,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,L),B(N,L),R(M,N)
      DO 30 J=1,N
        DO 20 I=1,M
          R(I,J)=0.
          DO 10 LL=1,L
            R(I,J)=R(I,J)+A(I,LL)*B(J,LL)
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
      END


      SUBROUTINE MTXMLT(A,B,R,M,L,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,L),B(L,N),R(M,N)
      DO 30 J=1,N
        DO 20 I=1,M
          R(I,J)=0.
          DO 10 LL=1,L
            R(I,J)=R(I,J)+A(I,LL)*B(LL,J)
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
      END


      SUBROUTINE MTXMAT(A,B,R,M,L,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(L,M),B(L,N),R(M,N)
      DO 30 J=1,N
        DO 20 I=1,M
          R(I,J)=0.
          DO 10 LL=1,L
            R(I,J)=R(I,J)+A(LL,I)*B(LL,J)
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
      END


      SUBROUTINE MTXSUB(A,B,R,M,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,N),B(M,N),R(M,N)
      DO 20 J=1,N
        DO 10 I=1,M
          R(I,J)=A(I,J)-B(I,J)
   10   CONTINUE
   20 CONTINUE
      END


      SUBROUTINE MTXUNT(R,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION R(N,N)
      DO 20 J=1,N
        DO 10 I=1,N
          IF(I.EQ.J) THEN
            R(I,J)=1.D0
          ELSE
            R(I,J)=0.D0
          END IF
   10   CONTINUE
   20 CONTINUE
      END


      SUBROUTINE MTXWRT(A,M,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,N)
      DO 20 I=1,M
        WRITE(2,10)(A(I,J),J=1,N)
   10 FORMAT(1X,4G20.13)
   20 CONTINUE
      END


      SUBROUTINE MTXGRW(A,R,M,N,I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,N), R(N)
      DO 10 J=1,N
        R(J)=A(I,J)
   10 CONTINUE
      END


      SUBROUTINE MTXHST(V,UP,B,C,N,LP,L)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(N),C(N)
      DUP=UP
      S=C(LP)*DUP
      WRITE (*,*) 'MTXHST V', S, V(1),V(2),V(3),V(4)
      WRITE (*,*) 'MTXHST C', C(1),C(2),C(3),C(4)
      DO 10 I=L,N
        S=S+C(I)*V(I)
        IF (V(2) .GT. 276.) WRITE(*,*) I,S,C(I),V(I)
   10 CONTINUE
      S=S*B
      C(LP)=C(LP)+S*DUP
      WRITE(*,*) 'MTXHST 1', C(1),C(2),C(3),C(4),S,B,DUP
      DO 20 I=L,N
        C(I)=C(I)+S*V(I)
   20 CONTINUE
      WRITE(*,*) 'MTXHST 2', C(1),C(2),C(3),C(4)
      END


      SUBROUTINE MTXPRW(A,R,M,N,I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,N), R(N)
      DO 10 J=1,N
        A(I,J)=R(J)
   10 CONTINUE
      END


      SUBROUTINE MTXSVD(A,B,X,R,M,N,NB,FRAC,OK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NMAX=1000)
      DIMENSION A(M,N),B(M,NB),X(N,NB),R(NB),D(NMAX),E(NMAX)
      COMMON /DASV01/ D,E
      LOGICAL OK
C STEP 1: Bidiagonalisation of A
      CALL MTXSV1(A,B,D,E,M,N,NB)
C STEP 2: Diagonalisation of bidiagonal matrix
      CALL MTXSV2(A,B,D,E,M,N,NB,OK)
C STEP 3: Order singular values and perform  permutations
      CALL MTXSV3(A,B,D,M,N,NB)
C STEP 4: Singular value analysis
      CALL MTXSV4(A,B,D,X,R,M,N,NB,FRAC)
      END


      SUBROUTINE MTXGCL(A,C,M,N,I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,N), C(M)
      DO 10 J=1,M
        C(J)=A(J,I)
   10 CONTINUE
      END


      SUBROUTINE MTXHSD(V,UP,B,N,LP,L)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(N)
      C=ABS(V(LP))
      DO 10 I=L,N
        C=MAX(ABS(V(I)),C)
   10 CONTINUE
      IF(C.LE.0.D0) GO TO 30
      C1=1./C
      SD=(V(LP)*C1)**2
      DO 20 I=L,N
        SD=SD+(V(I)*C1)**2
   20 CONTINUE
      VPPRIM=SD
      VPPRIM=C*SQRT(ABS(VPPRIM))
      IF(V(LP).GT.0.D0) VPPRIM=-VPPRIM
      UP=V(LP)-VPPRIM
      B=1./(VPPRIM*UP)
   30 CONTINUE
      END


      SUBROUTINE MTXPCL(A,C,M,N,I)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,N), C(M)
      DO 10 J=1,M
        A(J,I)=C(J)
   10 CONTINUE
      END


      SUBROUTINE MTXGVD(V1,V2,C,S)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      A1=ABS(V1)
      A2=ABS(V2)
      IF(A1.GT.A2) THEN
        W=V2/V1
        Q=SQRT(1.+W*W)
        C=1./Q
        IF(V1.LT.0.) C=-C
        S=C*W
      ELSE
        IF(V2.NE.0.) THEN
          W=V1/V2
          Q=SQRT(1.+W*W)
          S=1./Q
          IF(V2.LT.0.) S=-S
          C=S*W
        ELSE
          C=1.
          S=0.
        END IF
      END IF
      END


      SUBROUTINE MTXGVA(V1,V2,C,S)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      A1=ABS(V1)
      A2=ABS(V2)
      IF(A1.GT.A2) THEN
        W=V2/V1
        Q=SQRT(1.+W*W)
        C=1./Q
        IF(V1.LT.0.) C=-C
        S=C*W
        V1=A1*Q
        V2=0.
      ELSE
        IF(V2.NE.0.) THEN
          W=V1/V2
          Q=SQRT(1.+W*W)
          S=1./Q
          IF(V2.LT.0.) S=-S
          C=S*W
          V1=A2*Q
          V2=0.
        ELSE
          C=1.
          S=0.
        END IF
      END IF
      END


      SUBROUTINE MTXGVT(Z1,Z2,C,S)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      W=Z1*C+Z2*S
      Z2=-Z1*S+Z2*C
      Z1=W
      END


      SUBROUTINE MTXSV1(A,B,D,E,M,N,NB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NMAX=1000)
      DIMENSION A(M,N),B(M,NB),V(NMAX),S(NMAX),UPS(NMAX),D(N),E(N)
      DIMENSION BBS(NMAX)
      COMMON /DASV00/ V,S,UPS,BBS
      DO 50 I=1,N
C set up Householder Transformation Q(I)
        IF(I.LT.N .OR. M.GT.N) THEN
          CALL MTXGCL(A,V,M,N,I)
          CALL MTXHSD(V,UP,BB,M,I,I+1)
C apply Q(I) to A
          DO 10 J=I,N
            CALL MTXGCL(A,S,M,N,J)
            CALL MTXHST(V,UP,BB,S,M,I,I+1)
            CALL MTXPCL(A,S,M,N,J)
   10     CONTINUE
C apply Q(I) to B
          DO 20 K=1,NB
            CALL MTXGCL(B,S,M,NB,K)
            CALL MTXHST(V,UP,BB,S,M,I,I+1)
            CALL MTXPCL(B,S,M,NB,K)
   20     CONTINUE
        END IF
        IF(I.LT.N-1) THEN
C set up Householder Transformation H(I)
          CALL MTXGRW(A,V,M,N,I)
          CALL MTXHSD(V,UP,BB,N,I+1,I+2)
C save H(I)
          UPS(I)=UP
          BBS(I)=BB
C apply H(I) to A
          DO 40 J=I,M
            CALL MTXGRW(A,S,M,N,J)
            CALL MTXHST(V,UP,BB,S,N,I+1,I+2)
C save elements I+2,... in row J of matrix A
            IF (J.EQ.I) THEN
              DO 30 K=I+2,N
                S(K)=V(K)
   30         CONTINUE
            END IF
            CALL MTXPRW(A,S,M,N,J)
   40     CONTINUE
        END IF
   50 CONTINUE
C copy diagonal of transformed matrix A to D
C and upper parallel A to E
      IF(N.GT.1) THEN
        DO 60 I=2,N
          D(I)=A(I,I)
          E(I)=A(I-1,I)
   60   CONTINUE
      END IF
      D(1)=A(1,1)
      E(1)=0.
C construct product matrix H=H(2)*H(3)*...*H(N), H(N)=I
      DO 90 I=N,1,-1
        IF(I.LE.N-1) CALL MTXGRW(A,V,M,N,I)
        DO 70 K=1,N
          A(I,K)=0.
   70   CONTINUE
        A(I,I)=1.
        IF(I.LT.N-1) THEN
          DO 80 K=I,N
            CALL MTXGCL(A,S,M,N,K)
            CALL MTXHST(V,UPS(I),BBS(I),S,N,I+1,I+2)
            CALL MTXPCL(A,S,M,N,K)
   80     CONTINUE
        END IF
   90 CONTINUE
      END


      SUBROUTINE MTXSV2(A,B,D,E,M,N,NB,OK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,N),B(M,NB),D(N),E(N)
      LOGICAL OK,ELZERO
      PARAMETER(ZERO=0.D0)
      OK=.TRUE.
      NITERM=10*N
      NITER=0
      BMX=D(1)
      IF(N.GT.1) THEN
        DO 10 I =  2, N
          BMX=MAX(ABS(D(I))+ABS(E(I)),BMX)
   10   CONTINUE
      END IF
      DO 60 K=N,1,-1
   20   CONTINUE
        IF(K.NE.1) THEN
          IF((BMX+D(K))-BMX.EQ.ZERO) THEN
C Since D(K).EQ.0. perform Givens transform with result E(K)=0.
            CALL MTXS21(A,D,E,M,N,NB,K)
          END IF
C Find L (2. LE. L .LE. K) so that either E(L)=0. or D(L-1)=0.
C In the latter case transform E(L) to zero. In both cases the
C matrix splits and the bottom right minor begins with row L.
C If no such L is found set L=1
          DO 30 LL=K,1,-1
            L=LL
            IF(L.EQ.1) THEN
              ELZERO=.FALSE.
              GO TO 40
            ELSE IF((BMX-E(L))-BMX.EQ.ZERO) THEN
              ELZERO=.TRUE.
              GO TO 40
            ELSE IF((BMX+D(L-1))-BMX.EQ.ZERO) THEN
              ELZERO=.FALSE.
            END IF
   30     CONTINUE
   40     IF (L.GT. 1.AND. .NOT.ELZERO) THEN
            CALL MTXS22(B,D,E,M,N,NB,K,L)
          END IF
          IF(L.NE.K) THEN
C one more QR pass with order K
            CALL MTXS23(A,B,D,E,M,N,NB,K,L)
            NITER=NITER+1
            IF(NITER.LE.NITERM) GO TO 20
C set flag indicating non-convergence
            OK=.FALSE.
          END IF
        END IF
        IF(D(K).LT.0.) THEN
C for negative singular values perform change of sign
          D(K)=-D(K)
          DO 50 J=1,N
            A(J,K)=-A(J,K)
   50     CONTINUE
        END IF
C order is decreased by one in next pass
   60 CONTINUE
      END
C-----------------------------------------------------------------
      SUBROUTINE MTXS21(A,D,E,M,N,NB,K)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,N),D(N),E(N)
      DO 20 I = K- 1,  1, -1
        IF(I.EQ.K-1) THEN
          CALL MTXGVA(D(I),E(I+1),CS,SN)
        ELSE
          CALL MTXGVA(D(I),H,CS,SN)
        END IF
        IF(I.GT.1) THEN
          H=0.
          CALL MTXGVT(E(I),H,CS,SN)
        END IF
        DO 10 J =  1,N
          CALL MTXGVT(A(J,I),A(J,K),CS,SN)
   10   CONTINUE
   20 CONTINUE
      END
C-----------------------------------------------------------------
      SUBROUTINE MTXS22(B,D,E,M,N,NB,K,L)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION B(M,NB),D(N),E(N)
      DO 20 I=L,K
        IF(I.EQ.L) THEN
          CALL MTXGVA(D(I),E(I),CS,SN)
        ELSE
          CALL MTXGVA(D(I),H,CS,SN)
        END IF
        IF(I.LT.K) THEN
          H=0.
          CALL MTXGVT(E(I+1),H,CS,SN)
        END IF
        DO 10 J=1,NB
          CALL MTXGVT(CS,SN,B(I,J),B(L-1,J))
   10   CONTINUE
   20 CONTINUE
      END
C-----------------------------------------------------------------
      SUBROUTINE MTXS23(A,B,D,E,M,N,NB,K,L)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,N),B(M,NB),D(N),E(N)
C Determine shift parameter
      F=((D(K-1)-D(K))*(D(K-1)+D(K))+(E(K-1)-E(K))*(E(K-1)+E(K)))/
     +  (2.*E(K)*D(K-1))
      IF(ABS(F).GT.1.E10) THEN
        G=ABS(F)
      ELSE
        G=SQRT(1.+F*F)
      END IF
      IF(F.GE.0.) THEN
        T=F+G
      ELSE
        T=F-G
      END IF
      F=((D(L)-D(K))*(D(L)+D(K))+E(K)*(D(K-1)/T-E(K)))/D(L)
      DO 30 I = L , K-1
        IF(I.EQ.L) THEN
C Define R(L)
          CALL MTXGVD(F,E(I+1),CS,SN)
        ELSE
C Define R(I) , I.NE.L
          CALL MTXGVA(E(I),H,CS,SN)
        END IF
        CALL MTXGVT(D(I),E(I+1),CS,SN)
        H=0.
        CALL MTXGVT(H,D(I+1),CS,SN)
        DO 10 J =  1, N
          CALL MTXGVT(A(J,I),A(J,I+1),CS,SN)
   10   CONTINUE
C Define T(I)
        CALL MTXGVA(D(I),H,CS,SN)
        CALL MTXGVT(E(I+1),D(I+1),CS,SN)
        IF(I.LT.K-1) THEN
          H=0.
          CALL MTXGVT(H,E(I+2),CS,SN)
        END IF
        DO 20 J =  1, NB
          CALL MTXGVT(B(I,J),B(I+1,J),CS,SN)
   20   CONTINUE
   30 CONTINUE
      END


      SUBROUTINE MTXSV3(A,B,D,M,N,NB)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,N),B(M,NB),D(N)
C Order singular values
      IF(N.GT.1) THEN
   10   DO 20 I=2,N
          IF(D(I).GT.D(I-1)) GO TO 30
   20   CONTINUE
        RETURN
   30   CONTINUE
        DO 70 I =  2, N
          T=D(I-1)
          K=I-1
          DO 40 J = I , N
            IF(T.LT.D(J)) THEN
              T=D(J)
              K=J
            END IF
   40     CONTINUE
          IF(K.NE.I-1) THEN
C perform permutation on singular values
            D(K)=D(I-1)
            D(I-1)=T
C perform permutation on matrix A
            DO 50 J =  1, N
              T=A(J,K)
              A(J,K)=A(J,I-1)
              A(J,I-1)=T
   50       CONTINUE
C perform permutation on matrix B
            DO 60 J =  1, NB
              T=B(K,J)
              B(K,J)=B(I-1,J)
              B(I-1,J)=T
   60       CONTINUE
          END IF
   70   CONTINUE
        GO TO 10
      END IF
      END


      SUBROUTINE MTXSV4(A,B,D,X,R,M,N,NB,FRAC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(M,N),B(M,NB),D(N),X(N,NB),R(NB)
      PARAMETER (EPSILN=1.D-15)
      FRACT=ABS(FRAC)
      IF(FRACT.LT.EPSILN) FRACT=EPSILN
      SINMAX=0.
      DO 10 I =  1, N
        SINMAX=MAX(SINMAX,D(I))
   10 CONTINUE
      SINMIN=SINMAX*FRACT
      KK=N
      DO 20 I=1,N
        IF(D(I).LE.SINMIN) THEN
          KK=I-1
          GO TO 30
        END IF
   20 CONTINUE
   30 CONTINUE
      DO 60 I =  1, M
        IF(I.LE.KK) THEN
          S1=1./D(I)
          DO 40 J=1,NB
            B(I,J)=B(I,J)*S1
   40     CONTINUE
        ELSE
          DO 50 J=1,NB
            IF(I.EQ.KK+1) THEN
              R(J)=B(I,J)**2
            ELSE
              R(J)=R(J)+B(I,J)**2
            END IF
            IF(I.LE.N) B(I,J)=0.
   50     CONTINUE
        END IF
   60 CONTINUE
      DO 90 I=1,N
        DO 80 J=1,NB
          X(I,J)=0.
          DO 70 K=1,N
            X(I,J)=X(I,J)+A(I,K)*B(K,J)
   70     CONTINUE
   80   CONTINUE
   90 CONTINUE
      N=KK
      END
