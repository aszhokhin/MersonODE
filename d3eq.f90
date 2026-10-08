!=======================================================================
! D3EQ -- fully modernized Fortran 90 version
!
! The differential equations and the Merson integration algorithm are
! preserved from the original Fortran 77 program.  The old COMMON block,
! GOTO statements, EXTERNAL procedure arguments and implicit typing have
! been removed.
!=======================================================================

MODULE D3EQ_MODULE
  IMPLICIT NONE

  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  INTEGER, PARAMETER :: NVAR = 3
  REAL(DP), PARAMETER :: MERSON_R = 1.0D-13

CONTAINS

  !---------------------------------------------------------------------
  ! Differential equations of the original D3EQ system.
  ! The argument T is retained because the original system has the form
  ! F = F(T,X), although T is not explicitly used in these equations.
  !---------------------------------------------------------------------
  SUBROUTINE D3EQ(T, X, F)
    IMPLICIT NONE

    REAL(DP), INTENT(IN)  :: T
    REAL(DP), INTENT(IN)  :: X(:)
    REAL(DP), INTENT(OUT) :: F(:)

    F(1) = -0.2D0 * X(1) + (X(3) - 1.0D0) * X(2)
    F(2) = -1.08D0 * X(2) + (X(3) + 1.0D0) * X(1)
    F(3) = -0.1D0 * X(3) - 2.0D0 * X(1) * X(2) + 5.0D0
  END SUBROUTINE D3EQ


  !---------------------------------------------------------------------
  ! Integrates the D3EQ system and stores M points between Z and V.
  !---------------------------------------------------------------------
  SUBROUTINE SOLEQ(X, U, Z, V)
    IMPLICIT NONE

    REAL(DP), INTENT(INOUT) :: X(:)
    REAL(DP), INTENT(OUT)   :: U(:,:)
    REAL(DP), INTENT(IN)    :: Z, V

    INTEGER :: N, M, I, L
    INTEGER :: J
    REAL(DP) :: A, H, W, S, T, T1
    LOGICAL :: OK

    N = SIZE(X)
    M = SIZE(U, 1)

    IF (SIZE(U, 2) /= N) THEN
       PRINT *, 'SOLEQ: incompatible dimensions of X and U'
       STOP 1
    END IF

    IF (M < 2) THEN
       PRINT *, 'SOLEQ: at least two output points are required'
       STOP 1
    END IF

    A = 0.00001D0
    H = 0.01D0
    W = 0.00001D0
    J = 0

    T = 0.0D0
    CALL MERSON(T, Z, X, A, H, W, J, OK)

    IF (.NOT. OK) THEN
       PRINT *, 'SOLEQ: Merson integration failed'
          STOP 1
    END IF

    U(1,:) = X(:)

    S = (V - Z) / REAL(M - 1, DP)
    T = Z

    DO L = 2, M
       T1 = T + S
       CALL MERSON(T, T1, X, A, H, W, J, OK)

       IF (.NOT. OK) THEN
          PRINT *, 'SOLEQ: Merson integration failed'
          STOP 1
       END IF

       T = T1
       U(L,:) = X(:)
    END DO
  END SUBROUTINE SOLEQ


  !---------------------------------------------------------------------
  ! Merson adaptive-step integration method.
  !
  ! This is a structured version of the original MERSON routine.  The
  ! numerical operations and acceptance/rejection criteria are retained;
  ! only the control flow has been rewritten without GOTO statements.
  !---------------------------------------------------------------------
  SUBROUTINE MERSON(T, Q, Y, A, H, O, J, OK)
    IMPLICIT NONE

    REAL(DP), INTENT(INOUT) :: T
    REAL(DP), INTENT(IN)    :: Q, A, O
    REAL(DP), INTENT(INOUT) :: Y(:)
    REAL(DP), INTENT(INOUT) :: H
    INTEGER, INTENT(IN)     :: J
    LOGICAL, INTENT(OUT)    :: OK

    INTEGER :: N, K
    INTEGER :: N2, N3, N4, N31, N41
    INTEGER :: KN, KN2, KN3, KN4
    INTEGER :: IS
    REAL(DP) :: D(500)
    REAL(DP) :: R
    REAL(DP) :: Z, S, P, C, F
    LOGICAL :: reject_step, reduce_step

    N = SIZE(Y)

    IF (4 * N > SIZE(D)) THEN
       PRINT *, 'MERSON: state vector is too large'
       STOP 1
    END IF

    R = MERSON_R
    OK = .TRUE.

    N2  = 2 * N
    N3  = 3 * N
    N4  = 4 * N
    N31 = N3 + 1
    N41 = N4 + 1

    DO K = 1, N
       D(K + N4) = Y(K)
    END DO

    Z  = T
    S  = H
    IS = 0

    ! Adaptive integration loop (original label 2).
    DO
       P = S
       C = Q - Z

       IF (ABS(S) >= ABS(C)) THEN
          S = C
          IF (ABS(C / P) < R) THEN
             H = P
             T = Z
             Y(:) = D(N4+1:N4+N)
             RETURN
          END IF
          IS = 1
       END IF

       DO K = 1, N
          D(K) = D(K + N4)
       END DO

       ! Merson stage 1
       F = S / 3.0D0
       CALL D3EQ(Z, D(N41:N41+N-1), D(N31:N31+N-1))
       Z = Z + F

       DO K = 1, N
          KN  = K + N
          KN3 = K + N3
          KN4 = K + N4
          D(KN)  = F * D(KN3)
          D(KN4) = D(KN) + D(K)
       END DO

       ! Merson stage 2
       CALL D3EQ(Z, D(N41:N41+N-1), D(N31:N31+N-1))

       DO K = 1, N
          KN  = K + N
          KN3 = K + N3
          KN4 = K + N4
          D(KN)  = 0.5D0 * D(KN)
          D(KN4) = 0.5D0 * F * D(KN3) + D(KN) + D(K)
       END DO

       ! Merson stage 3
       CALL D3EQ(Z, D(N41:N41+N-1), D(N31:N31+N-1))

       Z = Z + 0.5D0 * F

       DO K = 1, N
          KN  = K + N
          KN2 = K + N2
          KN3 = K + N3
          KN4 = K + N4
          D(KN2) = 4.5D0 * F * D(KN3)
          D(KN4) = 0.25D0 * D(KN2) + 0.75D0 * D(KN) + D(K)
       END DO

       ! Merson stage 4
       CALL D3EQ(Z, D(N41:N41+N-1), D(N31:N31+N-1))

       Z = Z + 0.5D0 * S

       DO K = 1, N
          KN  = K + N
          KN2 = K + N2
          KN3 = K + N3
          KN4 = K + N4
          D(KN)  = 2.0D0 * F * D(KN3) + D(KN)
          D(KN4) = 3.0D0 * D(KN) - D(KN2) + D(K)
       END DO

       ! Merson stage 5
       CALL D3EQ(Z, D(N41:N41+N-1), D(N31:N31+N-1))

       reject_step = .FALSE.

       DO K = 1, N
          KN  = K + N
          KN2 = K + N2
          KN3 = K + N3
          KN4 = K + N4

          D(KN2) = -0.5D0 * F * D(KN3) - D(KN2) + 2.0D0 * D(KN)
          D(KN4) = D(K + N4) - D(KN2)
          D(KN)  = ABS(0.5D0 * A * D(KN4))
          D(KN2) = ABS(D(KN2))

          IF (ABS(D(KN4)) > R .AND. D(KN2) > D(KN)) THEN
             reject_step = .TRUE.
             EXIT
          END IF
       END DO

       ! Original label 13: reduce the step after an error estimate.
       IF (reject_step) THEN
          C = 0.5D0 * S

          IF (ABS(C) >= O) THEN
             DO K = 1, N
                KN4 = K + N4
                D(KN4) = D(K)
             END DO

             Z = Z - S
             S = C
             IS = 0
          ELSE IF (J == 0) THEN
             OK = .FALSE.
             H = P
             T = Z
             Y(:) = D(N4+1:N4+N)
             RETURN
          ELSE
             S = O
             IF (P < 0.0D0) S = -S

             IF (IS == 1) THEN
                H = P
                T = Z
                Y(:) = D(N4+1:N4+N)
                RETURN
             END IF
          END IF

          CYCLE
       END IF

       ! Original label 11: the requested interval has been reached.
       IF (IS == 1) THEN
          H = P
          T = Z
          Y(:) = D(N4+1:N4+N)
          RETURN
       END IF

       ! Original label 10 / 2: determine whether the step should be
       ! retried with the current or a larger step.
       reduce_step = .FALSE.

       DO K = 1, N
          KN  = K + N
          KN2 = K + N2
          IF (D(KN2) > D(KN) / 32.0D0) THEN
             reduce_step = .TRUE.
             EXIT
          END IF
       END DO

       IF (reduce_step) THEN
          CYCLE
       END IF

       S = 2.0D0 * S
    END DO
  END SUBROUTINE MERSON

END MODULE D3EQ_MODULE


!=======================================================================
! Main program
!=======================================================================
PROGRAM INT4G
  USE D3EQ_MODULE
  IMPLICIT NONE

  INTEGER, PARAMETER :: KD = 8000
  REAL(DP) :: U(KD, NVAR)
  REAL(DP) :: X(NVAR)
  REAL(DP) :: X0(NVAR)
  REAL(DP) :: T1, T
  REAL(DP) :: TI
  INTEGER :: I

  X0 = (/ 0.8D0, 0.9D0, 0.8D0 /)
  X  = X0

  T1 = 130.0D0
  T  = 160.0D0

  OPEN(UNIT=6, FILE='D3EQ.RES', STATUS='UNKNOWN')
  OPEN(UNIT=3, FILE='D3EQ.DAT', STATUS='UNKNOWN')

  CALL SOLEQ(X, U, T1, T)

  DO I = 1, KD
     TI = REAL(I - 1, DP) * (T - T1) / REAL(KD - 1, DP) + T1
     WRITE(3, '(5E15.5)') TI, U(I,1), U(I,2), U(I,3)
  END DO

  CLOSE(UNIT=3)
  CLOSE(UNIT=6)

END PROGRAM INT4G
