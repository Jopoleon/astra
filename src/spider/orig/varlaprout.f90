!
!
SUBROUTINE DGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
  !
  !  -- LAPACK routine (version 2.0) --
  !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
  !     Courant Institute, Argonne National Lab, and Rice University
  !     February 29, 1992
  !
  !     .. Scalar Arguments ..
  USE prec_rkind
  implicit none
  INTEGER            INFO, KL, KU, LDAB, M, N
  !     ..
  !     .. Array Arguments ..
  INTEGER            IPIV( * )
  REAL(RKIND) ::               AB( LDAB, * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  SGBTRF computes an LU factorization of a real m-by-n band matrix A
  !  using partial pivoting with row interchanges.
  !
  !  This is the blocked version of the algorithm, calling Level 3 BLAS.
  !
  !  Arguments
  !  =========
  !
  !  M       (input) INTEGER
  !          The number of rows of the matrix A.  M >= 0.
  !
  !  N       (input) INTEGER
  !          The number of columns of the matrix A.  N >= 0.
  !
  !  KL      (input) INTEGER
  !          The number of subdiagonals within the band of A.  KL >= 0.
  !
  !  KU      (input) INTEGER
  !          The number of superdiagonals within the band of A.  KU >= 0.
  !
  !  AB      (input/output) REAL array, dimension (LDAB,N)
  !          On entry, the matrix A in band storage, in rows KL+1 to
  !          2*KL+KU+1; rows 1 to KL of the array need not be set.
  !          The j-th column of A is stored in the j-th column of the
  !          array AB as follows:
  !          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
  !
  !          On exit, details of the factorization: U is stored as an
  !          upper triangular band matrix with KL+KU superdiagonals in
  !          rows 1 to KL+KU+1, and the multipliers used during the
  !          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
  !          See below for further details.
  !
  !  LDAB    (input) INTEGER
  !          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
  !
  !  IPIV    (output) INTEGER array, dimension (min(M,N))
  !          The pivot indices; for 1 <= i <= min(M,N), row i of the
  !          matrix was interchanged with row IPIV(i).
  !
  !  INFO    (output) INTEGER
  !          = 0: successful exit
  !          < 0: if INFO = -i, the i-th argument had an illegal value
  !          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
  !               has been completed, but the factor U is exactly
  !               singular, and division by zero will occur if it is used
  !               to solve a system of equations.
  !
  !  Further Details
  !  ===============
  !
  !  The band storage scheme is illustrated by the following example, when
  !  M = N = 6, KL = 2, KU = 1:
  !
  !  On entry:                       On exit:
  !
  !      *    *    *    +    +    +       *    *    *   u14  u25  u36
  !      *    *    +    +    +    +       *    *   u13  u24  u35  u46
  !      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
  !     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
  !     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
  !     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
  !
  !  Array elements marked * are not used by the routine; elements marked
  !  + need not be set on entry, but are required by the routine to store
  !  elements of U because of fill-in resulting from the row interchanges.
  !
  !  =====================================================================
  !
  !     .. Parameters ..
  REAL(RKIND) ::               ONE, ZERO
  PARAMETER          ( ONE = 1.0E+0_RKIND, ZERO = 0.0E+0_RKIND )
  INTEGER            NBMAX, LDWORK
  PARAMETER          ( NBMAX = 64, LDWORK = NBMAX+1 )
  !     ..
  !     .. Local Scalars ..
  INTEGER            I, I2, I3, II, IP, J, J2, J3, JB, JJ, JM, JP, &
       &                   JU, K2, KM, KV, NB, NW
  REAL(RKIND) ::               TEMP
  !     ..
  !     .. Local Arrays ..
  REAL(RKIND) ::               WORK13( LDWORK, NBMAX ), &
       &                   WORK31( LDWORK, NBMAX )
  !     ..
  !     .. External Functions ..
  INTEGER            ILAENV, ISAMAX
  !OS      EXTERNAL           ILAENV, ISAMAX
  !     ..
  !     .. External Subroutines ..
  !OS      EXTERNAL           SCOPY, SGBTF2, SGEMM, SGER, SLASWP, SSCAL, &
  !OS     &                   SSWAP, STRSM, XERBLA
  !     ..
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX, MIN
  !     ..
  !     .. Executable Statements ..
  !
  !     KV is the number of superdiagonals in the factor U, allowing for
  !     fill-in
  !
  KV = KU + KL
  !
  !     Test the input parameters.
  !
  INFO = 0
  IF( M.LT.0 ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( KL.LT.0 ) THEN
     INFO = -3
  ELSE IF( KU.LT.0 ) THEN
     INFO = -4
  ELSE IF( LDAB.LT.KL+KV+1 ) THEN
     INFO = -6
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'SGBTRF', -INFO )
     RETURN
  END IF
  !
  !     Quick return if possible
  !
  IF( M.EQ.0 .OR. N.EQ.0 ) &
       &   RETURN
  !
  !     Determine the block size for this environment
  !
  NB = ILAENV( 1, 'SGBTRF', ' ', M, N, KL, KU )
  !
  !     The block size must not exceed the limit set by the size of the
  !     local arrays WORK13 and WORK31.
  !
  NB = MIN( NB, NBMAX )
  !
  IF( NB.LE.1 .OR. NB.GT.KL ) THEN
     !
     !        Use unblocked code
     !
     CALL SGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
  ELSE
     !
     !        Use blocked code
     !
     !        Zero the superdiagonal elements of the work array WORK13
     !
     DO 20 J = 1, NB
        DO 10 I = 1, J - 1
           WORK13( I, J ) = ZERO
10      END DO
20   END DO
     !
     !        Zero the subdiagonal elements of the work array WORK31
     !
     DO 40 J = 1, NB
        DO 30 I = J + 1, NB
           WORK31( I, J ) = ZERO
30      END DO
40   END DO
     !
     !        Gaussian elimination with partial pivoting
     !
     !        Set fill-in elements in columns KU+2 to KV to zero
     !
     DO 60 J = KU + 2, MIN( KV, N )
        DO 50 I = KV - J + 2, KL
           AB( I, J ) = ZERO
50      END DO
60   END DO
     !
     !        JU is the index of the last column affected by the current
     !        stage of the factorization
     !
     JU = 1
     !
     DO 180 J = 1, MIN( M, N ), NB
        JB = MIN( NB, MIN( M, N )-J+1 )
        !
        !           The active part of the matrix is partitioned
        !
        !              A11   A12   A13
        !              A21   A22   A23
        !              A31   A32   A33
        !
        !           Here A11, A21 and A31 denote the current block of JB columns
        !           which is about to be factorized. The number of rows in the
        !           partitioning are JB, I2, I3 respectively, and the numbers
        !           of columns are JB, J2, J3. The superdiagonal elements of A13
        !           and the subdiagonal elements of A31 lie outside the band.
        !
        I2 = MIN( KL-JB, M-J-JB+1 )
        I3 = MIN( JB, M-J-KL+1 )
        !
        !           J2 and J3 are computed after JU has been updated.
        !
        !           Factorize the current block of JB columns
        !
        DO 80 JJ = J, J + JB - 1
           !
           !              Set fill-in elements in column JJ+KV to zero
           !
           IF( JJ+KV.LE.N ) THEN
              DO 70 I = 1, KL
                 AB( I, JJ+KV ) = ZERO
70            END DO
           END IF
           !
           !              Find pivot and test for singularity. KM is the number of
           !              subdiagonal elements in the current column.
           !
           KM = MIN( KL, M-JJ )
           JP = ISAMAX( KM+1, AB( KV+1, JJ ), 1 )
           IPIV( JJ ) = JP + JJ - J
           IF( AB( KV+JP, JJ ).NE.ZERO ) THEN
              JU = MAX( JU, MIN( JJ+KU+JP-1, N ) )
              IF( JP.NE.1 ) THEN
                 !
                 !                    Apply interchange to columns J to J+JB-1
                 !
                 IF( JP+JJ-1.LT.J+KL ) THEN
                    !
                    CALL SSWAP( JB, AB( KV+1+JJ-J, J ), LDAB-1, &
                         &                              AB( KV+JP+JJ-J, J ), LDAB-1 )
                 ELSE
                    !
                    !                       The interchange affects columns J to JJ-1 of A31
                    !                       which are stored in the work array WORK31
                    !
                    CALL SSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
                         &                              WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                    CALL SSWAP( J+JB-JJ, AB( KV+1, JJ ), LDAB-1, &
                         &                              AB( KV+JP, JJ ), LDAB-1 )
                 END IF
              END IF
              !
              !                 Compute multipliers
              !
              CALL SSCAL( KM, ONE / AB( KV+1, JJ ), AB( KV+2, JJ ), &
                   &                        1 )
              !
              !                 Update trailing submatrix within the band and within
              !                 the current block. JM is the index of the last column
              !                 which needs to be updated.
              !
              JM = MIN( JU, J+JB-1 )
              IF( JM.GT.JJ ) &
                   &               CALL SGER( KM, JM-JJ, -ONE, AB( KV+2, JJ ), 1, &
                   &                          AB( KV, JJ+1 ), LDAB-1, &
                   &                          AB( KV+1, JJ+1 ), LDAB-1 )
           ELSE
              !
              !                 If pivot is zero, set INFO to the index of the pivot
              !                 unless a zero pivot has already been found.
              !
              IF( INFO.EQ.0 ) &
                   &               INFO = JJ
           END IF
           !
           !              Copy current column of A31 into the work array WORK31
           !
           NW = MIN( JJ-J+1, I3 )
           IF( NW.GT.0 ) &
                &            CALL SCOPY( NW, AB( KV+KL+1-JJ+J, JJ ), 1, &
                &                        WORK31( 1, JJ-J+1 ), 1 )
80      END DO
        IF( J+JB.LE.N ) THEN
           !
           !              Apply the row interchanges to the other blocks.
           !
           J2 = MIN( JU-J+1, KV ) - JB
           J3 = MAX( 0, JU-J-KV+1 )
           !
           !              Use SLASWP to apply the row interchanges to A12, A22, and
           !              A32.
           !
           CALL SLASWP( J2, AB( KV+1-JB, J+JB ), LDAB-1, 1, JB, &
                &                      IPIV( J ), 1 )
           !
           !              Adjust the pivot indices.
           !
           DO 90 I = J, J + JB - 1
              IPIV( I ) = IPIV( I ) + J - 1
90         END DO
           !
           !              Apply the row interchanges to A13, A23, and A33
           !              columnwise.
           !
           K2 = J - 1 + JB + J2
           DO 110 I = 1, J3
              JJ = K2 + I
              DO 100 II = J + I - 1, J + JB - 1
                 IP = IPIV( II )
                 IF( IP.NE.II ) THEN
                    TEMP = AB( KV+1+II-JJ, JJ )
                    AB( KV+1+II-JJ, JJ ) = AB( KV+1+IP-JJ, JJ )
                    AB( KV+1+IP-JJ, JJ ) = TEMP
                 END IF
100           END DO
110        END DO
           !
           !              Update the relevant part of the trailing submatrix
           !
           IF( J2.GT.0 ) THEN
              !
              !                 Update A12
              !
              CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit', &
                   &                        JB, J2, ONE, AB( KV+1, J ), LDAB-1, &
                   &                        AB( KV+1-JB, J+JB ), LDAB-1 )
              !
              IF( I2.GT.0 ) THEN
                 !
                 !                    Update A22
                 !
                 CALL SGEMM( 'No transpose', 'No transpose', I2, J2, &
                      &                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1, &
                      &                           AB( KV+1-JB, J+JB ), LDAB-1, ONE, &
                      &                           AB( KV+1, J+JB ), LDAB-1 )
              END IF
              !
              IF( I3.GT.0 ) THEN
                 !
                 !                    Update A32
                 !
                 CALL SGEMM( 'No transpose', 'No transpose', I3, J2, &
                      &                           JB, -ONE, WORK31, LDWORK, &
                      &                           AB( KV+1-JB, J+JB ), LDAB-1, ONE, &
                      &                           AB( KV+KL+1-JB, J+JB ), LDAB-1 )
              END IF
           END IF
           !
           IF( J3.GT.0 ) THEN
              !
              !                 Copy the lower triangle of A13 into the work array
              !                 WORK13
              !
              DO 130 JJ = 1, J3
                 DO 120 II = JJ, JB
                    WORK13( II, JJ ) = AB( II-JJ+1, JJ+J+KV-1 )
120              END DO
130           END DO
              !
              !                 Update A13 in the work array
              !
              CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit', &
                   &                        JB, J3, ONE, AB( KV+1, J ), LDAB-1, &
                   &                        WORK13, LDWORK )
              !
              IF( I2.GT.0 ) THEN
                 !
                 !                    Update A23
                 !
                 CALL SGEMM( 'No transpose', 'No transpose', I2, J3, &
                      &                           JB, -ONE, AB( KV+1+JB, J ), LDAB-1, &
                      &                           WORK13, LDWORK, ONE, AB( 1+JB, J+KV ), &
                      &                           LDAB-1 )
              END IF
              !
              IF( I3.GT.0 ) THEN
                 !
                 !                    Update A33
                 !
                 CALL SGEMM( 'No transpose', 'No transpose', I3, J3, &
                      &                           JB, -ONE, WORK31, LDWORK, WORK13, &
                      &                           LDWORK, ONE, AB( 1+KL, J+KV ), LDAB-1 )
              END IF
              !
              !                 Copy the lower triangle of A13 back into place
              !
              DO 150 JJ = 1, J3
                 DO 140 II = JJ, JB
                    AB( II-JJ+1, JJ+J+KV-1 ) = WORK13( II, JJ )
140              END DO
150           END DO
           END IF
        ELSE
           !
           !              Adjust the pivot indices.
           !
           DO 160 I = J, J + JB - 1
              IPIV( I ) = IPIV( I ) + J - 1
160        END DO
        END IF
        !
        !           Partially undo the interchanges in the current block to
        !           restore the upper triangular form of A31 and copy the upper
        !           triangle of A31 back into place
        !
        DO 170 JJ = J + JB - 1, J, -1
           JP = IPIV( JJ ) - JJ + 1
           IF( JP.NE.1 ) THEN
              !
              !                 Apply interchange to columns J to JJ-1
              !
              IF( JP+JJ-1.LT.J+KL ) THEN
                 !
                 !                    The interchange does not affect A31
                 !
                 CALL SSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
                      &                           AB( KV+JP+JJ-J, J ), LDAB-1 )
              ELSE
                 !
                 !                    The interchange does affect A31
                 !
                 CALL SSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
                      &                           WORK31( JP+JJ-J-KL, 1 ), LDWORK )
              END IF
           END IF
           !
           !              Copy the current column of A31 back into place
           !
           NW = MIN( I3, JJ-J+1 )
           IF( NW.GT.0 ) &
                &            CALL SCOPY( NW, WORK31( 1, JJ-J+1 ), 1, &
                &                        AB( KV+KL+1-JJ+J, JJ ), 1 )
170     END DO
180  END DO
  END IF
  !
  RETURN
  !
  !     End of SGBTRF
  !
END SUBROUTINE DGBTRF
!
SUBROUTINE DGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, &
     &                   INFO )
  !
  !  -- LAPACK routine (version 2.0) --
  !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
  !     Courant Institute, Argonne National Lab, and Rice University
  !     March 31, 1993 
  !
  !     .. Scalar Arguments ..
  USE prec_rkind
  implicit none
  CHARACTER          TRANS
  INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
  !     ..
  !     .. Array Arguments ..
  INTEGER            IPIV( * )
  REAL(RKIND) ::               AB( LDAB, * ), B( LDB, * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  SGBTRS solves a system of linear equations
  !     A * X = B  or  A' * X = B
  !  with a general band matrix A using the LU factorization computed
  !  by SGBTRF.
  !
  !  Arguments
  !  =========
  !
  !  TRANS   (input) CHARACTER*1
  !          Specifies the form of the system of equations.
  !          = 'N':  A * X = B  (No transpose)
  !          = 'T':  A'* X = B  (Transpose)
  !          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
  !
  !  N       (input) INTEGER
  !          The order of the matrix A.  N >= 0.
  !
  !  KL      (input) INTEGER
  !          The number of subdiagonals within the band of A.  KL >= 0.
  !
  !  KU      (input) INTEGER
  !          The number of superdiagonals within the band of A.  KU >= 0.
  !
  !  NRHS    (input) INTEGER
  !          The number of right hand sides, i.e., the number of columns
  !          of the matrix B.  NRHS >= 0.
  !
  !  AB      (input) REAL array, dimension (LDAB,N)
  !          Details of the LU factorization of the band matrix A, as
  !          computed by SGBTRF.  U is stored as an upper triangular band
  !          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
  !          the multipliers used during the factorization are stored in
  !          rows KL+KU+2 to 2*KL+KU+1.
  !
  !  LDAB    (input) INTEGER
  !          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
  !
  !  IPIV    (input) INTEGER array, dimension (N)
  !          The pivot indices; for 1 <= i <= N, row i of the matrix was
  !          interchanged with row IPIV(i).
  !
  !  B       (input/output) REAL array, dimension (LDB,NRHS)
  !          On entry, the right hand side matrix B.
  !          On exit, the solution matrix X.
  !
  !  LDB     (input) INTEGER
  !          The leading dimension of the array B.  LDB >= max(1,N).
  !
  !  INFO    (output) INTEGER
  !          = 0:  successful exit
  !          < 0: if INFO = -i, the i-th argument had an illegal value
  !
  !  =====================================================================
  !
  !     .. Parameters ..
  REAL(RKIND) ::               ONE
  PARAMETER          ( ONE = 1.0E+0_RKIND )
  !     ..
  !     .. Local Scalars ..
  LOGICAL            LNOTI, NOTRAN
  INTEGER            I, J, KD, L, LM
  !     ..
  !     .. External Functions ..
  LOGICAL            LSAME
  !OS      EXTERNAL           LSAME
  !     ..
  !     .. External Subroutines ..
  !OS      EXTERNAL           SGEMV, SGER, SSWAP, STBSV, XERBLA
  !     ..
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX, MIN
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  INFO = 0
  NOTRAN = LSAME( TRANS, 'N' )
  IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. &
       &    LSAME( TRANS, 'C' ) ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( KL.LT.0 ) THEN
     INFO = -3
  ELSE IF( KU.LT.0 ) THEN
     INFO = -4
  ELSE IF( NRHS.LT.0 ) THEN
     INFO = -5
  ELSE IF( LDAB.LT.( 2*KL+KU+1 ) ) THEN
     INFO = -7
  ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
     INFO = -10
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'SGBTRS', -INFO )
     RETURN
  END IF
  !
  !     Quick return if possible
  !
  IF( N.EQ.0 .OR. NRHS.EQ.0 ) &
       &   RETURN
  !
  KD = KU + KL + 1
  LNOTI = KL.GT.0
  !
  IF( NOTRAN ) THEN
     !
     !        Solve  A*X = B.
     !
     !        Solve L*X = B, overwriting B with X.
     !
     !        L is represented as a product of permutations and unit lower
     !        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
     !        where each transformation L(i) is a rank-one modification of
     !        the identity matrix.
     !
     IF( LNOTI ) THEN
        DO 10 J = 1, N - 1
           LM = MIN( KL, N-J )
           L = IPIV( J )
           IF( L.NE.J ) &
                &            CALL SSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
           CALL SGER( LM, NRHS, -ONE, AB( KD+1, J ), 1, B( J, 1 ), &
                &                    LDB, B( J+1, 1 ), LDB )
10      END DO
     END IF
     !
     DO 20 I = 1, NRHS
        !
        !           Solve U*X = B, overwriting B with X.
        !
        CALL STBSV( 'Upper', 'No transpose', 'Non-unit', N, KL+KU, &
             &                  AB, LDAB, B( 1, I ), 1 )
20   END DO
     !
  ELSE
     !
     !        Solve A'*X = B.
     !
     DO 30 I = 1, NRHS
        !
        !           Solve U'*X = B, overwriting B with X.
        !
        CALL STBSV( 'Upper', 'Transpose', 'Non-unit', N, KL+KU, AB, &
             &                  LDAB, B( 1, I ), 1 )
30   END DO
     !
     !        Solve L'*X = B, overwriting B with X.
     !
     IF( LNOTI ) THEN
        DO 40 J = N - 1, 1, -1
           LM = MIN( KL, N-J )
           CALL SGEMV( 'Transpose', LM, NRHS, -ONE, B( J+1, 1 ), &
                &                     LDB, AB( KD+1, J ), 1, ONE, B( J, 1 ), LDB )
           L = IPIV( J )
           IF( L.NE.J ) &
                &            CALL SSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
40      END DO
     END IF
  END IF
  RETURN
  !
  !     End of SGBTRS
  !
END SUBROUTINE DGBTRS
!
SUBROUTINE SGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
  !
  !  -- LAPACK routine (version 2.0) --
  !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
  !     Courant Institute, Argonne National Lab, and Rice University
  !     February 29, 1992
  !
  !     .. Scalar Arguments ..
  USE prec_rkind
  implicit none
  INTEGER            INFO, KL, KU, LDAB, M, N
  !     ..
  !     .. Array Arguments ..
  INTEGER            IPIV( * )
  REAL(RKIND) ::               AB( LDAB, * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  SGBTF2 computes an LU factorization of a real m-by-n band matrix A
  !  using partial pivoting with row interchanges.
  !
  !  This is the unblocked version of the algorithm, calling Level 2 BLAS.
  !
  !  Arguments
  !  =========
  !
  !  M       (input) INTEGER
  !          The number of rows of the matrix A.  M >= 0.
  !
  !  N       (input) INTEGER
  !          The number of columns of the matrix A.  N >= 0.
  !
  !  KL      (input) INTEGER
  !          The number of subdiagonals within the band of A.  KL >= 0.
  !
  !  KU      (input) INTEGER
  !          The number of superdiagonals within the band of A.  KU >= 0.
  !
  !  AB      (input/output) REAL array, dimension (LDAB,N)
  !          On entry, the matrix A in band storage, in rows KL+1 to
  !          2*KL+KU+1; rows 1 to KL of the array need not be set.
  !          The j-th column of A is stored in the j-th column of the
  !          array AB as follows:
  !          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
  !
  !          On exit, details of the factorization: U is stored as an
  !          upper triangular band matrix with KL+KU superdiagonals in
  !          rows 1 to KL+KU+1, and the multipliers used during the
  !          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
  !          See below for further details.
  !
  !  LDAB    (input) INTEGER
  !          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
  !
  !  IPIV    (output) INTEGER array, dimension (min(M,N))
  !          The pivot indices; for 1 <= i <= min(M,N), row i of the
  !          matrix was interchanged with row IPIV(i).
  !
  !  INFO    (output) INTEGER
  !          = 0: successful exit
  !          < 0: if INFO = -i, the i-th argument had an illegal value
  !          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
  !               has been completed, but the factor U is exactly
  !               singular, and division by zero will occur if it is used
  !               to solve a system of equations.
  !
  !  Further Details
  !  ===============
  !
  !  The band storage scheme is illustrated by the following example, when
  !  M = N = 6, KL = 2, KU = 1:
  !
  !  On entry:                       On exit:
  !
  !      *    *    *    +    +    +       *    *    *   u14  u25  u36
  !      *    *    +    +    +    +       *    *   u13  u24  u35  u46
  !      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
  !     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
  !     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
  !     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
  !
  !  Array elements marked * are not used by the routine; elements marked
  !  + need not be set on entry, but are required by the routine to store
  !  elements of U, because of fill-in resulting from the row
  !  interchanges.
  !
  !  =====================================================================
  !
  !     .. Parameters ..
  REAL(RKIND) ::               ONE, ZERO
  PARAMETER          ( ONE = 1.0E+0_RKIND, ZERO = 0.0E+0_RKIND )
  !     ..
  !     .. Local Scalars ..
  INTEGER            I, J, JP, JU, KM, KV
  !     ..
  !     .. External Functions ..
  INTEGER            ISAMAX
  !OS      EXTERNAL           ISAMAX
  !     ..
  !     .. External Subroutines ..
  !OS      EXTERNAL           SGER, SSCAL, SSWAP, XERBLA
  !     ..
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX, MIN
  !     ..
  !     .. Executable Statements ..
  !
  !     KV is the number of superdiagonals in the factor U, allowing for
  !     fill-in.
  !
  KV = KU + KL
  !
  !     Test the input parameters.
  !
  INFO = 0
  IF( M.LT.0 ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( KL.LT.0 ) THEN
     INFO = -3
  ELSE IF( KU.LT.0 ) THEN
     INFO = -4
  ELSE IF( LDAB.LT.KL+KV+1 ) THEN
     INFO = -6
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'SGBTF2', -INFO )
     RETURN
  END IF
  !
  !     Quick return if possible
  !
  IF( M.EQ.0 .OR. N.EQ.0 ) &
       &   RETURN
  !
  !     Gaussian elimination with partial pivoting
  !
  !     Set fill-in elements in columns KU+2 to KV to zero.
  !
  DO 20 J = KU + 2, MIN( KV, N )
     DO 10 I = KV - J + 2, KL
        AB( I, J ) = ZERO
10   END DO
20 END DO
  !
  !     JU is the index of the last column affected by the current stage
  !     of the factorization.
  !
  JU = 1
  !
  DO 40 J = 1, MIN( M, N )
     !
     !        Set fill-in elements in column J+KV to zero.
     !
     IF( J+KV.LE.N ) THEN
        DO 30 I = 1, KL
           AB( I, J+KV ) = ZERO
30      END DO
     END IF
     !
     !        Find pivot and test for singularity. KM is the number of
     !        subdiagonal elements in the current column.
     !
     KM = MIN( KL, M-J )
     JP = ISAMAX( KM+1, AB( KV+1, J ), 1 )
     IPIV( J ) = JP + J - 1
     IF( AB( KV+JP, J ).NE.ZERO ) THEN
        JU = MAX( JU, MIN( J+KU+JP-1, N ) )
        !
        !           Apply interchange to columns J to JU.
        !
        IF( JP.NE.1 ) &
             &         CALL SSWAP( JU-J+1, AB( KV+JP, J ), LDAB-1, &
             &                     AB( KV+1, J ), LDAB-1 )
        !
        IF( KM.GT.0 ) THEN
           !
           !              Compute multipliers.
           !
           CALL SSCAL( KM, ONE / AB( KV+1, J ), AB( KV+2, J ), 1 )
           !
           !              Update trailing submatrix within the band.
           !
           IF( JU.GT.J ) &
                &            CALL SGER( KM, JU-J, -ONE, AB( KV+2, J ), 1, &
                &                       AB( KV, J+1 ), LDAB-1, AB( KV+1, J+1 ), &
                &                       LDAB-1 )
        END IF
     ELSE
        !
        !           If pivot is zero, set INFO to the index of the pivot
        !           unless a zero pivot has already been found.
        !
        IF( INFO.EQ.0 ) &
             &         INFO = J
     END IF
40 END DO
  RETURN
  !
  !     End of SGBTF2
  !
END SUBROUTINE SGBTF2
!
SUBROUTINE SLASWP( N, A, LDA, K1, K2, IPIV, INCX )
  !
  !  -- LAPACK auxiliary routine (version 2.0) --
  !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
  !     Courant Institute, Argonne National Lab, and Rice University
  !     October 31, 1992
  !
  !     .. Scalar Arguments ..
  USE prec_rkind
  implicit none
  INTEGER            INCX, K1, K2, LDA, N
  !     ..
  !     .. Array Arguments ..
  INTEGER            IPIV( * )
  REAL(RKIND) ::               A( LDA, * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  SLASWP performs a series of row interchanges on the matrix A.
  !  One row interchange is initiated for each of rows K1 through K2 of A.
  !
  !  Arguments
  !  =========
  !
  !  N       (input) INTEGER
  !          The number of columns of the matrix A.
  !
  !  A       (input/output) REAL array, dimension (LDA,N)
  !          On entry, the matrix of column dimension N to which the row
  !          interchanges will be applied.
  !          On exit, the permuted matrix.
  !
  !  LDA     (input) INTEGER
  !          The leading dimension of the array A.
  !
  !  K1      (input) INTEGER
  !          The first element of IPIV for which a row interchange will
  !          be done.
  !
  !  K2      (input) INTEGER
  !          The last element of IPIV for which a row interchange will
  !          be done.
  !
  !  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
  !          The vector of pivot indices.  Only the elements in positions
  !          K1 through K2 of IPIV are accessed.
  !          IPIV(K) = L implies rows K and L are to be interchanged.
  !
  !  INCX    (input) INTEGER
  !          The increment between successive values of IPIV.  If IPIV
  !          is negative, the pivots are applied in reverse order.
  !
  ! =====================================================================
  !
  !     .. Local Scalars ..
  INTEGER            I, IP, IX
  !     ..
  !     .. External Subroutines ..
  !OS      EXTERNAL           SSWAP
  !     ..
  !     .. Executable Statements ..
  !
  !     Interchange row I with row IPIV(I) for each of rows K1 through K2.
  !
  IF( INCX.EQ.0 ) &
       &   RETURN
  IF( INCX.GT.0 ) THEN
     IX = K1
  ELSE
     IX = 1 + ( 1-K2 )*INCX
  END IF
  IF( INCX.EQ.1 ) THEN
     DO 10 I = K1, K2
        IP = IPIV( I )
        IF( IP.NE.I ) &
             &         CALL SSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
10   END DO
  ELSE IF( INCX.GT.1 ) THEN
     DO 20 I = K1, K2
        IP = IPIV( IX )
        IF( IP.NE.I ) &
             &         CALL SSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
        IX = IX + INCX
20   END DO
  ELSE IF( INCX.LT.0 ) THEN
     DO 30 I = K2, K1, -1
        IP = IPIV( IX )
        IF( IP.NE.I ) &
             &         CALL SSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
        IX = IX + INCX
30   END DO
  END IF
  !
  RETURN
  !
  !     End of SLASWP
  !
END SUBROUTINE SLASWP
SUBROUTINE XERBLA( SRNAME, INFO )
  !
  !  -- LAPACK auxiliary routine (version 2.0) --
  !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
  !     Courant Institute, Argonne National Lab, and Rice University
  !     September 30, 1994
  !
  !     .. Scalar Arguments ..
  USE prec_rkind
  implicit none
  CHARACTER*6        SRNAME
  INTEGER            INFO
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  XERBLA  is an error handler for the LAPACK routines.
  !  It is called by an LAPACK routine if an input parameter has an
  !  invalid value.  A message is printed and execution stops.
  !
  !  Installers may consider modifying the STOP statement in order to
  !  call system-specific exception-handling facilities.
  !
  !  Arguments
  !  =========
  !
  !  SRNAME  (input) CHARACTER*6
  !          The name of the routine which called XERBLA.
  !
  !  INFO    (input) INTEGER
  !          The position of the invalid parameter in the parameter list
  !          of the calling routine.
  !
  ! =====================================================================
  !
  !     .. Executable Statements ..
  !
  WRITE( *, FMT = 9999 )SRNAME, INFO
  !
  STOP
  !
9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ', &
       &      'an illegal value' )
  !
  !     End of XERBLA
  !
END SUBROUTINE XERBLA

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     all_sgbtrf_s.blas.f
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine scopy(n,sx,incx,sy,incy)
  !
  !     copies a vector, x, to a vector, y.
  !     uses unrolled loops for increments equal to 1.
  !     jack dongarra, linpack, 3/11/78.
  !     modified 12/3/93, array(1) declarations changed to array(*)
  !
  USE prec_rkind
  implicit none
  REAL(RKIND) :: sx(*),sy(*)
  integer i,incx,incy,ix,iy,m,mp1,n
  !
  if(n.le.0)return
  if(incx.eq.1.and.incy.eq.1)go to 20
  !
  !        code for unequal increments or equal increments
  !          not equal to 1
  !
  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do 10 i = 1,n
     sy(iy) = sx(ix)
     ix = ix + incx
     iy = iy + incy
10 end do
  return
  !
  !        code for both increments equal to 1
  !
  !
  !        clean-up loop
  !
20 m = mod(n,7)
  if( m .eq. 0 ) go to 40
  do 30 i = 1,m
     sy(i) = sx(i)
30 end do
  if( n .lt. 7 ) return
40 mp1 = m + 1
  do 50 i = mp1,n,7
     sy(i) = sx(i)
     sy(i + 1) = sx(i + 1)
     sy(i + 2) = sx(i + 2)
     sy(i + 3) = sx(i + 3)
     sy(i + 4) = sx(i + 4)
     sy(i + 5) = sx(i + 5)
     sy(i + 6) = sx(i + 6)
50 end do
  return
end subroutine scopy
SUBROUTINE SGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, &
     &                   BETA, C, LDC )
  !     .. Scalar Arguments ..
  USE prec_rkind
  implicit none
  CHARACTER*1        TRANSA, TRANSB
  INTEGER            M, N, K, LDA, LDB, LDC
  REAL(RKIND) ::               ALPHA, BETA
  !     .. Array Arguments ..
  REAL(RKIND) ::               A( LDA, * ), B( LDB, * ), C( LDC, * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  SGEMM  performs one of the matrix-matrix operations
  !
  !     C := alpha*op( A )*op( B ) + beta*C,
  !
  !  where  op( X ) is one of
  !
  !     op( X ) = X   or   op( X ) = X',
  !
  !  alpha and beta are scalars, and A, B and C are matrices, with op( A )
  !  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
  !
  !  Parameters
  !  ==========
  !
  !  TRANSA - CHARACTER*1.
  !           On entry, TRANSA specifies the form of op( A ) to be used in
  !           the matrix multiplication as follows:
  !
  !              TRANSA = 'N' or 'n',  op( A ) = A.
  !
  !              TRANSA = 'T' or 't',  op( A ) = A'.
  !
  !              TRANSA = 'C' or 'c',  op( A ) = A'.
  !
  !           Unchanged on exit.
  !
  !  TRANSB - CHARACTER*1.
  !           On entry, TRANSB specifies the form of op( B ) to be used in
  !           the matrix multiplication as follows:
  !
  !              TRANSB = 'N' or 'n',  op( B ) = B.
  !
  !              TRANSB = 'T' or 't',  op( B ) = B'.
  !
  !              TRANSB = 'C' or 'c',  op( B ) = B'.
  !
  !           Unchanged on exit.
  !
  !  M      - INTEGER.
  !           On entry,  M  specifies  the number  of rows  of the  matrix
  !           op( A )  and of the  matrix  C.  M  must  be at least  zero.
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry,  N  specifies the number  of columns of the matrix
  !           op( B ) and the number of columns of the matrix C. N must be
  !           at least zero.
  !           Unchanged on exit.
  !
  !  K      - INTEGER.
  !           On entry,  K  specifies  the number of columns of the matrix
  !           op( A ) and the number of rows of the matrix op( B ). K must
  !           be at least  zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - REAL            .
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
  !           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
  !           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
  !           part of the array  A  must contain the matrix  A,  otherwise
  !           the leading  k by m  part of the array  A  must contain  the
  !           matrix A.
  !           Unchanged on exit.
  !
  !  LDA    - INTEGER.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
  !           LDA must be at least  max( 1, m ), otherwise  LDA must be at
  !           least  max( 1, k ).
  !           Unchanged on exit.
  !
  !  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is
  !           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
  !           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
  !           part of the array  B  must contain the matrix  B,  otherwise
  !           the leading  n by k  part of the array  B  must contain  the
  !           matrix B.
  !           Unchanged on exit.
  !
  !  LDB    - INTEGER.
  !           On entry, LDB specifies the first dimension of B as declared
  !           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
  !           LDB must be at least  max( 1, k ), otherwise  LDB must be at
  !           least  max( 1, n ).
  !           Unchanged on exit.
  !
  !  BETA   - REAL            .
  !           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
  !           supplied as zero then C need not be set on input.
  !           Unchanged on exit.
  !
  !  C      - REAL             array of DIMENSION ( LDC, n ).
  !           Before entry, the leading  m by n  part of the array  C must
  !           contain the matrix  C,  except when  beta  is zero, in which
  !           case C need not be set on entry.
  !           On exit, the array  C  is overwritten by the  m by n  matrix
  !           ( alpha*op( A )*op( B ) + beta*C ).
  !
  !  LDC    - INTEGER.
  !           On entry, LDC specifies the first dimension of C as declared
  !           in  the  calling  (sub)  program.   LDC  must  be  at  least
  !           max( 1, m ).
  !           Unchanged on exit.
  !
  !
  !  Level 3 Blas routine.
  !
  !  -- Written on 8-February-1989.
  !     Jack Dongarra, Argonne National Laboratory.
  !     Iain Duff, AERE Harwell.
  !     Jeremy Du Croz, Numerical Algorithms Group Ltd.
  !     Sven Hammarling, Numerical Algorithms Group Ltd.
  !
  !
  !     .. External Functions ..
  LOGICAL            LSAME
  !OS      EXTERNAL           LSAME
  !     .. External Subroutines ..
  !OS      EXTERNAL           XERBLA
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX
  !     .. Local Scalars ..
  LOGICAL            NOTA, NOTB
  INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
  REAL(RKIND) ::               TEMP
  !     .. Parameters ..
  REAL(RKIND) ::               ONE         , ZERO
  PARAMETER        ( ONE = 1.0E+0_RKIND, ZERO = 0.0E+0_RKIND )
  !     ..
  !     .. Executable Statements ..
  !
  !     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
  !     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
  !     and  columns of  A  and the  number of  rows  of  B  respectively.
  !
  NOTA  = LSAME( TRANSA, 'N' )
  NOTB  = LSAME( TRANSB, 'N' )
  IF( NOTA )THEN
     NROWA = M
     NCOLA = K
  ELSE
     NROWA = K
     NCOLA = M
  END IF
  IF( NOTB )THEN
     NROWB = K
  ELSE
     NROWB = N
  END IF
  !
  !     Test the input parameters.
  !
  INFO = 0
  IF(      ( .NOT.NOTA                 ).AND. &
       &         ( .NOT.LSAME( TRANSA, 'C' ) ).AND. &
       &         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
     INFO = 1
  ELSE IF( ( .NOT.NOTB                 ).AND. &
       &         ( .NOT.LSAME( TRANSB, 'C' ) ).AND. &
       &         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
     INFO = 2
  ELSE IF( M  .LT.0               )THEN
     INFO = 3
  ELSE IF( N  .LT.0               )THEN
     INFO = 4
  ELSE IF( K  .LT.0               )THEN
     INFO = 5
  ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
     INFO = 8
  ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
     INFO = 10
  ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
     INFO = 13
  END IF
  IF( INFO.NE.0 )THEN
     CALL XERBLA( 'SGEMM ', INFO )
     RETURN
  END IF
  !
  !     Quick return if possible.
  !
  IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR. &
       &    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) ) &
       &   RETURN
  !
  !     And if  alpha.eq.zero.
  !
  IF( ALPHA.EQ.ZERO )THEN
     IF( BETA.EQ.ZERO )THEN
        DO 20, J = 1, N
           DO 10, I = 1, M
              C( I, J ) = ZERO
10         END DO
20      END DO
     ELSE
        DO 40, J = 1, N
           DO 30, I = 1, M
              C( I, J ) = BETA*C( I, J )
30         END DO
40      END DO
     END IF
     RETURN
  END IF
  !
  !     Start the operations.
  !
  IF( NOTB )THEN
     IF( NOTA )THEN
        !
        !           Form  C := alpha*A*B + beta*C.
        !
        DO 90, J = 1, N
           IF( BETA.EQ.ZERO )THEN
              DO 50, I = 1, M
                 C( I, J ) = ZERO
50            END DO
           ELSE IF( BETA.NE.ONE )THEN
              DO 60, I = 1, M
                 C( I, J ) = BETA*C( I, J )
60            END DO
           END IF
           DO 80, L = 1, K
              IF( B( L, J ).NE.ZERO )THEN
                 TEMP = ALPHA*B( L, J )
                 DO 70, I = 1, M
                    C( I, J ) = C( I, J ) + TEMP*A( I, L )
70               END DO
              END IF
80         END DO
90      END DO
     ELSE
        !
        !           Form  C := alpha*A'*B + beta*C
        !
        DO 120, J = 1, N
           DO 110, I = 1, M
              TEMP = ZERO
              DO 100, L = 1, K
                 TEMP = TEMP + A( L, I )*B( L, J )
100           END DO
              IF( BETA.EQ.ZERO )THEN
                 C( I, J ) = ALPHA*TEMP
              ELSE
                 C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
              END IF
110        END DO
120     END DO
     END IF
  ELSE
     IF( NOTA )THEN
        !
        !           Form  C := alpha*A*B' + beta*C
        !
        DO 170, J = 1, N
           IF( BETA.EQ.ZERO )THEN
              DO 130, I = 1, M
                 C( I, J ) = ZERO
130           END DO
           ELSE IF( BETA.NE.ONE )THEN
              DO 140, I = 1, M
                 C( I, J ) = BETA*C( I, J )
140           END DO
           END IF
           DO 160, L = 1, K
              IF( B( J, L ).NE.ZERO )THEN
                 TEMP = ALPHA*B( J, L )
                 DO 150, I = 1, M
                    C( I, J ) = C( I, J ) + TEMP*A( I, L )
150              END DO
              END IF
160        END DO
170     END DO
     ELSE
        !
        !           Form  C := alpha*A'*B' + beta*C
        !
        DO 200, J = 1, N
           DO 190, I = 1, M
              TEMP = ZERO
              DO 180, L = 1, K
                 TEMP = TEMP + A( L, I )*B( J, L )
180           END DO
              IF( BETA.EQ.ZERO )THEN
                 C( I, J ) = ALPHA*TEMP
              ELSE
                 C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
              END IF
190        END DO
200     END DO
     END IF
  END IF
  !
  RETURN
  !
  !     End of SGEMM .
  !
END SUBROUTINE SGEMM
SUBROUTINE SGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, &
     &                   BETA, Y, INCY )
  !     .. Scalar Arguments ..
  USE prec_rkind
  implicit none
  REAL(RKIND) ::               ALPHA, BETA
  INTEGER            INCX, INCY, LDA, M, N
  CHARACTER*1        TRANS
  !     .. Array Arguments ..
  REAL(RKIND) ::               A( LDA, * ), X( * ), Y( * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  SGEMV  performs one of the matrix-vector operations
  !
  !     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
  !
  !  where alpha and beta are scalars, x and y are vectors and A is an
  !  m by n matrix.
  !
  !  Parameters
  !  ==========
  !
  !  TRANS  - CHARACTER*1.
  !           On entry, TRANS specifies the operation to be performed as
  !           follows:
  !
  !              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
  !
  !              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
  !
  !              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
  !
  !           Unchanged on exit.
  !
  !  M      - INTEGER.
  !           On entry, M specifies the number of rows of the matrix A.
  !           M must be at least zero.
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the number of columns of the matrix A.
  !           N must be at least zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - REAL            .
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  A      - REAL             array of DIMENSION ( LDA, n ).
  !           Before entry, the leading m by n part of the array A must
  !           contain the matrix of coefficients.
  !           Unchanged on exit.
  !
  !  LDA    - INTEGER.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program. LDA must be at least
  !           max( 1, m ).
  !           Unchanged on exit.
  !
  !  X      - REAL             array of DIMENSION at least
  !           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
  !           and at least
  !           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
  !           Before entry, the incremented array X must contain the
  !           vector x.
  !           Unchanged on exit.
  !
  !  INCX   - INTEGER.
  !           On entry, INCX specifies the increment for the elements of
  !           X. INCX must not be zero.
  !           Unchanged on exit.
  !
  !  BETA   - REAL            .
  !           On entry, BETA specifies the scalar beta. When BETA is
  !           supplied as zero then Y need not be set on input.
  !           Unchanged on exit.
  !
  !  Y      - REAL             array of DIMENSION at least
  !           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
  !           and at least
  !           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
  !           Before entry with BETA non-zero, the incremented array Y
  !           must contain the vector y. On exit, Y is overwritten by the
  !           updated vector y.
  !
  !  INCY   - INTEGER.
  !           On entry, INCY specifies the increment for the elements of
  !           Y. INCY must not be zero.
  !           Unchanged on exit.
  !
  !
  !  Level 2 Blas routine.
  !
  !  -- Written on 22-October-1986.
  !     Jack Dongarra, Argonne National Lab.
  !     Jeremy Du Croz, Nag Central Office.
  !     Sven Hammarling, Nag Central Office.
  !     Richard Hanson, Sandia National Labs.
  !
  !
  !     .. Parameters ..
  REAL(RKIND) ::               ONE         , ZERO
  PARAMETER        ( ONE = 1.0E+0_RKIND, ZERO = 0.0E+0_RKIND )
  !     .. Local Scalars ..
  REAL(RKIND) ::               TEMP
  INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
  !     .. External Functions ..
  LOGICAL            LSAME
  !OS      EXTERNAL           LSAME
  !     .. External Subroutines ..
  !OS      EXTERNAL           XERBLA
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  INFO = 0
  IF     ( .NOT.LSAME( TRANS, 'N' ).AND. &
       &         .NOT.LSAME( TRANS, 'T' ).AND. &
       &         .NOT.LSAME( TRANS, 'C' )      )THEN
     INFO = 1
  ELSE IF( M.LT.0 )THEN
     INFO = 2
  ELSE IF( N.LT.0 )THEN
     INFO = 3
  ELSE IF( LDA.LT.MAX( 1, M ) )THEN
     INFO = 6
  ELSE IF( INCX.EQ.0 )THEN
     INFO = 8
  ELSE IF( INCY.EQ.0 )THEN
     INFO = 11
  END IF
  IF( INFO.NE.0 )THEN
     CALL XERBLA( 'SGEMV ', INFO )
     RETURN
  END IF
  !
  !     Quick return if possible.
  !
  IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR. &
       &    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) ) &
       &   RETURN
  !
  !     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
  !     up the start points in  X  and  Y.
  !
  IF( LSAME( TRANS, 'N' ) )THEN
     LENX = N
     LENY = M
  ELSE
     LENX = M
     LENY = N
  END IF
  IF( INCX.GT.0 )THEN
     KX = 1
  ELSE
     KX = 1 - ( LENX - 1 )*INCX
  END IF
  IF( INCY.GT.0 )THEN
     KY = 1
  ELSE
     KY = 1 - ( LENY - 1 )*INCY
  END IF
  !
  !     Start the operations. In this version the elements of A are
  !     accessed sequentially with one pass through A.
  !
  !     First form  y := beta*y.
  !
  IF( BETA.NE.ONE )THEN
     IF( INCY.EQ.1 )THEN
        IF( BETA.EQ.ZERO )THEN
           DO 10, I = 1, LENY
              Y( I ) = ZERO
10         END DO
        ELSE
           DO 20, I = 1, LENY
              Y( I ) = BETA*Y( I )
20         END DO
        END IF
     ELSE
        IY = KY
        IF( BETA.EQ.ZERO )THEN
           DO 30, I = 1, LENY
              Y( IY ) = ZERO
              IY      = IY   + INCY
30         END DO
        ELSE
           DO 40, I = 1, LENY
              Y( IY ) = BETA*Y( IY )
              IY      = IY           + INCY
40         END DO
        END IF
     END IF
  END IF
  IF( ALPHA.EQ.ZERO ) &
       &   RETURN
  IF( LSAME( TRANS, 'N' ) )THEN
     !
     !        Form  y := alpha*A*x + y.
     !
     JX = KX
     IF( INCY.EQ.1 )THEN
        DO 60, J = 1, N
           IF( X( JX ).NE.ZERO )THEN
              TEMP = ALPHA*X( JX )
              DO 50, I = 1, M
                 Y( I ) = Y( I ) + TEMP*A( I, J )
50            END DO
           END IF
           JX = JX + INCX
60      END DO
     ELSE
        DO 80, J = 1, N
           IF( X( JX ).NE.ZERO )THEN
              TEMP = ALPHA*X( JX )
              IY   = KY
              DO 70, I = 1, M
                 Y( IY ) = Y( IY ) + TEMP*A( I, J )
                 IY      = IY      + INCY
70            END DO
           END IF
           JX = JX + INCX
80      END DO
     END IF
  ELSE
     !
     !        Form  y := alpha*A'*x + y.
     !
     JY = KY
     IF( INCX.EQ.1 )THEN
        DO 100, J = 1, N
           TEMP = ZERO
           DO 90, I = 1, M
              TEMP = TEMP + A( I, J )*X( I )
90         END DO
           Y( JY ) = Y( JY ) + ALPHA*TEMP
           JY      = JY      + INCY
100     END DO
     ELSE
        DO 120, J = 1, N
           TEMP = ZERO
           IX   = KX
           DO 110, I = 1, M
              TEMP = TEMP + A( I, J )*X( IX )
              IX   = IX   + INCX
110        END DO
           Y( JY ) = Y( JY ) + ALPHA*TEMP
           JY      = JY      + INCY
120     END DO
     END IF
  END IF
  !
  RETURN
  !
  !     End of SGEMV .
  !
END SUBROUTINE SGEMV
SUBROUTINE SGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
  !     .. Scalar Arguments ..
  USE prec_rkind
  implicit none
  REAL(RKIND) ::               ALPHA
  INTEGER            INCX, INCY, LDA, M, N
  !     .. Array Arguments ..
  REAL(RKIND) ::               A( LDA, * ), X( * ), Y( * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  SGER   performs the rank 1 operation
  !
  !     A := alpha*x*y' + A,
  !
  !  where alpha is a scalar, x is an m element vector, y is an n element
  !  vector and A is an m by n matrix.
  !
  !  Parameters
  !  ==========
  !
  !  M      - INTEGER.
  !           On entry, M specifies the number of rows of the matrix A.
  !           M must be at least zero.
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the number of columns of the matrix A.
  !           N must be at least zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - REAL            .
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  X      - REAL             array of dimension at least
  !           ( 1 + ( m - 1 )*abs( INCX ) ).
  !           Before entry, the incremented array X must contain the m
  !           element vector x.
  !           Unchanged on exit.
  !
  !  INCX   - INTEGER.
  !           On entry, INCX specifies the increment for the elements of
  !           X. INCX must not be zero.
  !           Unchanged on exit.
  !
  !  Y      - REAL             array of dimension at least
  !           ( 1 + ( n - 1 )*abs( INCY ) ).
  !           Before entry, the incremented array Y must contain the n
  !           element vector y.
  !           Unchanged on exit.
  !
  !  INCY   - INTEGER.
  !           On entry, INCY specifies the increment for the elements of
  !           Y. INCY must not be zero.
  !           Unchanged on exit.
  !
  !  A      - REAL             array of DIMENSION ( LDA, n ).
  !           Before entry, the leading m by n part of the array A must
  !           contain the matrix of coefficients. On exit, A is
  !           overwritten by the updated matrix.
  !
  !  LDA    - INTEGER.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program. LDA must be at least
  !           max( 1, m ).
  !           Unchanged on exit.
  !
  !
  !  Level 2 Blas routine.
  !
  !  -- Written on 22-October-1986.
  !     Jack Dongarra, Argonne National Lab.
  !     Jeremy Du Croz, Nag Central Office.
  !     Sven Hammarling, Nag Central Office.
  !     Richard Hanson, Sandia National Labs.
  !
  !
  !     .. Parameters ..
  REAL(RKIND) ::               ZERO
  PARAMETER        ( ZERO = 0.0E+0_RKIND )
  !     .. Local Scalars ..
  REAL(RKIND) ::               TEMP
  INTEGER            I, INFO, IX, J, JY, KX
  !     .. External Subroutines ..
  !OS      EXTERNAL           XERBLA
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  INFO = 0
  IF     ( M.LT.0 )THEN
     INFO = 1
  ELSE IF( N.LT.0 )THEN
     INFO = 2
  ELSE IF( INCX.EQ.0 )THEN
     INFO = 5
  ELSE IF( INCY.EQ.0 )THEN
     INFO = 7
  ELSE IF( LDA.LT.MAX( 1, M ) )THEN
     INFO = 9
  END IF
  IF( INFO.NE.0 )THEN
     CALL XERBLA( 'SGER  ', INFO )
     RETURN
  END IF
  !
  !     Quick return if possible.
  !
  IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) ) &
       &   RETURN
  !
  !     Start the operations. In this version the elements of A are
  !     accessed sequentially with one pass through A.
  !
  IF( INCY.GT.0 )THEN
     JY = 1
  ELSE
     JY = 1 - ( N - 1 )*INCY
  END IF
  IF( INCX.EQ.1 )THEN
     DO 20, J = 1, N
        IF( Y( JY ).NE.ZERO )THEN
           TEMP = ALPHA*Y( JY )
           DO 10, I = 1, M
              A( I, J ) = A( I, J ) + X( I )*TEMP
10         END DO
        END IF
        JY = JY + INCY
20   END DO
  ELSE
     IF( INCX.GT.0 )THEN
        KX = 1
     ELSE
        KX = 1 - ( M - 1 )*INCX
     END IF
     DO 40, J = 1, N
        IF( Y( JY ).NE.ZERO )THEN
           TEMP = ALPHA*Y( JY )
           IX   = KX
           DO 30, I = 1, M
              A( I, J ) = A( I, J ) + X( IX )*TEMP
              IX        = IX        + INCX
30         END DO
        END IF
        JY = JY + INCY
40   END DO
  END IF
  !
  RETURN
  !
  !     End of SGER  .
  !
END SUBROUTINE SGER
subroutine sscal(n,sa,sx,incx)
  !
  !     scales a vector by a constant.
  !     uses unrolled loops for increment equal to 1.
  !     jack dongarra, linpack, 3/11/78.
  !     modified 3/93 to return if incx .le. 0.
  !     modified 12/3/93, array(1) declarations changed to array(*)
  !
  USE prec_rkind
  implicit none
  REAL(RKIND) :: sa,sx(*)
  integer i,incx,m,mp1,n,nincx
  !
  if( n.le.0 .or. incx.le.0 )return
  if(incx.eq.1)go to 20
  !
  !        code for increment not equal to 1
  !
  nincx = n*incx
  do 10 i = 1,nincx,incx
     sx(i) = sa*sx(i)
10 end do
  return
  !
  !        code for increment equal to 1
  !
  !
  !        clean-up loop
  !
20 m = mod(n,5)
  if( m .eq. 0 ) go to 40
  do 30 i = 1,m
     sx(i) = sa*sx(i)
30 end do
  if( n .lt. 5 ) return
40 mp1 = m + 1
  do 50 i = mp1,n,5
     sx(i) = sa*sx(i)
     sx(i + 1) = sa*sx(i + 1)
     sx(i + 2) = sa*sx(i + 2)
     sx(i + 3) = sa*sx(i + 3)
     sx(i + 4) = sa*sx(i + 4)
50 end do
  return
end subroutine sscal
subroutine sswap (n,sx,incx,sy,incy)
  !
  !     interchanges two vectors.
  !     uses unrolled loops for increments equal to 1.
  !     jack dongarra, linpack, 3/11/78.
  !     modified 12/3/93, array(1) declarations changed to array(*)
  !
  USE prec_rkind
  implicit none
  REAL(RKIND) :: sx(*),sy(*),stemp
  integer i,incx,incy,ix,iy,m,mp1,n
  !
  if(n.le.0)return
  if(incx.eq.1.and.incy.eq.1)go to 20
  !
  !       code for unequal increments or equal increments not equal
  !         to 1
  !
  ix = 1
  iy = 1
  if(incx.lt.0)ix = (-n+1)*incx + 1
  if(incy.lt.0)iy = (-n+1)*incy + 1
  do 10 i = 1,n
     stemp = sx(ix)
     sx(ix) = sy(iy)
     sy(iy) = stemp
     ix = ix + incx
     iy = iy + incy
10 end do
  return
  !
  !       code for both increments equal to 1
  !
  !
  !       clean-up loop
  !
20 m = mod(n,3)
  if( m .eq. 0 ) go to 40
  do 30 i = 1,m
     stemp = sx(i)
     sx(i) = sy(i)
     sy(i) = stemp
30 end do
  if( n .lt. 3 ) return
40 mp1 = m + 1
  do 50 i = mp1,n,3
     stemp = sx(i)
     sx(i) = sy(i)
     sy(i) = stemp
     stemp = sx(i + 1)
     sx(i + 1) = sy(i + 1)
     sy(i + 1) = stemp
     stemp = sx(i + 2)
     sx(i + 2) = sy(i + 2)
     sy(i + 2) = stemp
50 end do
  return
end subroutine sswap
SUBROUTINE STBSV ( UPLO, TRANS, DIAG, N, K, A, LDA, X, INCX )
  !     .. Scalar Arguments ..
  USE prec_rkind
  implicit none
  INTEGER            INCX, K, LDA, N
  CHARACTER*1        DIAG, TRANS, UPLO
  !     .. Array Arguments ..
  REAL(RKIND) ::               A( LDA, * ), X( * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  STBSV  solves one of the systems of equations
  !
  !     A*x = b,   or   A'*x = b,
  !
  !  where b and x are n element vectors and A is an n by n unit, or
  !  non-unit, upper or lower triangular band matrix, with ( k + 1 )
  !  diagonals.
  !
  !  No test for singularity or near-singularity is included in this
  !  routine. Such tests must be performed before calling this routine.
  !
  !  Parameters
  !  ==========
  !
  !  UPLO   - CHARACTER*1.
  !           On entry, UPLO specifies whether the matrix is an upper or
  !           lower triangular matrix as follows:
  !
  !              UPLO = 'U' or 'u'   A is an upper triangular matrix.
  !
  !              UPLO = 'L' or 'l'   A is a lower triangular matrix.
  !
  !           Unchanged on exit.
  !
  !  TRANS  - CHARACTER*1.
  !           On entry, TRANS specifies the equations to be solved as
  !           follows:
  !
  !              TRANS = 'N' or 'n'   A*x = b.
  !
  !              TRANS = 'T' or 't'   A'*x = b.
  !
  !              TRANS = 'C' or 'c'   A'*x = b.
  !
  !           Unchanged on exit.
  !
  !  DIAG   - CHARACTER*1.
  !           On entry, DIAG specifies whether or not A is unit
  !           triangular as follows:
  !
  !              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
  !
  !              DIAG = 'N' or 'n'   A is not assumed to be unit
  !                                  triangular.
  !
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the order of the matrix A.
  !           N must be at least zero.
  !           Unchanged on exit.
  !
  !  K      - INTEGER.
  !           On entry with UPLO = 'U' or 'u', K specifies the number of
  !           super-diagonals of the matrix A.
  !           On entry with UPLO = 'L' or 'l', K specifies the number of
  !           sub-diagonals of the matrix A.
  !           K must satisfy  0 .le. K.
  !           Unchanged on exit.
  !
  !  A      - REAL             array of DIMENSION ( LDA, n ).
  !           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
  !           by n part of the array A must contain the upper triangular
  !           band part of the matrix of coefficients, supplied column by
  !           column, with the leading diagonal of the matrix in row
  !           ( k + 1 ) of the array, the first super-diagonal starting at
  !           position 2 in row k, and so on. The top left k by k triangle
  !           of the array A is not referenced.
  !           The following program segment will transfer an upper
  !           triangular band matrix from conventional full matrix storage
  !           to band storage:
  !
  !                 DO 20, J = 1, N
  !                    M = K + 1 - J
  !                    DO 10, I = MAX( 1, J - K ), J
  !                       A( M + I, J ) = matrix( I, J )
  !              10    END DO
  !              20 END DO
  !
  !           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
  !           by n part of the array A must contain the lower triangular
  !           band part of the matrix of coefficients, supplied column by
  !           column, with the leading diagonal of the matrix in row 1 of
  !           the array, the first sub-diagonal starting at position 1 in
  !           row 2, and so on. The bottom right k by k triangle of the
  !           array A is not referenced.
  !           The following program segment will transfer a lower
  !           triangular band matrix from conventional full matrix storage
  !           to band storage:
  !
  !                 DO 20, J = 1, N
  !                    M = 1 - J
  !                    DO 10, I = J, MIN( N, J + K )
  !                       A( M + I, J ) = matrix( I, J )
  !              10    END DO
  !              20 END DO
  !
  !           Note that when DIAG = 'U' or 'u' the elements of the array A
  !           corresponding to the diagonal elements of the matrix are not
  !           referenced, but are assumed to be unity.
  !           Unchanged on exit.
  !
  !  LDA    - INTEGER.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program. LDA must be at least
  !           ( k + 1 ).
  !           Unchanged on exit.
  !
  !  X      - REAL             array of dimension at least
  !           ( 1 + ( n - 1 )*abs( INCX ) ).
  !           Before entry, the incremented array X must contain the n
  !           element right-hand side vector b. On exit, X is overwritten
  !           with the solution vector x.
  !
  !  INCX   - INTEGER.
  !           On entry, INCX specifies the increment for the elements of
  !           X. INCX must not be zero.
  !           Unchanged on exit.
  !
  !
  !  Level 2 Blas routine.
  !
  !  -- Written on 22-October-1986.
  !     Jack Dongarra, Argonne National Lab.
  !     Jeremy Du Croz, Nag Central Office.
  !     Sven Hammarling, Nag Central Office.
  !     Richard Hanson, Sandia National Labs.
  !
  !
  !     .. Parameters ..
  REAL(RKIND) ::               ZERO
  PARAMETER        ( ZERO = 0.0E+0_RKIND )
  !     .. Local Scalars ..
  REAL(RKIND) ::               TEMP
  INTEGER            I, INFO, IX, J, JX, KPLUS1, KX, L
  LOGICAL            NOUNIT
  !     .. External Functions ..
  LOGICAL            LSAME
  !OS      EXTERNAL           LSAME
  !     .. External Subroutines ..
  !OS      EXTERNAL           XERBLA
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX, MIN
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  INFO = 0
  IF     ( .NOT.LSAME( UPLO , 'U' ).AND. &
       &         .NOT.LSAME( UPLO , 'L' )      )THEN
     INFO = 1
  ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND. &
       &         .NOT.LSAME( TRANS, 'T' ).AND. &
       &         .NOT.LSAME( TRANS, 'C' )      )THEN
     INFO = 2
  ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND. &
       &         .NOT.LSAME( DIAG , 'N' )      )THEN
     INFO = 3
  ELSE IF( N.LT.0 )THEN
     INFO = 4
  ELSE IF( K.LT.0 )THEN
     INFO = 5
  ELSE IF( LDA.LT.( K + 1 ) )THEN
     INFO = 7
  ELSE IF( INCX.EQ.0 )THEN
     INFO = 9
  END IF
  IF( INFO.NE.0 )THEN
     CALL XERBLA( 'STBSV ', INFO )
     RETURN
  END IF
  !
  !     Quick return if possible.
  !
  IF( N.EQ.0 ) &
       &   RETURN
  !
  NOUNIT = LSAME( DIAG, 'N' )
  !
  !     Set up the start point in X if the increment is not unity. This
  !     will be  ( N - 1 )*INCX  too small for descending loops.
  !
  IF( INCX.LE.0 )THEN
     KX = 1 - ( N - 1 )*INCX
  ELSE IF( INCX.NE.1 )THEN
     KX = 1
  END IF
  !
  !     Start the operations. In this version the elements of A are
  !     accessed by sequentially with one pass through A.
  !
  IF( LSAME( TRANS, 'N' ) )THEN
     !
     !        Form  x := inv( A )*x.
     !
     IF( LSAME( UPLO, 'U' ) )THEN
        KPLUS1 = K + 1
        IF( INCX.EQ.1 )THEN
           DO 20, J = N, 1, -1
              IF( X( J ).NE.ZERO )THEN
                 L = KPLUS1 - J
                 IF( NOUNIT ) &
                      &                  X( J ) = X( J )/A( KPLUS1, J )
                 TEMP = X( J )
                 DO 10, I = J - 1, MAX( 1, J - K ), -1
                    X( I ) = X( I ) - TEMP*A( L + I, J )
10               END DO
              END IF
20         END DO
        ELSE
           KX = KX + ( N - 1 )*INCX
           JX = KX
           DO 40, J = N, 1, -1
              KX = KX - INCX
              IF( X( JX ).NE.ZERO )THEN
                 IX = KX
                 L  = KPLUS1 - J
                 IF( NOUNIT ) &
                      &                  X( JX ) = X( JX )/A( KPLUS1, J )
                 TEMP = X( JX )
                 DO 30, I = J - 1, MAX( 1, J - K ), -1
                    X( IX ) = X( IX ) - TEMP*A( L + I, J )
                    IX      = IX      - INCX
30               END DO
              END IF
              JX = JX - INCX
40         END DO
        END IF
     ELSE
        IF( INCX.EQ.1 )THEN
           DO 60, J = 1, N
              IF( X( J ).NE.ZERO )THEN
                 L = 1 - J
                 IF( NOUNIT ) &
                      &                  X( J ) = X( J )/A( 1, J )
                 TEMP = X( J )
                 DO 50, I = J + 1, MIN( N, J + K )
                    X( I ) = X( I ) - TEMP*A( L + I, J )
50               END DO
              END IF
60         END DO
        ELSE
           JX = KX
           DO 80, J = 1, N
              KX = KX + INCX
              IF( X( JX ).NE.ZERO )THEN
                 IX = KX
                 L  = 1  - J
                 IF( NOUNIT ) &
                      &                  X( JX ) = X( JX )/A( 1, J )
                 TEMP = X( JX )
                 DO 70, I = J + 1, MIN( N, J + K )
                    X( IX ) = X( IX ) - TEMP*A( L + I, J )
                    IX      = IX      + INCX
70               END DO
              END IF
              JX = JX + INCX
80         END DO
        END IF
     END IF
  ELSE
     !
     !        Form  x := inv( A')*x.
     !
     IF( LSAME( UPLO, 'U' ) )THEN
        KPLUS1 = K + 1
        IF( INCX.EQ.1 )THEN
           DO 100, J = 1, N
              TEMP = X( J )
              L    = KPLUS1 - J
              DO 90, I = MAX( 1, J - K ), J - 1
                 TEMP = TEMP - A( L + I, J )*X( I )
90            END DO
              IF( NOUNIT ) &
                   &               TEMP = TEMP/A( KPLUS1, J )
              X( J ) = TEMP
100        END DO
        ELSE
           JX = KX
           DO 120, J = 1, N
              TEMP = X( JX )
              IX   = KX
              L    = KPLUS1  - J
              DO 110, I = MAX( 1, J - K ), J - 1
                 TEMP = TEMP - A( L + I, J )*X( IX )
                 IX   = IX   + INCX
110           END DO
              IF( NOUNIT ) &
                   &               TEMP = TEMP/A( KPLUS1, J )
              X( JX ) = TEMP
              JX      = JX   + INCX
              IF( J.GT.K ) &
                   &               KX = KX + INCX
120        END DO
        END IF
     ELSE
        IF( INCX.EQ.1 )THEN
           DO 140, J = N, 1, -1
              TEMP = X( J )
              L    = 1      - J
              DO 130, I = MIN( N, J + K ), J + 1, -1
                 TEMP = TEMP - A( L + I, J )*X( I )
130           END DO
              IF( NOUNIT ) &
                   &               TEMP = TEMP/A( 1, J )
              X( J ) = TEMP
140        END DO
        ELSE
           KX = KX + ( N - 1 )*INCX
           JX = KX
           DO 160, J = N, 1, -1
              TEMP = X( JX )
              IX   = KX
              L    = 1       - J
              DO 150, I = MIN( N, J + K ), J + 1, -1
                 TEMP = TEMP - A( L + I, J )*X( IX )
                 IX   = IX   - INCX
150           END DO
              IF( NOUNIT ) &
                   &               TEMP = TEMP/A( 1, J )
              X( JX ) = TEMP
              JX      = JX   - INCX
              IF( ( N - J ).GE.K ) &
                   &               KX = KX - INCX
160        END DO
        END IF
     END IF
  END IF
  !
  RETURN
  !
  !     End of STBSV .
  !
END SUBROUTINE STBSV
SUBROUTINE STRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, &
     &                   B, LDB )
  !     .. Scalar Arguments ..
  USE prec_rkind
  implicit none
  CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
  INTEGER            M, N, LDA, LDB
  REAL(RKIND) ::               ALPHA
  !     .. Array Arguments ..
  REAL(RKIND) ::               A( LDA, * ), B( LDB, * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  STRSM  solves one of the matrix equations
  !
  !     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
  !
  !  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
  !  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
  !
  !     op( A ) = A   or   op( A ) = A'.
  !
  !  The matrix X is overwritten on B.
  !
  !  Parameters
  !  ==========
  !
  !  SIDE   - CHARACTER*1.
  !           On entry, SIDE specifies whether op( A ) appears on the left
  !           or right of X as follows:
  !
  !              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
  !
  !              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
  !
  !           Unchanged on exit.
  !
  !  UPLO   - CHARACTER*1.
  !           On entry, UPLO specifies whether the matrix A is an upper or
  !           lower triangular matrix as follows:
  !
  !              UPLO = 'U' or 'u'   A is an upper triangular matrix.
  !
  !              UPLO = 'L' or 'l'   A is a lower triangular matrix.
  !
  !           Unchanged on exit.
  !
  !  TRANSA - CHARACTER*1.
  !           On entry, TRANSA specifies the form of op( A ) to be used in
  !           the matrix multiplication as follows:
  !
  !              TRANSA = 'N' or 'n'   op( A ) = A.
  !
  !              TRANSA = 'T' or 't'   op( A ) = A'.
  !
  !              TRANSA = 'C' or 'c'   op( A ) = A'.
  !
  !           Unchanged on exit.
  !
  !  DIAG   - CHARACTER*1.
  !           On entry, DIAG specifies whether or not A is unit triangular
  !           as follows:
  !
  !              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
  !
  !              DIAG = 'N' or 'n'   A is not assumed to be unit
  !                                  triangular.
  !
  !           Unchanged on exit.
  !
  !  M      - INTEGER.
  !           On entry, M specifies the number of rows of B. M must be at
  !           least zero.
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the number of columns of B.  N must be
  !           at least zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - REAL            .
  !           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
  !           zero then  A is not referenced and  B need not be set before
  !           entry.
  !           Unchanged on exit.
  !
  !  A      - REAL             array of DIMENSION ( LDA, k ), where k is m
  !           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
  !           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
  !           upper triangular part of the array  A must contain the upper
  !           triangular matrix  and the strictly lower triangular part of
  !           A is not referenced.
  !           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
  !           lower triangular part of the array  A must contain the lower
  !           triangular matrix  and the strictly upper triangular part of
  !           A is not referenced.
  !           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
  !           A  are not referenced either,  but are assumed to be  unity.
  !           Unchanged on exit.
  !
  !  LDA    - INTEGER.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
  !           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
  !           then LDA must be at least max( 1, n ).
  !           Unchanged on exit.
  !
  !  B      - REAL             array of DIMENSION ( LDB, n ).
  !           Before entry,  the leading  m by n part of the array  B must
  !           contain  the  right-hand  side  matrix  B,  and  on exit  is
  !           overwritten by the solution matrix  X.
  !
  !  LDB    - INTEGER.
  !           On entry, LDB specifies the first dimension of B as declared
  !           in  the  calling  (sub)  program.   LDB  must  be  at  least
  !           max( 1, m ).
  !           Unchanged on exit.
  !
  !
  !  Level 3 Blas routine.
  !
  !
  !  -- Written on 8-February-1989.
  !     Jack Dongarra, Argonne National Laboratory.
  !     Iain Duff, AERE Harwell.
  !     Jeremy Du Croz, Numerical Algorithms Group Ltd.
  !     Sven Hammarling, Numerical Algorithms Group Ltd.
  !
  !
  !     .. External Functions ..
  LOGICAL            LSAME
  !OS      EXTERNAL           LSAME
  !     .. External Subroutines ..
  !OS      EXTERNAL           XERBLA
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX
  !     .. Local Scalars ..
  LOGICAL            LSIDE, NOUNIT, UPPER
  INTEGER            I, INFO, J, K, NROWA
  REAL(RKIND) ::               TEMP
  !     .. Parameters ..
  REAL(RKIND) ::               ONE         , ZERO
  PARAMETER        ( ONE = 1.0E+0_RKIND, ZERO = 0.0E+0_RKIND )
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  LSIDE  = LSAME( SIDE  , 'L' )
  IF( LSIDE )THEN
     NROWA = M
  ELSE
     NROWA = N
  END IF
  NOUNIT = LSAME( DIAG  , 'N' )
  UPPER  = LSAME( UPLO  , 'U' )
  !
  INFO   = 0
  IF(      ( .NOT.LSIDE                ).AND. &
       &         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
     INFO = 1
  ELSE IF( ( .NOT.UPPER                ).AND. &
       &         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
     INFO = 2
  ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND. &
       &         ( .NOT.LSAME( TRANSA, 'T' ) ).AND. &
       &         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
     INFO = 3
  ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND. &
       &         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
     INFO = 4
  ELSE IF( M  .LT.0               )THEN
     INFO = 5
  ELSE IF( N  .LT.0               )THEN
     INFO = 6
  ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
     INFO = 9
  ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
     INFO = 11
  END IF
  IF( INFO.NE.0 )THEN
     CALL XERBLA( 'STRSM ', INFO )
     RETURN
  END IF
  !
  !     Quick return if possible.
  !
  IF( N.EQ.0 ) &
       &   RETURN
  !
  !     And when  alpha.eq.zero.
  !
  IF( ALPHA.EQ.ZERO )THEN
     DO 20, J = 1, N
        DO 10, I = 1, M
           B( I, J ) = ZERO
10      END DO
20   END DO
     RETURN
  END IF
  !
  !     Start the operations.
  !
  IF( LSIDE )THEN
     IF( LSAME( TRANSA, 'N' ) )THEN
        !
        !           Form  B := alpha*inv( A )*B.
        !
        IF( UPPER )THEN
           DO 60, J = 1, N
              IF( ALPHA.NE.ONE )THEN
                 DO 30, I = 1, M
                    B( I, J ) = ALPHA*B( I, J )
30               END DO
              END IF
              DO 50, K = M, 1, -1
                 IF( B( K, J ).NE.ZERO )THEN
                    IF( NOUNIT ) &
                         &                     B( K, J ) = B( K, J )/A( K, K )
                    DO 40, I = 1, K - 1
                       B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
40                  END DO
                 END IF
50            END DO
60         END DO
        ELSE
           DO 100, J = 1, N
              IF( ALPHA.NE.ONE )THEN
                 DO 70, I = 1, M
                    B( I, J ) = ALPHA*B( I, J )
70               END DO
              END IF
              DO 90 K = 1, M
                 IF( B( K, J ).NE.ZERO )THEN
                    IF( NOUNIT ) &
                         &                     B( K, J ) = B( K, J )/A( K, K )
                    DO 80, I = K + 1, M
                       B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
80                  END DO
                 END IF
90            END DO
100        END DO
        END IF
     ELSE
        !
        !           Form  B := alpha*inv( A' )*B.
        !
        IF( UPPER )THEN
           DO 130, J = 1, N
              DO 120, I = 1, M
                 TEMP = ALPHA*B( I, J )
                 DO 110, K = 1, I - 1
                    TEMP = TEMP - A( K, I )*B( K, J )
110              END DO
                 IF( NOUNIT ) &
                      &                  TEMP = TEMP/A( I, I )
                 B( I, J ) = TEMP
120           END DO
130        END DO
        ELSE
           DO 160, J = 1, N
              DO 150, I = M, 1, -1
                 TEMP = ALPHA*B( I, J )
                 DO 140, K = I + 1, M
                    TEMP = TEMP - A( K, I )*B( K, J )
140              END DO
                 IF( NOUNIT ) &
                      &                  TEMP = TEMP/A( I, I )
                 B( I, J ) = TEMP
150           END DO
160        END DO
        END IF
     END IF
  ELSE
     IF( LSAME( TRANSA, 'N' ) )THEN
        !
        !           Form  B := alpha*B*inv( A ).
        !
        IF( UPPER )THEN
           DO 210, J = 1, N
              IF( ALPHA.NE.ONE )THEN
                 DO 170, I = 1, M
                    B( I, J ) = ALPHA*B( I, J )
170              END DO
              END IF
              DO 190, K = 1, J - 1
                 IF( A( K, J ).NE.ZERO )THEN
                    DO 180, I = 1, M
                       B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
180                 END DO
                 END IF
190           END DO
              IF( NOUNIT )THEN
                 TEMP = ONE/A( J, J )
                 DO 200, I = 1, M
                    B( I, J ) = TEMP*B( I, J )
200              END DO
              END IF
210        END DO
        ELSE
           DO 260, J = N, 1, -1
              IF( ALPHA.NE.ONE )THEN
                 DO 220, I = 1, M
                    B( I, J ) = ALPHA*B( I, J )
220              END DO
              END IF
              DO 240, K = J + 1, N
                 IF( A( K, J ).NE.ZERO )THEN
                    DO 230, I = 1, M
                       B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
230                 END DO
                 END IF
240           END DO
              IF( NOUNIT )THEN
                 TEMP = ONE/A( J, J )
                 DO 250, I = 1, M
                    B( I, J ) = TEMP*B( I, J )
250              END DO
              END IF
260        END DO
        END IF
     ELSE
        !
        !           Form  B := alpha*B*inv( A' ).
        !
        IF( UPPER )THEN
           DO 310, K = N, 1, -1
              IF( NOUNIT )THEN
                 TEMP = ONE/A( K, K )
                 DO 270, I = 1, M
                    B( I, K ) = TEMP*B( I, K )
270              END DO
              END IF
              DO 290, J = 1, K - 1
                 IF( A( J, K ).NE.ZERO )THEN
                    TEMP = A( J, K )
                    DO 280, I = 1, M
                       B( I, J ) = B( I, J ) - TEMP*B( I, K )
280                 END DO
                 END IF
290           END DO
              IF( ALPHA.NE.ONE )THEN
                 DO 300, I = 1, M
                    B( I, K ) = ALPHA*B( I, K )
300              END DO
              END IF
310        END DO
        ELSE
           DO 360, K = 1, N
              IF( NOUNIT )THEN
                 TEMP = ONE/A( K, K )
                 DO 320, I = 1, M
                    B( I, K ) = TEMP*B( I, K )
320              END DO
              END IF
              DO 340, J = K + 1, N
                 IF( A( J, K ).NE.ZERO )THEN
                    TEMP = A( J, K )
                    DO 330, I = 1, M
                       B( I, J ) = B( I, J ) - TEMP*B( I, K )
330                 END DO
                 END IF
340           END DO
              IF( ALPHA.NE.ONE )THEN
                 DO 350, I = 1, M
                    B( I, K ) = ALPHA*B( I, K )
350              END DO
              END IF
360        END DO
        END IF
     END IF
  END IF
  !
  RETURN
  !
  !     End of STRSM .
  !
END SUBROUTINE STRSM
SUBROUTINE DPBTRF( UPLO, N, KD, AB, LDAB, INFO )
  !
  !  -- LAPACK routine (version 2.0) --
  !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
  !     Courant Institute, Argonne National Lab, and Rice University
  !     March 31, 1993 
  !
  !     .. Scalar Arguments ..
  USE prec_rkind
  implicit none
  CHARACTER          UPLO
  INTEGER            INFO, KD, LDAB, N
  !     ..
  !     .. Array Arguments ..
  REAL(RKIND) ::               AB( LDAB, * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  SPBTRF computes the Cholesky factorization of a real symmetric
  !  positive definite band matrix A.
  !
  !  The factorization has the form
  !     A = U**T * U,  if UPLO = 'U', or
  !     A = L  * L**T,  if UPLO = 'L',
  !  where U is an upper triangular matrix and L is lower triangular.
  !
  !  Arguments
  !  =========
  !
  !  UPLO    (input) CHARACTER*1
  !          = 'U':  Upper triangle of A is stored;
  !          = 'L':  Lower triangle of A is stored.
  !
  !  N       (input) INTEGER
  !          The order of the matrix A.  N >= 0.
  !
  !  KD      (input) INTEGER
  !          The number of superdiagonals of the matrix A if UPLO = 'U',
  !          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
  !
  !  AB      (input/output) REAL array, dimension (LDAB,N)
  !          On entry, the upper or lower triangle of the symmetric band
  !          matrix A, stored in the first KD+1 rows of the array.  The
  !          j-th column of A is stored in the j-th column of the array AB
  !          as follows:
  !          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
  !          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
  !
  !          On exit, if INFO = 0, the triangular factor U or L from the
  !          Cholesky factorization A = U**T*U or A = L*L**T of the band
  !          matrix A, in the same storage format as A.
  !
  !  LDAB    (input) INTEGER
  !          The leading dimension of the array AB.  LDAB >= KD+1.
  !
  !  INFO    (output) INTEGER
  !          = 0:  successful exit
  !          < 0:  if INFO = -i, the i-th argument had an illegal value
  !          > 0:  if INFO = i, the leading minor of order i is not
  !                positive definite, and the factorization could not be
  !                completed.
  !
  !  Further Details
  !  ===============
  !
  !  The band storage scheme is illustrated by the following example, when
  !  N = 6, KD = 2, and UPLO = 'U':
  !
  !  On entry:                       On exit:
  !
  !      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
  !      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
  !     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
  !
  !  Similarly, if UPLO = 'L' the format of A is as follows:
  !
  !  On entry:                       On exit:
  !
  !     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66
  !     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *
  !     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *
  !
  !  Array elements marked * are not used by the routine.
  !
  !  Contributed by
  !  Peter Mayes and Giuseppe Radicati, IBM ECSEC, Rome, March 23, 1989
  !
  !  =====================================================================
  !
  !     .. Parameters ..
  REAL(RKIND) ::               ONE, ZERO
  PARAMETER          ( ONE = 1.0E+0_RKIND, ZERO = 0.0E+0_RKIND )
  INTEGER            NBMAX, LDWORK
  PARAMETER          ( NBMAX = 32, LDWORK = NBMAX+1 )
  !     ..
  !     .. Local Scalars ..
  INTEGER            I, I2, I3, IB, II, J, JJ, NB
  !     ..
  !     .. Local Arrays ..
  REAL(RKIND) ::               WORK( LDWORK, NBMAX )
  !     ..
  !     .. External Functions ..
  LOGICAL            LSAME
  INTEGER            ILAENV
  !OS      EXTERNAL           LSAME, ILAENV
  !     ..
  !     .. External Subroutines ..
  !OS      EXTERNAL           SGEMM, SPBTF2, SPOTF2, SSYRK, STRSM, XERBLA
  !     ..
  !     .. Intrinsic Functions ..
  INTRINSIC          MIN
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  INFO = 0
  IF( ( .NOT.LSAME( UPLO, 'U' ) ) .AND. &
       &    ( .NOT.LSAME( UPLO, 'L' ) ) ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( KD.LT.0 ) THEN
     INFO = -3
  ELSE IF( LDAB.LT.KD+1 ) THEN
     INFO = -5
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'SPBTRF', -INFO )
     RETURN
  END IF
  !
  !     Quick return if possible
  !
  IF( N.EQ.0 ) &
       &   RETURN
  !
  !     Determine the block size for this environment
  !
  NB = ILAENV( 1, 'SPBTRF', UPLO, N, KD, -1, -1 )
  !
  !     The block size must not exceed the semi-bandwidth KD, and must not
  !     exceed the limit set by the size of the local array WORK.
  !
  NB = MIN( NB, NBMAX )
  !
  IF( NB.LE.1 .OR. NB.GT.KD ) THEN
     !
     !        Use unblocked code
     !
     CALL SPBTF2( UPLO, N, KD, AB, LDAB, INFO )
  ELSE
     !
     !        Use blocked code
     !
     IF( LSAME( UPLO, 'U' ) ) THEN
        !
        !           Compute the Cholesky factorization of a symmetric band
        !           matrix, given the upper triangle of the matrix in band
        !           storage.
        !
        !           Zero the upper triangle of the work array.
        !
        DO 20 J = 1, NB
           DO 10 I = 1, J - 1
              WORK( I, J ) = ZERO
10         END DO
20      END DO
        !
        !           Process the band matrix one diagonal block at a time.
        !
        DO 70 I = 1, N, NB
           IB = MIN( NB, N-I+1 )
           !
           !              Factorize the diagonal block
           !
           CALL SPOTF2( UPLO, IB, AB( KD+1, I ), LDAB-1, II )
           IF( II.NE.0 ) THEN
              INFO = I + II - 1
              GO TO 150
           END IF
           IF( I+IB.LE.N ) THEN
              !
              !                 Update the relevant part of the trailing submatrix.
              !                 If A11 denotes the diagonal block which has just been
              !                 factorized, then we need to update the remaining
              !                 blocks in the diagram:
              !
              !                    A11   A12   A13
              !                          A22   A23
              !                                A33
              !
              !                 The numbers of rows and columns in the partitioning
              !                 are IB, I2, I3 respectively. The blocks A12, A22 and
              !                 A23 are empty if IB = KD. The upper triangle of A13
              !                 lies outside the band.
              !
              I2 = MIN( KD-IB, N-I-IB+1 )
              I3 = MIN( IB, N-I-KD+1 )
              !
              IF( I2.GT.0 ) THEN
                 !
                 !                    Update A12
                 !
                 CALL STRSM( 'Left', 'Upper', 'Transpose', &
                      &                           'Non-unit', IB, I2, ONE, AB( KD+1, I ), &
                      &                           LDAB-1, AB( KD+1-IB, I+IB ), LDAB-1 )
                 !
                 !                    Update A22
                 !
                 CALL SSYRK( 'Upper', 'Transpose', I2, IB, -ONE, &
                      &                           AB( KD+1-IB, I+IB ), LDAB-1, ONE, &
                      &                           AB( KD+1, I+IB ), LDAB-1 )
              END IF
              !
              IF( I3.GT.0 ) THEN
                 !
                 !                    Copy the lower triangle of A13 into the work array.
                 !
                 DO 40 JJ = 1, I3
                    DO 30 II = JJ, IB
                       WORK( II, JJ ) = AB( II-JJ+1, JJ+I+KD-1 )
30                  END DO
40               END DO
                 !
                 !                    Update A13 (in the work array).
                 !
                 CALL STRSM( 'Left', 'Upper', 'Transpose', &
                      &                           'Non-unit', IB, I3, ONE, AB( KD+1, I ), &
                      &                           LDAB-1, WORK, LDWORK )
                 !
                 !                    Update A23
                 !
                 IF( I2.GT.0 ) &
                      &                  CALL SGEMM( 'Transpose', 'No Transpose', I2, I3, &
                      &                              IB, -ONE, AB( KD+1-IB, I+IB ), &
                      &                              LDAB-1, WORK, LDWORK, ONE, &
                      &                              AB( 1+IB, I+KD ), LDAB-1 )
                 !
                 !                    Update A33
                 !
                 CALL SSYRK( 'Upper', 'Transpose', I3, IB, -ONE, &
                      &                           WORK, LDWORK, ONE, AB( KD+1, I+KD ), &
                      &                           LDAB-1 )
                 !
                 !                    Copy the lower triangle of A13 back into place.
                 !
                 DO 60 JJ = 1, I3
                    DO 50 II = JJ, IB
                       AB( II-JJ+1, JJ+I+KD-1 ) = WORK( II, JJ )
50                  END DO
60               END DO
              END IF
           END IF
70      END DO
     ELSE
        !
        !           Compute the Cholesky factorization of a symmetric band
        !           matrix, given the lower triangle of the matrix in band
        !           storage.
        !
        !           Zero the lower triangle of the work array.
        !
        DO 90 J = 1, NB
           DO 80 I = J + 1, NB
              WORK( I, J ) = ZERO
80         END DO
90      END DO
        !
        !           Process the band matrix one diagonal block at a time.
        !
        DO 140 I = 1, N, NB
           IB = MIN( NB, N-I+1 )
           !
           !              Factorize the diagonal block
           !
           CALL SPOTF2( UPLO, IB, AB( 1, I ), LDAB-1, II )
           IF( II.NE.0 ) THEN
              INFO = I + II - 1
              GO TO 150
           END IF
           IF( I+IB.LE.N ) THEN
              !
              !                 Update the relevant part of the trailing submatrix.
              !                 If A11 denotes the diagonal block which has just been
              !                 factorized, then we need to update the remaining
              !                 blocks in the diagram:
              !
              !                    A11
              !                    A21   A22
              !                    A31   A32   A33
              !
              !                 The numbers of rows and columns in the partitioning
              !                 are IB, I2, I3 respectively. The blocks A21, A22 and
              !                 A32 are empty if IB = KD. The lower triangle of A31
              !                 lies outside the band.
              !
              I2 = MIN( KD-IB, N-I-IB+1 )
              I3 = MIN( IB, N-I-KD+1 )
              !
              IF( I2.GT.0 ) THEN
                 !
                 !                    Update A21
                 !
                 CALL STRSM( 'Right', 'Lower', 'Transpose', &
                      &                           'Non-unit', I2, IB, ONE, AB( 1, I ), &
                      &                           LDAB-1, AB( 1+IB, I ), LDAB-1 )
                 !
                 !                    Update A22
                 !
                 CALL SSYRK( 'Lower', 'No Transpose', I2, IB, -ONE, &
                      &                           AB( 1+IB, I ), LDAB-1, ONE, &
                      &                           AB( 1, I+IB ), LDAB-1 )
              END IF
              !
              IF( I3.GT.0 ) THEN
                 !
                 !                    Copy the upper triangle of A31 into the work array.
                 !
                 DO 110 JJ = 1, IB
                    DO 100 II = 1, MIN( JJ, I3 )
                       WORK( II, JJ ) = AB( KD+1-JJ+II, JJ+I-1 )
100                 END DO
110              END DO
                 !
                 !                    Update A31 (in the work array).
                 !
                 CALL STRSM( 'Right', 'Lower', 'Transpose', &
                      &                           'Non-unit', I3, IB, ONE, AB( 1, I ), &
                      &                           LDAB-1, WORK, LDWORK )
                 !
                 !                    Update A32
                 !
                 IF( I2.GT.0 ) &
                      &                  CALL SGEMM( 'No transpose', 'Transpose', I3, I2, &
                      &                              IB, -ONE, WORK, LDWORK, &
                      &                              AB( 1+IB, I ), LDAB-1, ONE, &
                      &                              AB( 1+KD-IB, I+IB ), LDAB-1 )
                 !
                 !                    Update A33
                 !
                 CALL SSYRK( 'Lower', 'No Transpose', I3, IB, -ONE, &
                      &                           WORK, LDWORK, ONE, AB( 1, I+KD ), &
                      &                           LDAB-1 )
                 !
                 !                    Copy the upper triangle of A31 back into place.
                 !
                 DO 130 JJ = 1, IB
                    DO 120 II = 1, MIN( JJ, I3 )
                       AB( KD+1-JJ+II, JJ+I-1 ) = WORK( II, JJ )
120                 END DO
130              END DO
              END IF
           END IF
140     END DO
     END IF
  END IF
  RETURN
  !
150 CONTINUE
  RETURN
  !
  !     End of SPBTRF
  !
END SUBROUTINE DPBTRF
SUBROUTINE DPBTRS( UPLO, N, KD, NRHS, AB, LDAB, B, LDB, INFO )
  !
  !  -- LAPACK routine (version 2.0) --
  !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
  !     Courant Institute, Argonne National Lab, and Rice University
  !     September 30, 1994
  !
  !     .. Scalar Arguments ..
  USE prec_rkind
  implicit none
  CHARACTER          UPLO
  INTEGER            INFO, KD, LDAB, LDB, N, NRHS
  !     ..
  !     .. Array Arguments ..
  REAL(RKIND) ::               AB( LDAB, * ), B( LDB, * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  SPBTRS solves a system of linear equations A*X = B with a symmetric
  !  positive definite band matrix A using the Cholesky factorization
  !  A = U**T*U or A = L*L**T computed by SPBTRF.
  !
  !  Arguments
  !  =========
  !
  !  UPLO    (input) CHARACTER*1
  !          = 'U':  Upper triangular factor stored in AB;
  !          = 'L':  Lower triangular factor stored in AB.
  !
  !  N       (input) INTEGER
  !          The order of the matrix A.  N >= 0.
  !
  !  KD      (input) INTEGER
  !          The number of superdiagonals of the matrix A if UPLO = 'U',
  !          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
  !
  !  NRHS    (input) INTEGER
  !          The number of right hand sides, i.e., the number of columns
  !          of the matrix B.  NRHS >= 0.
  !
  !  AB      (input) REAL array, dimension (LDAB,N)
  !          The triangular factor U or L from the Cholesky factorization
  !          A = U**T*U or A = L*L**T of the band matrix A, stored in the
  !          first KD+1 rows of the array.  The j-th column of U or L is
  !          stored in the j-th column of the array AB as follows:
  !          if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for max(1,j-kd)<=i<=j;
  !          if UPLO ='L', AB(1+i-j,j)    = L(i,j) for j<=i<=min(n,j+kd).
  !
  !  LDAB    (input) INTEGER
  !          The leading dimension of the array AB.  LDAB >= KD+1.
  !
  !  B       (input/output) REAL array, dimension (LDB,NRHS)
  !          On entry, the right hand side matrix B.
  !          On exit, the solution matrix X.
  !
  !  LDB     (input) INTEGER
  !          The leading dimension of the array B.  LDB >= max(1,N).
  !
  !  INFO    (output) INTEGER
  !          = 0:  successful exit
  !          < 0:  if INFO = -i, the i-th argument had an illegal value
  !
  !  =====================================================================
  !
  !     .. Local Scalars ..
  LOGICAL            UPPER
  INTEGER            J
  !     ..
  !     .. External Functions ..
  LOGICAL            LSAME
  !OS      EXTERNAL           LSAME
  !     ..
  !     .. External Subroutines ..
  !OS      EXTERNAL           STBSV, XERBLA
  !     ..
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  INFO = 0
  UPPER = LSAME( UPLO, 'U' )
  IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( KD.LT.0 ) THEN
     INFO = -3
  ELSE IF( NRHS.LT.0 ) THEN
     INFO = -4
  ELSE IF( LDAB.LT.KD+1 ) THEN
     INFO = -6
  ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
     INFO = -8
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'SPBTRS', -INFO )
     RETURN
  END IF
  !
  !     Quick return if possible
  !
  IF( N.EQ.0 .OR. NRHS.EQ.0 ) &
       &   RETURN
  !
  IF( UPPER ) THEN
     !
     !        Solve A*X = B where A = U'*U.
     !
     DO 10 J = 1, NRHS
        !
        !           Solve U'*X = B, overwriting B with X.
        !
        CALL STBSV( 'Upper', 'Transpose', 'Non-unit', N, KD, AB, &
             &                  LDAB, B( 1, J ), 1 )
        !
        !           Solve U*X = B, overwriting B with X.
        !
        CALL STBSV( 'Upper', 'No transpose', 'Non-unit', N, KD, AB, &
             &                  LDAB, B( 1, J ), 1 )
10   END DO
  ELSE
     !
     !        Solve A*X = B where A = L*L'.
     !
     DO 20 J = 1, NRHS
        !
        !           Solve L*X = B, overwriting B with X.
        !
        CALL STBSV( 'Lower', 'No transpose', 'Non-unit', N, KD, AB, &
             &                  LDAB, B( 1, J ), 1 )
        !
        !           Solve L'*X = B, overwriting B with X.
        !
        CALL STBSV( 'Lower', 'Transpose', 'Non-unit', N, KD, AB, &
             &                  LDAB, B( 1, J ), 1 )
20   END DO
  END IF
  !
  RETURN
  !
  !     End of SPBTRS
  !
END SUBROUTINE DPBTRS
SUBROUTINE SPBTF2( UPLO, N, KD, AB, LDAB, INFO )
  !
  !  -- LAPACK routine (version 2.0) --
  !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
  !     Courant Institute, Argonne National Lab, and Rice University
  !     February 29, 1992
  !
  !     .. Scalar Arguments ..
  USE prec_rkind
  implicit none
  CHARACTER          UPLO
  INTEGER            INFO, KD, LDAB, N
  !     ..
  !     .. Array Arguments ..
  REAL(RKIND) ::               AB( LDAB, * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  SPBTF2 computes the Cholesky factorization of a real symmetric
  !  positive definite band matrix A.
  !
  !  The factorization has the form
  !     A = U' * U ,  if UPLO = 'U', or
  !     A = L  * L',  if UPLO = 'L',
  !  where U is an upper triangular matrix, U' is the transpose of U, and
  !  L is lower triangular.
  !
  !  This is the unblocked version of the algorithm, calling Level 2 BLAS.
  !
  !  Arguments
  !  =========
  !
  !  UPLO    (input) CHARACTER*1
  !          Specifies whether the upper or lower triangular part of the
  !          symmetric matrix A is stored:
  !          = 'U':  Upper triangular
  !          = 'L':  Lower triangular
  !
  !  N       (input) INTEGER
  !          The order of the matrix A.  N >= 0.
  !
  !  KD      (input) INTEGER
  !          The number of super-diagonals of the matrix A if UPLO = 'U',
  !          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.
  !
  !  AB      (input/output) REAL array, dimension (LDAB,N)
  !          On entry, the upper or lower triangle of the symmetric band
  !          matrix A, stored in the first KD+1 rows of the array.  The
  !          j-th column of A is stored in the j-th column of the array AB
  !          as follows:
  !          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
  !          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
  !
  !          On exit, if INFO = 0, the triangular factor U or L from the
  !          Cholesky factorization A = U'*U or A = L*L' of the band
  !          matrix A, in the same storage format as A.
  !
  !  LDAB    (input) INTEGER
  !          The leading dimension of the array AB.  LDAB >= KD+1.
  !
  !  INFO    (output) INTEGER
  !          = 0: successful exit
  !          < 0: if INFO = -k, the k-th argument had an illegal value
  !          > 0: if INFO = k, the leading minor of order k is not
  !               positive definite, and the factorization could not be
  !               completed.
  !
  !  Further Details
  !  ===============
  !
  !  The band storage scheme is illustrated by the following example, when
  !  N = 6, KD = 2, and UPLO = 'U':
  !
  !  On entry:                       On exit:
  !
  !      *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46
  !      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
  !     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
  !
  !  Similarly, if UPLO = 'L' the format of A is as follows:
  !
  !  On entry:                       On exit:
  !
  !     a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66
  !     a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *
  !     a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *
  !
  !  Array elements marked * are not used by the routine.
  !
  !  =====================================================================
  !
  !     .. Parameters ..
  REAL(RKIND) ::               ONE, ZERO
  PARAMETER          ( ONE = 1.0E+0_RKIND, ZERO = 0.0E+0_RKIND )
  !     ..
  !     .. Local Scalars ..
  LOGICAL            UPPER
  INTEGER            J, KLD, KN
  REAL(RKIND) ::               AJJ
  !     ..
  !     .. External Functions ..
  LOGICAL            LSAME
  !OS      EXTERNAL           LSAME
  !     ..
  !     .. External Subroutines ..
  !OS      EXTERNAL           SSCAL, SSYR, XERBLA
  !     ..
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX, MIN, SQRT
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  INFO = 0
  UPPER = LSAME( UPLO, 'U' )
  IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( KD.LT.0 ) THEN
     INFO = -3
  ELSE IF( LDAB.LT.KD+1 ) THEN
     INFO = -5
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'SPBTF2', -INFO )
     RETURN
  END IF
  !
  !     Quick return if possible
  !
  IF( N.EQ.0 ) &
       &   RETURN
  !
  KLD = MAX( 1, LDAB-1 )
  !
  IF( UPPER ) THEN
     !
     !        Compute the Cholesky factorization A = U'*U.
     !
     DO 10 J = 1, N
        !
        !           Compute U(J,J) and test for non-positive-definiteness.
        !
        AJJ = AB( KD+1, J )
        IF( AJJ.LE.ZERO ) &
             &         GO TO 30
        AJJ = SQRT( AJJ )
        AB( KD+1, J ) = AJJ
        !
        !           Compute elements J+1:J+KN of row J and update the
        !           trailing submatrix within the band.
        !
        KN = MIN( KD, N-J )
        IF( KN.GT.0 ) THEN
           CALL SSCAL( KN, ONE / AJJ, AB( KD, J+1 ), KLD )
           CALL SSYR( 'Upper', KN, -ONE, AB( KD, J+1 ), KLD, &
                &                    AB( KD+1, J+1 ), KLD )
        END IF
10   END DO
  ELSE
     !
     !        Compute the Cholesky factorization A = L*L'.
     !
     DO 20 J = 1, N
        !
        !           Compute L(J,J) and test for non-positive-definiteness.
        !
        AJJ = AB( 1, J )
        IF( AJJ.LE.ZERO ) &
             &         GO TO 30
        AJJ = SQRT( AJJ )
        AB( 1, J ) = AJJ
        !
        !           Compute elements J+1:J+KN of column J and update the
        !           trailing submatrix within the band.
        !
        KN = MIN( KD, N-J )
        IF( KN.GT.0 ) THEN
           CALL SSCAL( KN, ONE / AJJ, AB( 2, J ), 1 )
           CALL SSYR( 'Lower', KN, -ONE, AB( 2, J ), 1, &
                &                    AB( 1, J+1 ), KLD )
        END IF
20   END DO
  END IF
  RETURN
  !
30 CONTINUE
  INFO = J
  RETURN
  !
  !     End of SPBTF2
  !
END SUBROUTINE SPBTF2
SUBROUTINE SPOTF2( UPLO, N, A, LDA, INFO )
  !
  !  -- LAPACK routine (version 2.0) --
  !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
  !     Courant Institute, Argonne National Lab, and Rice University
  !     February 29, 1992
  !
  !     .. Scalar Arguments ..
  USE prec_rkind
  implicit none
  CHARACTER          UPLO
  INTEGER            INFO, LDA, N
  !     ..
  !     .. Array Arguments ..
  REAL(RKIND) ::               A( LDA, * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  SPOTF2 computes the Cholesky factorization of a real symmetric
  !  positive definite matrix A.
  !
  !  The factorization has the form
  !     A = U' * U ,  if UPLO = 'U', or
  !     A = L  * L',  if UPLO = 'L',
  !  where U is an upper triangular matrix and L is lower triangular.
  !
  !  This is the unblocked version of the algorithm, calling Level 2 BLAS.
  !
  !  Arguments
  !  =========
  !
  !  UPLO    (input) CHARACTER*1
  !          Specifies whether the upper or lower triangular part of the
  !          symmetric matrix A is stored.
  !          = 'U':  Upper triangular
  !          = 'L':  Lower triangular
  !
  !  N       (input) INTEGER
  !          The order of the matrix A.  N >= 0.
  !
  !  A       (input/output) REAL array, dimension (LDA,N)
  !          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
  !          n by n upper triangular part of A contains the upper
  !          triangular part of the matrix A, and the strictly lower
  !          triangular part of A is not referenced.  If UPLO = 'L', the
  !          leading n by n lower triangular part of A contains the lower
  !          triangular part of the matrix A, and the strictly upper
  !          triangular part of A is not referenced.
  !
  !          On exit, if INFO = 0, the factor U or L from the Cholesky
  !          factorization A = U'*U  or A = L*L'.
  !
  !  LDA     (input) INTEGER
  !          The leading dimension of the array A.  LDA >= max(1,N).
  !
  !  INFO    (output) INTEGER
  !          = 0: successful exit
  !          < 0: if INFO = -k, the k-th argument had an illegal value
  !          > 0: if INFO = k, the leading minor of order k is not
  !               positive definite, and the factorization could not be
  !               completed.
  !
  !  =====================================================================
  !
  !     .. Parameters ..
  REAL(RKIND) ::               ONE, ZERO
  PARAMETER          ( ONE = 1.0E+0_RKIND, ZERO = 0.0E+0_RKIND )
  !     ..
  !     .. Local Scalars ..
  LOGICAL            UPPER
  INTEGER            J
  REAL(RKIND) ::               AJJ
  !     ..
  !     .. External Functions ..
  LOGICAL            LSAME
  REAL(RKIND) ::               SDOT
  !OS      EXTERNAL           LSAME, SDOT
  !     ..
  !     .. External Subroutines ..
  !OS      EXTERNAL           SGEMV, SSCAL, XERBLA
  !     ..
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX, SQRT
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  INFO = 0
  UPPER = LSAME( UPLO, 'U' )
  IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
     INFO = -1
  ELSE IF( N.LT.0 ) THEN
     INFO = -2
  ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
     INFO = -4
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'SPOTF2', -INFO )
     RETURN
  END IF
  !
  !     Quick return if possible
  !
  IF( N.EQ.0 ) &
       &   RETURN
  !
  IF( UPPER ) THEN
     !
     !        Compute the Cholesky factorization A = U'*U.
     !
     DO 10 J = 1, N
        !
        !           Compute U(J,J) and test for non-positive-definiteness.
        !
        AJJ = A( J, J ) - SDOT( J-1, A( 1, J ), 1, A( 1, J ), 1 )
        IF( AJJ.LE.ZERO ) THEN
           A( J, J ) = AJJ
           GO TO 30
        END IF
        AJJ = SQRT( AJJ )
        A( J, J ) = AJJ
        !
        !           Compute elements J+1:N of row J.
        !
        IF( J.LT.N ) THEN
           CALL SGEMV( 'Transpose', J-1, N-J, -ONE, A( 1, J+1 ), &
                &                     LDA, A( 1, J ), 1, ONE, A( J, J+1 ), LDA )
           CALL SSCAL( N-J, ONE / AJJ, A( J, J+1 ), LDA )
        END IF
10   END DO
  ELSE
     !
     !        Compute the Cholesky factorization A = L*L'.
     !
     DO 20 J = 1, N
        !
        !           Compute L(J,J) and test for non-positive-definiteness.
        !
        AJJ = A( J, J ) - SDOT( J-1, A( J, 1 ), LDA, A( J, 1 ), &
             &            LDA )
        IF( AJJ.LE.ZERO ) THEN
           A( J, J ) = AJJ
           GO TO 30
        END IF
        AJJ = SQRT( AJJ )
        A( J, J ) = AJJ
        !
        !           Compute elements J+1:N of column J.
        !
        IF( J.LT.N ) THEN
           CALL SGEMV( 'No transpose', N-J, J-1, -ONE, A( J+1, 1 ), &
                &                     LDA, A( J, 1 ), LDA, ONE, A( J+1, J ), 1 )
           CALL SSCAL( N-J, ONE / AJJ, A( J+1, J ), 1 )
        END IF
20   END DO
  END IF
  GO TO 40
  !
30 CONTINUE
  INFO = J
  !
40 CONTINUE
  RETURN
  !
  !     End of SPOTF2
  !
END SUBROUTINE SPOTF2


SUBROUTINE SSYR  ( UPLO, N, ALPHA, X, INCX, A, LDA )
  !     .. Scalar Arguments ..
  USE prec_rkind
  implicit none
  REAL(RKIND) ::               ALPHA
  INTEGER            INCX, LDA, N
  CHARACTER*1        UPLO
  !     .. Array Arguments ..
  REAL(RKIND) ::               A( LDA, * ), X( * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  SSYR   performs the symmetric rank 1 operation
  !
  !     A := alpha*x*x' + A,
  !
  !  where alpha is a real scalar, x is an n element vector and A is an
  !  n by n symmetric matrix.
  !
  !  Parameters
  !  ==========
  !
  !  UPLO   - CHARACTER*1.
  !           On entry, UPLO specifies whether the upper or lower
  !           triangular part of the array A is to be referenced as
  !           follows:
  !
  !              UPLO = 'U' or 'u'   Only the upper triangular part of A
  !                                  is to be referenced.
  !
  !              UPLO = 'L' or 'l'   Only the lower triangular part of A
  !                                  is to be referenced.
  !
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the order of the matrix A.
  !           N must be at least zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - REAL            .
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  X      - REAL             array of dimension at least
  !           ( 1 + ( n - 1 )*abs( INCX ) ).
  !           Before entry, the incremented array X must contain the n
  !           element vector x.
  !           Unchanged on exit.
  !
  !  INCX   - INTEGER.
  !           On entry, INCX specifies the increment for the elements of
  !           X. INCX must not be zero.
  !           Unchanged on exit.
  !
  !  A      - REAL             array of DIMENSION ( LDA, n ).
  !           Before entry with  UPLO = 'U' or 'u', the leading n by n
  !           upper triangular part of the array A must contain the upper
  !           triangular part of the symmetric matrix and the strictly
  !           lower triangular part of A is not referenced. On exit, the
  !           upper triangular part of the array A is overwritten by the
  !           upper triangular part of the updated matrix.
  !           Before entry with UPLO = 'L' or 'l', the leading n by n
  !           lower triangular part of the array A must contain the lower
  !           triangular part of the symmetric matrix and the strictly
  !           upper triangular part of A is not referenced. On exit, the
  !           lower triangular part of the array A is overwritten by the
  !           lower triangular part of the updated matrix.
  !
  !  LDA    - INTEGER.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program. LDA must be at least
  !           max( 1, n ).
  !           Unchanged on exit.
  !
  !
  !  Level 2 Blas routine.
  !
  !  -- Written on 22-October-1986.
  !     Jack Dongarra, Argonne National Lab.
  !     Jeremy Du Croz, Nag Central Office.
  !     Sven Hammarling, Nag Central Office.
  !     Richard Hanson, Sandia National Labs.
  !
  !
  !     .. Parameters ..
  REAL(RKIND) ::               ZERO
  PARAMETER        ( ZERO = 0.0E+0_RKIND )
  !     .. Local Scalars ..
  REAL(RKIND) ::               TEMP
  INTEGER            I, INFO, IX, J, JX, KX
  !     .. External Functions ..
  LOGICAL            LSAME
  !OS      EXTERNAL           LSAME
  !     .. External Subroutines ..
  !OS      EXTERNAL           XERBLA
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  INFO = 0
  IF     ( .NOT.LSAME( UPLO, 'U' ).AND. &
       &         .NOT.LSAME( UPLO, 'L' )      )THEN
     INFO = 1
  ELSE IF( N.LT.0 )THEN
     INFO = 2
  ELSE IF( INCX.EQ.0 )THEN
     INFO = 5
  ELSE IF( LDA.LT.MAX( 1, N ) )THEN
     INFO = 7
  END IF
  IF( INFO.NE.0 )THEN
     CALL XERBLA( 'SSYR  ', INFO )
     RETURN
  END IF
  !
  !     Quick return if possible.
  !
  IF( ( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) ) &
       &   RETURN
  !
  !     Set the start point in X if the increment is not unity.
  !
  IF( INCX.LE.0 )THEN
     KX = 1 - ( N - 1 )*INCX
  ELSE IF( INCX.NE.1 )THEN
     KX = 1
  END IF
  !
  !     Start the operations. In this version the elements of A are
  !     accessed sequentially with one pass through the triangular part
  !     of A.
  !
  IF( LSAME( UPLO, 'U' ) )THEN
     !
     !        Form  A  when A is stored in upper triangle.
     !
     IF( INCX.EQ.1 )THEN
        DO 20, J = 1, N
           IF( X( J ).NE.ZERO )THEN
              TEMP = ALPHA*X( J )
              DO 10, I = 1, J
                 A( I, J ) = A( I, J ) + X( I )*TEMP
10            END DO
           END IF
20      END DO
     ELSE
        JX = KX
        DO 40, J = 1, N
           IF( X( JX ).NE.ZERO )THEN
              TEMP = ALPHA*X( JX )
              IX   = KX
              DO 30, I = 1, J
                 A( I, J ) = A( I, J ) + X( IX )*TEMP
                 IX        = IX        + INCX
30            END DO
           END IF
           JX = JX + INCX
40      END DO
     END IF
  ELSE
     !
     !        Form  A  when A is stored in lower triangle.
     !
     IF( INCX.EQ.1 )THEN
        DO 60, J = 1, N
           IF( X( J ).NE.ZERO )THEN
              TEMP = ALPHA*X( J )
              DO 50, I = J, N
                 A( I, J ) = A( I, J ) + X( I )*TEMP
50            END DO
           END IF
60      END DO
     ELSE
        JX = KX
        DO 80, J = 1, N
           IF( X( JX ).NE.ZERO )THEN
              TEMP = ALPHA*X( JX )
              IX   = JX
              DO 70, I = J, N
                 A( I, J ) = A( I, J ) + X( IX )*TEMP
                 IX        = IX        + INCX
70            END DO
           END IF
           JX = JX + INCX
80      END DO
     END IF
  END IF
  !
  RETURN
  !
  !     End of SSYR  .
  !
END SUBROUTINE SSYR

SUBROUTINE SSYRK ( UPLO, TRANS, N, K, ALPHA, A, LDA, &
     &                   BETA, C, LDC )
  !     .. Scalar Arguments ..
  USE prec_rkind
  implicit none
  CHARACTER*1        UPLO, TRANS
  INTEGER            N, K, LDA, LDC
  REAL(RKIND) ::               ALPHA, BETA
  !     .. Array Arguments ..
  REAL(RKIND) ::               A( LDA, * ), C( LDC, * )
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  SSYRK  performs one of the symmetric rank k operations
  !
  !     C := alpha*A*A' + beta*C,
  !
  !  or
  !
  !     C := alpha*A'*A + beta*C,
  !
  !  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
  !  and  A  is an  n by k  matrix in the first case and a  k by n  matrix
  !  in the second case.
  !
  !  Parameters
  !  ==========
  !
  !  UPLO   - CHARACTER*1.
  !           On  entry,   UPLO  specifies  whether  the  upper  or  lower
  !           triangular  part  of the  array  C  is to be  referenced  as
  !           follows:
  !
  !              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
  !                                  is to be referenced.
  !
  !              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
  !                                  is to be referenced.
  !
  !           Unchanged on exit.
  !
  !  TRANS  - CHARACTER*1.
  !           On entry,  TRANS  specifies the operation to be performed as
  !           follows:
  !
  !              TRANS = 'N' or 'n'   C := alpha*A*A' + beta*C.
  !
  !              TRANS = 'T' or 't'   C := alpha*A'*A + beta*C.
  !
  !              TRANS = 'C' or 'c'   C := alpha*A'*A + beta*C.
  !
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry,  N specifies the order of the matrix C.  N must be
  !           at least zero.
  !           Unchanged on exit.
  !
  !  K      - INTEGER.
  !           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
  !           of  columns   of  the   matrix   A,   and  on   entry   with
  !           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
  !           of rows of the matrix  A.  K must be at least zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - REAL            .
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
  !           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
  !           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
  !           part of the array  A  must contain the matrix  A,  otherwise
  !           the leading  k by n  part of the array  A  must contain  the
  !           matrix A.
  !           Unchanged on exit.
  !
  !  LDA    - INTEGER.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
  !           then  LDA must be at least  max( 1, n ), otherwise  LDA must
  !           be at least  max( 1, k ).
  !           Unchanged on exit.
  !
  !  BETA   - REAL            .
  !           On entry, BETA specifies the scalar beta.
  !           Unchanged on exit.
  !
  !  C      - REAL             array of DIMENSION ( LDC, n ).
  !           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
  !           upper triangular part of the array C must contain the upper
  !           triangular part  of the  symmetric matrix  and the strictly
  !           lower triangular part of C is not referenced.  On exit, the
  !           upper triangular part of the array  C is overwritten by the
  !           upper triangular part of the updated matrix.
  !           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
  !           lower triangular part of the array C must contain the lower
  !           triangular part  of the  symmetric matrix  and the strictly
  !           upper triangular part of C is not referenced.  On exit, the
  !           lower triangular part of the array  C is overwritten by the
  !           lower triangular part of the updated matrix.
  !
  !  LDC    - INTEGER.
  !           On entry, LDC specifies the first dimension of C as declared
  !           in  the  calling  (sub)  program.   LDC  must  be  at  least
  !           max( 1, n ).
  !           Unchanged on exit.
  !
  !
  !  Level 3 Blas routine.
  !
  !  -- Written on 8-February-1989.
  !     Jack Dongarra, Argonne National Laboratory.
  !     Iain Duff, AERE Harwell.
  !     Jeremy Du Croz, Numerical Algorithms Group Ltd.
  !     Sven Hammarling, Numerical Algorithms Group Ltd.
  !
  !
  !     .. External Functions ..
  LOGICAL            LSAME
  !OS      EXTERNAL           LSAME
  !     .. External Subroutines ..
  !OS      EXTERNAL           XERBLA
  !     .. Intrinsic Functions ..
  INTRINSIC          MAX
  !     .. Local Scalars ..
  LOGICAL            UPPER
  INTEGER            I, INFO, J, L, NROWA
  REAL(RKIND) ::               TEMP
  !     .. Parameters ..
  REAL(RKIND) ::               ONE ,         ZERO
  PARAMETER        ( ONE = 1.0E+0_RKIND, ZERO = 0.0E+0_RKIND )
  !     ..
  !     .. Executable Statements ..
  !
  !     Test the input parameters.
  !
  IF( LSAME( TRANS, 'N' ) )THEN
     NROWA = N
  ELSE
     NROWA = K
  END IF
  UPPER = LSAME( UPLO, 'U' )
  !
  INFO = 0
  IF(      ( .NOT.UPPER               ).AND. &
       &         ( .NOT.LSAME( UPLO , 'L' ) )      )THEN
     INFO = 1
  ELSE IF( ( .NOT.LSAME( TRANS, 'N' ) ).AND. &
       &         ( .NOT.LSAME( TRANS, 'T' ) ).AND. &
       &         ( .NOT.LSAME( TRANS, 'C' ) )      )THEN
     INFO = 2
  ELSE IF( N  .LT.0               )THEN
     INFO = 3
  ELSE IF( K  .LT.0               )THEN
     INFO = 4
  ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
     INFO = 7
  ELSE IF( LDC.LT.MAX( 1, N     ) )THEN
     INFO = 10
  END IF
  IF( INFO.NE.0 )THEN
     CALL XERBLA( 'SSYRK ', INFO )
     RETURN
  END IF
  !
  !     Quick return if possible.
  !
  IF( ( N.EQ.0 ).OR. &
       &    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) ) &
       &   RETURN
  !
  !     And when  alpha.eq.zero.
  !
  IF( ALPHA.EQ.ZERO )THEN
     IF( UPPER )THEN
        IF( BETA.EQ.ZERO )THEN
           DO 20, J = 1, N
              DO 10, I = 1, J
                 C( I, J ) = ZERO
10            END DO
20         END DO
        ELSE
           DO 40, J = 1, N
              DO 30, I = 1, J
                 C( I, J ) = BETA*C( I, J )
30            END DO
40         END DO
        END IF
     ELSE
        IF( BETA.EQ.ZERO )THEN
           DO 60, J = 1, N
              DO 50, I = J, N
                 C( I, J ) = ZERO
50            END DO
60         END DO
        ELSE
           DO 80, J = 1, N
              DO 70, I = J, N
                 C( I, J ) = BETA*C( I, J )
70            END DO
80         END DO
        END IF
     END IF
     RETURN
  END IF
  !
  !     Start the operations.
  !
  IF( LSAME( TRANS, 'N' ) )THEN
     !
     !        Form  C := alpha*A*A' + beta*C.
     !
     IF( UPPER )THEN
        DO 130, J = 1, N
           IF( BETA.EQ.ZERO )THEN
              DO 90, I = 1, J
                 C( I, J ) = ZERO
90            END DO
           ELSE IF( BETA.NE.ONE )THEN
              DO 100, I = 1, J
                 C( I, J ) = BETA*C( I, J )
100           END DO
           END IF
           DO 120, L = 1, K
              IF( A( J, L ).NE.ZERO )THEN
                 TEMP = ALPHA*A( J, L )
                 DO 110, I = 1, J
                    C( I, J ) = C( I, J ) + TEMP*A( I, L )
110              END DO
              END IF
120        END DO
130     END DO
     ELSE
        DO 180, J = 1, N
           IF( BETA.EQ.ZERO )THEN
              DO 140, I = J, N
                 C( I, J ) = ZERO
140           END DO
           ELSE IF( BETA.NE.ONE )THEN
              DO 150, I = J, N
                 C( I, J ) = BETA*C( I, J )
150           END DO
           END IF
           DO 170, L = 1, K
              IF( A( J, L ).NE.ZERO )THEN
                 TEMP      = ALPHA*A( J, L )
                 DO 160, I = J, N
                    C( I, J ) = C( I, J ) + TEMP*A( I, L )
160              END DO
              END IF
170        END DO
180     END DO
     END IF
  ELSE
     !
     !        Form  C := alpha*A'*A + beta*C.
     !
     IF( UPPER )THEN
        DO 210, J = 1, N
           DO 200, I = 1, J
              TEMP = ZERO
              DO 190, L = 1, K
                 TEMP = TEMP + A( L, I )*A( L, J )
190           END DO
              IF( BETA.EQ.ZERO )THEN
                 C( I, J ) = ALPHA*TEMP
              ELSE
                 C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
              END IF
200        END DO
210     END DO
     ELSE
        DO 240, J = 1, N
           DO 230, I = J, N
              TEMP = ZERO
              DO 220, L = 1, K
                 TEMP = TEMP + A( L, I )*A( L, J )
220           END DO
              IF( BETA.EQ.ZERO )THEN
                 C( I, J ) = ALPHA*TEMP
              ELSE
                 C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
              END IF
230        END DO
240     END DO
     END IF
  END IF
  !
  RETURN
  !
  !     End of SSYRK .
  !
END SUBROUTINE SSYRK
!
!.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
!
subroutine copy8tor(px8,px,kn)
  !
  !   copy REAL(RKIND) :: array to real array
  !
  USE prec_rkind
  implicit none
  integer kn, i
  REAL(RKIND) :: px8(kn)
  REAL(RKIND) ::   px (kn)
  !
  do i=1,kn
     px(i) = px8(i)
  end do
  !
  return
end subroutine copy8tor
!     
subroutine copyrto8(px,px8,kn)
  !
  !   copy REAL(RKIND) :: array to real array
  !
  USE prec_rkind
  implicit none
  integer kn, i
  REAL(RKIND) :: px8(kn)
  REAL(RKIND) :: px (kn)
  !
  do i=1,kn
     px8(i) = px(i)
  end do
  !
  return
end subroutine copyrto8
