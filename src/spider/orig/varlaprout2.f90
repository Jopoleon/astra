integer function isamax(n,sx,incx)
  !
  !     finds the index of element having max. absolute value.
  !     jack dongarra, linpack, 3/11/78.
  !     modified 3/93 to return if incx .le. 0.
  !     modified 12/3/93, array(1) declarations changed to array(*)
  !
  USE prec_rkind
  implicit none
  REAL(RKIND) :: sx(*),smax
  integer i,incx,ix,n
  !
  isamax = 0
  if( n.lt.1 .or. incx.le.0 ) return
  isamax = 1
  if(n.eq.1)return
  if(incx.eq.1)go to 20
  !
  !        code for increment not equal to 1
  !
  ix = 1
  smax = abs(sx(1))
  ix = ix + incx
  do 10 i = 2,n
     if(abs(sx(ix)).le.smax) go to 5
     isamax = i
     smax = abs(sx(ix))
5    ix = ix + incx
10 END DO
  return
  !
  !        code for increment equal to 1
  !
20 smax = abs(sx(1))
  do 30 i = 2,n
     if(abs(sx(i)).le.smax) go to 30
     isamax = i
     smax = abs(sx(i))
30 END DO
  return
end function isamax
FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, &
     &                 N4 )
  !
  !  -- LAPACK auxiliary routine (version 2.0) --
  !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
  !     Courant Institute, Argonne National Lab, and Rice University
  !     September 30, 1994
  !
  !     .. Scalar Arguments ..
  USE prec_rkind
  implicit none
  CHARACTER*( * )    NAME, OPTS
  INTEGER            ISPEC, N1, N2, N3, N4, ILAENV
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  ILAENV is called from the LAPACK routines to choose problem-dependent
  !  parameters for the local environment.  See ISPEC for a description of
  !  the parameters.
  !
  !  This version provides a set of parameters which should give good,
  !  but not optimal, performance on many of the currently available
  !  computers.  Users are encouraged to modify this subroutine to set
  !  the tuning parameters for their particular machine using the option
  !  and problem size information in the arguments.
  !
  !  This routine will not function correctly if it is converted to all
  !  lower case.  Converting it to all upper case is allowed.
  !
  !  Arguments
  !  =========
  !
  !  ISPEC   (input) INTEGER
  !          Specifies the parameter to be returned as the value of
  !          ILAENV.
  !          = 1: the optimal blocksize; if this value is 1, an unblocked
  !               algorithm will give the best performance.
  !          = 2: the minimum block size for which the block routine
  !               should be used; if the usable block size is less than
  !               this value, an unblocked routine should be used.
  !          = 3: the crossover point (in a block routine, for N less
  !               than this value, an unblocked routine should be used)
  !          = 4: the number of shifts, used in the nonsymmetric
  !               eigenvalue routines
  !          = 5: the minimum column dimension for blocking to be used;
  !               rectangular blocks must have dimension at least k by m,
  !               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
  !          = 6: the crossover point for the SVD (when reducing an m by n
  !               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
  !               this value, a QR factorization is used first to reduce
  !               the matrix to a triangular form.)
  !          = 7: the number of processors
  !          = 8: the crossover point for the multishift QR and QZ methods
  !               for nonsymmetric eigenvalue problems.
  !
  !  NAME    (input) CHARACTER*(*)
  !          The name of the calling subroutine, in either upper case or
  !          lower case.
  !
  !  OPTS    (input) CHARACTER*(*)
  !          The character options to the subroutine NAME, concatenated
  !          into a single character string.  For example, UPLO = 'U',
  !          TRANS = 'T', and DIAG = 'N' for a triangular routine would
  !          be specified as OPTS = 'UTN'.
  !
  !  N1      (input) INTEGER
  !  N2      (input) INTEGER
  !  N3      (input) INTEGER
  !  N4      (input) INTEGER
  !          Problem dimensions for the subroutine NAME; these may not all
  !          be required.
  !
  ! (ILAENV) (output) INTEGER
  !          >= 0: the value of the parameter specified by ISPEC
  !          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
  !
  !  Further Details
  !  ===============
  !
  !  The following conventions have been used when calling ILAENV from the
  !  LAPACK routines:
  !  1)  OPTS is a concatenation of all of the character options to
  !      subroutine NAME, in the same order that they appear in the
  !      argument list for NAME, even if they are not used in determining
  !      the value of the parameter specified by ISPEC.
  !  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
  !      that they appear in the argument list for NAME.  N1 is used
  !      first, N2 second, and so on, and unused problem dimensions are
  !      passed a value of -1.
  !  3)  The parameter value returned by ILAENV is checked for validity in
  !      the calling subroutine.  For example, ILAENV is used to retrieve
  !      the optimal blocksize for STRTRI as follows:
  !
  !      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
  !      IF( NB.LE.1 ) NB = MAX( 1, N )
  !
  !  =====================================================================
  !
  !     .. Local Scalars ..
  LOGICAL            CNAME, SNAME
  CHARACTER*1        C1
  CHARACTER*2        C2, C4
  CHARACTER*3        C3
  CHARACTER*6        SUBNAM
  INTEGER            I, IC, IZ, NB, NBMIN, NX
  !     ..
  !     .. Intrinsic Functions ..
  INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
  !     ..
  !     .. Executable Statements ..
  !
  GO TO ( 100, 100, 100, 400, 500, 600, 700, 800 ) ISPEC
  !
  !     Invalid value for ISPEC
  !
  ILAENV = -1
  RETURN
  !
100 CONTINUE
  !
  !     Convert NAME to upper case if the first character is lower case.
  !
  ILAENV = 1
  SUBNAM = NAME
  IC = ICHAR( SUBNAM( 1:1 ) )
  IZ = ICHAR( 'Z' )
  IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
     !
     !        ASCII character set
     !
     IF( IC.GE.97 .AND. IC.LE.122 ) THEN
        SUBNAM( 1:1 ) = CHAR( IC-32 )
        DO 10 I = 2, 6
           IC = ICHAR( SUBNAM( I:I ) )
           IF( IC.GE.97 .AND. IC.LE.122 ) &
                &            SUBNAM( I:I ) = CHAR( IC-32 )
10      END DO
     END IF
     !
  ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
     !
     !        EBCDIC character set
     !
     IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
          &       ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
          &       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
        SUBNAM( 1:1 ) = CHAR( IC+64 )
        DO 20 I = 2, 6
           IC = ICHAR( SUBNAM( I:I ) )
           IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR. &
                &             ( IC.GE.145 .AND. IC.LE.153 ) .OR. &
                &             ( IC.GE.162 .AND. IC.LE.169 ) ) &
                &            SUBNAM( I:I ) = CHAR( IC+64 )
20      END DO
     END IF
     !
  ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
     !
     !        Prime machines:  ASCII+128
     !
     IF( IC.GE.225 .AND. IC.LE.250 ) THEN
        SUBNAM( 1:1 ) = CHAR( IC-32 )
        DO 30 I = 2, 6
           IC = ICHAR( SUBNAM( I:I ) )
           IF( IC.GE.225 .AND. IC.LE.250 ) &
                &            SUBNAM( I:I ) = CHAR( IC-32 )
30      END DO
     END IF
  END IF
  !
  C1 = SUBNAM( 1:1 )
  SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
  CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
  IF( .NOT.( CNAME .OR. SNAME ) ) &
       &   RETURN
  C2 = SUBNAM( 2:3 )
  C3 = SUBNAM( 4:6 )
  C4 = C3( 2:3 )
  !
  GO TO ( 110, 200, 300 ) ISPEC
  !
110 CONTINUE
  !
  !     ISPEC = 1:  block size
  !
  !     In these examples, separate code is provided for setting NB for
  !     real and complex.  We assume that NB will take the same value in
  !     single or double precision.
  !
  NB = 1
  !
  IF( C2.EQ.'GE' ) THEN
     IF( C3.EQ.'TRF' ) THEN
        IF( SNAME ) THEN
           NB = 64
        ELSE
           NB = 64
        END IF
     ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
          &            C3.EQ.'QLF' ) THEN
        IF( SNAME ) THEN
           NB = 32
        ELSE
           NB = 32
        END IF
     ELSE IF( C3.EQ.'HRD' ) THEN
        IF( SNAME ) THEN
           NB = 32
        ELSE
           NB = 32
        END IF
     ELSE IF( C3.EQ.'BRD' ) THEN
        IF( SNAME ) THEN
           NB = 32
        ELSE
           NB = 32
        END IF
     ELSE IF( C3.EQ.'TRI' ) THEN
        IF( SNAME ) THEN
           NB = 64
        ELSE
           NB = 64
        END IF
     END IF
  ELSE IF( C2.EQ.'PO' ) THEN
     IF( C3.EQ.'TRF' ) THEN
        IF( SNAME ) THEN
           NB = 64
        ELSE
           NB = 64
        END IF
     END IF
  ELSE IF( C2.EQ.'SY' ) THEN
     IF( C3.EQ.'TRF' ) THEN
        IF( SNAME ) THEN
           NB = 64
        ELSE
           NB = 64
        END IF
     ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
        NB = 1
     ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
        NB = 64
     END IF
  ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
     IF( C3.EQ.'TRF' ) THEN
        NB = 64
     ELSE IF( C3.EQ.'TRD' ) THEN
        NB = 1
     ELSE IF( C3.EQ.'GST' ) THEN
        NB = 64
     END IF
  ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
     IF( C3( 1:1 ).EQ.'G' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
             &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
             &          C4.EQ.'BR' ) THEN
           NB = 32
        END IF
     ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
             &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
             &          C4.EQ.'BR' ) THEN
           NB = 32
        END IF
     END IF
  ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
     IF( C3( 1:1 ).EQ.'G' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
             &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
             &          C4.EQ.'BR' ) THEN
           NB = 32
        END IF
     ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
             &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
             &          C4.EQ.'BR' ) THEN
           NB = 32
        END IF
     END IF
  ELSE IF( C2.EQ.'GB' ) THEN
     IF( C3.EQ.'TRF' ) THEN
        IF( SNAME ) THEN
           IF( N4.LE.64 ) THEN
              NB = 1
           ELSE
              NB = 32
           END IF
        ELSE
           IF( N4.LE.64 ) THEN
              NB = 1
           ELSE
              NB = 32
           END IF
        END IF
     END IF
  ELSE IF( C2.EQ.'PB' ) THEN
     IF( C3.EQ.'TRF' ) THEN
        IF( SNAME ) THEN
           IF( N2.LE.64 ) THEN
              NB = 1
           ELSE
              NB = 32
           END IF
        ELSE
           IF( N2.LE.64 ) THEN
              NB = 1
           ELSE
              NB = 32
           END IF
        END IF
     END IF
  ELSE IF( C2.EQ.'TR' ) THEN
     IF( C3.EQ.'TRI' ) THEN
        IF( SNAME ) THEN
           NB = 64
        ELSE
           NB = 64
        END IF
     END IF
  ELSE IF( C2.EQ.'LA' ) THEN
     IF( C3.EQ.'UUM' ) THEN
        IF( SNAME ) THEN
           NB = 64
        ELSE
           NB = 64
        END IF
     END IF
  ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
     IF( C3.EQ.'EBZ' ) THEN
        NB = 1
     END IF
  END IF
  ILAENV = NB
  RETURN
  !
200 CONTINUE
  !
  !     ISPEC = 2:  minimum block size
  !
  NBMIN = 2
  IF( C2.EQ.'GE' ) THEN
     IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
          &       C3.EQ.'QLF' ) THEN
        IF( SNAME ) THEN
           NBMIN = 2
        ELSE
           NBMIN = 2
        END IF
     ELSE IF( C3.EQ.'HRD' ) THEN
        IF( SNAME ) THEN
           NBMIN = 2
        ELSE
           NBMIN = 2
        END IF
     ELSE IF( C3.EQ.'BRD' ) THEN
        IF( SNAME ) THEN
           NBMIN = 2
        ELSE
           NBMIN = 2
        END IF
     ELSE IF( C3.EQ.'TRI' ) THEN
        IF( SNAME ) THEN
           NBMIN = 2
        ELSE
           NBMIN = 2
        END IF
     END IF
  ELSE IF( C2.EQ.'SY' ) THEN
     IF( C3.EQ.'TRF' ) THEN
        IF( SNAME ) THEN
           NBMIN = 8
        ELSE
           NBMIN = 8
        END IF
     ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
        NBMIN = 2
     END IF
  ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
     IF( C3.EQ.'TRD' ) THEN
        NBMIN = 2
     END IF
  ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
     IF( C3( 1:1 ).EQ.'G' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
             &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
             &          C4.EQ.'BR' ) THEN
           NBMIN = 2
        END IF
     ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
             &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
             &          C4.EQ.'BR' ) THEN
           NBMIN = 2
        END IF
     END IF
  ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
     IF( C3( 1:1 ).EQ.'G' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
             &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
             &          C4.EQ.'BR' ) THEN
           NBMIN = 2
        END IF
     ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
             &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
             &          C4.EQ.'BR' ) THEN
           NBMIN = 2
        END IF
     END IF
  END IF
  ILAENV = NBMIN
  RETURN
  !
300 CONTINUE
  !
  !     ISPEC = 3:  crossover point
  !
  NX = 0
  IF( C2.EQ.'GE' ) THEN
     IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR. &
          &       C3.EQ.'QLF' ) THEN
        IF( SNAME ) THEN
           NX = 128
        ELSE
           NX = 128
        END IF
     ELSE IF( C3.EQ.'HRD' ) THEN
        IF( SNAME ) THEN
           NX = 128
        ELSE
           NX = 128
        END IF
     ELSE IF( C3.EQ.'BRD' ) THEN
        IF( SNAME ) THEN
           NX = 128
        ELSE
           NX = 128
        END IF
     END IF
  ELSE IF( C2.EQ.'SY' ) THEN
     IF( SNAME .AND. C3.EQ.'TRD' ) THEN
        NX = 1
     END IF
  ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
     IF( C3.EQ.'TRD' ) THEN
        NX = 1
     END IF
  ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
     IF( C3( 1:1 ).EQ.'G' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
             &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
             &          C4.EQ.'BR' ) THEN
           NX = 128
        END IF
     END IF
  ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
     IF( C3( 1:1 ).EQ.'G' ) THEN
        IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR. &
             &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR. &
             &          C4.EQ.'BR' ) THEN
           NX = 128
        END IF
     END IF
  END IF
  ILAENV = NX
  RETURN
  !
400 CONTINUE
  !
  !     ISPEC = 4:  number of shifts (used by xHSEQR)
  !
  ILAENV = 6
  RETURN
  !
500 CONTINUE
  !
  !     ISPEC = 5:  minimum column dimension (not used)
  !
  ILAENV = 2
  RETURN
  !
600 CONTINUE 
  !
  !     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
  !
  ILAENV = INT( REAL( MIN( N1, N2 ),RKIND)*1.6_RKIND )
  RETURN
  !
700 CONTINUE
  !
  !     ISPEC = 7:  number of processors (not used)
  !
  ILAENV = 1
  RETURN
  !
800 CONTINUE
  !
  !     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
  !
  ILAENV = 50
  RETURN
  !
  !     End of ILAENV
  !
END FUNCTION ILAENV
function sdot(n,sx,incx,sy,incy)
  !
  !     forms the dot product of two vectors.
  !     uses unrolled loops for increments equal to one.
  !     jack dongarra, linpack, 3/11/78.
  !     modified 12/3/93, array(1) declarations changed to array(*)
  !
  USE prec_rkind
  implicit none
  REAL(RKIND) :: sx(*),sy(*),stemp, sdot
  integer i,incx,incy,ix,iy,m,mp1,n
  !
  stemp = 0.0e0_RKIND
  sdot = 0.0e0_RKIND
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
     stemp = stemp + sx(ix)*sy(iy)
     ix = ix + incx
     iy = iy + incy
10 END DO
  sdot = stemp
  return
  !
  !        code for both increments equal to 1
  !
  !
  !        clean-up loop
  !
20 m = mod(n,5)
  if( m .eq. 0 ) go to 40
  do 30 i = 1,m
     stemp = stemp + sx(i)*sy(i)
30 END DO
  if( n .lt. 5 ) go to 60
40 mp1 = m + 1
  do 50 i = mp1,n,5
     stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) + &
          &   sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
50 END DO
60 sdot = stemp
  return
end function sdot
LOGICAL          FUNCTION LSAME( CA, CB )
  !
  !  -- LAPACK auxiliary routine (version 2.0) --
  !     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
  !     Courant Institute, Argonne National Lab, and Rice University
  !     September 30, 1994
  !
  !     .. Scalar Arguments ..
  USE prec_rkind
  implicit none
  CHARACTER          CA, CB
  !     ..
  !
  !  Purpose
  !  =======
  !
  !  LSAME returns .TRUE. if CA is the same letter as CB regardless of
  !  case.
  !
  !  Arguments
  !  =========
  !
  !  CA      (input) CHARACTER*1
  !  CB      (input) CHARACTER*1
  !          CA and CB specify the single characters to be compared.
  !
  ! =====================================================================
  !
  !     .. Intrinsic Functions ..
  INTRINSIC          ICHAR
  !     ..
  !     .. Local Scalars ..
  INTEGER            INTA, INTB, ZCODE
  !     ..
  !     .. Executable Statements ..
  !
  !     Test if the characters are equal
  !
  LSAME = CA.EQ.CB
  IF( LSAME ) &
       &   RETURN
  !
  !     Now test for equivalence if both characters are alphabetic.
  !
  ZCODE = ICHAR( 'Z' )
  !
  !     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
  !     machines, on which ICHAR returns a value with bit 8 set.
  !     ICHAR('A') on Prime machines returns 193 which is the same as
  !     ICHAR('A') on an EBCDIC machine.
  !
  INTA = ICHAR( CA )
  INTB = ICHAR( CB )
  !
  IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
     !
     !        ASCII is assumed - ZCODE is the ASCII code of either lower or
     !        upper case 'Z'.
     !
     IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
     IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
     !
  ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
     !
     !        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
     !        upper case 'Z'.
     !
     IF( INTA.GE.129 .AND. INTA.LE.137 .OR. &
          &       INTA.GE.145 .AND. INTA.LE.153 .OR. &
          &       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
     IF( INTB.GE.129 .AND. INTB.LE.137 .OR. &
          &       INTB.GE.145 .AND. INTB.LE.153 .OR. &
          &       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
     !
  ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
     !
     !        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
     !        plus 128 of either lower or upper case 'Z'.
     !
     IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
     IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
  END IF
  LSAME = INTA.EQ.INTB
  !
  !     RETURN
  !
  !     End of LSAME
  !
END FUNCTION LSAME
