*> \brief \b ZGGBAK
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZGGBAK + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zggbak.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zggbak.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zggbak.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGGBAK( JOB, SIDE, N, ILO, IHI, LSCALE, RSCALE, M, V,
*                          LDV, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          JOB, SIDE
*       INTEGER            IHI, ILO, INFO, LDV, M, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   LSCALE( * ), RSCALE( * )
*       COMPLEX*16         V( LDV, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGGBAK forms the right or left eigenvectors of a complex generalized
*> eigenvalue problem A*x = lambda*B*x, by backward transformation on
*> the computed eigenvectors of the balanced pair of matrices output by
*> ZGGBAL.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOB
*> \verbatim
*>          JOB is CHARACTER*1
*>          Specifies the type of backward transformation required:
*>          = 'N':  do nothing, return immediately;
*>          = 'P':  do backward transformation for permutation only;
*>          = 'S':  do backward transformation for scaling only;
*>          = 'B':  do backward transformations for both permutation and
*>                  scaling.
*>          JOB must be the same as the argument JOB supplied to ZGGBAL.
*> \endverbatim
*>
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'R':  V contains right eigenvectors;
*>          = 'L':  V contains left eigenvectors.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of rows of the matrix V.  N >= 0.
*> \endverbatim
*>
*> \param[in] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[in] IHI
*> \verbatim
*>          IHI is INTEGER
*>          The integers ILO and IHI determined by ZGGBAL.
*>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
*> \endverbatim
*>
*> \param[in] LSCALE
*> \verbatim
*>          LSCALE is DOUBLE PRECISION array, dimension (N)
*>          Details of the permutations and/or scaling factors applied
*>          to the left side of A and B, as returned by ZGGBAL.
*> \endverbatim
*>
*> \param[in] RSCALE
*> \verbatim
*>          RSCALE is DOUBLE PRECISION array, dimension (N)
*>          Details of the permutations and/or scaling factors applied
*>          to the right side of A and B, as returned by ZGGBAL.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of columns of the matrix V.  M >= 0.
*> \endverbatim
*>
*> \param[in,out] V
*> \verbatim
*>          V is COMPLEX*16 array, dimension (LDV,M)
*>          On entry, the matrix of right or left eigenvectors to be
*>          transformed, as returned by ZTGEVC.
*>          On exit, V is overwritten by the transformed eigenvectors.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>          The leading dimension of the matrix V. LDV >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2015
*
*> \ingroup complex16GBcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  See R.C. Ward, Balancing the generalized eigenvalue problem,
*>                 SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZGGBAK( JOB, SIDE, N, ILO, IHI, LSCALE, RSCALE, M, V,
     $                   LDV, INFO )
*
*  -- LAPACK computational routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      CHARACTER          JOB, SIDE
      INTEGER            IHI, ILO, INFO, LDV, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   LSCALE( * ), RSCALE( * )
      COMPLEX*16         V( LDV, * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LEFTV, RIGHTV
      INTEGER            I, K
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZDSCAL, ZSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, INT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      RIGHTV = LSAME( SIDE, 'R' )
      LEFTV = LSAME( SIDE, 'L' )
*
      INFO = 0
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND.
     $    .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.RIGHTV .AND. .NOT.LEFTV ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 ) THEN
         INFO = -4
      ELSE IF( N.EQ.0 .AND. IHI.EQ.0 .AND. ILO.NE.1 ) THEN
         INFO = -4
      ELSE IF( N.GT.0 .AND. ( IHI.LT.ILO .OR. IHI.GT.MAX( 1, N ) ) )
     $   THEN
         INFO = -5
      ELSE IF( N.EQ.0 .AND. ILO.EQ.1 .AND. IHI.NE.0 ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -8
      ELSE IF( LDV.LT.MAX( 1, N ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGGBAK', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
      IF( M.EQ.0 )
     $   RETURN
      IF( LSAME( JOB, 'N' ) )
     $   RETURN
*
      IF( ILO.EQ.IHI )
     $   GO TO 30
*
*     Backward balance
*
      IF( LSAME( JOB, 'S' ) .OR. LSAME( JOB, 'B' ) ) THEN
*
*        Backward transformation on right eigenvectors
*
         IF( RIGHTV ) THEN
            DO 10 I = ILO, IHI
               CALL ZDSCAL( M, RSCALE( I ), V( I, 1 ), LDV )
   10       CONTINUE
         END IF
*
*        Backward transformation on left eigenvectors
*
         IF( LEFTV ) THEN
            DO 20 I = ILO, IHI
               CALL ZDSCAL( M, LSCALE( I ), V( I, 1 ), LDV )
   20       CONTINUE
         END IF
      END IF
*
*     Backward permutation
*
   30 CONTINUE
      IF( LSAME( JOB, 'P' ) .OR. LSAME( JOB, 'B' ) ) THEN
*
*        Backward permutation on right eigenvectors
*
         IF( RIGHTV ) THEN
            IF( ILO.EQ.1 )
     $         GO TO 50
            DO 40 I = ILO - 1, 1, -1
               K = INT(RSCALE( I ))
               IF( K.EQ.I )
     $            GO TO 40
               CALL ZSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
   40       CONTINUE
*
   50       CONTINUE
            IF( IHI.EQ.N )
     $         GO TO 70
            DO 60 I = IHI + 1, N
               K = INT(RSCALE( I ))
               IF( K.EQ.I )
     $            GO TO 60
               CALL ZSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
   60       CONTINUE
         END IF
*
*        Backward permutation on left eigenvectors
*
   70    CONTINUE
         IF( LEFTV ) THEN
            IF( ILO.EQ.1 )
     $         GO TO 90
            DO 80 I = ILO - 1, 1, -1
               K = INT(LSCALE( I ))
               IF( K.EQ.I )
     $            GO TO 80
               CALL ZSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
   80       CONTINUE
*
   90       CONTINUE
            IF( IHI.EQ.N )
     $         GO TO 110
            DO 100 I = IHI + 1, N
               K = INT(LSCALE( I ))
               IF( K.EQ.I )
     $            GO TO 100
               CALL ZSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
  100       CONTINUE
         END IF
      END IF
*
  110 CONTINUE
*
      RETURN
*
*     End of ZGGBAK
*
      END
*> \brief \b ZGGBAL
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZGGBAL + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zggbal.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zggbal.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zggbal.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGGBAL( JOB, N, A, LDA, B, LDB, ILO, IHI, LSCALE,
*                          RSCALE, WORK, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          JOB
*       INTEGER            IHI, ILO, INFO, LDA, LDB, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   LSCALE( * ), RSCALE( * ), WORK( * )
*       COMPLEX*16         A( LDA, * ), B( LDB, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGGBAL balances a pair of general complex matrices (A,B).  This
*> involves, first, permuting A and B by similarity transformations to
*> isolate eigenvalues in the first 1 to ILO$-$1 and last IHI+1 to N
*> elements on the diagonal; and second, applying a diagonal similarity
*> transformation to rows and columns ILO to IHI to make the rows
*> and columns as close in norm as possible. Both steps are optional.
*>
*> Balancing may reduce the 1-norm of the matrices, and improve the
*> accuracy of the computed eigenvalues and/or eigenvectors in the
*> generalized eigenvalue problem A*x = lambda*B*x.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOB
*> \verbatim
*>          JOB is CHARACTER*1
*>          Specifies the operations to be performed on A and B:
*>          = 'N':  none:  simply set ILO = 1, IHI = N, LSCALE(I) = 1.0
*>                  and RSCALE(I) = 1.0 for i=1,...,N;
*>          = 'P':  permute only;
*>          = 'S':  scale only;
*>          = 'B':  both permute and scale.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices A and B.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the input matrix A.
*>          On exit, A is overwritten by the balanced matrix.
*>          If JOB = 'N', A is not referenced.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A. LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX*16 array, dimension (LDB,N)
*>          On entry, the input matrix B.
*>          On exit, B is overwritten by the balanced matrix.
*>          If JOB = 'N', B is not referenced.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B. LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[out] IHI
*> \verbatim
*>          IHI is INTEGER
*>          ILO and IHI are set to integers such that on exit
*>          A(i,j) = 0 and B(i,j) = 0 if i > j and
*>          j = 1,...,ILO-1 or i = IHI+1,...,N.
*>          If JOB = 'N' or 'S', ILO = 1 and IHI = N.
*> \endverbatim
*>
*> \param[out] LSCALE
*> \verbatim
*>          LSCALE is DOUBLE PRECISION array, dimension (N)
*>          Details of the permutations and scaling factors applied
*>          to the left side of A and B.  If P(j) is the index of the
*>          row interchanged with row j, and D(j) is the scaling factor
*>          applied to row j, then
*>            LSCALE(j) = P(j)    for J = 1,...,ILO-1
*>                      = D(j)    for J = ILO,...,IHI
*>                      = P(j)    for J = IHI+1,...,N.
*>          The order in which the interchanges are made is N to IHI+1,
*>          then 1 to ILO-1.
*> \endverbatim
*>
*> \param[out] RSCALE
*> \verbatim
*>          RSCALE is DOUBLE PRECISION array, dimension (N)
*>          Details of the permutations and scaling factors applied
*>          to the right side of A and B.  If P(j) is the index of the
*>          column interchanged with column j, and D(j) is the scaling
*>          factor applied to column j, then
*>            RSCALE(j) = P(j)    for J = 1,...,ILO-1
*>                      = D(j)    for J = ILO,...,IHI
*>                      = P(j)    for J = IHI+1,...,N.
*>          The order in which the interchanges are made is N to IHI+1,
*>          then 1 to ILO-1.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (lwork)
*>          lwork must be at least max(1,6*N) when JOB = 'S' or 'B', and
*>          at least 1 when JOB = 'N' or 'P'.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2015
*
*> \ingroup complex16GBcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  See R.C. WARD, Balancing the generalized eigenvalue problem,
*>                 SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZGGBAL( JOB, N, A, LDA, B, LDB, ILO, IHI, LSCALE,
     $                   RSCALE, WORK, INFO )
*
*  -- LAPACK computational routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            IHI, ILO, INFO, LDA, LDB, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   LSCALE( * ), RSCALE( * ), WORK( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE
      PARAMETER          ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   THREE, SCLFAC
      PARAMETER          ( THREE = 3.0D+0, SCLFAC = 1.0D+1 )
      COMPLEX*16         CZERO
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ICAB, IFLOW, IP1, IR, IRAB, IT, J, JC, JP1,
     $                   K, KOUNT, L, LCAB, LM1, LRAB, LSFMAX, LSFMIN,
     $                   M, NR, NRP2
      DOUBLE PRECISION   ALPHA, BASL, BETA, CAB, CMAX, COEF, COEF2,
     $                   COEF5, COR, EW, EWC, GAMMA, PGAMMA, RAB, SFMAX,
     $                   SFMIN, SUM, T, TA, TB, TC
      COMPLEX*16         CDUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IZAMAX
      DOUBLE PRECISION   DDOT, DLAMCH
      EXTERNAL           LSAME, IZAMAX, DDOT, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DSCAL, XERBLA, ZDSCAL, ZSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG, INT, LOG10, MAX, MIN, SIGN
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND.
     $    .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGGBAL', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         ILO = 1
         IHI = N
         RETURN
      END IF
*
      IF( N.EQ.1 ) THEN
         ILO = 1
         IHI = N
         LSCALE( 1 ) = ONE
         RSCALE( 1 ) = ONE
         RETURN
      END IF
*
      IF( LSAME( JOB, 'N' ) ) THEN
         ILO = 1
         IHI = N
         DO 10 I = 1, N
            LSCALE( I ) = ONE
            RSCALE( I ) = ONE
   10    CONTINUE
         RETURN
      END IF
*
      K = 1
      L = N
      IF( LSAME( JOB, 'S' ) )
     $   GO TO 190
*
      GO TO 30
*
*     Permute the matrices A and B to isolate the eigenvalues.
*
*     Find row with one nonzero in columns 1 through L
*
   20 CONTINUE
      L = LM1
      IF( L.NE.1 )
     $   GO TO 30
*
      RSCALE( 1 ) = 1
      LSCALE( 1 ) = 1
      GO TO 190
*
   30 CONTINUE
      LM1 = L - 1
      DO 80 I = L, 1, -1
         DO 40 J = 1, LM1
            JP1 = J + 1
            IF( A( I, J ).NE.CZERO .OR. B( I, J ).NE.CZERO )
     $         GO TO 50
   40    CONTINUE
         J = L
         GO TO 70
*
   50    CONTINUE
         DO 60 J = JP1, L
            IF( A( I, J ).NE.CZERO .OR. B( I, J ).NE.CZERO )
     $         GO TO 80
   60    CONTINUE
         J = JP1 - 1
*
   70    CONTINUE
         M = L
         IFLOW = 1
         GO TO 160
   80 CONTINUE
      GO TO 100
*
*     Find column with one nonzero in rows K through N
*
   90 CONTINUE
      K = K + 1
*
  100 CONTINUE
      DO 150 J = K, L
         DO 110 I = K, LM1
            IP1 = I + 1
            IF( A( I, J ).NE.CZERO .OR. B( I, J ).NE.CZERO )
     $         GO TO 120
  110    CONTINUE
         I = L
         GO TO 140
  120    CONTINUE
         DO 130 I = IP1, L
            IF( A( I, J ).NE.CZERO .OR. B( I, J ).NE.CZERO )
     $         GO TO 150
  130    CONTINUE
         I = IP1 - 1
  140    CONTINUE
         M = K
         IFLOW = 2
         GO TO 160
  150 CONTINUE
      GO TO 190
*
*     Permute rows M and I
*
  160 CONTINUE
      LSCALE( M ) = I
      IF( I.EQ.M )
     $   GO TO 170
      CALL ZSWAP( N-K+1, A( I, K ), LDA, A( M, K ), LDA )
      CALL ZSWAP( N-K+1, B( I, K ), LDB, B( M, K ), LDB )
*
*     Permute columns M and J
*
  170 CONTINUE
      RSCALE( M ) = J
      IF( J.EQ.M )
     $   GO TO 180
      CALL ZSWAP( L, A( 1, J ), 1, A( 1, M ), 1 )
      CALL ZSWAP( L, B( 1, J ), 1, B( 1, M ), 1 )
*
  180 CONTINUE
      GO TO ( 20, 90 )IFLOW
*
  190 CONTINUE
      ILO = K
      IHI = L
*
      IF( LSAME( JOB, 'P' ) ) THEN
         DO 195 I = ILO, IHI
            LSCALE( I ) = ONE
            RSCALE( I ) = ONE
  195    CONTINUE
         RETURN
      END IF
*
      IF( ILO.EQ.IHI )
     $   RETURN
*
*     Balance the submatrix in rows ILO to IHI.
*
      NR = IHI - ILO + 1
      DO 200 I = ILO, IHI
         RSCALE( I ) = ZERO
         LSCALE( I ) = ZERO
*
         WORK( I ) = ZERO
         WORK( I+N ) = ZERO
         WORK( I+2*N ) = ZERO
         WORK( I+3*N ) = ZERO
         WORK( I+4*N ) = ZERO
         WORK( I+5*N ) = ZERO
  200 CONTINUE
*
*     Compute right side vector in resulting linear equations
*
      BASL = LOG10( SCLFAC )
      DO 240 I = ILO, IHI
         DO 230 J = ILO, IHI
            IF( A( I, J ).EQ.CZERO ) THEN
               TA = ZERO
               GO TO 210
            END IF
            TA = LOG10( CABS1( A( I, J ) ) ) / BASL
*
  210       CONTINUE
            IF( B( I, J ).EQ.CZERO ) THEN
               TB = ZERO
               GO TO 220
            END IF
            TB = LOG10( CABS1( B( I, J ) ) ) / BASL
*
  220       CONTINUE
            WORK( I+4*N ) = WORK( I+4*N ) - TA - TB
            WORK( J+5*N ) = WORK( J+5*N ) - TA - TB
  230    CONTINUE
  240 CONTINUE
*
      COEF = ONE / DBLE( 2*NR )
      COEF2 = COEF*COEF
      COEF5 = HALF*COEF2
      NRP2 = NR + 2
      BETA = ZERO
      IT = 1
*
*     Start generalized conjugate gradient iteration
*
  250 CONTINUE
*
      GAMMA = DDOT( NR, WORK( ILO+4*N ), 1, WORK( ILO+4*N ), 1 ) +
     $        DDOT( NR, WORK( ILO+5*N ), 1, WORK( ILO+5*N ), 1 )
*
      EW = ZERO
      EWC = ZERO
      DO 260 I = ILO, IHI
         EW = EW + WORK( I+4*N )
         EWC = EWC + WORK( I+5*N )
  260 CONTINUE
*
      GAMMA = COEF*GAMMA - COEF2*( EW**2+EWC**2 ) - COEF5*( EW-EWC )**2
      IF( GAMMA.EQ.ZERO )
     $   GO TO 350
      IF( IT.NE.1 )
     $   BETA = GAMMA / PGAMMA
      T = COEF5*( EWC-THREE*EW )
      TC = COEF5*( EW-THREE*EWC )
*
      CALL DSCAL( NR, BETA, WORK( ILO ), 1 )
      CALL DSCAL( NR, BETA, WORK( ILO+N ), 1 )
*
      CALL DAXPY( NR, COEF, WORK( ILO+4*N ), 1, WORK( ILO+N ), 1 )
      CALL DAXPY( NR, COEF, WORK( ILO+5*N ), 1, WORK( ILO ), 1 )
*
      DO 270 I = ILO, IHI
         WORK( I ) = WORK( I ) + TC
         WORK( I+N ) = WORK( I+N ) + T
  270 CONTINUE
*
*     Apply matrix to vector
*
      DO 300 I = ILO, IHI
         KOUNT = 0
         SUM = ZERO
         DO 290 J = ILO, IHI
            IF( A( I, J ).EQ.CZERO )
     $         GO TO 280
            KOUNT = KOUNT + 1
            SUM = SUM + WORK( J )
  280       CONTINUE
            IF( B( I, J ).EQ.CZERO )
     $         GO TO 290
            KOUNT = KOUNT + 1
            SUM = SUM + WORK( J )
  290    CONTINUE
         WORK( I+2*N ) = DBLE( KOUNT )*WORK( I+N ) + SUM
  300 CONTINUE
*
      DO 330 J = ILO, IHI
         KOUNT = 0
         SUM = ZERO
         DO 320 I = ILO, IHI
            IF( A( I, J ).EQ.CZERO )
     $         GO TO 310
            KOUNT = KOUNT + 1
            SUM = SUM + WORK( I+N )
  310       CONTINUE
            IF( B( I, J ).EQ.CZERO )
     $         GO TO 320
            KOUNT = KOUNT + 1
            SUM = SUM + WORK( I+N )
  320    CONTINUE
         WORK( J+3*N ) = DBLE( KOUNT )*WORK( J ) + SUM
  330 CONTINUE
*
      SUM = DDOT( NR, WORK( ILO+N ), 1, WORK( ILO+2*N ), 1 ) +
     $      DDOT( NR, WORK( ILO ), 1, WORK( ILO+3*N ), 1 )
      ALPHA = GAMMA / SUM
*
*     Determine correction to current iteration
*
      CMAX = ZERO
      DO 340 I = ILO, IHI
         COR = ALPHA*WORK( I+N )
         IF( ABS( COR ).GT.CMAX )
     $      CMAX = ABS( COR )
         LSCALE( I ) = LSCALE( I ) + COR
         COR = ALPHA*WORK( I )
         IF( ABS( COR ).GT.CMAX )
     $      CMAX = ABS( COR )
         RSCALE( I ) = RSCALE( I ) + COR
  340 CONTINUE
      IF( CMAX.LT.HALF )
     $   GO TO 350
*
      CALL DAXPY( NR, -ALPHA, WORK( ILO+2*N ), 1, WORK( ILO+4*N ), 1 )
      CALL DAXPY( NR, -ALPHA, WORK( ILO+3*N ), 1, WORK( ILO+5*N ), 1 )
*
      PGAMMA = GAMMA
      IT = IT + 1
      IF( IT.LE.NRP2 )
     $   GO TO 250
*
*     End generalized conjugate gradient iteration
*
  350 CONTINUE
      SFMIN = DLAMCH( 'S' )
      SFMAX = ONE / SFMIN
      LSFMIN = INT( LOG10( SFMIN ) / BASL+ONE )
      LSFMAX = INT( LOG10( SFMAX ) / BASL )
      DO 360 I = ILO, IHI
         IRAB = IZAMAX( N-ILO+1, A( I, ILO ), LDA )
         RAB = ABS( A( I, IRAB+ILO-1 ) )
         IRAB = IZAMAX( N-ILO+1, B( I, ILO ), LDB )
         RAB = MAX( RAB, ABS( B( I, IRAB+ILO-1 ) ) )
         LRAB = INT( LOG10( RAB+SFMIN ) / BASL+ONE )
         IR = INT(LSCALE( I ) + SIGN( HALF, LSCALE( I ) ))
         IR = MIN( MAX( IR, LSFMIN ), LSFMAX, LSFMAX-LRAB )
         LSCALE( I ) = SCLFAC**IR
         ICAB = IZAMAX( IHI, A( 1, I ), 1 )
         CAB = ABS( A( ICAB, I ) )
         ICAB = IZAMAX( IHI, B( 1, I ), 1 )
         CAB = MAX( CAB, ABS( B( ICAB, I ) ) )
         LCAB = INT( LOG10( CAB+SFMIN ) / BASL+ONE )
         JC = INT(RSCALE( I ) + SIGN( HALF, RSCALE( I ) ))
         JC = MIN( MAX( JC, LSFMIN ), LSFMAX, LSFMAX-LCAB )
         RSCALE( I ) = SCLFAC**JC
  360 CONTINUE
*
*     Row scaling of matrices A and B
*
      DO 370 I = ILO, IHI
         CALL ZDSCAL( N-ILO+1, LSCALE( I ), A( I, ILO ), LDA )
         CALL ZDSCAL( N-ILO+1, LSCALE( I ), B( I, ILO ), LDB )
  370 CONTINUE
*
*     Column scaling of matrices A and B
*
      DO 380 J = ILO, IHI
         CALL ZDSCAL( IHI, RSCALE( J ), A( 1, J ), 1 )
         CALL ZDSCAL( IHI, RSCALE( J ), B( 1, J ), 1 )
  380 CONTINUE
*
      RETURN
*
*     End of ZGGBAL
*
      END
*> \brief \b ZGGHRD
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZGGHRD + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgghrd.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgghrd.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgghrd.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGGHRD( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q,
*                          LDQ, Z, LDZ, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          COMPQ, COMPZ
*       INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
*      $                   Z( LDZ, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGGHRD reduces a pair of complex matrices (A,B) to generalized upper
*> Hessenberg form using unitary transformations, where A is a
*> general matrix and B is upper triangular.  The form of the
*> generalized eigenvalue problem is
*>    A*x = lambda*B*x,
*> and B is typically made upper triangular by computing its QR
*> factorization and moving the unitary matrix Q to the left side
*> of the equation.
*>
*> This subroutine simultaneously reduces A to a Hessenberg matrix H:
*>    Q**H*A*Z = H
*> and transforms B to another upper triangular matrix T:
*>    Q**H*B*Z = T
*> in order to reduce the problem to its standard form
*>    H*y = lambda*T*y
*> where y = Z**H*x.
*>
*> The unitary matrices Q and Z are determined as products of Givens
*> rotations.  They may either be formed explicitly, or they may be
*> postmultiplied into input matrices Q1 and Z1, so that
*>      Q1 * A * Z1**H = (Q1*Q) * H * (Z1*Z)**H
*>      Q1 * B * Z1**H = (Q1*Q) * T * (Z1*Z)**H
*> If Q1 is the unitary matrix from the QR factorization of B in the
*> original equation A*x = lambda*B*x, then ZGGHRD reduces the original
*> problem to generalized Hessenberg form.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] COMPQ
*> \verbatim
*>          COMPQ is CHARACTER*1
*>          = 'N': do not compute Q;
*>          = 'I': Q is initialized to the unit matrix, and the
*>                 unitary matrix Q is returned;
*>          = 'V': Q must contain a unitary matrix Q1 on entry,
*>                 and the product Q1*Q is returned.
*> \endverbatim
*>
*> \param[in] COMPZ
*> \verbatim
*>          COMPZ is CHARACTER*1
*>          = 'N': do not compute Z;
*>          = 'I': Z is initialized to the unit matrix, and the
*>                 unitary matrix Z is returned;
*>          = 'V': Z must contain a unitary matrix Z1 on entry,
*>                 and the product Z1*Z is returned.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices A and B.  N >= 0.
*> \endverbatim
*>
*> \param[in] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[in] IHI
*> \verbatim
*>          IHI is INTEGER
*>
*>          ILO and IHI mark the rows and columns of A which are to be
*>          reduced.  It is assumed that A is already upper triangular
*>          in rows and columns 1:ILO-1 and IHI+1:N.  ILO and IHI are
*>          normally set by a previous call to ZGGBAL; otherwise they
*>          should be set to 1 and N respectively.
*>          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA, N)
*>          On entry, the N-by-N general matrix to be reduced.
*>          On exit, the upper triangle and the first subdiagonal of A
*>          are overwritten with the upper Hessenberg matrix H, and the
*>          rest is set to zero.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX*16 array, dimension (LDB, N)
*>          On entry, the N-by-N upper triangular matrix B.
*>          On exit, the upper triangular matrix T = Q**H B Z.  The
*>          elements below the diagonal are set to zero.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] Q
*> \verbatim
*>          Q is COMPLEX*16 array, dimension (LDQ, N)
*>          On entry, if COMPQ = 'V', the unitary matrix Q1, typically
*>          from the QR factorization of B.
*>          On exit, if COMPQ='I', the unitary matrix Q, and if
*>          COMPQ = 'V', the product Q1*Q.
*>          Not referenced if COMPQ='N'.
*> \endverbatim
*>
*> \param[in] LDQ
*> \verbatim
*>          LDQ is INTEGER
*>          The leading dimension of the array Q.
*>          LDQ >= N if COMPQ='V' or 'I'; LDQ >= 1 otherwise.
*> \endverbatim
*>
*> \param[in,out] Z
*> \verbatim
*>          Z is COMPLEX*16 array, dimension (LDZ, N)
*>          On entry, if COMPZ = 'V', the unitary matrix Z1.
*>          On exit, if COMPZ='I', the unitary matrix Z, and if
*>          COMPZ = 'V', the product Z1*Z.
*>          Not referenced if COMPZ='N'.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z.
*>          LDZ >= N if COMPZ='V' or 'I'; LDZ >= 1 otherwise.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2015
*
*> \ingroup complex16OTHERcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  This routine reduces A to Hessenberg and B to triangular form by
*>  an unblocked reduction, as described in _Matrix_Computations_,
*>  by Golub and van Loan (Johns Hopkins Press).
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZGGHRD( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q,
     $                   LDQ, Z, LDZ, INFO )
*
*  -- LAPACK computational routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      CHARACTER          COMPQ, COMPZ
      INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
     $                   Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         CONE, CZERO
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ),
     $                   CZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILQ, ILZ
      INTEGER            ICOMPQ, ICOMPZ, JCOL, JROW
      DOUBLE PRECISION   C
      COMPLEX*16         CTEMP, S
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARTG, ZLASET, ZROT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
*     ..
*     .. Executable Statements ..
*
*     Decode COMPQ
*
      IF( LSAME( COMPQ, 'N' ) ) THEN
         ILQ = .FALSE.
         ICOMPQ = 1
      ELSE IF( LSAME( COMPQ, 'V' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 2
      ELSE IF( LSAME( COMPQ, 'I' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 3
      ELSE
         ICOMPQ = 0
      END IF
*
*     Decode COMPZ
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ILZ = .FALSE.
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 2
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 3
      ELSE
         ICOMPZ = 0
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF( ICOMPQ.LE.0 ) THEN
         INFO = -1
      ELSE IF( ICOMPZ.LE.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 ) THEN
         INFO = -4
      ELSE IF( IHI.GT.N .OR. IHI.LT.ILO-1 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( ( ILQ .AND. LDQ.LT.N ) .OR. LDQ.LT.1 ) THEN
         INFO = -11
      ELSE IF( ( ILZ .AND. LDZ.LT.N ) .OR. LDZ.LT.1 ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGGHRD', -INFO )
         RETURN
      END IF
*
*     Initialize Q and Z if desired.
*
      IF( ICOMPQ.EQ.3 )
     $   CALL ZLASET( 'Full', N, N, CZERO, CONE, Q, LDQ )
      IF( ICOMPZ.EQ.3 )
     $   CALL ZLASET( 'Full', N, N, CZERO, CONE, Z, LDZ )
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
*
*     Zero out lower triangle of B
*
      DO 20 JCOL = 1, N - 1
         DO 10 JROW = JCOL + 1, N
            B( JROW, JCOL ) = CZERO
   10    CONTINUE
   20 CONTINUE
*
*     Reduce A and B
*
      DO 40 JCOL = ILO, IHI - 2
*
         DO 30 JROW = IHI, JCOL + 2, -1
*
*           Step 1: rotate rows JROW-1, JROW to kill A(JROW,JCOL)
*
            CTEMP = A( JROW-1, JCOL )
            CALL ZLARTG( CTEMP, A( JROW, JCOL ), C, S,
     $                   A( JROW-1, JCOL ) )
            A( JROW, JCOL ) = CZERO
            CALL ZROT( N-JCOL, A( JROW-1, JCOL+1 ), LDA,
     $                 A( JROW, JCOL+1 ), LDA, C, S )
            CALL ZROT( N+2-JROW, B( JROW-1, JROW-1 ), LDB,
     $                 B( JROW, JROW-1 ), LDB, C, S )
            IF( ILQ )
     $         CALL ZROT( N, Q( 1, JROW-1 ), 1, Q( 1, JROW ), 1, C,
     $                    DCONJG( S ) )
*
*           Step 2: rotate columns JROW, JROW-1 to kill B(JROW,JROW-1)
*
            CTEMP = B( JROW, JROW )
            CALL ZLARTG( CTEMP, B( JROW, JROW-1 ), C, S,
     $                   B( JROW, JROW ) )
            B( JROW, JROW-1 ) = CZERO
            CALL ZROT( IHI, A( 1, JROW ), 1, A( 1, JROW-1 ), 1, C, S )
            CALL ZROT( JROW-1, B( 1, JROW ), 1, B( 1, JROW-1 ), 1, C,
     $                 S )
            IF( ILZ )
     $         CALL ZROT( N, Z( 1, JROW ), 1, Z( 1, JROW-1 ), 1, C, S )
   30    CONTINUE
   40 CONTINUE
*
      RETURN
*
*     End of ZGGHRD
*
      END
*> \brief \b ZHGEQZ
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZHGEQZ + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhgeqz.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhgeqz.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhgeqz.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZHGEQZ( JOB, COMPQ, COMPZ, N, ILO, IHI, H, LDH, T, LDT,
*                          ALPHA, BETA, Q, LDQ, Z, LDZ, WORK, LWORK,
*                          RWORK, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          COMPQ, COMPZ, JOB
*       INTEGER            IHI, ILO, INFO, LDH, LDQ, LDT, LDZ, LWORK, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   RWORK( * )
*       COMPLEX*16         ALPHA( * ), BETA( * ), H( LDH, * ),
*      $                   Q( LDQ, * ), T( LDT, * ), WORK( * ),
*      $                   Z( LDZ, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZHGEQZ computes the eigenvalues of a complex matrix pair (H,T),
*> where H is an upper Hessenberg matrix and T is upper triangular,
*> using the single-shift QZ method.
*> Matrix pairs of this type are produced by the reduction to
*> generalized upper Hessenberg form of a complex matrix pair (A,B):
*> 
*>    A = Q1*H*Z1**H,  B = Q1*T*Z1**H,
*> 
*> as computed by ZGGHRD.
*> 
*> If JOB='S', then the Hessenberg-triangular pair (H,T) is
*> also reduced to generalized Schur form,
*> 
*>    H = Q*S*Z**H,  T = Q*P*Z**H,
*> 
*> where Q and Z are unitary matrices and S and P are upper triangular.
*> 
*> Optionally, the unitary matrix Q from the generalized Schur
*> factorization may be postmultiplied into an input matrix Q1, and the
*> unitary matrix Z may be postmultiplied into an input matrix Z1.
*> If Q1 and Z1 are the unitary matrices from ZGGHRD that reduced
*> the matrix pair (A,B) to generalized Hessenberg form, then the output
*> matrices Q1*Q and Z1*Z are the unitary factors from the generalized
*> Schur factorization of (A,B):
*> 
*>    A = (Q1*Q)*S*(Z1*Z)**H,  B = (Q1*Q)*P*(Z1*Z)**H.
*> 
*> To avoid overflow, eigenvalues of the matrix pair (H,T)
*> (equivalently, of (A,B)) are computed as a pair of complex values
*> (alpha,beta).  If beta is nonzero, lambda = alpha / beta is an
*> eigenvalue of the generalized nonsymmetric eigenvalue problem (GNEP)
*>    A*x = lambda*B*x
*> and if alpha is nonzero, mu = beta / alpha is an eigenvalue of the
*> alternate form of the GNEP
*>    mu*A*y = B*y.
*> The values of alpha and beta for the i-th eigenvalue can be read
*> directly from the generalized Schur form:  alpha = S(i,i),
*> beta = P(i,i).
*>
*> Ref: C.B. Moler & G.W. Stewart, "An Algorithm for Generalized Matrix
*>      Eigenvalue Problems", SIAM J. Numer. Anal., 10(1973),
*>      pp. 241--256.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOB
*> \verbatim
*>          JOB is CHARACTER*1
*>          = 'E': Compute eigenvalues only;
*>          = 'S': Computer eigenvalues and the Schur form.
*> \endverbatim
*>
*> \param[in] COMPQ
*> \verbatim
*>          COMPQ is CHARACTER*1
*>          = 'N': Left Schur vectors (Q) are not computed;
*>          = 'I': Q is initialized to the unit matrix and the matrix Q
*>                 of left Schur vectors of (H,T) is returned;
*>          = 'V': Q must contain a unitary matrix Q1 on entry and
*>                 the product Q1*Q is returned.
*> \endverbatim
*>
*> \param[in] COMPZ
*> \verbatim
*>          COMPZ is CHARACTER*1
*>          = 'N': Right Schur vectors (Z) are not computed;
*>          = 'I': Q is initialized to the unit matrix and the matrix Z
*>                 of right Schur vectors of (H,T) is returned;
*>          = 'V': Z must contain a unitary matrix Z1 on entry and
*>                 the product Z1*Z is returned.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices H, T, Q, and Z.  N >= 0.
*> \endverbatim
*>
*> \param[in] ILO
*> \verbatim
*>          ILO is INTEGER
*> \endverbatim
*>
*> \param[in] IHI
*> \verbatim
*>          IHI is INTEGER
*>          ILO and IHI mark the rows and columns of H which are in
*>          Hessenberg form.  It is assumed that A is already upper
*>          triangular in rows and columns 1:ILO-1 and IHI+1:N.
*>          If N > 0, 1 <= ILO <= IHI <= N; if N = 0, ILO=1 and IHI=0.
*> \endverbatim
*>
*> \param[in,out] H
*> \verbatim
*>          H is COMPLEX*16 array, dimension (LDH, N)
*>          On entry, the N-by-N upper Hessenberg matrix H.
*>          On exit, if JOB = 'S', H contains the upper triangular
*>          matrix S from the generalized Schur factorization.
*>          If JOB = 'E', the diagonal of H matches that of S, but
*>          the rest of H is unspecified.
*> \endverbatim
*>
*> \param[in] LDH
*> \verbatim
*>          LDH is INTEGER
*>          The leading dimension of the array H.  LDH >= max( 1, N ).
*> \endverbatim
*>
*> \param[in,out] T
*> \verbatim
*>          T is COMPLEX*16 array, dimension (LDT, N)
*>          On entry, the N-by-N upper triangular matrix T.
*>          On exit, if JOB = 'S', T contains the upper triangular
*>          matrix P from the generalized Schur factorization.
*>          If JOB = 'E', the diagonal of T matches that of P, but
*>          the rest of T is unspecified.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T.  LDT >= max( 1, N ).
*> \endverbatim
*>
*> \param[out] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX*16 array, dimension (N)
*>          The complex scalars alpha that define the eigenvalues of
*>          GNEP.  ALPHA(i) = S(i,i) in the generalized Schur
*>          factorization.
*> \endverbatim
*>
*> \param[out] BETA
*> \verbatim
*>          BETA is COMPLEX*16 array, dimension (N)
*>          The real non-negative scalars beta that define the
*>          eigenvalues of GNEP.  BETA(i) = P(i,i) in the generalized
*>          Schur factorization.
*>
*>          Together, the quantities alpha = ALPHA(j) and beta = BETA(j)
*>          represent the j-th eigenvalue of the matrix pair (A,B), in
*>          one of the forms lambda = alpha/beta or mu = beta/alpha.
*>          Since either lambda or mu may overflow, they should not,
*>          in general, be computed.
*> \endverbatim
*>
*> \param[in,out] Q
*> \verbatim
*>          Q is COMPLEX*16 array, dimension (LDQ, N)
*>          On entry, if COMPZ = 'V', the unitary matrix Q1 used in the
*>          reduction of (A,B) to generalized Hessenberg form.
*>          On exit, if COMPZ = 'I', the unitary matrix of left Schur
*>          vectors of (H,T), and if COMPZ = 'V', the unitary matrix of
*>          left Schur vectors of (A,B).
*>          Not referenced if COMPZ = 'N'.
*> \endverbatim
*>
*> \param[in] LDQ
*> \verbatim
*>          LDQ is INTEGER
*>          The leading dimension of the array Q.  LDQ >= 1.
*>          If COMPQ='V' or 'I', then LDQ >= N.
*> \endverbatim
*>
*> \param[in,out] Z
*> \verbatim
*>          Z is COMPLEX*16 array, dimension (LDZ, N)
*>          On entry, if COMPZ = 'V', the unitary matrix Z1 used in the
*>          reduction of (A,B) to generalized Hessenberg form.
*>          On exit, if COMPZ = 'I', the unitary matrix of right Schur
*>          vectors of (H,T), and if COMPZ = 'V', the unitary matrix of
*>          right Schur vectors of (A,B).
*>          Not referenced if COMPZ = 'N'.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z.  LDZ >= 1.
*>          If COMPZ='V' or 'I', then LDZ >= N.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
*>          On exit, if INFO >= 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.  LWORK >= max(1,N).
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is DOUBLE PRECISION array, dimension (N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*>          = 1,...,N: the QZ iteration did not converge.  (H,T) is not
*>                     in Schur form, but ALPHA(i) and BETA(i),
*>                     i=INFO+1,...,N should be correct.
*>          = N+1,...,2*N: the shift calculation failed.  (H,T) is not
*>                     in Schur form, but ALPHA(i) and BETA(i),
*>                     i=INFO-N+1,...,N should be correct.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date April 2012
*
*> \ingroup complex16GEcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  We assume that complex ABS works as long as its value is less than
*>  overflow.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZHGEQZ( JOB, COMPQ, COMPZ, N, ILO, IHI, H, LDH, T, LDT,
     $                   ALPHA, BETA, Q, LDQ, Z, LDZ, WORK, LWORK,
     $                   RWORK, INFO )
*
*  -- LAPACK computational routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     April 2012
*
*     .. Scalar Arguments ..
      CHARACTER          COMPQ, COMPZ, JOB
      INTEGER            IHI, ILO, INFO, LDH, LDQ, LDT, LDZ, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         ALPHA( * ), BETA( * ), H( LDH, * ),
     $                   Q( LDQ, * ), T( LDT, * ), WORK( * ),
     $                   Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILAZR2, ILAZRO, ILQ, ILSCHR, ILZ, LQUERY
      INTEGER            ICOMPQ, ICOMPZ, IFIRST, IFRSTM, IITER, ILAST,
     $                   ILASTM, IN, ISCHUR, ISTART, J, JC, JCH, JITER,
     $                   JR, MAXIT
      DOUBLE PRECISION   ABSB, ANORM, ASCALE, ATOL, BNORM, BSCALE, BTOL,
     $                   C, SAFMIN, TEMP, TEMP2, TEMPR, ULP
      COMPLEX*16         ABI22, AD11, AD12, AD21, AD22, CTEMP, CTEMP2,
     $                   CTEMP3, ESHIFT, RTDISC, S, SHIFT, SIGNBC, T1,
     $                   U12, X
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, ZLANHS
      EXTERNAL           LSAME, DLAMCH, ZLANHS
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARTG, ZLASET, ZROT, ZSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DCONJG, DIMAG, MAX, MIN,
     $                   SQRT
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   ABS1
*     ..
*     .. Statement Function definitions ..
      ABS1( X ) = ABS( DBLE( X ) ) + ABS( DIMAG( X ) )
*     ..
*     .. Executable Statements ..
*
*     Decode JOB, COMPQ, COMPZ
*
      IF( LSAME( JOB, 'E' ) ) THEN
         ILSCHR = .FALSE.
         ISCHUR = 1
      ELSE IF( LSAME( JOB, 'S' ) ) THEN
         ILSCHR = .TRUE.
         ISCHUR = 2
      ELSE
         ISCHUR = 0
      END IF
*
      IF( LSAME( COMPQ, 'N' ) ) THEN
         ILQ = .FALSE.
         ICOMPQ = 1
      ELSE IF( LSAME( COMPQ, 'V' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 2
      ELSE IF( LSAME( COMPQ, 'I' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 3
      ELSE
         ICOMPQ = 0
      END IF
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ILZ = .FALSE.
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 2
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 3
      ELSE
         ICOMPZ = 0
      END IF
*
*     Check Argument Values
*
      INFO = 0
      WORK( 1 ) = MAX( 1, N )
      LQUERY = ( LWORK.EQ.-1 )
      IF( ISCHUR.EQ.0 ) THEN
         INFO = -1
      ELSE IF( ICOMPQ.EQ.0 ) THEN
         INFO = -2
      ELSE IF( ICOMPZ.EQ.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( ILO.LT.1 ) THEN
         INFO = -5
      ELSE IF( IHI.GT.N .OR. IHI.LT.ILO-1 ) THEN
         INFO = -6
      ELSE IF( LDH.LT.N ) THEN
         INFO = -8
      ELSE IF( LDT.LT.N ) THEN
         INFO = -10
      ELSE IF( LDQ.LT.1 .OR. ( ILQ .AND. LDQ.LT.N ) ) THEN
         INFO = -14
      ELSE IF( LDZ.LT.1 .OR. ( ILZ .AND. LDZ.LT.N ) ) THEN
         INFO = -16
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -18
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHGEQZ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
*     WORK( 1 ) = CMPLX( 1 )
      IF( N.LE.0 ) THEN
         WORK( 1 ) = DCMPLX( 1 )
         RETURN
      END IF
*
*     Initialize Q and Z
*
      IF( ICOMPQ.EQ.3 )
     $   CALL ZLASET( 'Full', N, N, CZERO, CONE, Q, LDQ )
      IF( ICOMPZ.EQ.3 )
     $   CALL ZLASET( 'Full', N, N, CZERO, CONE, Z, LDZ )
*
*     Machine Constants
*
      IN = IHI + 1 - ILO
      SAFMIN = DLAMCH( 'S' )
      ULP = DLAMCH( 'E' )*DLAMCH( 'B' )
      ANORM = ZLANHS( 'F', IN, H( ILO, ILO ), LDH, RWORK )
      BNORM = ZLANHS( 'F', IN, T( ILO, ILO ), LDT, RWORK )
      ATOL = MAX( SAFMIN, ULP*ANORM )
      BTOL = MAX( SAFMIN, ULP*BNORM )
      ASCALE = ONE / MAX( SAFMIN, ANORM )
      BSCALE = ONE / MAX( SAFMIN, BNORM )
*
*
*     Set Eigenvalues IHI+1:N
*
      DO 10 J = IHI + 1, N
         ABSB = ABS( T( J, J ) )
         IF( ABSB.GT.SAFMIN ) THEN
            SIGNBC = DCONJG( T( J, J ) / ABSB )
            T( J, J ) = ABSB
            IF( ILSCHR ) THEN
               CALL ZSCAL( J-1, SIGNBC, T( 1, J ), 1 )
               CALL ZSCAL( J, SIGNBC, H( 1, J ), 1 )
            ELSE
               CALL ZSCAL( 1, SIGNBC, H( J, J ), 1 )
            END IF
            IF( ILZ )
     $         CALL ZSCAL( N, SIGNBC, Z( 1, J ), 1 )
         ELSE
            T( J, J ) = CZERO
         END IF
         ALPHA( J ) = H( J, J )
         BETA( J ) = T( J, J )
   10 CONTINUE
*
*     If IHI < ILO, skip QZ steps
*
      IF( IHI.LT.ILO )
     $   GO TO 190
*
*     MAIN QZ ITERATION LOOP
*
*     Initialize dynamic indices
*
*     Eigenvalues ILAST+1:N have been found.
*        Column operations modify rows IFRSTM:whatever
*        Row operations modify columns whatever:ILASTM
*
*     If only eigenvalues are being computed, then
*        IFRSTM is the row of the last splitting row above row ILAST;
*        this is always at least ILO.
*     IITER counts iterations since the last eigenvalue was found,
*        to tell when to use an extraordinary shift.
*     MAXIT is the maximum number of QZ sweeps allowed.
*
      ILAST = IHI
      IF( ILSCHR ) THEN
         IFRSTM = 1
         ILASTM = N
      ELSE
         IFRSTM = ILO
         ILASTM = IHI
      END IF
      IITER = 0
      ESHIFT = CZERO
      MAXIT = 30*( IHI-ILO+1 )
*
      DO 170 JITER = 1, MAXIT
*
*        Check for too many iterations.
*
         IF( JITER.GT.MAXIT )
     $      GO TO 180
*
*        Split the matrix if possible.
*
*        Two tests:
*           1: H(j,j-1)=0  or  j=ILO
*           2: T(j,j)=0
*
*        Special case: j=ILAST
*
         IF( ILAST.EQ.ILO ) THEN
            GO TO 60
         ELSE
            IF( ABS1( H( ILAST, ILAST-1 ) ).LE.ATOL ) THEN
               H( ILAST, ILAST-1 ) = CZERO
               GO TO 60
            END IF
         END IF
*
         IF( ABS( T( ILAST, ILAST ) ).LE.BTOL ) THEN
            T( ILAST, ILAST ) = CZERO
            GO TO 50
         END IF
*
*        General case: j<ILAST
*
         DO 40 J = ILAST - 1, ILO, -1
*
*           Test 1: for H(j,j-1)=0 or j=ILO
*
            IF( J.EQ.ILO ) THEN
               ILAZRO = .TRUE.
            ELSE
               IF( ABS1( H( J, J-1 ) ).LE.ATOL ) THEN
                  H( J, J-1 ) = CZERO
                  ILAZRO = .TRUE.
               ELSE
                  ILAZRO = .FALSE.
               END IF
            END IF
*
*           Test 2: for T(j,j)=0
*
            IF( ABS( T( J, J ) ).LT.BTOL ) THEN
               T( J, J ) = CZERO
*
*              Test 1a: Check for 2 consecutive small subdiagonals in A
*
               ILAZR2 = .FALSE.
               IF( .NOT.ILAZRO ) THEN
                  IF( ABS1( H( J, J-1 ) )*( ASCALE*ABS1( H( J+1,
     $                J ) ) ).LE.ABS1( H( J, J ) )*( ASCALE*ATOL ) )
     $                ILAZR2 = .TRUE.
               END IF
*
*              If both tests pass (1 & 2), i.e., the leading diagonal
*              element of B in the block is zero, split a 1x1 block off
*              at the top. (I.e., at the J-th row/column) The leading
*              diagonal element of the remainder can also be zero, so
*              this may have to be done repeatedly.
*
               IF( ILAZRO .OR. ILAZR2 ) THEN
                  DO 20 JCH = J, ILAST - 1
                     CTEMP = H( JCH, JCH )
                     CALL ZLARTG( CTEMP, H( JCH+1, JCH ), C, S,
     $                            H( JCH, JCH ) )
                     H( JCH+1, JCH ) = CZERO
                     CALL ZROT( ILASTM-JCH, H( JCH, JCH+1 ), LDH,
     $                          H( JCH+1, JCH+1 ), LDH, C, S )
                     CALL ZROT( ILASTM-JCH, T( JCH, JCH+1 ), LDT,
     $                          T( JCH+1, JCH+1 ), LDT, C, S )
                     IF( ILQ )
     $                  CALL ZROT( N, Q( 1, JCH ), 1, Q( 1, JCH+1 ), 1,
     $                             C, DCONJG( S ) )
                     IF( ILAZR2 )
     $                  H( JCH, JCH-1 ) = H( JCH, JCH-1 )*C
                     ILAZR2 = .FALSE.
                     IF( ABS1( T( JCH+1, JCH+1 ) ).GE.BTOL ) THEN
                        IF( JCH+1.GE.ILAST ) THEN
                           GO TO 60
                        ELSE
                           IFIRST = JCH + 1
                           GO TO 70
                        END IF
                     END IF
                     T( JCH+1, JCH+1 ) = CZERO
   20             CONTINUE
                  GO TO 50
               ELSE
*
*                 Only test 2 passed -- chase the zero to T(ILAST,ILAST)
*                 Then process as in the case T(ILAST,ILAST)=0
*
                  DO 30 JCH = J, ILAST - 1
                     CTEMP = T( JCH, JCH+1 )
                     CALL ZLARTG( CTEMP, T( JCH+1, JCH+1 ), C, S,
     $                            T( JCH, JCH+1 ) )
                     T( JCH+1, JCH+1 ) = CZERO
                     IF( JCH.LT.ILASTM-1 )
     $                  CALL ZROT( ILASTM-JCH-1, T( JCH, JCH+2 ), LDT,
     $                             T( JCH+1, JCH+2 ), LDT, C, S )
                     CALL ZROT( ILASTM-JCH+2, H( JCH, JCH-1 ), LDH,
     $                          H( JCH+1, JCH-1 ), LDH, C, S )
                     IF( ILQ )
     $                  CALL ZROT( N, Q( 1, JCH ), 1, Q( 1, JCH+1 ), 1,
     $                             C, DCONJG( S ) )
                     CTEMP = H( JCH+1, JCH )
                     CALL ZLARTG( CTEMP, H( JCH+1, JCH-1 ), C, S,
     $                            H( JCH+1, JCH ) )
                     H( JCH+1, JCH-1 ) = CZERO
                     CALL ZROT( JCH+1-IFRSTM, H( IFRSTM, JCH ), 1,
     $                          H( IFRSTM, JCH-1 ), 1, C, S )
                     CALL ZROT( JCH-IFRSTM, T( IFRSTM, JCH ), 1,
     $                          T( IFRSTM, JCH-1 ), 1, C, S )
                     IF( ILZ )
     $                  CALL ZROT( N, Z( 1, JCH ), 1, Z( 1, JCH-1 ), 1,
     $                             C, S )
   30             CONTINUE
                  GO TO 50
               END IF
            ELSE IF( ILAZRO ) THEN
*
*              Only test 1 passed -- work on J:ILAST
*
               IFIRST = J
               GO TO 70
            END IF
*
*           Neither test passed -- try next J
*
   40    CONTINUE
*
*        (Drop-through is "impossible")
*
         INFO = 2*N + 1
         GO TO 210
*
*        T(ILAST,ILAST)=0 -- clear H(ILAST,ILAST-1) to split off a
*        1x1 block.
*
   50    CONTINUE
         CTEMP = H( ILAST, ILAST )
         CALL ZLARTG( CTEMP, H( ILAST, ILAST-1 ), C, S,
     $                H( ILAST, ILAST ) )
         H( ILAST, ILAST-1 ) = CZERO
         CALL ZROT( ILAST-IFRSTM, H( IFRSTM, ILAST ), 1,
     $              H( IFRSTM, ILAST-1 ), 1, C, S )
         CALL ZROT( ILAST-IFRSTM, T( IFRSTM, ILAST ), 1,
     $              T( IFRSTM, ILAST-1 ), 1, C, S )
         IF( ILZ )
     $      CALL ZROT( N, Z( 1, ILAST ), 1, Z( 1, ILAST-1 ), 1, C, S )
*
*        H(ILAST,ILAST-1)=0 -- Standardize B, set ALPHA and BETA
*
   60    CONTINUE
         ABSB = ABS( T( ILAST, ILAST ) )
         IF( ABSB.GT.SAFMIN ) THEN
            SIGNBC = DCONJG( T( ILAST, ILAST ) / ABSB )
            T( ILAST, ILAST ) = ABSB
            IF( ILSCHR ) THEN
               CALL ZSCAL( ILAST-IFRSTM, SIGNBC, T( IFRSTM, ILAST ), 1 )
               CALL ZSCAL( ILAST+1-IFRSTM, SIGNBC, H( IFRSTM, ILAST ),
     $                     1 )
            ELSE
               CALL ZSCAL( 1, SIGNBC, H( ILAST, ILAST ), 1 )
            END IF
            IF( ILZ )
     $         CALL ZSCAL( N, SIGNBC, Z( 1, ILAST ), 1 )
         ELSE
            T( ILAST, ILAST ) = CZERO
         END IF
         ALPHA( ILAST ) = H( ILAST, ILAST )
         BETA( ILAST ) = T( ILAST, ILAST )
*
*        Go to next block -- exit if finished.
*
         ILAST = ILAST - 1
         IF( ILAST.LT.ILO )
     $      GO TO 190
*
*        Reset counters
*
         IITER = 0
         ESHIFT = CZERO
         IF( .NOT.ILSCHR ) THEN
            ILASTM = ILAST
            IF( IFRSTM.GT.ILAST )
     $         IFRSTM = ILO
         END IF
         GO TO 160
*
*        QZ step
*
*        This iteration only involves rows/columns IFIRST:ILAST.  We
*        assume IFIRST < ILAST, and that the diagonal of B is non-zero.
*
   70    CONTINUE
         IITER = IITER + 1
         IF( .NOT.ILSCHR ) THEN
            IFRSTM = IFIRST
         END IF
*
*        Compute the Shift.
*
*        At this point, IFIRST < ILAST, and the diagonal elements of
*        T(IFIRST:ILAST,IFIRST,ILAST) are larger than BTOL (in
*        magnitude)
*
         IF( ( IITER / 10 )*10.NE.IITER ) THEN
*
*           The Wilkinson shift (AEP p.512), i.e., the eigenvalue of
*           the bottom-right 2x2 block of A inv(B) which is nearest to
*           the bottom-right element.
*
*           We factor B as U*D, where U has unit diagonals, and
*           compute (A*inv(D))*inv(U).
*
            U12 = ( BSCALE*T( ILAST-1, ILAST ) ) /
     $            ( BSCALE*T( ILAST, ILAST ) )
            AD11 = ( ASCALE*H( ILAST-1, ILAST-1 ) ) /
     $             ( BSCALE*T( ILAST-1, ILAST-1 ) )
            AD21 = ( ASCALE*H( ILAST, ILAST-1 ) ) /
     $             ( BSCALE*T( ILAST-1, ILAST-1 ) )
            AD12 = ( ASCALE*H( ILAST-1, ILAST ) ) /
     $             ( BSCALE*T( ILAST, ILAST ) )
            AD22 = ( ASCALE*H( ILAST, ILAST ) ) /
     $             ( BSCALE*T( ILAST, ILAST ) )
            ABI22 = AD22 - U12*AD21
*
            T1 = HALF*( AD11+ABI22 )
            RTDISC = SQRT( T1**2+AD12*AD21-AD11*AD22 )
            TEMP = DBLE( T1-ABI22 )*DBLE( RTDISC ) +
     $             DIMAG( T1-ABI22 )*DIMAG( RTDISC )
            IF( TEMP.LE.ZERO ) THEN
               SHIFT = T1 + RTDISC
            ELSE
               SHIFT = T1 - RTDISC
            END IF
         ELSE
*
*           Exceptional shift.  Chosen for no particularly good reason.
*
            ESHIFT = ESHIFT + (ASCALE*H(ILAST,ILAST-1))/
     $                        (BSCALE*T(ILAST-1,ILAST-1))
            SHIFT = ESHIFT
         END IF
*
*        Now check for two consecutive small subdiagonals.
*
         DO 80 J = ILAST - 1, IFIRST + 1, -1
            ISTART = J
            CTEMP = ASCALE*H( J, J ) - SHIFT*( BSCALE*T( J, J ) )
            TEMP = ABS1( CTEMP )
            TEMP2 = ASCALE*ABS1( H( J+1, J ) )
            TEMPR = MAX( TEMP, TEMP2 )
            IF( TEMPR.LT.ONE .AND. TEMPR.NE.ZERO ) THEN
               TEMP = TEMP / TEMPR
               TEMP2 = TEMP2 / TEMPR
            END IF
            IF( ABS1( H( J, J-1 ) )*TEMP2.LE.TEMP*ATOL )
     $         GO TO 90
   80    CONTINUE
*
         ISTART = IFIRST
         CTEMP = ASCALE*H( IFIRST, IFIRST ) -
     $           SHIFT*( BSCALE*T( IFIRST, IFIRST ) )
   90    CONTINUE
*
*        Do an implicit-shift QZ sweep.
*
*        Initial Q
*
         CTEMP2 = ASCALE*H( ISTART+1, ISTART )
         CALL ZLARTG( CTEMP, CTEMP2, C, S, CTEMP3 )
*
*        Sweep
*
         DO 150 J = ISTART, ILAST - 1
            IF( J.GT.ISTART ) THEN
               CTEMP = H( J, J-1 )
               CALL ZLARTG( CTEMP, H( J+1, J-1 ), C, S, H( J, J-1 ) )
               H( J+1, J-1 ) = CZERO
            END IF
*
            DO 100 JC = J, ILASTM
               CTEMP = C*H( J, JC ) + S*H( J+1, JC )
               H( J+1, JC ) = -DCONJG( S )*H( J, JC ) + C*H( J+1, JC )
               H( J, JC ) = CTEMP
               CTEMP2 = C*T( J, JC ) + S*T( J+1, JC )
               T( J+1, JC ) = -DCONJG( S )*T( J, JC ) + C*T( J+1, JC )
               T( J, JC ) = CTEMP2
  100       CONTINUE
            IF( ILQ ) THEN
               DO 110 JR = 1, N
                  CTEMP = C*Q( JR, J ) + DCONJG( S )*Q( JR, J+1 )
                  Q( JR, J+1 ) = -S*Q( JR, J ) + C*Q( JR, J+1 )
                  Q( JR, J ) = CTEMP
  110          CONTINUE
            END IF
*
            CTEMP = T( J+1, J+1 )
            CALL ZLARTG( CTEMP, T( J+1, J ), C, S, T( J+1, J+1 ) )
            T( J+1, J ) = CZERO
*
            DO 120 JR = IFRSTM, MIN( J+2, ILAST )
               CTEMP = C*H( JR, J+1 ) + S*H( JR, J )
               H( JR, J ) = -DCONJG( S )*H( JR, J+1 ) + C*H( JR, J )
               H( JR, J+1 ) = CTEMP
  120       CONTINUE
            DO 130 JR = IFRSTM, J
               CTEMP = C*T( JR, J+1 ) + S*T( JR, J )
               T( JR, J ) = -DCONJG( S )*T( JR, J+1 ) + C*T( JR, J )
               T( JR, J+1 ) = CTEMP
  130       CONTINUE
            IF( ILZ ) THEN
               DO 140 JR = 1, N
                  CTEMP = C*Z( JR, J+1 ) + S*Z( JR, J )
                  Z( JR, J ) = -DCONJG( S )*Z( JR, J+1 ) + C*Z( JR, J )
                  Z( JR, J+1 ) = CTEMP
  140          CONTINUE
            END IF
  150    CONTINUE
*
  160    CONTINUE
*
  170 CONTINUE
*
*     Drop-through = non-convergence
*
  180 CONTINUE
      INFO = ILAST
      GO TO 210
*
*     Successful completion of all QZ steps
*
  190 CONTINUE
*
*     Set Eigenvalues 1:ILO-1
*
      DO 200 J = 1, ILO - 1
         ABSB = ABS( T( J, J ) )
         IF( ABSB.GT.SAFMIN ) THEN
            SIGNBC = DCONJG( T( J, J ) / ABSB )
            T( J, J ) = ABSB
            IF( ILSCHR ) THEN
               CALL ZSCAL( J-1, SIGNBC, T( 1, J ), 1 )
               CALL ZSCAL( J, SIGNBC, H( 1, J ), 1 )
            ELSE
               CALL ZSCAL( 1, SIGNBC, H( J, J ), 1 )
            END IF
            IF( ILZ )
     $         CALL ZSCAL( N, SIGNBC, Z( 1, J ), 1 )
         ELSE
            T( J, J ) = CZERO
         END IF
         ALPHA( J ) = H( J, J )
         BETA( J ) = T( J, J )
  200 CONTINUE
*
*     Normal Termination
*
      INFO = 0
*
*     Exit (other than argument error) -- return optimal workspace size
*
  210 CONTINUE
      WORK( 1 ) = DCMPLX( N )
      RETURN
*
*     End of ZHGEQZ
*
      END
