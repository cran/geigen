*> \brief <b> ZGGSVD3 computes the singular value decomposition (SVD) for OTHER matrices</b>
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZGGSVD3 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zggsvd3.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zggsvd3.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zggsvd3.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGGSVD3( JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, B,
*                           LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK,
*                           LWORK, RWORK, IWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBQ, JOBU, JOBV
*       INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P, LWORK
*       ..
*       .. Array Arguments ..
*       INTEGER            IWORK( * )
*       DOUBLE PRECISION   ALPHA( * ), BETA( * ), RWORK( * )
*       COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
*      $                   U( LDU, * ), V( LDV, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGGSVD3 computes the generalized singular value decomposition (GSVD)
*> of an M-by-N complex matrix A and P-by-N complex matrix B:
*>
*>       U**H*A*Q = D1*( 0 R ),    V**H*B*Q = D2*( 0 R )
*>
*> where U, V and Q are unitary matrices.
*> Let K+L = the effective numerical rank of the
*> matrix (A**H,B**H)**H, then R is a (K+L)-by-(K+L) nonsingular upper
*> triangular matrix, D1 and D2 are M-by-(K+L) and P-by-(K+L) "diagonal"
*> matrices and of the following structures, respectively:
*>
*> If M-K-L >= 0,
*>
*>                     K  L
*>        D1 =     K ( I  0 )
*>                 L ( 0  C )
*>             M-K-L ( 0  0 )
*>
*>                   K  L
*>        D2 =   L ( 0  S )
*>             P-L ( 0  0 )
*>
*>                 N-K-L  K    L
*>   ( 0 R ) = K (  0   R11  R12 )
*>             L (  0    0   R22 )
*> where
*>
*>   C = diag( ALPHA(K+1), ... , ALPHA(K+L) ),
*>   S = diag( BETA(K+1),  ... , BETA(K+L) ),
*>   C**2 + S**2 = I.
*>
*>   R is stored in A(1:K+L,N-K-L+1:N) on exit.
*>
*> If M-K-L < 0,
*>
*>                   K M-K K+L-M
*>        D1 =   K ( I  0    0   )
*>             M-K ( 0  C    0   )
*>
*>                     K M-K K+L-M
*>        D2 =   M-K ( 0  S    0  )
*>             K+L-M ( 0  0    I  )
*>               P-L ( 0  0    0  )
*>
*>                    N-K-L  K   M-K  K+L-M
*>   ( 0 R ) =     K ( 0    R11  R12  R13  )
*>               M-K ( 0     0   R22  R23  )
*>             K+L-M ( 0     0    0   R33  )
*>
*> where
*>
*>   C = diag( ALPHA(K+1), ... , ALPHA(M) ),
*>   S = diag( BETA(K+1),  ... , BETA(M) ),
*>   C**2 + S**2 = I.
*>
*>   (R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N), and R33 is stored
*>   ( 0  R22 R23 )
*>   in B(M-K+1:L,N+M-K-L+1:N) on exit.
*>
*> The routine computes C, S, R, and optionally the unitary
*> transformation matrices U, V and Q.
*>
*> In particular, if B is an N-by-N nonsingular matrix, then the GSVD of
*> A and B implicitly gives the SVD of A*inv(B):
*>                      A*inv(B) = U*(D1*inv(D2))*V**H.
*> If ( A**H,B**H)**H has orthonormal columns, then the GSVD of A and B is also
*> equal to the CS decomposition of A and B. Furthermore, the GSVD can
*> be used to derive the solution of the eigenvalue problem:
*>                      A**H*A x = lambda* B**H*B x.
*> In some literature, the GSVD of A and B is presented in the form
*>                  U**H*A*X = ( 0 D1 ),   V**H*B*X = ( 0 D2 )
*> where U and V are orthogonal and X is nonsingular, and D1 and D2 are
*> ``diagonal''.  The former GSVD form can be converted to the latter
*> form by taking the nonsingular matrix X as
*>
*>                       X = Q*(  I   0    )
*>                             (  0 inv(R) )
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOBU
*> \verbatim
*>          JOBU is CHARACTER*1
*>          = 'U':  Unitary matrix U is computed;
*>          = 'N':  U is not computed.
*> \endverbatim
*>
*> \param[in] JOBV
*> \verbatim
*>          JOBV is CHARACTER*1
*>          = 'V':  Unitary matrix V is computed;
*>          = 'N':  V is not computed.
*> \endverbatim
*>
*> \param[in] JOBQ
*> \verbatim
*>          JOBQ is CHARACTER*1
*>          = 'Q':  Unitary matrix Q is computed;
*>          = 'N':  Q is not computed.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrices A and B.  N >= 0.
*> \endverbatim
*>
*> \param[in] P
*> \verbatim
*>          P is INTEGER
*>          The number of rows of the matrix B.  P >= 0.
*> \endverbatim
*>
*> \param[out] K
*> \verbatim
*>          K is INTEGER
*> \endverbatim
*>
*> \param[out] L
*> \verbatim
*>          L is INTEGER
*>
*>          On exit, K and L specify the dimension of the subblocks
*>          described in Purpose.
*>          K + L = effective numerical rank of (A**H,B**H)**H.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the M-by-N matrix A.
*>          On exit, A contains the triangular matrix R, or part of R.
*>          See Purpose for details.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A. LDA >= max(1,M).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX*16 array, dimension (LDB,N)
*>          On entry, the P-by-N matrix B.
*>          On exit, B contains part of the triangular matrix R if
*>          M-K-L < 0.  See Purpose for details.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B. LDB >= max(1,P).
*> \endverbatim
*>
*> \param[out] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION array, dimension (N)
*> \endverbatim
*>
*> \param[out] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION array, dimension (N)
*>
*>          On exit, ALPHA and BETA contain the generalized singular
*>          value pairs of A and B;
*>            ALPHA(1:K) = 1,
*>            BETA(1:K)  = 0,
*>          and if M-K-L >= 0,
*>            ALPHA(K+1:K+L) = C,
*>            BETA(K+1:K+L)  = S,
*>          or if M-K-L < 0,
*>            ALPHA(K+1:M)=C, ALPHA(M+1:K+L)=0
*>            BETA(K+1:M) =S, BETA(M+1:K+L) =1
*>          and
*>            ALPHA(K+L+1:N) = 0
*>            BETA(K+L+1:N)  = 0
*> \endverbatim
*>
*> \param[out] U
*> \verbatim
*>          U is COMPLEX*16 array, dimension (LDU,M)
*>          If JOBU = 'U', U contains the M-by-M unitary matrix U.
*>          If JOBU = 'N', U is not referenced.
*> \endverbatim
*>
*> \param[in] LDU
*> \verbatim
*>          LDU is INTEGER
*>          The leading dimension of the array U. LDU >= max(1,M) if
*>          JOBU = 'U'; LDU >= 1 otherwise.
*> \endverbatim
*>
*> \param[out] V
*> \verbatim
*>          V is COMPLEX*16 array, dimension (LDV,P)
*>          If JOBV = 'V', V contains the P-by-P unitary matrix V.
*>          If JOBV = 'N', V is not referenced.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>          The leading dimension of the array V. LDV >= max(1,P) if
*>          JOBV = 'V'; LDV >= 1 otherwise.
*> \endverbatim
*>
*> \param[out] Q
*> \verbatim
*>          Q is COMPLEX*16 array, dimension (LDQ,N)
*>          If JOBQ = 'Q', Q contains the N-by-N unitary matrix Q.
*>          If JOBQ = 'N', Q is not referenced.
*> \endverbatim
*>
*> \param[in] LDQ
*> \verbatim
*>          LDQ is INTEGER
*>          The leading dimension of the array Q. LDQ >= max(1,N) if
*>          JOBQ = 'Q'; LDQ >= 1 otherwise.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is DOUBLE PRECISION array, dimension (2*N)
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (N)
*>          On exit, IWORK stores the sorting information. More
*>          precisely, the following loop will sort ALPHA
*>             for I = K+1, min(M,K+L)
*>                 swap ALPHA(I) and ALPHA(IWORK(I))
*>             endfor
*>          such that ALPHA(1) >= ALPHA(2) >= ... >= ALPHA(N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          > 0:  if INFO = 1, the Jacobi-type procedure failed to
*>                converge.  For further details, see subroutine ZTGSJA.
*> \endverbatim
*
*> \par Internal Parameters:
*  =========================
*>
*> \verbatim
*>  TOLA    DOUBLE PRECISION
*>  TOLB    DOUBLE PRECISION
*>          TOLA and TOLB are the thresholds to determine the effective
*>          rank of (A**H,B**H)**H. Generally, they are set to
*>                   TOLA = MAX(M,N)*norm(A)*MACHEPS,
*>                   TOLB = MAX(P,N)*norm(B)*MACHEPS.
*>          The size of TOLA and TOLB may affect the size of backward
*>          errors of the decomposition.
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
*> \date August 2015
*
*> \ingroup complex16OTHERsing
*
*> \par Contributors:
*  ==================
*>
*>     Ming Gu and Huan Ren, Computer Science Division, University of
*>     California at Berkeley, USA
*>
*
*> \par Further Details:
*  =====================
*>
*>  ZGGSVD3 replaces the deprecated subroutine ZGGSVD.
*>
*  =====================================================================
      SUBROUTINE ZGGSVD3( JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, B,
     $                    LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ,
     $                    WORK, LWORK, RWORK, IWORK, INFO )
*
*  -- LAPACK driver routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     August 2015
*
*     .. Scalar Arguments ..
      CHARACTER          JOBQ, JOBU, JOBV
      INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P,
     $                   LWORK
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   ALPHA( * ), BETA( * ), RWORK( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
     $                   U( LDU, * ), V( LDV, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            WANTQ, WANTU, WANTV, LQUERY
      INTEGER            I, IBND, ISUB, J, NCYCLE, LWKOPT
      DOUBLE PRECISION   ANORM, BNORM, SMAX, TEMP, TOLA, TOLB, ULP, UNFL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, ZLANGE
      EXTERNAL           LSAME, DLAMCH, ZLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, XERBLA, ZGGSVP3, ZTGSJA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input parameters
*
      WANTU = LSAME( JOBU, 'U' )
      WANTV = LSAME( JOBV, 'V' )
      WANTQ = LSAME( JOBQ, 'Q' )
      LQUERY = ( LWORK.EQ.-1 )
      LWKOPT = 1
*
*     Test the input arguments
*
      INFO = 0
      IF( .NOT.( WANTU .OR. LSAME( JOBU, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( WANTV .OR. LSAME( JOBV, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( WANTQ .OR. LSAME( JOBQ, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( P.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LDB.LT.MAX( 1, P ) ) THEN
         INFO = -12
      ELSE IF( LDU.LT.1 .OR. ( WANTU .AND. LDU.LT.M ) ) THEN
         INFO = -16
      ELSE IF( LDV.LT.1 .OR. ( WANTV .AND. LDV.LT.P ) ) THEN
         INFO = -18
      ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.N ) ) THEN
         INFO = -20
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -24
      END IF
*
*     Compute workspace
*
      IF( INFO.EQ.0 ) THEN
         CALL ZGGSVP3( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, TOLA,
     $                 TOLB, K, L, U, LDU, V, LDV, Q, LDQ, IWORK, RWORK,
     $                 WORK, WORK, -1, INFO )
         LWKOPT = N + INT( WORK( 1 ) )
         LWKOPT = MAX( 2*N, LWKOPT )
         LWKOPT = MAX( 1, LWKOPT )
         WORK( 1 ) = DCMPLX( LWKOPT )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGGSVD3', -INFO )
         RETURN
      END IF
      IF( LQUERY ) THEN
         RETURN
      ENDIF
*
*     Compute the Frobenius norm of matrices A and B
*
      ANORM = ZLANGE( '1', M, N, A, LDA, RWORK )
      BNORM = ZLANGE( '1', P, N, B, LDB, RWORK )
*
*     Get machine precision and set up threshold for determining
*     the effective numerical rank of the matrices A and B.
*
      ULP = DLAMCH( 'Precision' )
      UNFL = DLAMCH( 'Safe Minimum' )
      TOLA = MAX( M, N )*MAX( ANORM, UNFL )*ULP
      TOLB = MAX( P, N )*MAX( BNORM, UNFL )*ULP
*
      CALL ZGGSVP3( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, TOLA,
     $              TOLB, K, L, U, LDU, V, LDV, Q, LDQ, IWORK, RWORK,
     $              WORK, WORK( N+1 ), LWORK-N, INFO )
*
*     Compute the GSVD of two upper "triangular" matrices
*
      CALL ZTGSJA( JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B, LDB,
     $             TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ,
     $             WORK, NCYCLE, INFO )
*
*     Sort the singular values and store the pivot indices in IWORK
*     Copy ALPHA to RWORK, then sort ALPHA in RWORK
*
      CALL DCOPY( N, ALPHA, 1, RWORK, 1 )
      IBND = MIN( L, M-K )
      DO 20 I = 1, IBND
*
*        Scan for largest ALPHA(K+I)
*
         ISUB = I
         SMAX = RWORK( K+I )
         DO 10 J = I + 1, IBND
            TEMP = RWORK( K+J )
            IF( TEMP.GT.SMAX ) THEN
               ISUB = J
               SMAX = TEMP
            END IF
   10    CONTINUE
         IF( ISUB.NE.I ) THEN
            RWORK( K+ISUB ) = RWORK( K+I )
            RWORK( K+I ) = SMAX
            IWORK( K+I ) = K + ISUB
         ELSE
            IWORK( K+I ) = K + I
         END IF
   20 CONTINUE
*
      WORK( 1 ) = DCMPLX( LWKOPT )
      RETURN
*
*     End of ZGGSVD3
*
      END
*> \brief \b ZGGSVP3
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZGGSVP3 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zggsvp3.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zggsvp3.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zggsvp3.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGGSVP3( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB,
*                           TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ,
*                           IWORK, RWORK, TAU, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBQ, JOBU, JOBV
*       INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P, LWORK
*       DOUBLE PRECISION   TOLA, TOLB
*       ..
*       .. Array Arguments ..
*       INTEGER            IWORK( * )
*       DOUBLE PRECISION   RWORK( * )
*       COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
*      $                   TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGGSVP3 computes unitary matrices U, V and Q such that
*>
*>                    N-K-L  K    L
*>  U**H*A*Q =     K ( 0    A12  A13 )  if M-K-L >= 0;
*>                 L ( 0     0   A23 )
*>             M-K-L ( 0     0    0  )
*>
*>                  N-K-L  K    L
*>         =     K ( 0    A12  A13 )  if M-K-L < 0;
*>             M-K ( 0     0   A23 )
*>
*>                  N-K-L  K    L
*>  V**H*B*Q =   L ( 0     0   B13 )
*>             P-L ( 0     0    0  )
*>
*> where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular
*> upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0,
*> otherwise A23 is (M-K)-by-L upper trapezoidal.  K+L = the effective
*> numerical rank of the (M+P)-by-N matrix (A**H,B**H)**H.
*>
*> This decomposition is the preprocessing step for computing the
*> Generalized Singular Value Decomposition (GSVD), see subroutine
*> ZGGSVD3.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOBU
*> \verbatim
*>          JOBU is CHARACTER*1
*>          = 'U':  Unitary matrix U is computed;
*>          = 'N':  U is not computed.
*> \endverbatim
*>
*> \param[in] JOBV
*> \verbatim
*>          JOBV is CHARACTER*1
*>          = 'V':  Unitary matrix V is computed;
*>          = 'N':  V is not computed.
*> \endverbatim
*>
*> \param[in] JOBQ
*> \verbatim
*>          JOBQ is CHARACTER*1
*>          = 'Q':  Unitary matrix Q is computed;
*>          = 'N':  Q is not computed.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] P
*> \verbatim
*>          P is INTEGER
*>          The number of rows of the matrix B.  P >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrices A and B.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the M-by-N matrix A.
*>          On exit, A contains the triangular (or trapezoidal) matrix
*>          described in the Purpose section.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A. LDA >= max(1,M).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX*16 array, dimension (LDB,N)
*>          On entry, the P-by-N matrix B.
*>          On exit, B contains the triangular matrix described in
*>          the Purpose section.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B. LDB >= max(1,P).
*> \endverbatim
*>
*> \param[in] TOLA
*> \verbatim
*>          TOLA is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] TOLB
*> \verbatim
*>          TOLB is DOUBLE PRECISION
*>
*>          TOLA and TOLB are the thresholds to determine the effective
*>          numerical rank of matrix B and a subblock of A. Generally,
*>          they are set to
*>             TOLA = MAX(M,N)*norm(A)*MAZHEPS,
*>             TOLB = MAX(P,N)*norm(B)*MAZHEPS.
*>          The size of TOLA and TOLB may affect the size of backward
*>          errors of the decomposition.
*> \endverbatim
*>
*> \param[out] K
*> \verbatim
*>          K is INTEGER
*> \endverbatim
*>
*> \param[out] L
*> \verbatim
*>          L is INTEGER
*>
*>          On exit, K and L specify the dimension of the subblocks
*>          described in Purpose section.
*>          K + L = effective numerical rank of (A**H,B**H)**H.
*> \endverbatim
*>
*> \param[out] U
*> \verbatim
*>          U is COMPLEX*16 array, dimension (LDU,M)
*>          If JOBU = 'U', U contains the unitary matrix U.
*>          If JOBU = 'N', U is not referenced.
*> \endverbatim
*>
*> \param[in] LDU
*> \verbatim
*>          LDU is INTEGER
*>          The leading dimension of the array U. LDU >= max(1,M) if
*>          JOBU = 'U'; LDU >= 1 otherwise.
*> \endverbatim
*>
*> \param[out] V
*> \verbatim
*>          V is COMPLEX*16 array, dimension (LDV,P)
*>          If JOBV = 'V', V contains the unitary matrix V.
*>          If JOBV = 'N', V is not referenced.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>          The leading dimension of the array V. LDV >= max(1,P) if
*>          JOBV = 'V'; LDV >= 1 otherwise.
*> \endverbatim
*>
*> \param[out] Q
*> \verbatim
*>          Q is COMPLEX*16 array, dimension (LDQ,N)
*>          If JOBQ = 'Q', Q contains the unitary matrix Q.
*>          If JOBQ = 'N', Q is not referenced.
*> \endverbatim
*>
*> \param[in] LDQ
*> \verbatim
*>          LDQ is INTEGER
*>          The leading dimension of the array Q. LDQ >= max(1,N) if
*>          JOBQ = 'Q'; LDQ >= 1 otherwise.
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (N)
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is DOUBLE PRECISION array, dimension (2*N)
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is COMPLEX*16 array, dimension (N)
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
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
*> \date August 2015
*
*> \ingroup complex16OTHERcomputational
*
*> \par Further Details:
*  =====================
*
*> \verbatim
*>
*>  The subroutine uses LAPACK subroutine ZGEQP3 for the QR factorization
*>  with column pivoting to detect the effective numerical rank of the
*>  a matrix. It may be replaced by a better rank determination strategy.
*>
*>  ZGGSVP3 replaces the deprecated subroutine ZGGSVP.
*>
*> \endverbatim
*>
* =====================================================================
      SUBROUTINE ZGGSVP3( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB,
     $                    TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ,
     $                    IWORK, RWORK, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     August 2015
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      CHARACTER          JOBQ, JOBU, JOBV
      INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P,
     $                   LWORK
      DOUBLE PRECISION   TOLA, TOLB
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
     $                   TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            FORWRD, WANTQ, WANTU, WANTV, LQUERY
      INTEGER            I, J, LWKOPT
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEQP3, ZGEQR2, ZGERQ2, ZLACPY, ZLAPMT,
     $                   ZLASET, ZUNG2R, ZUNM2R, ZUNMR2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      WANTU = LSAME( JOBU, 'U' )
      WANTV = LSAME( JOBV, 'V' )
      WANTQ = LSAME( JOBQ, 'Q' )
      FORWRD = .TRUE.
      LQUERY = ( LWORK.EQ.-1 )
      LWKOPT = 1
*
*     Test the input arguments
*
      INFO = 0
      IF( .NOT.( WANTU .OR. LSAME( JOBU, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( WANTV .OR. LSAME( JOBV, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( WANTQ .OR. LSAME( JOBQ, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( P.LT.0 ) THEN
         INFO = -5
      ELSE IF( N.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -8
      ELSE IF( LDB.LT.MAX( 1, P ) ) THEN
         INFO = -10
      ELSE IF( LDU.LT.1 .OR. ( WANTU .AND. LDU.LT.M ) ) THEN
         INFO = -16
      ELSE IF( LDV.LT.1 .OR. ( WANTV .AND. LDV.LT.P ) ) THEN
         INFO = -18
      ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.N ) ) THEN
         INFO = -20
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -24
      END IF
*
*     Compute workspace
*
      IF( INFO.EQ.0 ) THEN
         CALL ZGEQP3( P, N, B, LDB, IWORK, TAU, WORK, -1, RWORK, INFO )
         LWKOPT = INT( WORK ( 1 ) )
         IF( WANTV ) THEN
            LWKOPT = MAX( LWKOPT, P )
         END IF
         LWKOPT = MAX( LWKOPT, MIN( N, P ) )
         LWKOPT = MAX( LWKOPT, M )
         IF( WANTQ ) THEN
            LWKOPT = MAX( LWKOPT, N )
         END IF
         CALL ZGEQP3( M, N, A, LDA, IWORK, TAU, WORK, -1, RWORK, INFO )
         LWKOPT = MAX( LWKOPT, INT( WORK ( 1 ) ) )
         LWKOPT = MAX( 1, LWKOPT )
         WORK( 1 ) = DCMPLX( LWKOPT )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGGSVP3', -INFO )
         RETURN
      END IF
      IF( LQUERY ) THEN
         RETURN
      ENDIF
*
*     QR with column pivoting of B: B*P = V*( S11 S12 )
*                                           (  0   0  )
*
      DO 10 I = 1, N
         IWORK( I ) = 0
   10 CONTINUE
      CALL ZGEQP3( P, N, B, LDB, IWORK, TAU, WORK, LWORK, RWORK, INFO )
*
*     Update A := A*P
*
      CALL ZLAPMT( FORWRD, M, N, A, LDA, IWORK )
*
*     Determine the effective rank of matrix B.
*
      L = 0
      DO 20 I = 1, MIN( P, N )
         IF( ABS( B( I, I ) ).GT.TOLB )
     $      L = L + 1
   20 CONTINUE
*
      IF( WANTV ) THEN
*
*        Copy the details of V, and form V.
*
         CALL ZLASET( 'Full', P, P, CZERO, CZERO, V, LDV )
         IF( P.GT.1 )
     $      CALL ZLACPY( 'Lower', P-1, N, B( 2, 1 ), LDB, V( 2, 1 ),
     $                   LDV )
         CALL ZUNG2R( P, P, MIN( P, N ), V, LDV, TAU, WORK, INFO )
      END IF
*
*     Clean up B
*
      DO 40 J = 1, L - 1
         DO 30 I = J + 1, L
            B( I, J ) = CZERO
   30    CONTINUE
   40 CONTINUE
      IF( P.GT.L )
     $   CALL ZLASET( 'Full', P-L, N, CZERO, CZERO, B( L+1, 1 ), LDB )
*
      IF( WANTQ ) THEN
*
*        Set Q = I and Update Q := Q*P
*
         CALL ZLASET( 'Full', N, N, CZERO, CONE, Q, LDQ )
         CALL ZLAPMT( FORWRD, N, N, Q, LDQ, IWORK )
      END IF
*
      IF( P.GE.L .AND. N.NE.L ) THEN
*
*        RQ factorization of ( S11 S12 ) = ( 0 S12 )*Z
*
         CALL ZGERQ2( L, N, B, LDB, TAU, WORK, INFO )
*
*        Update A := A*Z**H
*
         CALL ZUNMR2( 'Right', 'Conjugate transpose', M, N, L, B, LDB,
     $                TAU, A, LDA, WORK, INFO )
         IF( WANTQ ) THEN
*
*           Update Q := Q*Z**H
*
            CALL ZUNMR2( 'Right', 'Conjugate transpose', N, N, L, B,
     $                   LDB, TAU, Q, LDQ, WORK, INFO )
         END IF
*
*        Clean up B
*
         CALL ZLASET( 'Full', L, N-L, CZERO, CZERO, B, LDB )
         DO 60 J = N - L + 1, N
            DO 50 I = J - N + L + 1, L
               B( I, J ) = CZERO
   50       CONTINUE
   60    CONTINUE
*
      END IF
*
*     Let              N-L     L
*                A = ( A11    A12 ) M,
*
*     then the following does the complete QR decomposition of A11:
*
*              A11 = U*(  0  T12 )*P1**H
*                      (  0   0  )
*
      DO 70 I = 1, N - L
         IWORK( I ) = 0
   70 CONTINUE
      CALL ZGEQP3( M, N-L, A, LDA, IWORK, TAU, WORK, LWORK, RWORK,
     $             INFO )
*
*     Determine the effective rank of A11
*
      K = 0
      DO 80 I = 1, MIN( M, N-L )
         IF( ABS( A( I, I ) ).GT.TOLA )
     $      K = K + 1
   80 CONTINUE
*
*     Update A12 := U**H*A12, where A12 = A( 1:M, N-L+1:N )
*
      CALL ZUNM2R( 'Left', 'Conjugate transpose', M, L, MIN( M, N-L ),
     $             A, LDA, TAU, A( 1, N-L+1 ), LDA, WORK, INFO )
*
      IF( WANTU ) THEN
*
*        Copy the details of U, and form U
*
         CALL ZLASET( 'Full', M, M, CZERO, CZERO, U, LDU )
         IF( M.GT.1 )
     $      CALL ZLACPY( 'Lower', M-1, N-L, A( 2, 1 ), LDA, U( 2, 1 ),
     $                   LDU )
         CALL ZUNG2R( M, M, MIN( M, N-L ), U, LDU, TAU, WORK, INFO )
      END IF
*
      IF( WANTQ ) THEN
*
*        Update Q( 1:N, 1:N-L )  = Q( 1:N, 1:N-L )*P1
*
         CALL ZLAPMT( FORWRD, N, N-L, Q, LDQ, IWORK )
      END IF
*
*     Clean up A: set the strictly lower triangular part of
*     A(1:K, 1:K) = 0, and A( K+1:M, 1:N-L ) = 0.
*
      DO 100 J = 1, K - 1
         DO 90 I = J + 1, K
            A( I, J ) = CZERO
   90    CONTINUE
  100 CONTINUE
      IF( M.GT.K )
     $   CALL ZLASET( 'Full', M-K, N-L, CZERO, CZERO, A( K+1, 1 ), LDA )
*
      IF( N-L.GT.K ) THEN
*
*        RQ factorization of ( T11 T12 ) = ( 0 T12 )*Z1
*
         CALL ZGERQ2( K, N-L, A, LDA, TAU, WORK, INFO )
*
         IF( WANTQ ) THEN
*
*           Update Q( 1:N,1:N-L ) = Q( 1:N,1:N-L )*Z1**H
*
            CALL ZUNMR2( 'Right', 'Conjugate transpose', N, N-L, K, A,
     $                   LDA, TAU, Q, LDQ, WORK, INFO )
         END IF
*
*        Clean up A
*
         CALL ZLASET( 'Full', K, N-L-K, CZERO, CZERO, A, LDA )
         DO 120 J = N - L - K + 1, N - L
            DO 110 I = J - N + L + K + 1, K
               A( I, J ) = CZERO
  110       CONTINUE
  120    CONTINUE
*
      END IF
*
      IF( M.GT.K ) THEN
*
*        QR factorization of A( K+1:M,N-L+1:N )
*
         CALL ZGEQR2( M-K, L, A( K+1, N-L+1 ), LDA, TAU, WORK, INFO )
*
         IF( WANTU ) THEN
*
*           Update U(:,K+1:M) := U(:,K+1:M)*U1
*
            CALL ZUNM2R( 'Right', 'No transpose', M, M-K, MIN( M-K, L ),
     $                   A( K+1, N-L+1 ), LDA, TAU, U( 1, K+1 ), LDU,
     $                   WORK, INFO )
         END IF
*
*        Clean up
*
         DO 140 J = N - L + 1, N
            DO 130 I = J - N + K + L + 1, M
               A( I, J ) = CZERO
  130       CONTINUE
  140    CONTINUE
*
      END IF
*
      WORK( 1 ) = DCMPLX( LWKOPT )
      RETURN
*
*     End of ZGGSVP3
*
      END
*> \brief \b ZTGSJA
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZTGSJA + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztgsja.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztgsja.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztgsja.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZTGSJA( JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B,
*                          LDB, TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV,
*                          Q, LDQ, WORK, NCYCLE, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          JOBQ, JOBU, JOBV
*       INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N,
*      $                   NCYCLE, P
*       DOUBLE PRECISION   TOLA, TOLB
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   ALPHA( * ), BETA( * )
*       COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
*      $                   U( LDU, * ), V( LDV, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZTGSJA computes the generalized singular value decomposition (GSVD)
*> of two complex upper triangular (or trapezoidal) matrices A and B.
*>
*> On entry, it is assumed that matrices A and B have the following
*> forms, which may be obtained by the preprocessing subroutine ZGGSVP
*> from a general M-by-N matrix A and P-by-N matrix B:
*>
*>              N-K-L  K    L
*>    A =    K ( 0    A12  A13 ) if M-K-L >= 0;
*>           L ( 0     0   A23 )
*>       M-K-L ( 0     0    0  )
*>
*>            N-K-L  K    L
*>    A =  K ( 0    A12  A13 ) if M-K-L < 0;
*>       M-K ( 0     0   A23 )
*>
*>            N-K-L  K    L
*>    B =  L ( 0     0   B13 )
*>       P-L ( 0     0    0  )
*>
*> where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular
*> upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0,
*> otherwise A23 is (M-K)-by-L upper trapezoidal.
*>
*> On exit,
*>
*>        U**H *A*Q = D1*( 0 R ),    V**H *B*Q = D2*( 0 R ),
*>
*> where U, V and Q are unitary matrices.
*> R is a nonsingular upper triangular matrix, and D1
*> and D2 are ``diagonal'' matrices, which are of the following
*> structures:
*>
*> If M-K-L >= 0,
*>
*>                     K  L
*>        D1 =     K ( I  0 )
*>                 L ( 0  C )
*>             M-K-L ( 0  0 )
*>
*>                    K  L
*>        D2 = L   ( 0  S )
*>             P-L ( 0  0 )
*>
*>                N-K-L  K    L
*>   ( 0 R ) = K (  0   R11  R12 ) K
*>             L (  0    0   R22 ) L
*>
*> where
*>
*>   C = diag( ALPHA(K+1), ... , ALPHA(K+L) ),
*>   S = diag( BETA(K+1),  ... , BETA(K+L) ),
*>   C**2 + S**2 = I.
*>
*>   R is stored in A(1:K+L,N-K-L+1:N) on exit.
*>
*> If M-K-L < 0,
*>
*>                K M-K K+L-M
*>     D1 =   K ( I  0    0   )
*>          M-K ( 0  C    0   )
*>
*>                  K M-K K+L-M
*>     D2 =   M-K ( 0  S    0   )
*>          K+L-M ( 0  0    I   )
*>            P-L ( 0  0    0   )
*>
*>                N-K-L  K   M-K  K+L-M
*> ( 0 R ) =    K ( 0    R11  R12  R13  )
*>           M-K ( 0     0   R22  R23  )
*>         K+L-M ( 0     0    0   R33  )
*>
*> where
*> C = diag( ALPHA(K+1), ... , ALPHA(M) ),
*> S = diag( BETA(K+1),  ... , BETA(M) ),
*> C**2 + S**2 = I.
*>
*> R = ( R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N) and R33 is stored
*>     (  0  R22 R23 )
*> in B(M-K+1:L,N+M-K-L+1:N) on exit.
*>
*> The computation of the unitary transformation matrices U, V or Q
*> is optional.  These matrices may either be formed explicitly, or they
*> may be postmultiplied into input matrices U1, V1, or Q1.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOBU
*> \verbatim
*>          JOBU is CHARACTER*1
*>          = 'U':  U must contain a unitary matrix U1 on entry, and
*>                  the product U1*U is returned;
*>          = 'I':  U is initialized to the unit matrix, and the
*>                  unitary matrix U is returned;
*>          = 'N':  U is not computed.
*> \endverbatim
*>
*> \param[in] JOBV
*> \verbatim
*>          JOBV is CHARACTER*1
*>          = 'V':  V must contain a unitary matrix V1 on entry, and
*>                  the product V1*V is returned;
*>          = 'I':  V is initialized to the unit matrix, and the
*>                  unitary matrix V is returned;
*>          = 'N':  V is not computed.
*> \endverbatim
*>
*> \param[in] JOBQ
*> \verbatim
*>          JOBQ is CHARACTER*1
*>          = 'Q':  Q must contain a unitary matrix Q1 on entry, and
*>                  the product Q1*Q is returned;
*>          = 'I':  Q is initialized to the unit matrix, and the
*>                  unitary matrix Q is returned;
*>          = 'N':  Q is not computed.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] P
*> \verbatim
*>          P is INTEGER
*>          The number of rows of the matrix B.  P >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrices A and B.  N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*> \endverbatim
*>
*> \param[in] L
*> \verbatim
*>          L is INTEGER
*>
*>          K and L specify the subblocks in the input matrices A and B:
*>          A23 = A(K+1:MIN(K+L,M),N-L+1:N) and B13 = B(1:L,,N-L+1:N)
*>          of A and B, whose GSVD is going to be computed by ZTGSJA.
*>          See Further Details.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the M-by-N matrix A.
*>          On exit, A(N-K+1:N,1:MIN(K+L,M) ) contains the triangular
*>          matrix R or part of R.  See Purpose for details.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A. LDA >= max(1,M).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX*16 array, dimension (LDB,N)
*>          On entry, the P-by-N matrix B.
*>          On exit, if necessary, B(M-K+1:L,N+M-K-L+1:N) contains
*>          a part of R.  See Purpose for details.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B. LDB >= max(1,P).
*> \endverbatim
*>
*> \param[in] TOLA
*> \verbatim
*>          TOLA is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] TOLB
*> \verbatim
*>          TOLB is DOUBLE PRECISION
*>
*>          TOLA and TOLB are the convergence criteria for the Jacobi-
*>          Kogbetliantz iteration procedure. Generally, they are the
*>          same as used in the preprocessing step, say
*>              TOLA = MAX(M,N)*norm(A)*MAZHEPS,
*>              TOLB = MAX(P,N)*norm(B)*MAZHEPS.
*> \endverbatim
*>
*> \param[out] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION array, dimension (N)
*> \endverbatim
*>
*> \param[out] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION array, dimension (N)
*>
*>          On exit, ALPHA and BETA contain the generalized singular
*>          value pairs of A and B;
*>            ALPHA(1:K) = 1,
*>            BETA(1:K)  = 0,
*>          and if M-K-L >= 0,
*>            ALPHA(K+1:K+L) = diag(C),
*>            BETA(K+1:K+L)  = diag(S),
*>          or if M-K-L < 0,
*>            ALPHA(K+1:M)= C, ALPHA(M+1:K+L)= 0
*>            BETA(K+1:M) = S, BETA(M+1:K+L) = 1.
*>          Furthermore, if K+L < N,
*>            ALPHA(K+L+1:N) = 0 and
*>            BETA(K+L+1:N)  = 0.
*> \endverbatim
*>
*> \param[in,out] U
*> \verbatim
*>          U is COMPLEX*16 array, dimension (LDU,M)
*>          On entry, if JOBU = 'U', U must contain a matrix U1 (usually
*>          the unitary matrix returned by ZGGSVP).
*>          On exit,
*>          if JOBU = 'I', U contains the unitary matrix U;
*>          if JOBU = 'U', U contains the product U1*U.
*>          If JOBU = 'N', U is not referenced.
*> \endverbatim
*>
*> \param[in] LDU
*> \verbatim
*>          LDU is INTEGER
*>          The leading dimension of the array U. LDU >= max(1,M) if
*>          JOBU = 'U'; LDU >= 1 otherwise.
*> \endverbatim
*>
*> \param[in,out] V
*> \verbatim
*>          V is COMPLEX*16 array, dimension (LDV,P)
*>          On entry, if JOBV = 'V', V must contain a matrix V1 (usually
*>          the unitary matrix returned by ZGGSVP).
*>          On exit,
*>          if JOBV = 'I', V contains the unitary matrix V;
*>          if JOBV = 'V', V contains the product V1*V.
*>          If JOBV = 'N', V is not referenced.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>          The leading dimension of the array V. LDV >= max(1,P) if
*>          JOBV = 'V'; LDV >= 1 otherwise.
*> \endverbatim
*>
*> \param[in,out] Q
*> \verbatim
*>          Q is COMPLEX*16 array, dimension (LDQ,N)
*>          On entry, if JOBQ = 'Q', Q must contain a matrix Q1 (usually
*>          the unitary matrix returned by ZGGSVP).
*>          On exit,
*>          if JOBQ = 'I', Q contains the unitary matrix Q;
*>          if JOBQ = 'Q', Q contains the product Q1*Q.
*>          If JOBQ = 'N', Q is not referenced.
*> \endverbatim
*>
*> \param[in] LDQ
*> \verbatim
*>          LDQ is INTEGER
*>          The leading dimension of the array Q. LDQ >= max(1,N) if
*>          JOBQ = 'Q'; LDQ >= 1 otherwise.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (2*N)
*> \endverbatim
*>
*> \param[out] NCYCLE
*> \verbatim
*>          NCYCLE is INTEGER
*>          The number of cycles required for convergence.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          = 1:  the procedure does not converge after MAXIT cycles.
*> \endverbatim
*
*> \par Internal Parameters:
*  =========================
*>
*> \verbatim
*>  MAXIT   INTEGER
*>          MAXIT specifies the total loops that the iterative procedure
*>          may take. If after MAXIT cycles, the routine fails to
*>          converge, we return INFO = 1.
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
*> \date November 2011
*
*> \ingroup complex16OTHERcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  ZTGSJA essentially uses a variant of Kogbetliantz algorithm to reduce
*>  min(L,M-K)-by-L triangular (or trapezoidal) matrix A23 and L-by-L
*>  matrix B13 to the form:
*>
*>           U1**H *A13*Q1 = C1*R1; V1**H *B13*Q1 = S1*R1,
*>
*>  where U1, V1 and Q1 are unitary matrix.
*>  C1 and S1 are diagonal matrices satisfying
*>
*>                C1**2 + S1**2 = I,
*>
*>  and R1 is an L-by-L nonsingular upper triangular matrix.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZTGSJA( JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B,
     $                   LDB, TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV,
     $                   Q, LDQ, WORK, NCYCLE, INFO )
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER          JOBQ, JOBU, JOBV
      INTEGER            INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N,
     $                   NCYCLE, P
      DOUBLE PRECISION   TOLA, TOLB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   ALPHA( * ), BETA( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
     $                   U( LDU, * ), V( LDV, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 40 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
*
      LOGICAL            INITQ, INITU, INITV, UPPER, WANTQ, WANTU, WANTV
      INTEGER            I, J, KCYCLE
      DOUBLE PRECISION   A1, A3, B1, B3, CSQ, CSU, CSV, ERROR, GAMMA,
     $                   RWK, SSMIN
      COMPLEX*16         A2, B2, SNQ, SNU, SNV
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARTG, XERBLA, ZCOPY, ZDSCAL, ZLAGS2, ZLAPLL,
     $                   ZLASET, ZROT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCONJG, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input parameters
*
      INITU = LSAME( JOBU, 'I' )
      WANTU = INITU .OR. LSAME( JOBU, 'U' )
*
      INITV = LSAME( JOBV, 'I' )
      WANTV = INITV .OR. LSAME( JOBV, 'V' )
*
      INITQ = LSAME( JOBQ, 'I' )
      WANTQ = INITQ .OR. LSAME( JOBQ, 'Q' )
*
      INFO = 0
      IF( .NOT.( INITU .OR. WANTU .OR. LSAME( JOBU, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( INITV .OR. WANTV .OR. LSAME( JOBV, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( INITQ .OR. WANTQ .OR. LSAME( JOBQ, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( P.LT.0 ) THEN
         INFO = -5
      ELSE IF( N.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LDB.LT.MAX( 1, P ) ) THEN
         INFO = -12
      ELSE IF( LDU.LT.1 .OR. ( WANTU .AND. LDU.LT.M ) ) THEN
         INFO = -18
      ELSE IF( LDV.LT.1 .OR. ( WANTV .AND. LDV.LT.P ) ) THEN
         INFO = -20
      ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.N ) ) THEN
         INFO = -22
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTGSJA', -INFO )
         RETURN
      END IF
*
*     Initialize U, V and Q, if necessary
*
      IF( INITU )
     $   CALL ZLASET( 'Full', M, M, CZERO, CONE, U, LDU )
      IF( INITV )
     $   CALL ZLASET( 'Full', P, P, CZERO, CONE, V, LDV )
      IF( INITQ )
     $   CALL ZLASET( 'Full', N, N, CZERO, CONE, Q, LDQ )
*
*     Loop until convergence
*
      UPPER = .FALSE.
      DO 40 KCYCLE = 1, MAXIT
*
         UPPER = .NOT.UPPER
*
         DO 20 I = 1, L - 1
            DO 10 J = I + 1, L
*
               A1 = ZERO
               A2 = CZERO
               A3 = ZERO
               IF( K+I.LE.M )
     $            A1 = DBLE( A( K+I, N-L+I ) )
               IF( K+J.LE.M )
     $            A3 = DBLE( A( K+J, N-L+J ) )
*
               B1 = DBLE( B( I, N-L+I ) )
               B3 = DBLE( B( J, N-L+J ) )
*
               IF( UPPER ) THEN
                  IF( K+I.LE.M )
     $               A2 = A( K+I, N-L+J )
                  B2 = B( I, N-L+J )
               ELSE
                  IF( K+J.LE.M )
     $               A2 = A( K+J, N-L+I )
                  B2 = B( J, N-L+I )
               END IF
*
               CALL ZLAGS2( UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU,
     $                      CSV, SNV, CSQ, SNQ )
*
*              Update (K+I)-th and (K+J)-th rows of matrix A: U**H *A
*
               IF( K+J.LE.M )
     $            CALL ZROT( L, A( K+J, N-L+1 ), LDA, A( K+I, N-L+1 ),
     $                       LDA, CSU, DCONJG( SNU ) )
*
*              Update I-th and J-th rows of matrix B: V**H *B
*
               CALL ZROT( L, B( J, N-L+1 ), LDB, B( I, N-L+1 ), LDB,
     $                    CSV, DCONJG( SNV ) )
*
*              Update (N-L+I)-th and (N-L+J)-th columns of matrices
*              A and B: A*Q and B*Q
*
               CALL ZROT( MIN( K+L, M ), A( 1, N-L+J ), 1,
     $                    A( 1, N-L+I ), 1, CSQ, SNQ )
*
               CALL ZROT( L, B( 1, N-L+J ), 1, B( 1, N-L+I ), 1, CSQ,
     $                    SNQ )
*
               IF( UPPER ) THEN
                  IF( K+I.LE.M )
     $               A( K+I, N-L+J ) = CZERO
                  B( I, N-L+J ) = CZERO
               ELSE
                  IF( K+J.LE.M )
     $               A( K+J, N-L+I ) = CZERO
                  B( J, N-L+I ) = CZERO
               END IF
*
*              Ensure that the diagonal elements of A and B are real.
*
               IF( K+I.LE.M )
     $            A( K+I, N-L+I ) = DBLE( A( K+I, N-L+I ) )
               IF( K+J.LE.M )
     $            A( K+J, N-L+J ) = DBLE( A( K+J, N-L+J ) )
               B( I, N-L+I ) = DBLE( B( I, N-L+I ) )
               B( J, N-L+J ) = DBLE( B( J, N-L+J ) )
*
*              Update unitary matrices U, V, Q, if desired.
*
               IF( WANTU .AND. K+J.LE.M )
     $            CALL ZROT( M, U( 1, K+J ), 1, U( 1, K+I ), 1, CSU,
     $                       SNU )
*
               IF( WANTV )
     $            CALL ZROT( P, V( 1, J ), 1, V( 1, I ), 1, CSV, SNV )
*
               IF( WANTQ )
     $            CALL ZROT( N, Q( 1, N-L+J ), 1, Q( 1, N-L+I ), 1, CSQ,
     $                       SNQ )
*
   10       CONTINUE
   20    CONTINUE
*
         IF( .NOT.UPPER ) THEN
*
*           The matrices A13 and B13 were lower triangular at the start
*           of the cycle, and are now upper triangular.
*
*           Convergence test: test the parallelism of the corresponding
*           rows of A and B.
*
            ERROR = ZERO
            DO 30 I = 1, MIN( L, M-K )
               CALL ZCOPY( L-I+1, A( K+I, N-L+I ), LDA, WORK, 1 )
               CALL ZCOPY( L-I+1, B( I, N-L+I ), LDB, WORK( L+1 ), 1 )
               CALL ZLAPLL( L-I+1, WORK, 1, WORK( L+1 ), 1, SSMIN )
               ERROR = MAX( ERROR, SSMIN )
   30       CONTINUE
*
            IF( ABS( ERROR ).LE.MIN( TOLA, TOLB ) )
     $         GO TO 50
         END IF
*
*        End of cycle loop
*
   40 CONTINUE
*
*     The algorithm has not converged after MAXIT cycles.
*
      INFO = 1
      GO TO 100
*
   50 CONTINUE
*
*     If ERROR <= MIN(TOLA,TOLB), then the algorithm has converged.
*     Compute the generalized singular value pairs (ALPHA, BETA), and
*     set the triangular matrix R to array A.
*
      DO 60 I = 1, K
         ALPHA( I ) = ONE
         BETA( I ) = ZERO
   60 CONTINUE
*
      DO 70 I = 1, MIN( L, M-K )
*
         A1 = DBLE( A( K+I, N-L+I ) )
         B1 = DBLE( B( I, N-L+I ) )
*
         IF( A1.NE.ZERO ) THEN
            GAMMA = B1 / A1
*
            IF( GAMMA.LT.ZERO ) THEN
               CALL ZDSCAL( L-I+1, -ONE, B( I, N-L+I ), LDB )
               IF( WANTV )
     $            CALL ZDSCAL( P, -ONE, V( 1, I ), 1 )
            END IF
*
            CALL DLARTG( ABS( GAMMA ), ONE, BETA( K+I ), ALPHA( K+I ),
     $                   RWK )
*
            IF( ALPHA( K+I ).GE.BETA( K+I ) ) THEN
               CALL ZDSCAL( L-I+1, ONE / ALPHA( K+I ), A( K+I, N-L+I ),
     $                      LDA )
            ELSE
               CALL ZDSCAL( L-I+1, ONE / BETA( K+I ), B( I, N-L+I ),
     $                      LDB )
               CALL ZCOPY( L-I+1, B( I, N-L+I ), LDB, A( K+I, N-L+I ),
     $                     LDA )
            END IF
*
         ELSE
*
            ALPHA( K+I ) = ZERO
            BETA( K+I ) = ONE
            CALL ZCOPY( L-I+1, B( I, N-L+I ), LDB, A( K+I, N-L+I ),
     $                  LDA )
         END IF
   70 CONTINUE
*
*     Post-assignment
*
      DO 80 I = M + 1, K + L
         ALPHA( I ) = ZERO
         BETA( I ) = ONE
   80 CONTINUE
*
      IF( K+L.LT.N ) THEN
         DO 90 I = K + L + 1, N
            ALPHA( I ) = ZERO
            BETA( I ) = ZERO
   90    CONTINUE
      END IF
*
  100 CONTINUE
      NCYCLE = KCYCLE
*
      RETURN
*
*     End of ZTGSJA
*
      END
*> \brief \b ZLAGS2
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLAGS2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlags2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlags2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlags2.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLAGS2( UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU, CSV,
*                          SNV, CSQ, SNQ )
*
*       .. Scalar Arguments ..
*       LOGICAL            UPPER
*       DOUBLE PRECISION   A1, A3, B1, B3, CSQ, CSU, CSV
*       COMPLEX*16         A2, B2, SNQ, SNU, SNV
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLAGS2 computes 2-by-2 unitary matrices U, V and Q, such
*> that if ( UPPER ) then
*>
*>           U**H *A*Q = U**H *( A1 A2 )*Q = ( x  0  )
*>                             ( 0  A3 )     ( x  x  )
*> and
*>           V**H*B*Q = V**H *( B1 B2 )*Q = ( x  0  )
*>                            ( 0  B3 )     ( x  x  )
*>
*> or if ( .NOT.UPPER ) then
*>
*>           U**H *A*Q = U**H *( A1 0  )*Q = ( x  x  )
*>                             ( A2 A3 )     ( 0  x  )
*> and
*>           V**H *B*Q = V**H *( B1 0  )*Q = ( x  x  )
*>                             ( B2 B3 )     ( 0  x  )
*> where
*>
*>   U = (   CSU    SNU ), V = (  CSV    SNV ),
*>       ( -SNU**H  CSU )      ( -SNV**H CSV )
*>
*>   Q = (   CSQ    SNQ )
*>       ( -SNQ**H  CSQ )
*>
*> The rows of the transformed A and B are parallel. Moreover, if the
*> input 2-by-2 matrix A is not zero, then the transformed (1,1) entry
*> of A is not zero. If the input matrices A and B are both not zero,
*> then the transformed (2,2) element of B is not zero, except when the
*> first rows of input A and B are parallel and the second rows are
*> zero.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPPER
*> \verbatim
*>          UPPER is LOGICAL
*>          = .TRUE.: the input matrices A and B are upper triangular.
*>          = .FALSE.: the input matrices A and B are lower triangular.
*> \endverbatim
*>
*> \param[in] A1
*> \verbatim
*>          A1 is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] A2
*> \verbatim
*>          A2 is COMPLEX*16
*> \endverbatim
*>
*> \param[in] A3
*> \verbatim
*>          A3 is DOUBLE PRECISION
*>          On entry, A1, A2 and A3 are elements of the input 2-by-2
*>          upper (lower) triangular matrix A.
*> \endverbatim
*>
*> \param[in] B1
*> \verbatim
*>          B1 is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[in] B2
*> \verbatim
*>          B2 is COMPLEX*16
*> \endverbatim
*>
*> \param[in] B3
*> \verbatim
*>          B3 is DOUBLE PRECISION
*>          On entry, B1, B2 and B3 are elements of the input 2-by-2
*>          upper (lower) triangular matrix B.
*> \endverbatim
*>
*> \param[out] CSU
*> \verbatim
*>          CSU is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[out] SNU
*> \verbatim
*>          SNU is COMPLEX*16
*>          The desired unitary matrix U.
*> \endverbatim
*>
*> \param[out] CSV
*> \verbatim
*>          CSV is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[out] SNV
*> \verbatim
*>          SNV is COMPLEX*16
*>          The desired unitary matrix V.
*> \endverbatim
*>
*> \param[out] CSQ
*> \verbatim
*>          CSQ is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[out] SNQ
*> \verbatim
*>          SNQ is COMPLEX*16
*>          The desired unitary matrix Q.
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
*> \date November 2011
*
*> \ingroup complex16OTHERauxiliary
*
*  =====================================================================
      SUBROUTINE ZLAGS2( UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU, CSV,
     $                   SNV, CSQ, SNQ )
*
*  -- LAPACK auxiliary routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      LOGICAL            UPPER
      DOUBLE PRECISION   A1, A3, B1, B3, CSQ, CSU, CSV
      COMPLEX*16         A2, B2, SNQ, SNU, SNV
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   A, AUA11, AUA12, AUA21, AUA22, AVB12, AVB11,
     $                   AVB21, AVB22, CSL, CSR, D, FB, FC, S1, S2,
     $                   SNL, SNR, UA11R, UA22R, VB11R, VB22R
      COMPLEX*16         B, C, D1, R, T, UA11, UA12, UA21, UA22, VB11,
     $                   VB12, VB21, VB22
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASV2, ZLARTG
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DCONJG, DIMAG
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   ABS1
*     ..
*     .. Statement Function definitions ..
      ABS1( T ) = ABS( DBLE( T ) ) + ABS( DIMAG( T ) )
*     ..
*     .. Executable Statements ..
*
      IF( UPPER ) THEN
*
*        Input matrices A and B are upper triangular matrices
*
*        Form matrix C = A*adj(B) = ( a b )
*                                   ( 0 d )
*
         A = A1*B3
         D = A3*B1
         B = A2*B1 - A1*B2
         FB = ABS( B )
*
*        Transform complex 2-by-2 matrix C to real matrix by unitary
*        diagonal matrix diag(1,D1).
*
         D1 = ONE
         IF( FB.NE.ZERO )
     $      D1 = B / FB
*
*        The SVD of real 2 by 2 triangular C
*
*         ( CSL -SNL )*( A B )*(  CSR  SNR ) = ( R 0 )
*         ( SNL  CSL ) ( 0 D ) ( -SNR  CSR )   ( 0 T )
*
         CALL DLASV2( A, FB, D, S1, S2, SNR, CSR, SNL, CSL )
*
         IF( ABS( CSL ).GE.ABS( SNL ) .OR. ABS( CSR ).GE.ABS( SNR ) )
     $        THEN
*
*           Compute the (1,1) and (1,2) elements of U**H *A and V**H *B,
*           and (1,2) element of |U|**H *|A| and |V|**H *|B|.
*
            UA11R = CSL*A1
            UA12 = CSL*A2 + D1*SNL*A3
*
            VB11R = CSR*B1
            VB12 = CSR*B2 + D1*SNR*B3
*
            AUA12 = ABS( CSL )*ABS1( A2 ) + ABS( SNL )*ABS( A3 )
            AVB12 = ABS( CSR )*ABS1( B2 ) + ABS( SNR )*ABS( B3 )
*
*           zero (1,2) elements of U**H *A and V**H *B
*
            IF( ( ABS( UA11R )+ABS1( UA12 ) ).EQ.ZERO ) THEN
               CALL ZLARTG( -DCMPLX( VB11R ), DCONJG( VB12 ), CSQ, SNQ,
     $                      R )
            ELSE IF( ( ABS( VB11R )+ABS1( VB12 ) ).EQ.ZERO ) THEN
               CALL ZLARTG( -DCMPLX( UA11R ), DCONJG( UA12 ), CSQ, SNQ,
     $                      R )
            ELSE IF( AUA12 / ( ABS( UA11R )+ABS1( UA12 ) ).LE.AVB12 /
     $               ( ABS( VB11R )+ABS1( VB12 ) ) ) THEN
               CALL ZLARTG( -DCMPLX( UA11R ), DCONJG( UA12 ), CSQ, SNQ,
     $                      R )
            ELSE
               CALL ZLARTG( -DCMPLX( VB11R ), DCONJG( VB12 ), CSQ, SNQ,
     $                      R )
            END IF
*
            CSU = CSL
            SNU = -D1*SNL
            CSV = CSR
            SNV = -D1*SNR
*
         ELSE
*
*           Compute the (2,1) and (2,2) elements of U**H *A and V**H *B,
*           and (2,2) element of |U|**H *|A| and |V|**H *|B|.
*
            UA21 = -DCONJG( D1 )*SNL*A1
            UA22 = -DCONJG( D1 )*SNL*A2 + CSL*A3
*
            VB21 = -DCONJG( D1 )*SNR*B1
            VB22 = -DCONJG( D1 )*SNR*B2 + CSR*B3
*
            AUA22 = ABS( SNL )*ABS1( A2 ) + ABS( CSL )*ABS( A3 )
            AVB22 = ABS( SNR )*ABS1( B2 ) + ABS( CSR )*ABS( B3 )
*
*           zero (2,2) elements of U**H *A and V**H *B, and then swap.
*
            IF( ( ABS1( UA21 )+ABS1( UA22 ) ).EQ.ZERO ) THEN
               CALL ZLARTG( -DCONJG( VB21 ), DCONJG( VB22 ), CSQ, SNQ,
     $                      R )
            ELSE IF( ( ABS1( VB21 )+ABS( VB22 ) ).EQ.ZERO ) THEN
               CALL ZLARTG( -DCONJG( UA21 ), DCONJG( UA22 ), CSQ, SNQ,
     $                      R )
            ELSE IF( AUA22 / ( ABS1( UA21 )+ABS1( UA22 ) ).LE.AVB22 /
     $               ( ABS1( VB21 )+ABS1( VB22 ) ) ) THEN
               CALL ZLARTG( -DCONJG( UA21 ), DCONJG( UA22 ), CSQ, SNQ,
     $                      R )
            ELSE
               CALL ZLARTG( -DCONJG( VB21 ), DCONJG( VB22 ), CSQ, SNQ,
     $                      R )
            END IF
*
            CSU = SNL
            SNU = D1*CSL
            CSV = SNR
            SNV = D1*CSR
*
         END IF
*
      ELSE
*
*        Input matrices A and B are lower triangular matrices
*
*        Form matrix C = A*adj(B) = ( a 0 )
*                                   ( c d )
*
         A = A1*B3
         D = A3*B1
         C = A2*B3 - A3*B2
         FC = ABS( C )
*
*        Transform complex 2-by-2 matrix C to real matrix by unitary
*        diagonal matrix diag(d1,1).
*
         D1 = ONE
         IF( FC.NE.ZERO )
     $      D1 = C / FC
*
*        The SVD of real 2 by 2 triangular C
*
*         ( CSL -SNL )*( A 0 )*(  CSR  SNR ) = ( R 0 )
*         ( SNL  CSL ) ( C D ) ( -SNR  CSR )   ( 0 T )
*
         CALL DLASV2( A, FC, D, S1, S2, SNR, CSR, SNL, CSL )
*
         IF( ABS( CSR ).GE.ABS( SNR ) .OR. ABS( CSL ).GE.ABS( SNL ) )
     $        THEN
*
*           Compute the (2,1) and (2,2) elements of U**H *A and V**H *B,
*           and (2,1) element of |U|**H *|A| and |V|**H *|B|.
*
            UA21 = -D1*SNR*A1 + CSR*A2
            UA22R = CSR*A3
*
            VB21 = -D1*SNL*B1 + CSL*B2
            VB22R = CSL*B3
*
            AUA21 = ABS( SNR )*ABS( A1 ) + ABS( CSR )*ABS1( A2 )
            AVB21 = ABS( SNL )*ABS( B1 ) + ABS( CSL )*ABS1( B2 )
*
*           zero (2,1) elements of U**H *A and V**H *B.
*
            IF( ( ABS1( UA21 )+ABS( UA22R ) ).EQ.ZERO ) THEN
               CALL ZLARTG( DCMPLX( VB22R ), VB21, CSQ, SNQ, R )
            ELSE IF( ( ABS1( VB21 )+ABS( VB22R ) ).EQ.ZERO ) THEN
               CALL ZLARTG( DCMPLX( UA22R ), UA21, CSQ, SNQ, R )
            ELSE IF( AUA21 / ( ABS1( UA21 )+ABS( UA22R ) ).LE.AVB21 /
     $               ( ABS1( VB21 )+ABS( VB22R ) ) ) THEN
               CALL ZLARTG( DCMPLX( UA22R ), UA21, CSQ, SNQ, R )
            ELSE
               CALL ZLARTG( DCMPLX( VB22R ), VB21, CSQ, SNQ, R )
            END IF
*
            CSU = CSR
            SNU = -DCONJG( D1 )*SNR
            CSV = CSL
            SNV = -DCONJG( D1 )*SNL
*
         ELSE
*
*           Compute the (1,1) and (1,2) elements of U**H *A and V**H *B,
*           and (1,1) element of |U|**H *|A| and |V|**H *|B|.
*
            UA11 = CSR*A1 + DCONJG( D1 )*SNR*A2
            UA12 = DCONJG( D1 )*SNR*A3
*
            VB11 = CSL*B1 + DCONJG( D1 )*SNL*B2
            VB12 = DCONJG( D1 )*SNL*B3
*
            AUA11 = ABS( CSR )*ABS( A1 ) + ABS( SNR )*ABS1( A2 )
            AVB11 = ABS( CSL )*ABS( B1 ) + ABS( SNL )*ABS1( B2 )
*
*           zero (1,1) elements of U**H *A and V**H *B, and then swap.
*
            IF( ( ABS1( UA11 )+ABS1( UA12 ) ).EQ.ZERO ) THEN
               CALL ZLARTG( VB12, VB11, CSQ, SNQ, R )
            ELSE IF( ( ABS1( VB11 )+ABS1( VB12 ) ).EQ.ZERO ) THEN
               CALL ZLARTG( UA12, UA11, CSQ, SNQ, R )
            ELSE IF( AUA11 / ( ABS1( UA11 )+ABS1( UA12 ) ).LE.AVB11 /
     $               ( ABS1( VB11 )+ABS1( VB12 ) ) ) THEN
               CALL ZLARTG( UA12, UA11, CSQ, SNQ, R )
            ELSE
               CALL ZLARTG( VB12, VB11, CSQ, SNQ, R )
            END IF
*
            CSU = SNR
            SNU = DCONJG( D1 )*CSR
            CSV = SNL
            SNV = DCONJG( D1 )*CSL
*
         END IF
*
      END IF
*
      RETURN
*
*     End of ZLAGS2
*
      END
*> \brief \b ZLAPLL measures the linear dependence of two vectors.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLAPLL + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlapll.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlapll.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlapll.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLAPLL( N, X, INCX, Y, INCY, SSMIN )
*
*       .. Scalar Arguments ..
*       INTEGER            INCX, INCY, N
*       DOUBLE PRECISION   SSMIN
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         X( * ), Y( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> Given two column vectors X and Y, let
*>
*>                      A = ( X Y ).
*>
*> The subroutine first computes the QR factorization of A = Q*R,
*> and then computes the SVD of the 2-by-2 upper triangular matrix R.
*> The smaller singular value of R is returned in SSMIN, which is used
*> as the measurement of the linear dependency of the vectors X and Y.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The length of the vectors X and Y.
*> \endverbatim
*>
*> \param[in,out] X
*> \verbatim
*>          X is COMPLEX*16 array, dimension (1+(N-1)*INCX)
*>          On entry, X contains the N-vector X.
*>          On exit, X is overwritten.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>          The increment between successive elements of X. INCX > 0.
*> \endverbatim
*>
*> \param[in,out] Y
*> \verbatim
*>          Y is COMPLEX*16 array, dimension (1+(N-1)*INCY)
*>          On entry, Y contains the N-vector Y.
*>          On exit, Y is overwritten.
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>          The increment between successive elements of Y. INCY > 0.
*> \endverbatim
*>
*> \param[out] SSMIN
*> \verbatim
*>          SSMIN is DOUBLE PRECISION
*>          The smallest singular value of the N-by-2 matrix A = ( X Y ).
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
*> \date September 2012
*
*> \ingroup complex16OTHERauxiliary
*
*  =====================================================================
      SUBROUTINE ZLAPLL( N, X, INCX, Y, INCY, SSMIN )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
      DOUBLE PRECISION   SSMIN
*     ..
*     .. Array Arguments ..
      COMPLEX*16         X( * ), Y( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   SSMAX
      COMPLEX*16         A11, A12, A22, C, TAU
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DCONJG
*     ..
*     .. External Functions ..
      COMPLEX*16         ZDOTC
      EXTERNAL           ZDOTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAS2, ZAXPY, ZLARFG
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.1 ) THEN
         SSMIN = ZERO
         RETURN
      END IF
*
*     Compute the QR factorization of the N-by-2 matrix ( X Y )
*
      CALL ZLARFG( N, X( 1 ), X( 1+INCX ), INCX, TAU )
      A11 = X( 1 )
      X( 1 ) = CONE
*
      C = -DCONJG( TAU )*ZDOTC( N, X, INCX, Y, INCY )
      CALL ZAXPY( N, C, X, INCX, Y, INCY )
*
      CALL ZLARFG( N-1, Y( 1+INCY ), Y( 1+2*INCY ), INCY, TAU )
*
      A12 = Y( 1 )
      A22 = Y( 1+INCY )
*
*     Compute the SVD of 2-by-2 Upper triangular matrix.
*
      CALL DLAS2( ABS( A11 ), ABS( A12 ), ABS( A22 ), SSMIN, SSMAX )
*
      RETURN
*
*     End of ZLAPLL
*
      END
*> \brief \b ZGERQ2 computes the RQ factorization of a general rectangular matrix using an unblocked algorithm.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZGERQ2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgerq2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgerq2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgerq2.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGERQ2( M, N, A, LDA, TAU, WORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, M, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGERQ2 computes an RQ factorization of a complex m by n matrix A:
*> A = R * Q.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the m by n matrix A.
*>          On exit, if m <= n, the upper triangle of the subarray
*>          A(1:m,n-m+1:n) contains the m by m upper triangular matrix R;
*>          if m >= n, the elements on and above the (m-n)-th subdiagonal
*>          contain the m by n upper trapezoidal matrix R; the remaining
*>          elements, with the array TAU, represent the unitary matrix
*>          Q as a product of elementary reflectors (see Further
*>          Details).
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is COMPLEX*16 array, dimension (min(M,N))
*>          The scalar factors of the elementary reflectors (see Further
*>          Details).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (M)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
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
*> \date September 2012
*
*> \ingroup complex16GEcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The matrix Q is represented as a product of elementary reflectors
*>
*>     Q = H(1)**H H(2)**H . . . H(k)**H, where k = min(m,n).
*>
*>  Each H(i) has the form
*>
*>     H(i) = I - tau * v * v**H
*>
*>  where tau is a complex scalar, and v is a complex vector with
*>  v(n-k+i+1:n) = 0 and v(n-k+i) = 1; conjg(v(1:n-k+i-1)) is stored on
*>  exit in A(m-k+i,1:n-k+i-1), and tau in TAU(i).
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZGERQ2( M, N, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK computational routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, K
      COMPLEX*16         ALPHA
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLACGV, ZLARF, ZLARFG
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGERQ2', -INFO )
         RETURN
      END IF
*
      K = MIN( M, N )
*
      DO 10 I = K, 1, -1
*
*        Generate elementary reflector H(i) to annihilate
*        A(m-k+i,1:n-k+i-1)
*
         CALL ZLACGV( N-K+I, A( M-K+I, 1 ), LDA )
         ALPHA = A( M-K+I, N-K+I )
         CALL ZLARFG( N-K+I, ALPHA, A( M-K+I, 1 ), LDA, TAU( I ) )
*
*        Apply H(i) to A(1:m-k+i-1,1:n-k+i) from the right
*
         A( M-K+I, N-K+I ) = ONE
         CALL ZLARF( 'Right', M-K+I-1, N-K+I, A( M-K+I, 1 ), LDA,
     $               TAU( I ), A, LDA, WORK )
         A( M-K+I, N-K+I ) = ALPHA
         CALL ZLACGV( N-K+I-1, A( M-K+I, 1 ), LDA )
   10 CONTINUE
      RETURN
*
*     End of ZGERQ2
*
      END
*> \brief \b ZLAPMT performs a forward or backward permutation of the columns of a matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLAPMT + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlapmt.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlapmt.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlapmt.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLAPMT( FORWRD, M, N, X, LDX, K )
*
*       .. Scalar Arguments ..
*       LOGICAL            FORWRD
*       INTEGER            LDX, M, N
*       ..
*       .. Array Arguments ..
*       INTEGER            K( * )
*       COMPLEX*16         X( LDX, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLAPMT rearranges the columns of the M by N matrix X as specified
*> by the permutation K(1),K(2),...,K(N) of the integers 1,...,N.
*> If FORWRD = .TRUE.,  forward permutation:
*>
*>      X(*,K(J)) is moved X(*,J) for J = 1,2,...,N.
*>
*> If FORWRD = .FALSE., backward permutation:
*>
*>      X(*,J) is moved to X(*,K(J)) for J = 1,2,...,N.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] FORWRD
*> \verbatim
*>          FORWRD is LOGICAL
*>          = .TRUE., forward permutation
*>          = .FALSE., backward permutation
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix X. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix X. N >= 0.
*> \endverbatim
*>
*> \param[in,out] X
*> \verbatim
*>          X is COMPLEX*16 array, dimension (LDX,N)
*>          On entry, the M by N matrix X.
*>          On exit, X contains the permuted matrix X.
*> \endverbatim
*>
*> \param[in] LDX
*> \verbatim
*>          LDX is INTEGER
*>          The leading dimension of the array X, LDX >= MAX(1,M).
*> \endverbatim
*>
*> \param[in,out] K
*> \verbatim
*>          K is INTEGER array, dimension (N)
*>          On entry, K contains the permutation vector. K is used as
*>          internal workspace, but reset to its original value on
*>          output.
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
*> \date September 2012
*
*> \ingroup complex16OTHERauxiliary
*
*  =====================================================================
      SUBROUTINE ZLAPMT( FORWRD, M, N, X, LDX, K )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      LOGICAL            FORWRD
      INTEGER            LDX, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            K( * )
      COMPLEX*16         X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, II, IN, J
      COMPLEX*16         TEMP
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.1 )
     $   RETURN
*
      DO 10 I = 1, N
         K( I ) = -K( I )
   10 CONTINUE
*
      IF( FORWRD ) THEN
*
*        Forward permutation
*
         DO 50 I = 1, N
*
            IF( K( I ).GT.0 )
     $         GO TO 40
*
            J = I
            K( J ) = -K( J )
            IN = K( J )
*
   20       CONTINUE
            IF( K( IN ).GT.0 )
     $         GO TO 40
*
            DO 30 II = 1, M
               TEMP = X( II, J )
               X( II, J ) = X( II, IN )
               X( II, IN ) = TEMP
   30       CONTINUE
*
            K( IN ) = -K( IN )
            J = IN
            IN = K( IN )
            GO TO 20
*
   40       CONTINUE
*
   50    CONTINUE
*
      ELSE
*
*        Backward permutation
*
         DO 90 I = 1, N
*
            IF( K( I ).GT.0 )
     $         GO TO 80
*
            K( I ) = -K( I )
            J = K( I )
   60       CONTINUE
            IF( J.EQ.I )
     $         GO TO 80
*
            DO 70 II = 1, M
               TEMP = X( II, I )
               X( II, I ) = X( II, J )
               X( II, J ) = TEMP
   70       CONTINUE
*
            K( J ) = -K( J )
            J = K( J )
            GO TO 60
*
   80       CONTINUE
*
   90    CONTINUE
*
      END IF
*
      RETURN
*
*     End of ZLAPMT
*
      END
*> \brief \b ZUNMR2 multiplies a general matrix by the unitary matrix from a RQ factorization determined by cgerqf (unblocked algorithm).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZUNMR2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunmr2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunmr2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunmr2.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZUNMR2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
*                          WORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          SIDE, TRANS
*       INTEGER            INFO, K, LDA, LDC, M, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZUNMR2 overwrites the general complex m-by-n matrix C with
*>
*>       Q * C  if SIDE = 'L' and TRANS = 'N', or
*>
*>       Q**H* C  if SIDE = 'L' and TRANS = 'C', or
*>
*>       C * Q  if SIDE = 'R' and TRANS = 'N', or
*>
*>       C * Q**H if SIDE = 'R' and TRANS = 'C',
*>
*> where Q is a complex unitary matrix defined as the product of k
*> elementary reflectors
*>
*>       Q = H(1)**H H(2)**H . . . H(k)**H
*>
*> as returned by ZGERQF. Q is of order m if SIDE = 'L' and of order n
*> if SIDE = 'R'.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': apply Q or Q**H from the Left
*>          = 'R': apply Q or Q**H from the Right
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          = 'N': apply Q  (No transpose)
*>          = 'C': apply Q**H (Conjugate transpose)
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix C. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix C. N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The number of elementary reflectors whose product defines
*>          the matrix Q.
*>          If SIDE = 'L', M >= K >= 0;
*>          if SIDE = 'R', N >= K >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension
*>                               (LDA,M) if SIDE = 'L',
*>                               (LDA,N) if SIDE = 'R'
*>          The i-th row must contain the vector which defines the
*>          elementary reflector H(i), for i = 1,2,...,k, as returned by
*>          ZGERQF in the last k rows of its array argument A.
*>          A is modified by the routine but restored on exit.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A. LDA >= max(1,K).
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is COMPLEX*16 array, dimension (K)
*>          TAU(i) must contain the scalar factor of the elementary
*>          reflector H(i), as returned by ZGERQF.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is COMPLEX*16 array, dimension (LDC,N)
*>          On entry, the m-by-n matrix C.
*>          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C. LDC >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension
*>                                   (N) if SIDE = 'L',
*>                                   (M) if SIDE = 'R'
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
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
*> \date September 2012
*
*> \ingroup complex16OTHERcomputational
*
*  =====================================================================
      SUBROUTINE ZUNMR2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, INFO )
*
*  -- LAPACK computational routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, MI, NI, NQ
      COMPLEX*16         AII, TAUI
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLACGV, ZLARF
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
*
*     NQ is the order of Q
*
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, K ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNMR2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 )
     $   RETURN
*
      IF( ( LEFT .AND. .NOT.NOTRAN .OR. .NOT.LEFT .AND. NOTRAN ) ) THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
*
      IF( LEFT ) THEN
         NI = N
      ELSE
         MI = M
      END IF
*
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
*
*           H(i) or H(i)**H is applied to C(1:m-k+i,1:n)
*
            MI = M - K + I
         ELSE
*
*           H(i) or H(i)**H is applied to C(1:m,1:n-k+i)
*
            NI = N - K + I
         END IF
*
*        Apply H(i) or H(i)**H
*
         IF( NOTRAN ) THEN
            TAUI = DCONJG( TAU( I ) )
         ELSE
            TAUI = TAU( I )
         END IF
         CALL ZLACGV( NQ-K+I-1, A( I, 1 ), LDA )
         AII = A( I, NQ-K+I )
         A( I, NQ-K+I ) = ONE
         CALL ZLARF( SIDE, MI, NI, A( I, 1 ), LDA, TAUI, C, LDC, WORK )
         A( I, NQ-K+I ) = AII
         CALL ZLACGV( NQ-K+I-1, A( I, 1 ), LDA )
   10 CONTINUE
      RETURN
*
*     End of ZUNMR2
*
      END
