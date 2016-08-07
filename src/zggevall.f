*> \brief <b> ZGGEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices</b>
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZGGEV + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zggev.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zggev.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zggev.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA,
*                         VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          JOBVL, JOBVR
*       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   RWORK( * )
*       COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),
*      $                   BETA( * ), VL( LDVL, * ), VR( LDVR, * ),
*      $                   WORK( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGGEV computes for a pair of N-by-N complex nonsymmetric matrices
*> (A,B), the generalized eigenvalues, and optionally, the left and/or
*> right generalized eigenvectors.
*>
*> A generalized eigenvalue for a pair of matrices (A,B) is a scalar
*> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
*> singular. It is usually represented as the pair (alpha,beta), as
*> there is a reasonable interpretation for beta=0, and even for both
*> being zero.
*>
*> The right generalized eigenvector v(j) corresponding to the
*> generalized eigenvalue lambda(j) of (A,B) satisfies
*>
*>              A * v(j) = lambda(j) * B * v(j).
*>
*> The left generalized eigenvector u(j) corresponding to the
*> generalized eigenvalues lambda(j) of (A,B) satisfies
*>
*>              u(j)**H * A = lambda(j) * u(j)**H * B
*>
*> where u(j)**H is the conjugate-transpose of u(j).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOBVL
*> \verbatim
*>          JOBVL is CHARACTER*1
*>          = 'N':  do not compute the left generalized eigenvectors;
*>          = 'V':  compute the left generalized eigenvectors.
*> \endverbatim
*>
*> \param[in] JOBVR
*> \verbatim
*>          JOBVR is CHARACTER*1
*>          = 'N':  do not compute the right generalized eigenvectors;
*>          = 'V':  compute the right generalized eigenvectors.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices A, B, VL, and VR.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA, N)
*>          On entry, the matrix A in the pair (A,B).
*>          On exit, A has been overwritten.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX*16 array, dimension (LDB, N)
*>          On entry, the matrix B in the pair (A,B).
*>          On exit, B has been overwritten.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX*16 array, dimension (N)
*> \endverbatim
*>
*> \param[out] BETA
*> \verbatim
*>          BETA is COMPLEX*16 array, dimension (N)
*>          On exit, ALPHA(j)/BETA(j), j=1,...,N, will be the
*>          generalized eigenvalues.
*>
*>          Note: the quotients ALPHA(j)/BETA(j) may easily over- or
*>          underflow, and BETA(j) may even be zero.  Thus, the user
*>          should avoid naively computing the ratio alpha/beta.
*>          However, ALPHA will be always less than and usually
*>          comparable with norm(A) in magnitude, and BETA always less
*>          than and usually comparable with norm(B).
*> \endverbatim
*>
*> \param[out] VL
*> \verbatim
*>          VL is COMPLEX*16 array, dimension (LDVL,N)
*>          If JOBVL = 'V', the left generalized eigenvectors u(j) are
*>          stored one after another in the columns of VL, in the same
*>          order as their eigenvalues.
*>          Each eigenvector is scaled so the largest component has
*>          abs(real part) + abs(imag. part) = 1.
*>          Not referenced if JOBVL = 'N'.
*> \endverbatim
*>
*> \param[in] LDVL
*> \verbatim
*>          LDVL is INTEGER
*>          The leading dimension of the matrix VL. LDVL >= 1, and
*>          if JOBVL = 'V', LDVL >= N.
*> \endverbatim
*>
*> \param[out] VR
*> \verbatim
*>          VR is COMPLEX*16 array, dimension (LDVR,N)
*>          If JOBVR = 'V', the right generalized eigenvectors v(j) are
*>          stored one after another in the columns of VR, in the same
*>          order as their eigenvalues.
*>          Each eigenvector is scaled so the largest component has
*>          abs(real part) + abs(imag. part) = 1.
*>          Not referenced if JOBVR = 'N'.
*> \endverbatim
*>
*> \param[in] LDVR
*> \verbatim
*>          LDVR is INTEGER
*>          The leading dimension of the matrix VR. LDVR >= 1, and
*>          if JOBVR = 'V', LDVR >= N.
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
*>          The dimension of the array WORK.  LWORK >= max(1,2*N).
*>          For good performance, LWORK must generally be larger.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is DOUBLE PRECISION array, dimension (8*N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          =1,...,N:
*>                The QZ iteration failed.  No eigenvectors have been
*>                calculated, but ALPHA(j) and BETA(j) should be
*>                correct for j=INFO+1,...,N.
*>          > N:  =N+1: other then QZ iteration failed in DHGEQZ,
*>                =N+2: error return from DTGEVC.
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
*> \ingroup complex16GEeigen
*
*  =====================================================================
      SUBROUTINE ZGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA,
     $                  VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )
*
*  -- LAPACK driver routine (version 3.4.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     April 2012
*
*     .. Scalar Arguments ..
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),
     $                   BETA( * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D0, 0.0D0 ),
     $                   CONE = ( 1.0D0, 0.0D0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILASCL, ILBSCL, ILV, ILVL, ILVR, LQUERY
      CHARACTER          CHTEMP
      INTEGER            ICOLS, IERR, IHI, IJOBVL, IJOBVR, ILEFT, ILO,
     $                   IN, IRIGHT, IROWS, IRWRK, ITAU, IWRK, JC, JR,
     $                   LWKMIN, LWKOPT
      DOUBLE PRECISION   ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS,
     $                   SMLNUM, TEMP
      COMPLEX*16         X
*     ..
*     .. Local Arrays ..
      LOGICAL            LDUMMA( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLABAD, XERBLA, ZGEQRF, ZGGBAK, ZGGBAL, ZGGHRD,
     $                   ZHGEQZ, ZLACPY, ZLASCL, ZLASET, ZTGEVC, ZUNGQR,
     $                   ZUNMQR
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, ZLANGE
      EXTERNAL           LSAME, ILAENV, DLAMCH, ZLANGE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG, MAX, SQRT
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   ABS1
*     ..
*     .. Statement Function definitions ..
      ABS1( X ) = ABS( DBLE( X ) ) + ABS( DIMAG( X ) )
*     ..
*     .. Executable Statements ..
*
*     Decode the input arguments
*
      IF( LSAME( JOBVL, 'N' ) ) THEN
         IJOBVL = 1
         ILVL = .FALSE.
      ELSE IF( LSAME( JOBVL, 'V' ) ) THEN
         IJOBVL = 2
         ILVL = .TRUE.
      ELSE
         IJOBVL = -1
         ILVL = .FALSE.
      END IF
*
      IF( LSAME( JOBVR, 'N' ) ) THEN
         IJOBVR = 1
         ILVR = .FALSE.
      ELSE IF( LSAME( JOBVR, 'V' ) ) THEN
         IJOBVR = 2
         ILVR = .TRUE.
      ELSE
         IJOBVR = -1
         ILVR = .FALSE.
      END IF
      ILV = ILVL .OR. ILVR
*
*     Test the input arguments
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( IJOBVL.LE.0 ) THEN
         INFO = -1
      ELSE IF( IJOBVR.LE.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDVL.LT.1 .OR. ( ILVL .AND. LDVL.LT.N ) ) THEN
         INFO = -11
      ELSE IF( LDVR.LT.1 .OR. ( ILVR .AND. LDVR.LT.N ) ) THEN
         INFO = -13
      END IF
*
*     Compute workspace
*      (Note: Comments in the code beginning "Workspace:" describe the
*       minimal amount of workspace needed at that point in the code,
*       as well as the preferred amount for good performance.
*       NB refers to the optimal block size for the immediately
*       following subroutine, as returned by ILAENV. The workspace is
*       computed assuming ILO = 1 and IHI = N, the worst case.)
*
      IF( INFO.EQ.0 ) THEN
         LWKMIN = MAX( 1, 2*N )
         LWKOPT = MAX( 1, N + N*ILAENV( 1, 'ZGEQRF', ' ', N, 1, N, 0 ) )
         LWKOPT = MAX( LWKOPT, N +
     $                 N*ILAENV( 1, 'ZUNMQR', ' ', N, 1, N, 0 ) )
         IF( ILVL ) THEN
            LWKOPT = MAX( LWKOPT, N +
     $                    N*ILAENV( 1, 'ZUNGQR', ' ', N, 1, N, -1 ) )
         END IF
         WORK( 1 ) = LWKOPT
*
         IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY )
     $      INFO = -15
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGGEV ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Get machine constants
*
      EPS = DLAMCH( 'E' )*DLAMCH( 'B' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM
*
*     Scale A if max element outside range [SMLNUM,BIGNUM]
*
      ANRM = ZLANGE( 'M', N, N, A, LDA, RWORK )
      ILASCL = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         ANRMTO = SMLNUM
         ILASCL = .TRUE.
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         ANRMTO = BIGNUM
         ILASCL = .TRUE.
      END IF
      IF( ILASCL )
     $   CALL ZLASCL( 'G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR )
*
*     Scale B if max element outside range [SMLNUM,BIGNUM]
*
      BNRM = ZLANGE( 'M', N, N, B, LDB, RWORK )
      ILBSCL = .FALSE.
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
         BNRMTO = SMLNUM
         ILBSCL = .TRUE.
      ELSE IF( BNRM.GT.BIGNUM ) THEN
         BNRMTO = BIGNUM
         ILBSCL = .TRUE.
      END IF
      IF( ILBSCL )
     $   CALL ZLASCL( 'G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR )
*
*     Permute the matrices A, B to isolate eigenvalues if possible
*     (Real Workspace: need 6*N)
*
      ILEFT = 1
      IRIGHT = N + 1
      IRWRK = IRIGHT + N
      CALL ZGGBAL( 'P', N, A, LDA, B, LDB, ILO, IHI, RWORK( ILEFT ),
     $             RWORK( IRIGHT ), RWORK( IRWRK ), IERR )
*
*     Reduce B to triangular form (QR decomposition of B)
*     (Complex Workspace: need N, prefer N*NB)
*
      IROWS = IHI + 1 - ILO
      IF( ILV ) THEN
         ICOLS = N + 1 - ILO
      ELSE
         ICOLS = IROWS
      END IF
      ITAU = 1
      IWRK = ITAU + IROWS
      CALL ZGEQRF( IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ),
     $             WORK( IWRK ), LWORK+1-IWRK, IERR )
*
*     Apply the orthogonal transformation to matrix A
*     (Complex Workspace: need N, prefer N*NB)
*
      CALL ZUNMQR( 'L', 'C', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB,
     $             WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ),
     $             LWORK+1-IWRK, IERR )
*
*     Initialize VL
*     (Complex Workspace: need N, prefer N*NB)
*
      IF( ILVL ) THEN
         CALL ZLASET( 'Full', N, N, CZERO, CONE, VL, LDVL )
         IF( IROWS.GT.1 ) THEN
            CALL ZLACPY( 'L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB,
     $                   VL( ILO+1, ILO ), LDVL )
         END IF
         CALL ZUNGQR( IROWS, IROWS, IROWS, VL( ILO, ILO ), LDVL,
     $                WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR )
      END IF
*
*     Initialize VR
*
      IF( ILVR )
     $   CALL ZLASET( 'Full', N, N, CZERO, CONE, VR, LDVR )
*
*     Reduce to generalized Hessenberg form
*
      IF( ILV ) THEN
*
*        Eigenvectors requested -- work on whole matrix.
*
         CALL ZGGHRD( JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, VL,
     $                LDVL, VR, LDVR, IERR )
      ELSE
         CALL ZGGHRD( 'N', 'N', IROWS, 1, IROWS, A( ILO, ILO ), LDA,
     $                B( ILO, ILO ), LDB, VL, LDVL, VR, LDVR, IERR )
      END IF
*
*     Perform QZ algorithm (Compute eigenvalues, and optionally, the
*     Schur form and Schur vectors)
*     (Complex Workspace: need N)
*     (Real Workspace: need N)
*
      IWRK = ITAU
      IF( ILV ) THEN
         CHTEMP = 'S'
      ELSE
         CHTEMP = 'E'
      END IF
      CALL ZHGEQZ( CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB,
     $             ALPHA, BETA, VL, LDVL, VR, LDVR, WORK( IWRK ),
     $             LWORK+1-IWRK, RWORK( IRWRK ), IERR )
      IF( IERR.NE.0 ) THEN
         IF( IERR.GT.0 .AND. IERR.LE.N ) THEN
            INFO = IERR
         ELSE IF( IERR.GT.N .AND. IERR.LE.2*N ) THEN
            INFO = IERR - N
         ELSE
            INFO = N + 1
         END IF
         GO TO 70
      END IF
*
*     Compute Eigenvectors
*     (Real Workspace: need 2*N)
*     (Complex Workspace: need 2*N)
*
      IF( ILV ) THEN
         IF( ILVL ) THEN
            IF( ILVR ) THEN
               CHTEMP = 'B'
            ELSE
               CHTEMP = 'L'
            END IF
         ELSE
            CHTEMP = 'R'
         END IF
*
         CALL ZTGEVC( CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL,
     $                VR, LDVR, N, IN, WORK( IWRK ), RWORK( IRWRK ),
     $                IERR )
         IF( IERR.NE.0 ) THEN
            INFO = N + 2
            GO TO 70
         END IF
*
*        Undo balancing on VL and VR and normalization
*        (Workspace: none needed)
*
         IF( ILVL ) THEN
            CALL ZGGBAK( 'P', 'L', N, ILO, IHI, RWORK( ILEFT ),
     $                   RWORK( IRIGHT ), N, VL, LDVL, IERR )
            DO 30 JC = 1, N
               TEMP = ZERO
               DO 10 JR = 1, N
                  TEMP = MAX( TEMP, ABS1( VL( JR, JC ) ) )
   10          CONTINUE
               IF( TEMP.LT.SMLNUM )
     $            GO TO 30
               TEMP = ONE / TEMP
               DO 20 JR = 1, N
                  VL( JR, JC ) = VL( JR, JC )*TEMP
   20          CONTINUE
   30       CONTINUE
         END IF
         IF( ILVR ) THEN
            CALL ZGGBAK( 'P', 'R', N, ILO, IHI, RWORK( ILEFT ),
     $                   RWORK( IRIGHT ), N, VR, LDVR, IERR )
            DO 60 JC = 1, N
               TEMP = ZERO
               DO 40 JR = 1, N
                  TEMP = MAX( TEMP, ABS1( VR( JR, JC ) ) )
   40          CONTINUE
               IF( TEMP.LT.SMLNUM )
     $            GO TO 60
               TEMP = ONE / TEMP
               DO 50 JR = 1, N
                  VR( JR, JC ) = VR( JR, JC )*TEMP
   50          CONTINUE
   60       CONTINUE
         END IF
      END IF
*
*     Undo scaling if necessary
*
   70 CONTINUE
*
      IF( ILASCL )
     $   CALL ZLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHA, N, IERR )
*
      IF( ILBSCL )
     $   CALL ZLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR )
*
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of ZGGEV
*
      END
*> \brief \b ZTGEVC
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZTGEVC + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztgevc.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztgevc.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztgevc.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZTGEVC( SIDE, HOWMNY, SELECT, N, S, LDS, P, LDP, VL,
*                          LDVL, VR, LDVR, MM, M, WORK, RWORK, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          HOWMNY, SIDE
*       INTEGER            INFO, LDP, LDS, LDVL, LDVR, M, MM, N
*       ..
*       .. Array Arguments ..
*       LOGICAL            SELECT( * )
*       DOUBLE PRECISION   RWORK( * )
*       COMPLEX*16         P( LDP, * ), S( LDS, * ), VL( LDVL, * ),
*      $                   VR( LDVR, * ), WORK( * )
*       ..
*  
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZTGEVC computes some or all of the right and/or left eigenvectors of
*> a pair of complex matrices (S,P), where S and P are upper triangular.
*> Matrix pairs of this type are produced by the generalized Schur
*> factorization of a complex matrix pair (A,B):
*> 
*>    A = Q*S*Z**H,  B = Q*P*Z**H
*> 
*> as computed by ZGGHRD + ZHGEQZ.
*> 
*> The right eigenvector x and the left eigenvector y of (S,P)
*> corresponding to an eigenvalue w are defined by:
*> 
*>    S*x = w*P*x,  (y**H)*S = w*(y**H)*P,
*> 
*> where y**H denotes the conjugate tranpose of y.
*> The eigenvalues are not input to this routine, but are computed
*> directly from the diagonal elements of S and P.
*> 
*> This routine returns the matrices X and/or Y of right and left
*> eigenvectors of (S,P), or the products Z*X and/or Q*Y,
*> where Z and Q are input matrices.
*> If Q and Z are the unitary factors from the generalized Schur
*> factorization of a matrix pair (A,B), then Z*X and Q*Y
*> are the matrices of right and left eigenvectors of (A,B).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'R': compute right eigenvectors only;
*>          = 'L': compute left eigenvectors only;
*>          = 'B': compute both right and left eigenvectors.
*> \endverbatim
*>
*> \param[in] HOWMNY
*> \verbatim
*>          HOWMNY is CHARACTER*1
*>          = 'A': compute all right and/or left eigenvectors;
*>          = 'B': compute all right and/or left eigenvectors,
*>                 backtransformed by the matrices in VR and/or VL;
*>          = 'S': compute selected right and/or left eigenvectors,
*>                 specified by the logical array SELECT.
*> \endverbatim
*>
*> \param[in] SELECT
*> \verbatim
*>          SELECT is LOGICAL array, dimension (N)
*>          If HOWMNY='S', SELECT specifies the eigenvectors to be
*>          computed.  The eigenvector corresponding to the j-th
*>          eigenvalue is computed if SELECT(j) = .TRUE..
*>          Not referenced if HOWMNY = 'A' or 'B'.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices S and P.  N >= 0.
*> \endverbatim
*>
*> \param[in] S
*> \verbatim
*>          S is COMPLEX*16 array, dimension (LDS,N)
*>          The upper triangular matrix S from a generalized Schur
*>          factorization, as computed by ZHGEQZ.
*> \endverbatim
*>
*> \param[in] LDS
*> \verbatim
*>          LDS is INTEGER
*>          The leading dimension of array S.  LDS >= max(1,N).
*> \endverbatim
*>
*> \param[in] P
*> \verbatim
*>          P is COMPLEX*16 array, dimension (LDP,N)
*>          The upper triangular matrix P from a generalized Schur
*>          factorization, as computed by ZHGEQZ.  P must have real
*>          diagonal elements.
*> \endverbatim
*>
*> \param[in] LDP
*> \verbatim
*>          LDP is INTEGER
*>          The leading dimension of array P.  LDP >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] VL
*> \verbatim
*>          VL is COMPLEX*16 array, dimension (LDVL,MM)
*>          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
*>          contain an N-by-N matrix Q (usually the unitary matrix Q
*>          of left Schur vectors returned by ZHGEQZ).
*>          On exit, if SIDE = 'L' or 'B', VL contains:
*>          if HOWMNY = 'A', the matrix Y of left eigenvectors of (S,P);
*>          if HOWMNY = 'B', the matrix Q*Y;
*>          if HOWMNY = 'S', the left eigenvectors of (S,P) specified by
*>                      SELECT, stored consecutively in the columns of
*>                      VL, in the same order as their eigenvalues.
*>          Not referenced if SIDE = 'R'.
*> \endverbatim
*>
*> \param[in] LDVL
*> \verbatim
*>          LDVL is INTEGER
*>          The leading dimension of array VL.  LDVL >= 1, and if
*>          SIDE = 'L' or 'l' or 'B' or 'b', LDVL >= N.
*> \endverbatim
*>
*> \param[in,out] VR
*> \verbatim
*>          VR is COMPLEX*16 array, dimension (LDVR,MM)
*>          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
*>          contain an N-by-N matrix Q (usually the unitary matrix Z
*>          of right Schur vectors returned by ZHGEQZ).
*>          On exit, if SIDE = 'R' or 'B', VR contains:
*>          if HOWMNY = 'A', the matrix X of right eigenvectors of (S,P);
*>          if HOWMNY = 'B', the matrix Z*X;
*>          if HOWMNY = 'S', the right eigenvectors of (S,P) specified by
*>                      SELECT, stored consecutively in the columns of
*>                      VR, in the same order as their eigenvalues.
*>          Not referenced if SIDE = 'L'.
*> \endverbatim
*>
*> \param[in] LDVR
*> \verbatim
*>          LDVR is INTEGER
*>          The leading dimension of the array VR.  LDVR >= 1, and if
*>          SIDE = 'R' or 'B', LDVR >= N.
*> \endverbatim
*>
*> \param[in] MM
*> \verbatim
*>          MM is INTEGER
*>          The number of columns in the arrays VL and/or VR. MM >= M.
*> \endverbatim
*>
*> \param[out] M
*> \verbatim
*>          M is INTEGER
*>          The number of columns in the arrays VL and/or VR actually
*>          used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M
*>          is set to N.  Each selected eigenvector occupies one column.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (2*N)
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is DOUBLE PRECISION array, dimension (2*N)
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
*> \date November 2011
*
*> \ingroup complex16GEcomputational
*
*  =====================================================================
      SUBROUTINE ZTGEVC( SIDE, HOWMNY, SELECT, N, S, LDS, P, LDP, VL,
     $                   LDVL, VR, LDVR, MM, M, WORK, RWORK, INFO )
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER          HOWMNY, SIDE
      INTEGER            INFO, LDP, LDS, LDVL, LDVR, M, MM, N
*     ..
*     .. Array Arguments ..
      LOGICAL            SELECT( * )
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         P( LDP, * ), S( LDS, * ), VL( LDVL, * ),
     $                   VR( LDVR, * ), WORK( * )
*     ..
*
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            COMPL, COMPR, ILALL, ILBACK, ILBBAD, ILCOMP,
     $                   LSA, LSB
      INTEGER            I, IBEG, IEIG, IEND, IHWMNY, IM, ISIDE, ISRC,
     $                   J, JE, JR
      DOUBLE PRECISION   ACOEFA, ACOEFF, ANORM, ASCALE, BCOEFA, BIG,
     $                   BIGNUM, BNORM, BSCALE, DMIN, SAFMIN, SBETA,
     $                   SCALE, SMALL, TEMP, ULP, XMAX
      COMPLEX*16         BCOEFF, CA, CB, D, SALPHA, SUM, SUMA, SUMB, X
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      COMPLEX*16         ZLADIV
      EXTERNAL           LSAME, DLAMCH, ZLADIV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLABAD, XERBLA, ZGEMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DCONJG, DIMAG, MAX, MIN
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   ABS1
*     ..
*     .. Statement Function definitions ..
      ABS1( X ) = ABS( DBLE( X ) ) + ABS( DIMAG( X ) )
*     ..
*     .. Executable Statements ..
*
*     Decode and Test the input parameters
*
      IF( LSAME( HOWMNY, 'A' ) ) THEN
         IHWMNY = 1
         ILALL = .TRUE.
         ILBACK = .FALSE.
      ELSE IF( LSAME( HOWMNY, 'S' ) ) THEN
         IHWMNY = 2
         ILALL = .FALSE.
         ILBACK = .FALSE.
      ELSE IF( LSAME( HOWMNY, 'B' ) ) THEN
         IHWMNY = 3
         ILALL = .TRUE.
         ILBACK = .TRUE.
      ELSE
         IHWMNY = -1
      END IF
*
      IF( LSAME( SIDE, 'R' ) ) THEN
         ISIDE = 1
         COMPL = .FALSE.
         COMPR = .TRUE.
      ELSE IF( LSAME( SIDE, 'L' ) ) THEN
         ISIDE = 2
         COMPL = .TRUE.
         COMPR = .FALSE.
      ELSE IF( LSAME( SIDE, 'B' ) ) THEN
         ISIDE = 3
         COMPL = .TRUE.
         COMPR = .TRUE.
      ELSE
         ISIDE = -1
      END IF
*
      INFO = 0
      IF( ISIDE.LT.0 ) THEN
         INFO = -1
      ELSE IF( IHWMNY.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDS.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDP.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTGEVC', -INFO )
         RETURN
      END IF
*
*     Count the number of eigenvectors
*
      IF( .NOT.ILALL ) THEN
         IM = 0
         DO 10 J = 1, N
            IF( SELECT( J ) )
     $         IM = IM + 1
   10    CONTINUE
      ELSE
         IM = N
      END IF
*
*     Check diagonal of B
*
      ILBBAD = .FALSE.
      DO 20 J = 1, N
         IF( DIMAG( P( J, J ) ).NE.ZERO )
     $      ILBBAD = .TRUE.
   20 CONTINUE
*
      IF( ILBBAD ) THEN
         INFO = -7
      ELSE IF( COMPL .AND. LDVL.LT.N .OR. LDVL.LT.1 ) THEN
         INFO = -10
      ELSE IF( COMPR .AND. LDVR.LT.N .OR. LDVR.LT.1 ) THEN
         INFO = -12
      ELSE IF( MM.LT.IM ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTGEVC', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      M = IM
      IF( N.EQ.0 )
     $   RETURN
*
*     Machine Constants
*
      SAFMIN = DLAMCH( 'Safe minimum' )
      BIG = ONE / SAFMIN
      CALL DLABAD( SAFMIN, BIG )
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
      SMALL = SAFMIN*N / ULP
      BIG = ONE / SMALL
      BIGNUM = ONE / ( SAFMIN*N )
*
*     Compute the 1-norm of each column of the strictly upper triangular
*     part of A and B to check for possible overflow in the triangular
*     solver.
*
      ANORM = ABS1( S( 1, 1 ) )
      BNORM = ABS1( P( 1, 1 ) )
      RWORK( 1 ) = ZERO
      RWORK( N+1 ) = ZERO
      DO 40 J = 2, N
         RWORK( J ) = ZERO
         RWORK( N+J ) = ZERO
         DO 30 I = 1, J - 1
            RWORK( J ) = RWORK( J ) + ABS1( S( I, J ) )
            RWORK( N+J ) = RWORK( N+J ) + ABS1( P( I, J ) )
   30    CONTINUE
         ANORM = MAX( ANORM, RWORK( J )+ABS1( S( J, J ) ) )
         BNORM = MAX( BNORM, RWORK( N+J )+ABS1( P( J, J ) ) )
   40 CONTINUE
*
      ASCALE = ONE / MAX( ANORM, SAFMIN )
      BSCALE = ONE / MAX( BNORM, SAFMIN )
*
*     Left eigenvectors
*
      IF( COMPL ) THEN
         IEIG = 0
*
*        Main loop over eigenvalues
*
         DO 140 JE = 1, N
            IF( ILALL ) THEN
               ILCOMP = .TRUE.
            ELSE
               ILCOMP = SELECT( JE )
            END IF
            IF( ILCOMP ) THEN
               IEIG = IEIG + 1
*
               IF( ABS1( S( JE, JE ) ).LE.SAFMIN .AND.
     $             ABS( DBLE( P( JE, JE ) ) ).LE.SAFMIN ) THEN
*
*                 Singular matrix pencil -- return unit eigenvector
*
                  DO 50 JR = 1, N
                     VL( JR, IEIG ) = CZERO
   50             CONTINUE
                  VL( IEIG, IEIG ) = CONE
                  GO TO 140
               END IF
*
*              Non-singular eigenvalue:
*              Compute coefficients  a  and  b  in
*                   H
*                 y  ( a A - b B ) = 0
*
               TEMP = ONE / MAX( ABS1( S( JE, JE ) )*ASCALE,
     $                ABS( DBLE( P( JE, JE ) ) )*BSCALE, SAFMIN )
               SALPHA = ( TEMP*S( JE, JE ) )*ASCALE
               SBETA = ( TEMP*DBLE( P( JE, JE ) ) )*BSCALE
               ACOEFF = SBETA*ASCALE
               BCOEFF = SALPHA*BSCALE
*
*              Scale to avoid underflow
*
               LSA = ABS( SBETA ).GE.SAFMIN .AND. ABS( ACOEFF ).LT.SMALL
               LSB = ABS1( SALPHA ).GE.SAFMIN .AND. ABS1( BCOEFF ).LT.
     $               SMALL
*
               SCALE = ONE
               IF( LSA )
     $            SCALE = ( SMALL / ABS( SBETA ) )*MIN( ANORM, BIG )
               IF( LSB )
     $            SCALE = MAX( SCALE, ( SMALL / ABS1( SALPHA ) )*
     $                    MIN( BNORM, BIG ) )
               IF( LSA .OR. LSB ) THEN
                  SCALE = MIN( SCALE, ONE /
     $                    ( SAFMIN*MAX( ONE, ABS( ACOEFF ),
     $                    ABS1( BCOEFF ) ) ) )
                  IF( LSA ) THEN
                     ACOEFF = ASCALE*( SCALE*SBETA )
                  ELSE
                     ACOEFF = SCALE*ACOEFF
                  END IF
                  IF( LSB ) THEN
                     BCOEFF = BSCALE*( SCALE*SALPHA )
                  ELSE
                     BCOEFF = SCALE*BCOEFF
                  END IF
               END IF
*
               ACOEFA = ABS( ACOEFF )
               BCOEFA = ABS1( BCOEFF )
               XMAX = ONE
               DO 60 JR = 1, N
                  WORK( JR ) = CZERO
   60          CONTINUE
               WORK( JE ) = CONE
               DMIN = MAX( ULP*ACOEFA*ANORM, ULP*BCOEFA*BNORM, SAFMIN )
*
*                                              H
*              Triangular solve of  (a A - b B)  y = 0
*
*                                      H
*              (rowwise in  (a A - b B) , or columnwise in a A - b B)
*
               DO 100 J = JE + 1, N
*
*                 Compute
*                       j-1
*                 SUM = sum  conjg( a*S(k,j) - b*P(k,j) )*x(k)
*                       k=je
*                 (Scale if necessary)
*
                  TEMP = ONE / XMAX
                  IF( ACOEFA*RWORK( J )+BCOEFA*RWORK( N+J ).GT.BIGNUM*
     $                TEMP ) THEN
                     DO 70 JR = JE, J - 1
                        WORK( JR ) = TEMP*WORK( JR )
   70                CONTINUE
                     XMAX = ONE
                  END IF
                  SUMA = CZERO
                  SUMB = CZERO
*
                  DO 80 JR = JE, J - 1
                     SUMA = SUMA + DCONJG( S( JR, J ) )*WORK( JR )
                     SUMB = SUMB + DCONJG( P( JR, J ) )*WORK( JR )
   80             CONTINUE
                  SUM = ACOEFF*SUMA - DCONJG( BCOEFF )*SUMB
*
*                 Form x(j) = - SUM / conjg( a*S(j,j) - b*P(j,j) )
*
*                 with scaling and perturbation of the denominator
*
                  D = DCONJG( ACOEFF*S( J, J )-BCOEFF*P( J, J ) )
                  IF( ABS1( D ).LE.DMIN )
     $               D = DCMPLX( DMIN )
*
                  IF( ABS1( D ).LT.ONE ) THEN
                     IF( ABS1( SUM ).GE.BIGNUM*ABS1( D ) ) THEN
                        TEMP = ONE / ABS1( SUM )
                        DO 90 JR = JE, J - 1
                           WORK( JR ) = TEMP*WORK( JR )
   90                   CONTINUE
                        XMAX = TEMP*XMAX
                        SUM = TEMP*SUM
                     END IF
                  END IF
                  WORK( J ) = ZLADIV( -SUM, D )
                  XMAX = MAX( XMAX, ABS1( WORK( J ) ) )
  100          CONTINUE
*
*              Back transform eigenvector if HOWMNY='B'.
*
               IF( ILBACK ) THEN
                  CALL ZGEMV( 'N', N, N+1-JE, CONE, VL( 1, JE ), LDVL,
     $                        WORK( JE ), 1, CZERO, WORK( N+1 ), 1 )
                  ISRC = 2
                  IBEG = 1
               ELSE
                  ISRC = 1
                  IBEG = JE
               END IF
*
*              Copy and scale eigenvector into column of VL
*
               XMAX = ZERO
               DO 110 JR = IBEG, N
                  XMAX = MAX( XMAX, ABS1( WORK( ( ISRC-1 )*N+JR ) ) )
  110          CONTINUE
*
               IF( XMAX.GT.SAFMIN ) THEN
                  TEMP = ONE / XMAX
                  DO 120 JR = IBEG, N
                     VL( JR, IEIG ) = TEMP*WORK( ( ISRC-1 )*N+JR )
  120             CONTINUE
               ELSE
                  IBEG = N + 1
               END IF
*
               DO 130 JR = 1, IBEG - 1
                  VL( JR, IEIG ) = CZERO
  130          CONTINUE
*
            END IF
  140    CONTINUE
      END IF
*
*     Right eigenvectors
*
      IF( COMPR ) THEN
         IEIG = IM + 1
*
*        Main loop over eigenvalues
*
         DO 250 JE = N, 1, -1
            IF( ILALL ) THEN
               ILCOMP = .TRUE.
            ELSE
               ILCOMP = SELECT( JE )
            END IF
            IF( ILCOMP ) THEN
               IEIG = IEIG - 1
*
               IF( ABS1( S( JE, JE ) ).LE.SAFMIN .AND.
     $             ABS( DBLE( P( JE, JE ) ) ).LE.SAFMIN ) THEN
*
*                 Singular matrix pencil -- return unit eigenvector
*
                  DO 150 JR = 1, N
                     VR( JR, IEIG ) = CZERO
  150             CONTINUE
                  VR( IEIG, IEIG ) = CONE
                  GO TO 250
               END IF
*
*              Non-singular eigenvalue:
*              Compute coefficients  a  and  b  in
*
*              ( a A - b B ) x  = 0
*
               TEMP = ONE / MAX( ABS1( S( JE, JE ) )*ASCALE,
     $                ABS( DBLE( P( JE, JE ) ) )*BSCALE, SAFMIN )
               SALPHA = ( TEMP*S( JE, JE ) )*ASCALE
               SBETA = ( TEMP*DBLE( P( JE, JE ) ) )*BSCALE
               ACOEFF = SBETA*ASCALE
               BCOEFF = SALPHA*BSCALE
*
*              Scale to avoid underflow
*
               LSA = ABS( SBETA ).GE.SAFMIN .AND. ABS( ACOEFF ).LT.SMALL
               LSB = ABS1( SALPHA ).GE.SAFMIN .AND. ABS1( BCOEFF ).LT.
     $               SMALL
*
               SCALE = ONE
               IF( LSA )
     $            SCALE = ( SMALL / ABS( SBETA ) )*MIN( ANORM, BIG )
               IF( LSB )
     $            SCALE = MAX( SCALE, ( SMALL / ABS1( SALPHA ) )*
     $                    MIN( BNORM, BIG ) )
               IF( LSA .OR. LSB ) THEN
                  SCALE = MIN( SCALE, ONE /
     $                    ( SAFMIN*MAX( ONE, ABS( ACOEFF ),
     $                    ABS1( BCOEFF ) ) ) )
                  IF( LSA ) THEN
                     ACOEFF = ASCALE*( SCALE*SBETA )
                  ELSE
                     ACOEFF = SCALE*ACOEFF
                  END IF
                  IF( LSB ) THEN
                     BCOEFF = BSCALE*( SCALE*SALPHA )
                  ELSE
                     BCOEFF = SCALE*BCOEFF
                  END IF
               END IF
*
               ACOEFA = ABS( ACOEFF )
               BCOEFA = ABS1( BCOEFF )
               XMAX = ONE
               DO 160 JR = 1, N
                  WORK( JR ) = CZERO
  160          CONTINUE
               WORK( JE ) = CONE
               DMIN = MAX( ULP*ACOEFA*ANORM, ULP*BCOEFA*BNORM, SAFMIN )
*
*              Triangular solve of  (a A - b B) x = 0  (columnwise)
*
*              WORK(1:j-1) contains sums w,
*              WORK(j+1:JE) contains x
*
               DO 170 JR = 1, JE - 1
                  WORK( JR ) = ACOEFF*S( JR, JE ) - BCOEFF*P( JR, JE )
  170          CONTINUE
               WORK( JE ) = CONE
*
               DO 210 J = JE - 1, 1, -1
*
*                 Form x(j) := - w(j) / d
*                 with scaling and perturbation of the denominator
*
                  D = ACOEFF*S( J, J ) - BCOEFF*P( J, J )
                  IF( ABS1( D ).LE.DMIN )
     $               D = DCMPLX( DMIN )
*
                  IF( ABS1( D ).LT.ONE ) THEN
                     IF( ABS1( WORK( J ) ).GE.BIGNUM*ABS1( D ) ) THEN
                        TEMP = ONE / ABS1( WORK( J ) )
                        DO 180 JR = 1, JE
                           WORK( JR ) = TEMP*WORK( JR )
  180                   CONTINUE
                     END IF
                  END IF
*
                  WORK( J ) = ZLADIV( -WORK( J ), D )
*
                  IF( J.GT.1 ) THEN
*
*                    w = w + x(j)*(a S(*,j) - b P(*,j) ) with scaling
*
                     IF( ABS1( WORK( J ) ).GT.ONE ) THEN
                        TEMP = ONE / ABS1( WORK( J ) )
                        IF( ACOEFA*RWORK( J )+BCOEFA*RWORK( N+J ).GE.
     $                      BIGNUM*TEMP ) THEN
                           DO 190 JR = 1, JE
                              WORK( JR ) = TEMP*WORK( JR )
  190                      CONTINUE
                        END IF
                     END IF
*
                     CA = ACOEFF*WORK( J )
                     CB = BCOEFF*WORK( J )
                     DO 200 JR = 1, J - 1
                        WORK( JR ) = WORK( JR ) + CA*S( JR, J ) -
     $                               CB*P( JR, J )
  200                CONTINUE
                  END IF
  210          CONTINUE
*
*              Back transform eigenvector if HOWMNY='B'.
*
               IF( ILBACK ) THEN
                  CALL ZGEMV( 'N', N, JE, CONE, VR, LDVR, WORK, 1,
     $                        CZERO, WORK( N+1 ), 1 )
                  ISRC = 2
                  IEND = N
               ELSE
                  ISRC = 1
                  IEND = JE
               END IF
*
*              Copy and scale eigenvector into column of VR
*
               XMAX = ZERO
               DO 220 JR = 1, IEND
                  XMAX = MAX( XMAX, ABS1( WORK( ( ISRC-1 )*N+JR ) ) )
  220          CONTINUE
*
               IF( XMAX.GT.SAFMIN ) THEN
                  TEMP = ONE / XMAX
                  DO 230 JR = 1, IEND
                     VR( JR, IEIG ) = TEMP*WORK( ( ISRC-1 )*N+JR )
  230             CONTINUE
               ELSE
                  IEND = 0
               END IF
*
               DO 240 JR = IEND + 1, N
                  VR( JR, IEIG ) = CZERO
  240          CONTINUE
*
            END IF
  250    CONTINUE
      END IF
*
      RETURN
*
*     End of ZTGEVC
*
      END
