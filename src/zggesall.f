*> \brief <b> ZGGES computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors for GE matrices</b>
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZGGES + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgges.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgges.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgges.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGGES( JOBVSL, JOBVSR, SORT, SELCTG, N, A, LDA, B, LDB,
*                         SDIM, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK,
*                         LWORK, RWORK, BWORK, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          JOBVSL, JOBVSR, SORT
*       INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N, SDIM
*       ..
*       .. Array Arguments ..
*       LOGICAL            BWORK( * )
*       DOUBLE PRECISION   RWORK( * )
*       COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),
*      $                   BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ),
*      $                   WORK( * )
*       ..
*       .. Function Arguments ..
*       LOGICAL            SELCTG
*       EXTERNAL           SELCTG
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGGES computes for a pair of N-by-N complex nonsymmetric matrices
*> (A,B), the generalized eigenvalues, the generalized complex Schur
*> form (S, T), and optionally left and/or right Schur vectors (VSL
*> and VSR). This gives the generalized Schur factorization
*>
*>         (A,B) = ( (VSL)*S*(VSR)**H, (VSL)*T*(VSR)**H )
*>
*> where (VSR)**H is the conjugate-transpose of VSR.
*>
*> Optionally, it also orders the eigenvalues so that a selected cluster
*> of eigenvalues appears in the leading diagonal blocks of the upper
*> triangular matrix S and the upper triangular matrix T. The leading
*> columns of VSL and VSR then form an unitary basis for the
*> corresponding left and right eigenspaces (deflating subspaces).
*>
*> (If only the generalized eigenvalues are needed, use the driver
*> ZGGEV instead, which is faster.)
*>
*> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
*> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
*> usually represented as the pair (alpha,beta), as there is a
*> reasonable interpretation for beta=0, and even for both being zero.
*>
*> A pair of matrices (S,T) is in generalized complex Schur form if S
*> and T are upper triangular and, in addition, the diagonal elements
*> of T are non-negative real numbers.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] JOBVSL
*> \verbatim
*>          JOBVSL is CHARACTER*1
*>          = 'N':  do not compute the left Schur vectors;
*>          = 'V':  compute the left Schur vectors.
*> \endverbatim
*>
*> \param[in] JOBVSR
*> \verbatim
*>          JOBVSR is CHARACTER*1
*>          = 'N':  do not compute the right Schur vectors;
*>          = 'V':  compute the right Schur vectors.
*> \endverbatim
*>
*> \param[in] SORT
*> \verbatim
*>          SORT is CHARACTER*1
*>          Specifies whether or not to order the eigenvalues on the
*>          diagonal of the generalized Schur form.
*>          = 'N':  Eigenvalues are not ordered;
*>          = 'S':  Eigenvalues are ordered (see SELCTG).
*> \endverbatim
*>
*> \param[in] SELCTG
*> \verbatim
*>          SELCTG is a LOGICAL FUNCTION of two COMPLEX*16 arguments
*>          SELCTG must be declared EXTERNAL in the calling subroutine.
*>          If SORT = 'N', SELCTG is not referenced.
*>          If SORT = 'S', SELCTG is used to select eigenvalues to sort
*>          to the top left of the Schur form.
*>          An eigenvalue ALPHA(j)/BETA(j) is selected if
*>          SELCTG(ALPHA(j),BETA(j)) is true.
*>
*>          Note that a selected complex eigenvalue may no longer satisfy
*>          SELCTG(ALPHA(j),BETA(j)) = .TRUE. after ordering, since
*>          ordering may change the value of complex eigenvalues
*>          (especially if the eigenvalue is ill-conditioned), in this
*>          case INFO is set to N+2 (See INFO below).
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices A, B, VSL, and VSR.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA, N)
*>          On entry, the first of the pair of matrices.
*>          On exit, A has been overwritten by its generalized Schur
*>          form S.
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
*>          On entry, the second of the pair of matrices.
*>          On exit, B has been overwritten by its generalized Schur
*>          form T.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] SDIM
*> \verbatim
*>          SDIM is INTEGER
*>          If SORT = 'N', SDIM = 0.
*>          If SORT = 'S', SDIM = number of eigenvalues (after sorting)
*>          for which SELCTG is true.
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
*>          On exit,  ALPHA(j)/BETA(j), j=1,...,N, will be the
*>          generalized eigenvalues.  ALPHA(j), j=1,...,N  and  BETA(j),
*>          j=1,...,N  are the diagonals of the complex Schur form (A,B)
*>          output by ZGGES. The  BETA(j) will be non-negative real.
*>
*>          Note: the quotients ALPHA(j)/BETA(j) may easily over- or
*>          underflow, and BETA(j) may even be zero.  Thus, the user
*>          should avoid naively computing the ratio alpha/beta.
*>          However, ALPHA will be always less than and usually
*>          comparable with norm(A) in magnitude, and BETA always less
*>          than and usually comparable with norm(B).
*> \endverbatim
*>
*> \param[out] VSL
*> \verbatim
*>          VSL is COMPLEX*16 array, dimension (LDVSL,N)
*>          If JOBVSL = 'V', VSL will contain the left Schur vectors.
*>          Not referenced if JOBVSL = 'N'.
*> \endverbatim
*>
*> \param[in] LDVSL
*> \verbatim
*>          LDVSL is INTEGER
*>          The leading dimension of the matrix VSL. LDVSL >= 1, and
*>          if JOBVSL = 'V', LDVSL >= N.
*> \endverbatim
*>
*> \param[out] VSR
*> \verbatim
*>          VSR is COMPLEX*16 array, dimension (LDVSR,N)
*>          If JOBVSR = 'V', VSR will contain the right Schur vectors.
*>          Not referenced if JOBVSR = 'N'.
*> \endverbatim
*>
*> \param[in] LDVSR
*> \verbatim
*>          LDVSR is INTEGER
*>          The leading dimension of the matrix VSR. LDVSR >= 1, and
*>          if JOBVSR = 'V', LDVSR >= N.
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
*> \param[out] BWORK
*> \verbatim
*>          BWORK is LOGICAL array, dimension (N)
*>          Not referenced if SORT = 'N'.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          =1,...,N:
*>                The QZ iteration failed.  (A,B) are not in Schur
*>                form, but ALPHA(j) and BETA(j) should be correct for
*>                j=INFO+1,...,N.
*>          > N:  =N+1: other than QZ iteration failed in ZHGEQZ
*>                =N+2: after reordering, roundoff changed values of
*>                      some complex eigenvalues so that leading
*>                      eigenvalues in the Generalized Schur form no
*>                      longer satisfy SELCTG=.TRUE.  This could also
*>                      be caused due to scaling.
*>                =N+3: reordering failed in ZTGSEN.
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
*> \ingroup complex16GEeigen
*
*  =====================================================================
      SUBROUTINE ZGGES( JOBVSL, JOBVSR, SORT, SELCTG, N, A, LDA, B, LDB,
     $                  SDIM, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK,
     $                  LWORK, RWORK, BWORK, INFO )
*
*  -- LAPACK driver routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      CHARACTER          JOBVSL, JOBVSR, SORT
      INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N, SDIM
*     ..
*     .. Array Arguments ..
      LOGICAL            BWORK( * )
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),
     $                   BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ),
     $                   WORK( * )
*     ..
*     .. Function Arguments ..
      LOGICAL            SELCTG
      EXTERNAL           SELCTG
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
      LOGICAL            CURSL, ILASCL, ILBSCL, ILVSL, ILVSR, LASTSL,
     $                   LQUERY, WANTST
      INTEGER            I, ICOLS, IERR, IHI, IJOBVL, IJOBVR, ILEFT,
     $                   ILO, IRIGHT, IROWS, IRWRK, ITAU, IWRK, LWKMIN,
     $                   LWKOPT
      DOUBLE PRECISION   ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS, PVSL,
     $                   PVSR, SMLNUM
*     ..
*     .. Local Arrays ..
      INTEGER            IDUM( 1 )
      DOUBLE PRECISION   DIF( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLABAD, XERBLA, ZGEQRF, ZGGBAK, ZGGBAL, ZGGHRD,
     $                   ZHGEQZ, ZLACPY, ZLASCL, ZLASET, ZTGSEN, ZUNGQR,
     $                   ZUNMQR
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, ZLANGE
      EXTERNAL           LSAME, ILAENV, DLAMCH, ZLANGE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Decode the input arguments
*
      IF( LSAME( JOBVSL, 'N' ) ) THEN
         IJOBVL = 1
         ILVSL = .FALSE.
      ELSE IF( LSAME( JOBVSL, 'V' ) ) THEN
         IJOBVL = 2
         ILVSL = .TRUE.
      ELSE
         IJOBVL = -1
         ILVSL = .FALSE.
      END IF
*
      IF( LSAME( JOBVSR, 'N' ) ) THEN
         IJOBVR = 1
         ILVSR = .FALSE.
      ELSE IF( LSAME( JOBVSR, 'V' ) ) THEN
         IJOBVR = 2
         ILVSR = .TRUE.
      ELSE
         IJOBVR = -1
         ILVSR = .FALSE.
      END IF
*
      WANTST = LSAME( SORT, 'S' )
*
*     Test the input arguments
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      IF( IJOBVL.LE.0 ) THEN
         INFO = -1
      ELSE IF( IJOBVR.LE.0 ) THEN
         INFO = -2
      ELSE IF( ( .NOT.WANTST ) .AND. ( .NOT.LSAME( SORT, 'N' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDVSL.LT.1 .OR. ( ILVSL .AND. LDVSL.LT.N ) ) THEN
         INFO = -14
      ELSE IF( LDVSR.LT.1 .OR. ( ILVSR .AND. LDVSR.LT.N ) ) THEN
         INFO = -16
      END IF
*
*     Compute workspace
*      (Note: Comments in the code beginning "Workspace:" describe the
*       minimal amount of workspace needed at that point in the code,
*       as well as the preferred amount for good performance.
*       NB refers to the optimal block size for the immediately
*       following subroutine, as returned by ILAENV.)
*
      IF( INFO.EQ.0 ) THEN
         LWKMIN = MAX( 1, 2*N )
         LWKOPT = MAX( 1, N + N*ILAENV( 1, 'ZGEQRF', ' ', N, 1, N, 0 ) )
         LWKOPT = MAX( LWKOPT, N +
     $                 N*ILAENV( 1, 'ZUNMQR', ' ', N, 1, N, -1 ) )
         IF( ILVSL ) THEN
            LWKOPT = MAX( LWKOPT, N +
     $                    N*ILAENV( 1, 'ZUNGQR', ' ', N, 1, N, -1 ) )
         END IF
         WORK( 1 ) = LWKOPT
*
         IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY )
     $      INFO = -18
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGGES ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         SDIM = 0
         RETURN
      END IF
*
*     Get machine constants
*
      EPS = DLAMCH( 'P' )
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
*
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
*
      IF( ILBSCL )
     $   CALL ZLASCL( 'G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR )
*
*     Permute the matrix to make it more nearly triangular
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
      ICOLS = N + 1 - ILO
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
*     Initialize VSL
*     (Complex Workspace: need N, prefer N*NB)
*
      IF( ILVSL ) THEN
         CALL ZLASET( 'Full', N, N, CZERO, CONE, VSL, LDVSL )
         IF( IROWS.GT.1 ) THEN
            CALL ZLACPY( 'L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB,
     $                   VSL( ILO+1, ILO ), LDVSL )
         END IF
         CALL ZUNGQR( IROWS, IROWS, IROWS, VSL( ILO, ILO ), LDVSL,
     $                WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR )
      END IF
*
*     Initialize VSR
*
      IF( ILVSR )
     $   CALL ZLASET( 'Full', N, N, CZERO, CONE, VSR, LDVSR )
*
*     Reduce to generalized Hessenberg form
*     (Workspace: none needed)
*
      CALL ZGGHRD( JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB, VSL,
     $             LDVSL, VSR, LDVSR, IERR )
*
      SDIM = 0
*
*     Perform QZ algorithm, computing Schur vectors if desired
*     (Complex Workspace: need N)
*     (Real Workspace: need N)
*
      IWRK = ITAU
      CALL ZHGEQZ( 'S', JOBVSL, JOBVSR, N, ILO, IHI, A, LDA, B, LDB,
     $             ALPHA, BETA, VSL, LDVSL, VSR, LDVSR, WORK( IWRK ),
     $             LWORK+1-IWRK, RWORK( IRWRK ), IERR )
      IF( IERR.NE.0 ) THEN
         IF( IERR.GT.0 .AND. IERR.LE.N ) THEN
            INFO = IERR
         ELSE IF( IERR.GT.N .AND. IERR.LE.2*N ) THEN
            INFO = IERR - N
         ELSE
            INFO = N + 1
         END IF
         GO TO 30
      END IF
*
*     Sort eigenvalues ALPHA/BETA if desired
*     (Workspace: none needed)
*
      IF( WANTST ) THEN
*
*        Undo scaling on eigenvalues before selecting
*
         IF( ILASCL )
     $      CALL ZLASCL( 'G', 0, 0, ANRM, ANRMTO, N, 1, ALPHA, N, IERR )
         IF( ILBSCL )
     $      CALL ZLASCL( 'G', 0, 0, BNRM, BNRMTO, N, 1, BETA, N, IERR )
*
*        Select eigenvalues
*
         DO 10 I = 1, N
            BWORK( I ) = SELCTG( ALPHA( I ), BETA( I ) )
   10    CONTINUE
*
         CALL ZTGSEN( 0, ILVSL, ILVSR, BWORK, N, A, LDA, B, LDB, ALPHA,
     $                BETA, VSL, LDVSL, VSR, LDVSR, SDIM, PVSL, PVSR,
     $                DIF, WORK( IWRK ), LWORK-IWRK+1, IDUM, 1, IERR )
         IF( IERR.EQ.1 )
     $      INFO = N + 3
*
      END IF
*
*     Apply back-permutation to VSL and VSR
*     (Workspace: none needed)
*
      IF( ILVSL )
     $   CALL ZGGBAK( 'P', 'L', N, ILO, IHI, RWORK( ILEFT ),
     $                RWORK( IRIGHT ), N, VSL, LDVSL, IERR )
      IF( ILVSR )
     $   CALL ZGGBAK( 'P', 'R', N, ILO, IHI, RWORK( ILEFT ),
     $                RWORK( IRIGHT ), N, VSR, LDVSR, IERR )
*
*     Undo scaling
*
      IF( ILASCL ) THEN
         CALL ZLASCL( 'U', 0, 0, ANRMTO, ANRM, N, N, A, LDA, IERR )
         CALL ZLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHA, N, IERR )
      END IF
*
      IF( ILBSCL ) THEN
         CALL ZLASCL( 'U', 0, 0, BNRMTO, BNRM, N, N, B, LDB, IERR )
         CALL ZLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR )
      END IF
*
      IF( WANTST ) THEN
*
*        Check if reordering is correct
*
         LASTSL = .TRUE.
         SDIM = 0
         DO 20 I = 1, N
            CURSL = SELCTG( ALPHA( I ), BETA( I ) )
            IF( CURSL )
     $         SDIM = SDIM + 1
            IF( CURSL .AND. .NOT.LASTSL )
     $         INFO = N + 2
            LASTSL = CURSL
   20    CONTINUE
*
      END IF
*
   30 CONTINUE
*
      WORK( 1 ) = LWKOPT
*
      RETURN
*
*     End of ZGGES
*
      END
*> \brief \b ZGESC2 solves a system of linear equations using the LU factorization with complete pivoting computed by sgetc2.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZGESC2 + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgesc2.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgesc2.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgesc2.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGESC2( N, A, LDA, RHS, IPIV, JPIV, SCALE )
* 
*       .. Scalar Arguments ..
*       INTEGER            LDA, N
*       DOUBLE PRECISION   SCALE
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * ), JPIV( * )
*       COMPLEX*16         A( LDA, * ), RHS( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGESC2 solves a system of linear equations
*>
*>           A * X = scale* RHS
*>
*> with a general N-by-N matrix A using the LU factorization with
*> complete pivoting computed by ZGETC2.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA, N)
*>          On entry, the  LU part of the factorization of the n-by-n
*>          matrix A computed by ZGETC2:  A = P * L * U * Q
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1, N).
*> \endverbatim
*>
*> \param[in,out] RHS
*> \verbatim
*>          RHS is COMPLEX*16 array, dimension N.
*>          On entry, the right hand side vector b.
*>          On exit, the solution vector X.
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N).
*>          The pivot indices; for 1 <= i <= N, row i of the
*>          matrix has been interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[in] JPIV
*> \verbatim
*>          JPIV is INTEGER array, dimension (N).
*>          The pivot indices; for 1 <= j <= N, column j of the
*>          matrix has been interchanged with column JPIV(j).
*> \endverbatim
*>
*> \param[out] SCALE
*> \verbatim
*>          SCALE is DOUBLE PRECISION
*>           On exit, SCALE contains the scale factor. SCALE is chosen
*>           0 <= SCALE <= 1 to prevent owerflow in the solution.
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
*> \ingroup complex16GEauxiliary
*
*> \par Contributors:
*  ==================
*>
*>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
*>     Umea University, S-901 87 Umea, Sweden.
*
*  =====================================================================
      SUBROUTINE ZGESC2( N, A, LDA, RHS, IPIV, JPIV, SCALE )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      INTEGER            LDA, N
      DOUBLE PRECISION   SCALE
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * ), JPIV( * )
      COMPLEX*16         A( LDA, * ), RHS( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, TWO = 2.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   BIGNUM, EPS, SMLNUM
      COMPLEX*16         TEMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZLASWP, ZSCAL
*     ..
*     .. External Functions ..
      INTEGER            IZAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           IZAMAX, DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX
*     ..
*     .. Executable Statements ..
*
*     Set constant to control overflow
*
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
*
*     Apply permutations IPIV to RHS
*
      CALL ZLASWP( 1, RHS, LDA, 1, N-1, IPIV, 1 )
*
*     Solve for L part
*
      DO 20 I = 1, N - 1
         DO 10 J = I + 1, N
            RHS( J ) = RHS( J ) - A( J, I )*RHS( I )
   10    CONTINUE
   20 CONTINUE
*
*     Solve for U part
*
      SCALE = ONE
*
*     Check for scaling
*
      I = IZAMAX( N, RHS, 1 )
      IF( TWO*SMLNUM*ABS( RHS( I ) ).GT.ABS( A( N, N ) ) ) THEN
         TEMP = DCMPLX( ONE / TWO, ZERO ) / ABS( RHS( I ) )
         CALL ZSCAL( N, TEMP, RHS( 1 ), 1 )
         SCALE = SCALE*DBLE( TEMP )
      END IF
      DO 40 I = N, 1, -1
         TEMP = DCMPLX( ONE, ZERO ) / A( I, I )
         RHS( I ) = RHS( I )*TEMP
         DO 30 J = I + 1, N
            RHS( I ) = RHS( I ) - RHS( J )*( A( I, J )*TEMP )
   30    CONTINUE
   40 CONTINUE
*
*     Apply permutations JPIV to the solution (RHS)
*
      CALL ZLASWP( 1, RHS, LDA, 1, N-1, JPIV, -1 )
      RETURN
*
*     End of ZGESC2
*
      END
*> \brief \b ZGETC2 computes the LU factorization with complete pivoting of the general n-by-n matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZGETC2 + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgetc2.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgetc2.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgetc2.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGETC2( N, A, LDA, IPIV, JPIV, INFO )
* 
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * ), JPIV( * )
*       COMPLEX*16         A( LDA, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGETC2 computes an LU factorization, using complete pivoting, of the
*> n-by-n matrix A. The factorization has the form A = P * L * U * Q,
*> where P and Q are permutation matrices, L is lower triangular with
*> unit diagonal elements and U is upper triangular.
*>
*> This is a level 1 BLAS version of the algorithm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A. N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA, N)
*>          On entry, the n-by-n matrix to be factored.
*>          On exit, the factors L and U from the factorization
*>          A = P*L*U*Q; the unit diagonal elements of L are not stored.
*>          If U(k, k) appears to be less than SMIN, U(k, k) is given the
*>          value of SMIN, giving a nonsingular perturbed system.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1, N).
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N).
*>          The pivot indices; for 1 <= i <= N, row i of the
*>          matrix has been interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[out] JPIV
*> \verbatim
*>          JPIV is INTEGER array, dimension (N).
*>          The pivot indices; for 1 <= j <= N, column j of the
*>          matrix has been interchanged with column JPIV(j).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>           = 0: successful exit
*>           > 0: if INFO = k, U(k, k) is likely to produce overflow if
*>                one tries to solve for x in Ax = b. So U is perturbed
*>                to avoid the overflow.
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
*> \date November 2013
*
*> \ingroup complex16GEauxiliary
*
*> \par Contributors:
*  ==================
*>
*>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
*>     Umea University, S-901 87 Umea, Sweden.
*
*  =====================================================================
      SUBROUTINE ZGETC2( N, A, LDA, IPIV, JPIV, INFO )
*
*  -- LAPACK auxiliary routine (version 3.5.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2013
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * ), JPIV( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IP, IPV, J, JP, JPV
      DOUBLE PRECISION   BIGNUM, EPS, SMIN, SMLNUM, XMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZGERU, ZSWAP
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DCMPLX, MAX
*     ..
*     .. Executable Statements ..
*
*     Set constants to control overflow
*
      INFO = 0
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
*
*     Factorize A using complete pivoting.
*     Set pivots less than SMIN to SMIN
*
      DO 40 I = 1, N - 1
*
*        Find max element in matrix A
*
         XMAX = ZERO
         DO 20 IP = I, N
            DO 10 JP = I, N
               IF( ABS( A( IP, JP ) ).GE.XMAX ) THEN
                  XMAX = ABS( A( IP, JP ) )
                  IPV = IP
                  JPV = JP
               END IF
   10       CONTINUE
   20    CONTINUE
         IF( I.EQ.1 )
     $      SMIN = MAX( EPS*XMAX, SMLNUM )
*
*        Swap rows
*
         IF( IPV.NE.I )
     $      CALL ZSWAP( N, A( IPV, 1 ), LDA, A( I, 1 ), LDA )
         IPIV( I ) = IPV
*
*        Swap columns
*
         IF( JPV.NE.I )
     $      CALL ZSWAP( N, A( 1, JPV ), 1, A( 1, I ), 1 )
         JPIV( I ) = JPV
*
*        Check for singularity
*
         IF( ABS( A( I, I ) ).LT.SMIN ) THEN
            INFO = I
            A( I, I ) = DCMPLX( SMIN, ZERO )
         END IF
         DO 30 J = I + 1, N
            A( J, I ) = A( J, I ) / A( I, I )
   30    CONTINUE
         CALL ZGERU( N-I, N-I, -DCMPLX( ONE ), A( I+1, I ), 1,
     $               A( I, I+1 ), LDA, A( I+1, I+1 ), LDA )
   40 CONTINUE
*
      IF( ABS( A( N, N ) ).LT.SMIN ) THEN
         INFO = N
         A( N, N ) = DCMPLX( SMIN, ZERO )
      END IF
*
*     Set last pivots to N
*
      IPIV( N ) = N
      JPIV( N ) = N
*
      RETURN
*
*     End of ZGETC2
*
      END
*> \brief \b ZLATDF uses the LU factorization of the n-by-n matrix computed by sgetc2 and computes a contribution to the reciprocal Dif-estimate.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZLATDF + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlatdf.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlatdf.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlatdf.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLATDF( IJOB, N, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV,
*                          JPIV )
* 
*       .. Scalar Arguments ..
*       INTEGER            IJOB, LDZ, N
*       DOUBLE PRECISION   RDSCAL, RDSUM
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * ), JPIV( * )
*       COMPLEX*16         RHS( * ), Z( LDZ, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLATDF computes the contribution to the reciprocal Dif-estimate
*> by solving for x in Z * x = b, where b is chosen such that the norm
*> of x is as large as possible. It is assumed that LU decomposition
*> of Z has been computed by ZGETC2. On entry RHS = f holds the
*> contribution from earlier solved sub-systems, and on return RHS = x.
*>
*> The factorization of Z returned by ZGETC2 has the form
*> Z = P * L * U * Q, where P and Q are permutation matrices. L is lower
*> triangular with unit diagonal elements and U is upper triangular.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] IJOB
*> \verbatim
*>          IJOB is INTEGER
*>          IJOB = 2: First compute an approximative null-vector e
*>              of Z using ZGECON, e is normalized and solve for
*>              Zx = +-e - f with the sign giving the greater value of
*>              2-norm(x).  About 5 times as expensive as Default.
*>          IJOB .ne. 2: Local look ahead strategy where
*>              all entries of the r.h.s. b is choosen as either +1 or
*>              -1.  Default.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix Z.
*> \endverbatim
*>
*> \param[in] Z
*> \verbatim
*>          Z is DOUBLE PRECISION array, dimension (LDZ, N)
*>          On entry, the LU part of the factorization of the n-by-n
*>          matrix Z computed by ZGETC2:  Z = P * L * U * Q
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z.  LDA >= max(1, N).
*> \endverbatim
*>
*> \param[in,out] RHS
*> \verbatim
*>          RHS is DOUBLE PRECISION array, dimension (N).
*>          On entry, RHS contains contributions from other subsystems.
*>          On exit, RHS contains the solution of the subsystem with
*>          entries according to the value of IJOB (see above).
*> \endverbatim
*>
*> \param[in,out] RDSUM
*> \verbatim
*>          RDSUM is DOUBLE PRECISION
*>          On entry, the sum of squares of computed contributions to
*>          the Dif-estimate under computation by ZTGSYL, where the
*>          scaling factor RDSCAL (see below) has been factored out.
*>          On exit, the corresponding sum of squares updated with the
*>          contributions from the current sub-system.
*>          If TRANS = 'T' RDSUM is not touched.
*>          NOTE: RDSUM only makes sense when ZTGSY2 is called by CTGSYL.
*> \endverbatim
*>
*> \param[in,out] RDSCAL
*> \verbatim
*>          RDSCAL is DOUBLE PRECISION
*>          On entry, scaling factor used to prevent overflow in RDSUM.
*>          On exit, RDSCAL is updated w.r.t. the current contributions
*>          in RDSUM.
*>          If TRANS = 'T', RDSCAL is not touched.
*>          NOTE: RDSCAL only makes sense when ZTGSY2 is called by
*>          ZTGSYL.
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N).
*>          The pivot indices; for 1 <= i <= N, row i of the
*>          matrix has been interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[in] JPIV
*> \verbatim
*>          JPIV is INTEGER array, dimension (N).
*>          The pivot indices; for 1 <= j <= N, column j of the
*>          matrix has been interchanged with column JPIV(j).
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
*> \par Further Details:
*  =====================
*>
*>  This routine is a further developed implementation of algorithm
*>  BSOLVE in [1] using complete pivoting in the LU factorization.
*
*> \par Contributors:
*  ==================
*>
*>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
*>     Umea University, S-901 87 Umea, Sweden.
*
*> \par References:
*  ================
*>
*>   [1]   Bo Kagstrom and Lars Westin,
*>         Generalized Schur Methods with Condition Estimators for
*>         Solving the Generalized Sylvester Equation, IEEE Transactions
*>         on Automatic Control, Vol. 34, No. 7, July 1989, pp 745-751.
*>\n
*>   [2]   Peter Poromaa,
*>         On Efficient and Robust Estimators for the Separation
*>         between two Regular Matrix Pairs with Applications in
*>         Condition Estimation. Report UMINF-95.05, Department of
*>         Computing Science, Umea University, S-901 87 Umea, Sweden,
*>         1995.
*
*  =====================================================================
      SUBROUTINE ZLATDF( IJOB, N, Z, LDZ, RHS, RDSUM, RDSCAL, IPIV,
     $                   JPIV )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      INTEGER            IJOB, LDZ, N
      DOUBLE PRECISION   RDSCAL, RDSUM
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * ), JPIV( * )
      COMPLEX*16         RHS( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            MAXDIM
      PARAMETER          ( MAXDIM = 2 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J, K
      DOUBLE PRECISION   RTEMP, SCALE, SMINU, SPLUS
      COMPLEX*16         BM, BP, PMONE, TEMP
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   RWORK( MAXDIM )
      COMPLEX*16         WORK( 4*MAXDIM ), XM( MAXDIM ), XP( MAXDIM )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZAXPY, ZCOPY, ZGECON, ZGESC2, ZLASSQ, ZLASWP,
     $                   ZSCAL
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DZASUM
      COMPLEX*16         ZDOTC
      EXTERNAL           DZASUM, ZDOTC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( IJOB.NE.2 ) THEN
*
*        Apply permutations IPIV to RHS
*
         CALL ZLASWP( 1, RHS, LDZ, 1, N-1, IPIV, 1 )
*
*        Solve for L-part choosing RHS either to +1 or -1.
*
         PMONE = -CONE
         DO 10 J = 1, N - 1
            BP = RHS( J ) + CONE
            BM = RHS( J ) - CONE
            SPLUS = ONE
*
*           Lockahead for L- part RHS(1:N-1) = +-1
*           SPLUS and SMIN computed more efficiently than in BSOLVE[1].
*
            SPLUS = SPLUS + DBLE( ZDOTC( N-J, Z( J+1, J ), 1, Z( J+1,
     $              J ), 1 ) )
            SMINU = DBLE( ZDOTC( N-J, Z( J+1, J ), 1, RHS( J+1 ), 1 ) )
            SPLUS = SPLUS*DBLE( RHS( J ) )
            IF( SPLUS.GT.SMINU ) THEN
               RHS( J ) = BP
            ELSE IF( SMINU.GT.SPLUS ) THEN
               RHS( J ) = BM
            ELSE
*
*              In this case the updating sums are equal and we can
*              choose RHS(J) +1 or -1. The first time this happens we
*              choose -1, thereafter +1. This is a simple way to get
*              good estimates of matrices like Byers well-known example
*              (see [1]). (Not done in BSOLVE.)
*
               RHS( J ) = RHS( J ) + PMONE
               PMONE = CONE
            END IF
*
*           Compute the remaining r.h.s.
*
            TEMP = -RHS( J )
            CALL ZAXPY( N-J, TEMP, Z( J+1, J ), 1, RHS( J+1 ), 1 )
   10    CONTINUE
*
*        Solve for U- part, lockahead for RHS(N) = +-1. This is not done
*        In BSOLVE and will hopefully give us a better estimate because
*        any ill-conditioning of the original matrix is transfered to U
*        and not to L. U(N, N) is an approximation to sigma_min(LU).
*
         CALL ZCOPY( N-1, RHS, 1, WORK, 1 )
         WORK( N ) = RHS( N ) + CONE
         RHS( N ) = RHS( N ) - CONE
         SPLUS = ZERO
         SMINU = ZERO
         DO 30 I = N, 1, -1
            TEMP = CONE / Z( I, I )
            WORK( I ) = WORK( I )*TEMP
            RHS( I ) = RHS( I )*TEMP
            DO 20 K = I + 1, N
               WORK( I ) = WORK( I ) - WORK( K )*( Z( I, K )*TEMP )
               RHS( I ) = RHS( I ) - RHS( K )*( Z( I, K )*TEMP )
   20       CONTINUE
            SPLUS = SPLUS + ABS( WORK( I ) )
            SMINU = SMINU + ABS( RHS( I ) )
   30    CONTINUE
         IF( SPLUS.GT.SMINU )
     $      CALL ZCOPY( N, WORK, 1, RHS, 1 )
*
*        Apply the permutations JPIV to the computed solution (RHS)
*
         CALL ZLASWP( 1, RHS, LDZ, 1, N-1, JPIV, -1 )
*
*        Compute the sum of squares
*
         CALL ZLASSQ( N, RHS, 1, RDSCAL, RDSUM )
         RETURN
      END IF
*
*     ENTRY IJOB = 2
*
*     Compute approximate nullvector XM of Z
*
      CALL ZGECON( 'I', N, Z, LDZ, ONE, RTEMP, WORK, RWORK, INFO )
      CALL ZCOPY( N, WORK( N+1 ), 1, XM, 1 )
*
*     Compute RHS
*
      CALL ZLASWP( 1, XM, LDZ, 1, N-1, IPIV, -1 )
      TEMP = CONE / SQRT( ZDOTC( N, XM, 1, XM, 1 ) )
      CALL ZSCAL( N, TEMP, XM, 1 )
      CALL ZCOPY( N, XM, 1, XP, 1 )
      CALL ZAXPY( N, CONE, RHS, 1, XP, 1 )
      CALL ZAXPY( N, -CONE, XM, 1, RHS, 1 )
      CALL ZGESC2( N, Z, LDZ, RHS, IPIV, JPIV, SCALE )
      CALL ZGESC2( N, Z, LDZ, XP, IPIV, JPIV, SCALE )
      IF( DZASUM( N, XP, 1 ).GT.DZASUM( N, RHS, 1 ) )
     $   CALL ZCOPY( N, XP, 1, RHS, 1 )
*
*     Compute the sum of squares
*
      CALL ZLASSQ( N, RHS, 1, RDSCAL, RDSUM )
      RETURN
*
*     End of ZLATDF
*
      END
*> \brief \b ZTGEX2 swaps adjacent diagonal blocks in an upper (quasi) triangular matrix pair by an unitary equivalence transformation.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZTGEX2 + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztgex2.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztgex2.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztgex2.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZTGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z,
*                          LDZ, J1, INFO )
* 
*       .. Scalar Arguments ..
*       LOGICAL            WANTQ, WANTZ
*       INTEGER            INFO, J1, LDA, LDB, LDQ, LDZ, N
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
*> ZTGEX2 swaps adjacent diagonal 1 by 1 blocks (A11,B11) and (A22,B22)
*> in an upper triangular matrix pair (A, B) by an unitary equivalence
*> transformation.
*>
*> (A, B) must be in generalized Schur canonical form, that is, A and
*> B are both upper triangular.
*>
*> Optionally, the matrices Q and Z of generalized Schur vectors are
*> updated.
*>
*>        Q(in) * A(in) * Z(in)**H = Q(out) * A(out) * Z(out)**H
*>        Q(in) * B(in) * Z(in)**H = Q(out) * B(out) * Z(out)**H
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] WANTQ
*> \verbatim
*>          WANTQ is LOGICAL
*>          .TRUE. : update the left transformation matrix Q;
*>          .FALSE.: do not update Q.
*> \endverbatim
*>
*> \param[in] WANTZ
*> \verbatim
*>          WANTZ is LOGICAL
*>          .TRUE. : update the right transformation matrix Z;
*>          .FALSE.: do not update Z.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices A and B. N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 arrays, dimensions (LDA,N)
*>          On entry, the matrix A in the pair (A, B).
*>          On exit, the updated matrix A.
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
*>          B is COMPLEX*16 arrays, dimensions (LDB,N)
*>          On entry, the matrix B in the pair (A, B).
*>          On exit, the updated matrix B.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B. LDB >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] Q
*> \verbatim
*>          Q is COMPLEX*16 array, dimension (LDZ,N)
*>          If WANTQ = .TRUE, on entry, the unitary matrix Q. On exit,
*>          the updated matrix Q.
*>          Not referenced if WANTQ = .FALSE..
*> \endverbatim
*>
*> \param[in] LDQ
*> \verbatim
*>          LDQ is INTEGER
*>          The leading dimension of the array Q. LDQ >= 1;
*>          If WANTQ = .TRUE., LDQ >= N.
*> \endverbatim
*>
*> \param[in,out] Z
*> \verbatim
*>          Z is COMPLEX*16 array, dimension (LDZ,N)
*>          If WANTZ = .TRUE, on entry, the unitary matrix Z. On exit,
*>          the updated matrix Z.
*>          Not referenced if WANTZ = .FALSE..
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z. LDZ >= 1;
*>          If WANTZ = .TRUE., LDZ >= N.
*> \endverbatim
*>
*> \param[in] J1
*> \verbatim
*>          J1 is INTEGER
*>          The index to the first block (A11, B11).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>           =0:  Successful exit.
*>           =1:  The transformed matrix pair (A, B) would be too far
*>                from generalized Schur form; the problem is ill-
*>                conditioned. 
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
*> \ingroup complex16GEauxiliary
*
*> \par Further Details:
*  =====================
*>
*>  In the current code both weak and strong stability tests are
*>  performed. The user can omit the strong stability test by changing
*>  the internal logical parameter WANDS to .FALSE.. See ref. [2] for
*>  details.
*
*> \par Contributors:
*  ==================
*>
*>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
*>     Umea University, S-901 87 Umea, Sweden.
*
*> \par References:
*  ================
*>
*>  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the
*>      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in
*>      M.S. Moonen et al (eds), Linear Algebra for Large Scale and
*>      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.
*> \n
*>  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified
*>      Eigenvalues of a Regular Matrix Pair (A, B) and Condition
*>      Estimation: Theory, Algorithms and Software, Report UMINF-94.04,
*>      Department of Computing Science, Umea University, S-901 87 Umea,
*>      Sweden, 1994. Also as LAPACK Working Note 87. To appear in
*>      Numerical Algorithms, 1996.
*>
*  =====================================================================
      SUBROUTINE ZTGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z,
     $                   LDZ, J1, INFO )
*
*  -- LAPACK auxiliary routine (version 3.4.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     September 2012
*
*     .. Scalar Arguments ..
      LOGICAL            WANTQ, WANTZ
      INTEGER            INFO, J1, LDA, LDB, LDQ, LDZ, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
     $                   Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
      DOUBLE PRECISION   TWENTY
      PARAMETER          ( TWENTY = 2.0D+1 )
      INTEGER            LDST
      PARAMETER          ( LDST = 2 )
      LOGICAL            WANDS
      PARAMETER          ( WANDS = .TRUE. )
*     ..
*     .. Local Scalars ..
      LOGICAL            DTRONG, WEAK
      INTEGER            I, M
      DOUBLE PRECISION   CQ, CZ, EPS, SA, SB, SCALE, SMLNUM, SS, SUM,
     $                   THRESH, WS
      COMPLEX*16         CDUM, F, G, SQ, SZ
*     ..
*     .. Local Arrays ..
      COMPLEX*16         S( LDST, LDST ), T( LDST, LDST ), WORK( 8 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZLACPY, ZLARTG, ZLASSQ, ZROT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCONJG, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
      INFO = 0
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
*
      M = LDST
      WEAK = .FALSE.
      DTRONG = .FALSE.
*
*     Make a local copy of selected block in (A, B)
*
      CALL ZLACPY( 'Full', M, M, A( J1, J1 ), LDA, S, LDST )
      CALL ZLACPY( 'Full', M, M, B( J1, J1 ), LDB, T, LDST )
*
*     Compute the threshold for testing the acceptance of swapping.
*
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' ) / EPS
      SCALE = DBLE( CZERO )
      SUM = DBLE( CONE )
      CALL ZLACPY( 'Full', M, M, S, LDST, WORK, M )
      CALL ZLACPY( 'Full', M, M, T, LDST, WORK( M*M+1 ), M )
      CALL ZLASSQ( 2*M*M, WORK, 1, SCALE, SUM )
      SA = SCALE*SQRT( SUM )
*
*     THRES has been changed from 
*        THRESH = MAX( TEN*EPS*SA, SMLNUM )
*     to
*        THRESH = MAX( TWENTY*EPS*SA, SMLNUM )
*     on 04/01/10.
*     "Bug" reported by Ondra Kamenik, confirmed by Julie Langou, fixed by
*     Jim Demmel and Guillaume Revy. See forum post 1783.
*
      THRESH = MAX( TWENTY*EPS*SA, SMLNUM )
*
*     Compute unitary QL and RQ that swap 1-by-1 and 1-by-1 blocks
*     using Givens rotations and perform the swap tentatively.
*
      F = S( 2, 2 )*T( 1, 1 ) - T( 2, 2 )*S( 1, 1 )
      G = S( 2, 2 )*T( 1, 2 ) - T( 2, 2 )*S( 1, 2 )
      SA = ABS( S( 2, 2 ) )
      SB = ABS( T( 2, 2 ) )
      CALL ZLARTG( G, F, CZ, SZ, CDUM )
      SZ = -SZ
      CALL ZROT( 2, S( 1, 1 ), 1, S( 1, 2 ), 1, CZ, DCONJG( SZ ) )
      CALL ZROT( 2, T( 1, 1 ), 1, T( 1, 2 ), 1, CZ, DCONJG( SZ ) )
      IF( SA.GE.SB ) THEN
         CALL ZLARTG( S( 1, 1 ), S( 2, 1 ), CQ, SQ, CDUM )
      ELSE
         CALL ZLARTG( T( 1, 1 ), T( 2, 1 ), CQ, SQ, CDUM )
      END IF
      CALL ZROT( 2, S( 1, 1 ), LDST, S( 2, 1 ), LDST, CQ, SQ )
      CALL ZROT( 2, T( 1, 1 ), LDST, T( 2, 1 ), LDST, CQ, SQ )
*
*     Weak stability test: |S21| + |T21| <= O(EPS F-norm((S, T)))
*
      WS = ABS( S( 2, 1 ) ) + ABS( T( 2, 1 ) )
      WEAK = WS.LE.THRESH
      IF( .NOT.WEAK )
     $   GO TO 20
*
      IF( WANDS ) THEN
*
*        Strong stability test:
*           F-norm((A-QL**H*S*QR, B-QL**H*T*QR)) <= O(EPS*F-norm((A, B)))
*
         CALL ZLACPY( 'Full', M, M, S, LDST, WORK, M )
         CALL ZLACPY( 'Full', M, M, T, LDST, WORK( M*M+1 ), M )
         CALL ZROT( 2, WORK, 1, WORK( 3 ), 1, CZ, -DCONJG( SZ ) )
         CALL ZROT( 2, WORK( 5 ), 1, WORK( 7 ), 1, CZ, -DCONJG( SZ ) )
         CALL ZROT( 2, WORK, 2, WORK( 2 ), 2, CQ, -SQ )
         CALL ZROT( 2, WORK( 5 ), 2, WORK( 6 ), 2, CQ, -SQ )
         DO 10 I = 1, 2
            WORK( I ) = WORK( I ) - A( J1+I-1, J1 )
            WORK( I+2 ) = WORK( I+2 ) - A( J1+I-1, J1+1 )
            WORK( I+4 ) = WORK( I+4 ) - B( J1+I-1, J1 )
            WORK( I+6 ) = WORK( I+6 ) - B( J1+I-1, J1+1 )
   10    CONTINUE
         SCALE = DBLE( CZERO )
         SUM = DBLE( CONE )
         CALL ZLASSQ( 2*M*M, WORK, 1, SCALE, SUM )
         SS = SCALE*SQRT( SUM )
         DTRONG = SS.LE.THRESH
         IF( .NOT.DTRONG )
     $      GO TO 20
      END IF
*
*     If the swap is accepted ("weakly" and "strongly"), apply the
*     equivalence transformations to the original matrix pair (A,B)
*
      CALL ZROT( J1+1, A( 1, J1 ), 1, A( 1, J1+1 ), 1, CZ,
     $           DCONJG( SZ ) )
      CALL ZROT( J1+1, B( 1, J1 ), 1, B( 1, J1+1 ), 1, CZ,
     $           DCONJG( SZ ) )
      CALL ZROT( N-J1+1, A( J1, J1 ), LDA, A( J1+1, J1 ), LDA, CQ, SQ )
      CALL ZROT( N-J1+1, B( J1, J1 ), LDB, B( J1+1, J1 ), LDB, CQ, SQ )
*
*     Set  N1 by N2 (2,1) blocks to 0
*
      A( J1+1, J1 ) = CZERO
      B( J1+1, J1 ) = CZERO
*
*     Accumulate transformations into Q and Z if requested.
*
      IF( WANTZ )
     $   CALL ZROT( N, Z( 1, J1 ), 1, Z( 1, J1+1 ), 1, CZ,
     $              DCONJG( SZ ) )
      IF( WANTQ )
     $   CALL ZROT( N, Q( 1, J1 ), 1, Q( 1, J1+1 ), 1, CQ,
     $              DCONJG( SQ ) )
*
*     Exit with INFO = 0 if swap was successfully performed.
*
      RETURN
*
*     Exit with INFO = 1 if swap was rejected.
*
   20 CONTINUE
      INFO = 1
      RETURN
*
*     End of ZTGEX2
*
      END
*> \brief \b ZTGEXC
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZTGEXC + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztgexc.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztgexc.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztgexc.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZTGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z,
*                          LDZ, IFST, ILST, INFO )
* 
*       .. Scalar Arguments ..
*       LOGICAL            WANTQ, WANTZ
*       INTEGER            IFST, ILST, INFO, LDA, LDB, LDQ, LDZ, N
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
*> ZTGEXC reorders the generalized Schur decomposition of a complex
*> matrix pair (A,B), using an unitary equivalence transformation
*> (A, B) := Q * (A, B) * Z**H, so that the diagonal block of (A, B) with
*> row index IFST is moved to row ILST.
*>
*> (A, B) must be in generalized Schur canonical form, that is, A and
*> B are both upper triangular.
*>
*> Optionally, the matrices Q and Z of generalized Schur vectors are
*> updated.
*>
*>        Q(in) * A(in) * Z(in)**H = Q(out) * A(out) * Z(out)**H
*>        Q(in) * B(in) * Z(in)**H = Q(out) * B(out) * Z(out)**H
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] WANTQ
*> \verbatim
*>          WANTQ is LOGICAL
*>          .TRUE. : update the left transformation matrix Q;
*>          .FALSE.: do not update Q.
*> \endverbatim
*>
*> \param[in] WANTZ
*> \verbatim
*>          WANTZ is LOGICAL
*>          .TRUE. : update the right transformation matrix Z;
*>          .FALSE.: do not update Z.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices A and B. N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the upper triangular matrix A in the pair (A, B).
*>          On exit, the updated matrix A.
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
*>          On entry, the upper triangular matrix B in the pair (A, B).
*>          On exit, the updated matrix B.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B. LDB >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] Q
*> \verbatim
*>          Q is COMPLEX*16 array, dimension (LDZ,N)
*>          On entry, if WANTQ = .TRUE., the unitary matrix Q.
*>          On exit, the updated matrix Q.
*>          If WANTQ = .FALSE., Q is not referenced.
*> \endverbatim
*>
*> \param[in] LDQ
*> \verbatim
*>          LDQ is INTEGER
*>          The leading dimension of the array Q. LDQ >= 1;
*>          If WANTQ = .TRUE., LDQ >= N.
*> \endverbatim
*>
*> \param[in,out] Z
*> \verbatim
*>          Z is COMPLEX*16 array, dimension (LDZ,N)
*>          On entry, if WANTZ = .TRUE., the unitary matrix Z.
*>          On exit, the updated matrix Z.
*>          If WANTZ = .FALSE., Z is not referenced.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z. LDZ >= 1;
*>          If WANTZ = .TRUE., LDZ >= N.
*> \endverbatim
*>
*> \param[in] IFST
*> \verbatim
*>          IFST is INTEGER
*> \endverbatim
*>
*> \param[in,out] ILST
*> \verbatim
*>          ILST is INTEGER
*>          Specify the reordering of the diagonal blocks of (A, B).
*>          The block with row index IFST is moved to row ILST, by a
*>          sequence of swapping between adjacent blocks.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>           =0:  Successful exit.
*>           <0:  if INFO = -i, the i-th argument had an illegal value.
*>           =1:  The transformed matrix pair (A, B) would be too far
*>                from generalized Schur form; the problem is ill-
*>                conditioned. (A, B) may have been partially reordered,
*>                and ILST points to the first row of the current
*>                position of the block being moved.
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
*> \par Contributors:
*  ==================
*>
*>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
*>     Umea University, S-901 87 Umea, Sweden.
*
*> \par References:
*  ================
*>
*>  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the
*>      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in
*>      M.S. Moonen et al (eds), Linear Algebra for Large Scale and
*>      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.
*> \n
*>  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified
*>      Eigenvalues of a Regular Matrix Pair (A, B) and Condition
*>      Estimation: Theory, Algorithms and Software, Report
*>      UMINF - 94.04, Department of Computing Science, Umea University,
*>      S-901 87 Umea, Sweden, 1994. Also as LAPACK Working Note 87.
*>      To appear in Numerical Algorithms, 1996.
*> \n
*>  [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software
*>      for Solving the Generalized Sylvester Equation and Estimating the
*>      Separation between Regular Matrix Pairs, Report UMINF - 93.23,
*>      Department of Computing Science, Umea University, S-901 87 Umea,
*>      Sweden, December 1993, Revised April 1994, Also as LAPACK working
*>      Note 75. To appear in ACM Trans. on Math. Software, Vol 22, No 1,
*>      1996.
*>
*  =====================================================================
      SUBROUTINE ZTGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z,
     $                   LDZ, IFST, ILST, INFO )
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      LOGICAL            WANTQ, WANTZ
      INTEGER            IFST, ILST, INFO, LDA, LDB, LDQ, LDZ, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
     $                   Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            HERE
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZTGEX2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Decode and test input arguments.
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDQ.LT.1 .OR. WANTQ .AND. ( LDQ.LT.MAX( 1, N ) ) ) THEN
         INFO = -9
      ELSE IF( LDZ.LT.1 .OR. WANTZ .AND. ( LDZ.LT.MAX( 1, N ) ) ) THEN
         INFO = -11
      ELSE IF( IFST.LT.1 .OR. IFST.GT.N ) THEN
         INFO = -12
      ELSE IF( ILST.LT.1 .OR. ILST.GT.N ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTGEXC', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
      IF( IFST.EQ.ILST )
     $   RETURN
*
      IF( IFST.LT.ILST ) THEN
*
         HERE = IFST
*
   10    CONTINUE
*
*        Swap with next one below
*
         CALL ZTGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ,
     $                HERE, INFO )
         IF( INFO.NE.0 ) THEN
            ILST = HERE
            RETURN
         END IF
         HERE = HERE + 1
         IF( HERE.LT.ILST )
     $      GO TO 10
         HERE = HERE - 1
      ELSE
         HERE = IFST - 1
*
   20    CONTINUE
*
*        Swap with next one above
*
         CALL ZTGEX2( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z, LDZ,
     $                HERE, INFO )
         IF( INFO.NE.0 ) THEN
            ILST = HERE
            RETURN
         END IF
         HERE = HERE - 1
         IF( HERE.GE.ILST )
     $      GO TO 20
         HERE = HERE + 1
      END IF
      ILST = HERE
      RETURN
*
*     End of ZTGEXC
*
      END
*> \brief \b ZTGSEN
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZTGSEN + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztgsen.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztgsen.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztgsen.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZTGSEN( IJOB, WANTQ, WANTZ, SELECT, N, A, LDA, B, LDB,
*                          ALPHA, BETA, Q, LDQ, Z, LDZ, M, PL, PR, DIF,
*                          WORK, LWORK, IWORK, LIWORK, INFO )
* 
*       .. Scalar Arguments ..
*       LOGICAL            WANTQ, WANTZ
*       INTEGER            IJOB, INFO, LDA, LDB, LDQ, LDZ, LIWORK, LWORK,
*      $                   M, N
*       DOUBLE PRECISION   PL, PR
*       ..
*       .. Array Arguments ..
*       LOGICAL            SELECT( * )
*       INTEGER            IWORK( * )
*       DOUBLE PRECISION   DIF( * )
*       COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),
*      $                   BETA( * ), Q( LDQ, * ), WORK( * ), Z( LDZ, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZTGSEN reorders the generalized Schur decomposition of a complex
*> matrix pair (A, B) (in terms of an unitary equivalence trans-
*> formation Q**H * (A, B) * Z), so that a selected cluster of eigenvalues
*> appears in the leading diagonal blocks of the pair (A,B). The leading
*> columns of Q and Z form unitary bases of the corresponding left and
*> right eigenspaces (deflating subspaces). (A, B) must be in
*> generalized Schur canonical form, that is, A and B are both upper
*> triangular.
*>
*> ZTGSEN also computes the generalized eigenvalues
*>
*>          w(j)= ALPHA(j) / BETA(j)
*>
*> of the reordered matrix pair (A, B).
*>
*> Optionally, the routine computes estimates of reciprocal condition
*> numbers for eigenvalues and eigenspaces. These are Difu[(A11,B11),
*> (A22,B22)] and Difl[(A11,B11), (A22,B22)], i.e. the separation(s)
*> between the matrix pairs (A11, B11) and (A22,B22) that correspond to
*> the selected cluster and the eigenvalues outside the cluster, resp.,
*> and norms of "projections" onto left and right eigenspaces w.r.t.
*> the selected cluster in the (1,1)-block.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] IJOB
*> \verbatim
*>          IJOB is integer
*>          Specifies whether condition numbers are required for the
*>          cluster of eigenvalues (PL and PR) or the deflating subspaces
*>          (Difu and Difl):
*>           =0: Only reorder w.r.t. SELECT. No extras.
*>           =1: Reciprocal of norms of "projections" onto left and right
*>               eigenspaces w.r.t. the selected cluster (PL and PR).
*>           =2: Upper bounds on Difu and Difl. F-norm-based estimate
*>               (DIF(1:2)).
*>           =3: Estimate of Difu and Difl. 1-norm-based estimate
*>               (DIF(1:2)).
*>               About 5 times as expensive as IJOB = 2.
*>           =4: Compute PL, PR and DIF (i.e. 0, 1 and 2 above): Economic
*>               version to get it all.
*>           =5: Compute PL, PR and DIF (i.e. 0, 1 and 3 above)
*> \endverbatim
*>
*> \param[in] WANTQ
*> \verbatim
*>          WANTQ is LOGICAL
*>          .TRUE. : update the left transformation matrix Q;
*>          .FALSE.: do not update Q.
*> \endverbatim
*>
*> \param[in] WANTZ
*> \verbatim
*>          WANTZ is LOGICAL
*>          .TRUE. : update the right transformation matrix Z;
*>          .FALSE.: do not update Z.
*> \endverbatim
*>
*> \param[in] SELECT
*> \verbatim
*>          SELECT is LOGICAL array, dimension (N)
*>          SELECT specifies the eigenvalues in the selected cluster. To
*>          select an eigenvalue w(j), SELECT(j) must be set to
*>          .TRUE..
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices A and B. N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension(LDA,N)
*>          On entry, the upper triangular matrix A, in generalized
*>          Schur canonical form.
*>          On exit, A is overwritten by the reordered matrix A.
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
*>          B is COMPLEX*16 array, dimension(LDB,N)
*>          On entry, the upper triangular matrix B, in generalized
*>          Schur canonical form.
*>          On exit, B is overwritten by the reordered matrix B.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B. LDB >= max(1,N).
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
*>
*>          The diagonal elements of A and B, respectively,
*>          when the pair (A,B) has been reduced to generalized Schur
*>          form.  ALPHA(i)/BETA(i) i=1,...,N are the generalized
*>          eigenvalues.
*> \endverbatim
*>
*> \param[in,out] Q
*> \verbatim
*>          Q is COMPLEX*16 array, dimension (LDQ,N)
*>          On entry, if WANTQ = .TRUE., Q is an N-by-N matrix.
*>          On exit, Q has been postmultiplied by the left unitary
*>          transformation matrix which reorder (A, B); The leading M
*>          columns of Q form orthonormal bases for the specified pair of
*>          left eigenspaces (deflating subspaces).
*>          If WANTQ = .FALSE., Q is not referenced.
*> \endverbatim
*>
*> \param[in] LDQ
*> \verbatim
*>          LDQ is INTEGER
*>          The leading dimension of the array Q. LDQ >= 1.
*>          If WANTQ = .TRUE., LDQ >= N.
*> \endverbatim
*>
*> \param[in,out] Z
*> \verbatim
*>          Z is COMPLEX*16 array, dimension (LDZ,N)
*>          On entry, if WANTZ = .TRUE., Z is an N-by-N matrix.
*>          On exit, Z has been postmultiplied by the left unitary
*>          transformation matrix which reorder (A, B); The leading M
*>          columns of Z form orthonormal bases for the specified pair of
*>          left eigenspaces (deflating subspaces).
*>          If WANTZ = .FALSE., Z is not referenced.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z. LDZ >= 1.
*>          If WANTZ = .TRUE., LDZ >= N.
*> \endverbatim
*>
*> \param[out] M
*> \verbatim
*>          M is INTEGER
*>          The dimension of the specified pair of left and right
*>          eigenspaces, (deflating subspaces) 0 <= M <= N.
*> \endverbatim
*>
*> \param[out] PL
*> \verbatim
*>          PL is DOUBLE PRECISION
*> \endverbatim
*>
*> \param[out] PR
*> \verbatim
*>          PR is DOUBLE PRECISION
*>
*>          If IJOB = 1, 4 or 5, PL, PR are lower bounds on the
*>          reciprocal  of the norm of "projections" onto left and right
*>          eigenspace with respect to the selected cluster.
*>          0 < PL, PR <= 1.
*>          If M = 0 or M = N, PL = PR  = 1.
*>          If IJOB = 0, 2 or 3 PL, PR are not referenced.
*> \endverbatim
*>
*> \param[out] DIF
*> \verbatim
*>          DIF is DOUBLE PRECISION array, dimension (2).
*>          If IJOB >= 2, DIF(1:2) store the estimates of Difu and Difl.
*>          If IJOB = 2 or 4, DIF(1:2) are F-norm-based upper bounds on
*>          Difu and Difl. If IJOB = 3 or 5, DIF(1:2) are 1-norm-based
*>          estimates of Difu and Difl, computed using reversed
*>          communication with ZLACN2.
*>          If M = 0 or N, DIF(1:2) = F-norm([A, B]).
*>          If IJOB = 0 or 1, DIF is not referenced.
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
*>          The dimension of the array WORK. LWORK >=  1
*>          If IJOB = 1, 2 or 4, LWORK >=  2*M*(N-M)
*>          If IJOB = 3 or 5, LWORK >=  4*M*(N-M)
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
*>          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*> \endverbatim
*>
*> \param[in] LIWORK
*> \verbatim
*>          LIWORK is INTEGER
*>          The dimension of the array IWORK. LIWORK >= 1.
*>          If IJOB = 1, 2 or 4, LIWORK >=  N+2;
*>          If IJOB = 3 or 5, LIWORK >= MAX(N+2, 2*M*(N-M));
*>
*>          If LIWORK = -1, then a workspace query is assumed; the
*>          routine only calculates the optimal size of the IWORK array,
*>          returns this value as the first entry of the IWORK array, and
*>          no error message related to LIWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>            =0: Successful exit.
*>            <0: If INFO = -i, the i-th argument had an illegal value.
*>            =1: Reordering of (A, B) failed because the transformed
*>                matrix pair (A, B) would be too far from generalized
*>                Schur form; the problem is very ill-conditioned.
*>                (A, B) may have been partially reordered.
*>                If requested, 0 is returned in DIF(*), PL and PR.
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
*>  ZTGSEN first collects the selected eigenvalues by computing unitary
*>  U and W that move them to the top left corner of (A, B). In other
*>  words, the selected eigenvalues are the eigenvalues of (A11, B11) in
*>
*>              U**H*(A, B)*W = (A11 A12) (B11 B12) n1
*>                              ( 0  A22),( 0  B22) n2
*>                                n1  n2    n1  n2
*>
*>  where N = n1+n2 and U**H means the conjugate transpose of U. The first
*>  n1 columns of U and W span the specified pair of left and right
*>  eigenspaces (deflating subspaces) of (A, B).
*>
*>  If (A, B) has been obtained from the generalized real Schur
*>  decomposition of a matrix pair (C, D) = Q*(A, B)*Z**H, then the
*>  reordered generalized Schur form of (C, D) is given by
*>
*>           (C, D) = (Q*U)*(U**H *(A, B)*W)*(Z*W)**H,
*>
*>  and the first n1 columns of Q*U and Z*W span the corresponding
*>  deflating subspaces of (C, D) (Q and Z store Q*U and Z*W, resp.).
*>
*>  Note that if the selected eigenvalue is sufficiently ill-conditioned,
*>  then its value may differ significantly from its value before
*>  reordering.
*>
*>  The reciprocal condition numbers of the left and right eigenspaces
*>  spanned by the first n1 columns of U and W (or Q*U and Z*W) may
*>  be returned in DIF(1:2), corresponding to Difu and Difl, resp.
*>
*>  The Difu and Difl are defined as:
*>
*>       Difu[(A11, B11), (A22, B22)] = sigma-min( Zu )
*>  and
*>       Difl[(A11, B11), (A22, B22)] = Difu[(A22, B22), (A11, B11)],
*>
*>  where sigma-min(Zu) is the smallest singular value of the
*>  (2*n1*n2)-by-(2*n1*n2) matrix
*>
*>       Zu = [ kron(In2, A11)  -kron(A22**H, In1) ]
*>            [ kron(In2, B11)  -kron(B22**H, In1) ].
*>
*>  Here, Inx is the identity matrix of size nx and A22**H is the
*>  conjugate transpose of A22. kron(X, Y) is the Kronecker product between
*>  the matrices X and Y.
*>
*>  When DIF(2) is small, small changes in (A, B) can cause large changes
*>  in the deflating subspace. An approximate (asymptotic) bound on the
*>  maximum angular error in the computed deflating subspaces is
*>
*>       EPS * norm((A, B)) / DIF(2),
*>
*>  where EPS is the machine precision.
*>
*>  The reciprocal norm of the projectors on the left and right
*>  eigenspaces associated with (A11, B11) may be returned in PL and PR.
*>  They are computed as follows. First we compute L and R so that
*>  P*(A, B)*Q is block diagonal, where
*>
*>       P = ( I -L ) n1           Q = ( I R ) n1
*>           ( 0  I ) n2    and        ( 0 I ) n2
*>             n1 n2                    n1 n2
*>
*>  and (L, R) is the solution to the generalized Sylvester equation
*>
*>       A11*R - L*A22 = -A12
*>       B11*R - L*B22 = -B12
*>
*>  Then PL = (F-norm(L)**2+1)**(-1/2) and PR = (F-norm(R)**2+1)**(-1/2).
*>  An approximate (asymptotic) bound on the average absolute error of
*>  the selected eigenvalues is
*>
*>       EPS * norm((A, B)) / PL.
*>
*>  There are also global error bounds which valid for perturbations up
*>  to a certain restriction:  A lower bound (x) on the smallest
*>  F-norm(E,F) for which an eigenvalue of (A11, B11) may move and
*>  coalesce with an eigenvalue of (A22, B22) under perturbation (E,F),
*>  (i.e. (A + E, B + F), is
*>
*>   x = min(Difu,Difl)/((1/(PL*PL)+1/(PR*PR))**(1/2)+2*max(1/PL,1/PR)).
*>
*>  An approximate bound on x can be computed from DIF(1:2), PL and PR.
*>
*>  If y = ( F-norm(E,F) / x) <= 1, the angles between the perturbed
*>  (L', R') and unperturbed (L, R) left and right deflating subspaces
*>  associated with the selected cluster in the (1,1)-blocks can be
*>  bounded as
*>
*>   max-angle(L, L') <= arctan( y * PL / (1 - y * (1 - PL * PL)**(1/2))
*>   max-angle(R, R') <= arctan( y * PR / (1 - y * (1 - PR * PR)**(1/2))
*>
*>  See LAPACK User's Guide section 4.11 or the following references
*>  for more information.
*>
*>  Note that if the default method for computing the Frobenius-norm-
*>  based estimate DIF is not wanted (see ZLATDF), then the parameter
*>  IDIFJB (see below) should be changed from 3 to 4 (routine ZLATDF
*>  (IJOB = 2 will be used)). See ZTGSYL for more details.
*> \endverbatim
*
*> \par Contributors:
*  ==================
*>
*>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
*>     Umea University, S-901 87 Umea, Sweden.
*
*> \par References:
*  ================
*>
*>  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the
*>      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in
*>      M.S. Moonen et al (eds), Linear Algebra for Large Scale and
*>      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.
*> \n
*>  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified
*>      Eigenvalues of a Regular Matrix Pair (A, B) and Condition
*>      Estimation: Theory, Algorithms and Software, Report
*>      UMINF - 94.04, Department of Computing Science, Umea University,
*>      S-901 87 Umea, Sweden, 1994. Also as LAPACK Working Note 87.
*>      To appear in Numerical Algorithms, 1996.
*> \n
*>  [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software
*>      for Solving the Generalized Sylvester Equation and Estimating the
*>      Separation between Regular Matrix Pairs, Report UMINF - 93.23,
*>      Department of Computing Science, Umea University, S-901 87 Umea,
*>      Sweden, December 1993, Revised April 1994, Also as LAPACK working
*>      Note 75. To appear in ACM Trans. on Math. Software, Vol 22, No 1,
*>      1996.
*>
*  =====================================================================
      SUBROUTINE ZTGSEN( IJOB, WANTQ, WANTZ, SELECT, N, A, LDA, B, LDB,
     $                   ALPHA, BETA, Q, LDQ, Z, LDZ, M, PL, PR, DIF,
     $                   WORK, LWORK, IWORK, LIWORK, INFO )
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      LOGICAL            WANTQ, WANTZ
      INTEGER            IJOB, INFO, LDA, LDB, LDQ, LDZ, LIWORK, LWORK,
     $                   M, N
      DOUBLE PRECISION   PL, PR
*     ..
*     .. Array Arguments ..
      LOGICAL            SELECT( * )
      INTEGER            IWORK( * )
      DOUBLE PRECISION   DIF( * )
      COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),
     $                   BETA( * ), Q( LDQ, * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            IDIFJB
      PARAMETER          ( IDIFJB = 3 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, SWAP, WANTD, WANTD1, WANTD2, WANTP
      INTEGER            I, IERR, IJB, K, KASE, KS, LIWMIN, LWMIN, MN2,
     $                   N1, N2
      DOUBLE PRECISION   DSCALE, DSUM, RDSCAL, SAFMIN
      COMPLEX*16         TEMP1, TEMP2
*     ..
*     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLACN2, ZLACPY, ZLASSQ, ZSCAL, ZTGEXC,
     $                   ZTGSYL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DCMPLX, DCONJG, MAX, SQRT
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input parameters
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
*
      IF( IJOB.LT.0 .OR. IJOB.GT.5 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.N ) ) THEN
         INFO = -13
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -15
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTGSEN', -INFO )
         RETURN
      END IF
*
      IERR = 0
*
      WANTP = IJOB.EQ.1 .OR. IJOB.GE.4
      WANTD1 = IJOB.EQ.2 .OR. IJOB.EQ.4
      WANTD2 = IJOB.EQ.3 .OR. IJOB.EQ.5
      WANTD = WANTD1 .OR. WANTD2
*
*     Set M to the dimension of the specified pair of deflating
*     subspaces.
*
      M = 0
      DO 10 K = 1, N
         ALPHA( K ) = A( K, K )
         BETA( K ) = B( K, K )
         IF( K.LT.N ) THEN
            IF( SELECT( K ) )
     $         M = M + 1
         ELSE
            IF( SELECT( N ) )
     $         M = M + 1
         END IF
   10 CONTINUE
*
      IF( IJOB.EQ.1 .OR. IJOB.EQ.2 .OR. IJOB.EQ.4 ) THEN
         LWMIN = MAX( 1, 2*M*( N-M ) )
         LIWMIN = MAX( 1, N+2 )
      ELSE IF( IJOB.EQ.3 .OR. IJOB.EQ.5 ) THEN
         LWMIN = MAX( 1, 4*M*( N-M ) )
         LIWMIN = MAX( 1, 2*M*( N-M ), N+2 )
      ELSE
         LWMIN = 1
         LIWMIN = 1
      END IF
*
      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN
*
      IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -21
      ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT.LQUERY ) THEN
         INFO = -23
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTGSEN', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( M.EQ.N .OR. M.EQ.0 ) THEN
         IF( WANTP ) THEN
            PL = ONE
            PR = ONE
         END IF
         IF( WANTD ) THEN
            DSCALE = ZERO
            DSUM = ONE
            DO 20 I = 1, N
               CALL ZLASSQ( N, A( 1, I ), 1, DSCALE, DSUM )
               CALL ZLASSQ( N, B( 1, I ), 1, DSCALE, DSUM )
   20       CONTINUE
            DIF( 1 ) = DSCALE*SQRT( DSUM )
            DIF( 2 ) = DIF( 1 )
         END IF
         GO TO 70
      END IF
*
*     Get machine constant
*
      SAFMIN = DLAMCH( 'S' )
*
*     Collect the selected blocks at the top-left corner of (A, B).
*
      KS = 0
      DO 30 K = 1, N
         SWAP = SELECT( K )
         IF( SWAP ) THEN
            KS = KS + 1
*
*           Swap the K-th block to position KS. Compute unitary Q
*           and Z that will swap adjacent diagonal blocks in (A, B).
*
            IF( K.NE.KS )
     $         CALL ZTGEXC( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ, Z,
     $                      LDZ, K, KS, IERR )
*
            IF( IERR.GT.0 ) THEN
*
*              Swap is rejected: exit.
*
               INFO = 1
               IF( WANTP ) THEN
                  PL = ZERO
                  PR = ZERO
               END IF
               IF( WANTD ) THEN
                  DIF( 1 ) = ZERO
                  DIF( 2 ) = ZERO
               END IF
               GO TO 70
            END IF
         END IF
   30 CONTINUE
      IF( WANTP ) THEN
*
*        Solve generalized Sylvester equation for R and L:
*                   A11 * R - L * A22 = A12
*                   B11 * R - L * B22 = B12
*
         N1 = M
         N2 = N - M
         I = N1 + 1
         CALL ZLACPY( 'Full', N1, N2, A( 1, I ), LDA, WORK, N1 )
         CALL ZLACPY( 'Full', N1, N2, B( 1, I ), LDB, WORK( N1*N2+1 ),
     $                N1 )
         IJB = 0
         CALL ZTGSYL( 'N', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK,
     $                N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ), N1,
     $                DSCALE, DIF( 1 ), WORK( N1*N2*2+1 ),
     $                LWORK-2*N1*N2, IWORK, IERR )
*
*        Estimate the reciprocal of norms of "projections" onto
*        left and right eigenspaces
*
         RDSCAL = ZERO
         DSUM = ONE
         CALL ZLASSQ( N1*N2, WORK, 1, RDSCAL, DSUM )
         PL = RDSCAL*SQRT( DSUM )
         IF( PL.EQ.ZERO ) THEN
            PL = ONE
         ELSE
            PL = DSCALE / ( SQRT( DSCALE*DSCALE / PL+PL )*SQRT( PL ) )
         END IF
         RDSCAL = ZERO
         DSUM = ONE
         CALL ZLASSQ( N1*N2, WORK( N1*N2+1 ), 1, RDSCAL, DSUM )
         PR = RDSCAL*SQRT( DSUM )
         IF( PR.EQ.ZERO ) THEN
            PR = ONE
         ELSE
            PR = DSCALE / ( SQRT( DSCALE*DSCALE / PR+PR )*SQRT( PR ) )
         END IF
      END IF
      IF( WANTD ) THEN
*
*        Compute estimates Difu and Difl.
*
         IF( WANTD1 ) THEN
            N1 = M
            N2 = N - M
            I = N1 + 1
            IJB = IDIFJB
*
*           Frobenius norm-based Difu estimate.
*
            CALL ZTGSYL( 'N', IJB, N1, N2, A, LDA, A( I, I ), LDA, WORK,
     $                   N1, B, LDB, B( I, I ), LDB, WORK( N1*N2+1 ),
     $                   N1, DSCALE, DIF( 1 ), WORK( N1*N2*2+1 ),
     $                   LWORK-2*N1*N2, IWORK, IERR )
*
*           Frobenius norm-based Difl estimate.
*
            CALL ZTGSYL( 'N', IJB, N2, N1, A( I, I ), LDA, A, LDA, WORK,
     $                   N2, B( I, I ), LDB, B, LDB, WORK( N1*N2+1 ),
     $                   N2, DSCALE, DIF( 2 ), WORK( N1*N2*2+1 ),
     $                   LWORK-2*N1*N2, IWORK, IERR )
         ELSE
*
*           Compute 1-norm-based estimates of Difu and Difl using
*           reversed communication with ZLACN2. In each step a
*           generalized Sylvester equation or a transposed variant
*           is solved.
*
            KASE = 0
            N1 = M
            N2 = N - M
            I = N1 + 1
            IJB = 0
            MN2 = 2*N1*N2
*
*           1-norm-based estimate of Difu.
*
   40       CONTINUE
            CALL ZLACN2( MN2, WORK( MN2+1 ), WORK, DIF( 1 ), KASE,
     $                   ISAVE )
            IF( KASE.NE.0 ) THEN
               IF( KASE.EQ.1 ) THEN
*
*                 Solve generalized Sylvester equation
*
                  CALL ZTGSYL( 'N', IJB, N1, N2, A, LDA, A( I, I ), LDA,
     $                         WORK, N1, B, LDB, B( I, I ), LDB,
     $                         WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ),
     $                         WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK,
     $                         IERR )
               ELSE
*
*                 Solve the transposed variant.
*
                  CALL ZTGSYL( 'C', IJB, N1, N2, A, LDA, A( I, I ), LDA,
     $                         WORK, N1, B, LDB, B( I, I ), LDB,
     $                         WORK( N1*N2+1 ), N1, DSCALE, DIF( 1 ),
     $                         WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK,
     $                         IERR )
               END IF
               GO TO 40
            END IF
            DIF( 1 ) = DSCALE / DIF( 1 )
*
*           1-norm-based estimate of Difl.
*
   50       CONTINUE
            CALL ZLACN2( MN2, WORK( MN2+1 ), WORK, DIF( 2 ), KASE,
     $                   ISAVE )
            IF( KASE.NE.0 ) THEN
               IF( KASE.EQ.1 ) THEN
*
*                 Solve generalized Sylvester equation
*
                  CALL ZTGSYL( 'N', IJB, N2, N1, A( I, I ), LDA, A, LDA,
     $                         WORK, N2, B( I, I ), LDB, B, LDB,
     $                         WORK( N1*N2+1 ), N2, DSCALE, DIF( 2 ),
     $                         WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK,
     $                         IERR )
               ELSE
*
*                 Solve the transposed variant.
*
                  CALL ZTGSYL( 'C', IJB, N2, N1, A( I, I ), LDA, A, LDA,
     $                         WORK, N2, B, LDB, B( I, I ), LDB,
     $                         WORK( N1*N2+1 ), N2, DSCALE, DIF( 2 ),
     $                         WORK( N1*N2*2+1 ), LWORK-2*N1*N2, IWORK,
     $                         IERR )
               END IF
               GO TO 50
            END IF
            DIF( 2 ) = DSCALE / DIF( 2 )
         END IF
      END IF
*
*     If B(K,K) is complex, make it real and positive (normalization
*     of the generalized Schur form) and Store the generalized
*     eigenvalues of reordered pair (A, B)
*
      DO 60 K = 1, N
         DSCALE = ABS( B( K, K ) )
         IF( DSCALE.GT.SAFMIN ) THEN
            TEMP1 = DCONJG( B( K, K ) / DSCALE )
            TEMP2 = B( K, K ) / DSCALE
            B( K, K ) = DSCALE
            CALL ZSCAL( N-K, TEMP1, B( K, K+1 ), LDB )
            CALL ZSCAL( N-K+1, TEMP1, A( K, K ), LDA )
            IF( WANTQ )
     $         CALL ZSCAL( N, TEMP2, Q( 1, K ), 1 )
         ELSE
            B( K, K ) = DCMPLX( ZERO, ZERO )
         END IF
*
         ALPHA( K ) = A( K, K )
         BETA( K ) = B( K, K )
*
   60 CONTINUE
*
   70 CONTINUE
*
      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN
*
      RETURN
*
*     End of ZTGSEN
*
      END
*> \brief \b ZTGSY2 solves the generalized Sylvester equation (unblocked algorithm).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZTGSY2 + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztgsy2.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztgsy2.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztgsy2.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZTGSY2( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D,
*                          LDD, E, LDE, F, LDF, SCALE, RDSUM, RDSCAL,
*                          INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          TRANS
*       INTEGER            IJOB, INFO, LDA, LDB, LDC, LDD, LDE, LDF, M, N
*       DOUBLE PRECISION   RDSCAL, RDSUM, SCALE
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * ),
*      $                   D( LDD, * ), E( LDE, * ), F( LDF, * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZTGSY2 solves the generalized Sylvester equation
*>
*>             A * R - L * B = scale * C               (1)
*>             D * R - L * E = scale * F
*>
*> using Level 1 and 2 BLAS, where R and L are unknown M-by-N matrices,
*> (A, D), (B, E) and (C, F) are given matrix pairs of size M-by-M,
*> N-by-N and M-by-N, respectively. A, B, D and E are upper triangular
*> (i.e., (A,D) and (B,E) in generalized Schur form).
*>
*> The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1 is an output
*> scaling factor chosen to avoid overflow.
*>
*> In matrix notation solving equation (1) corresponds to solve
*> Zx = scale * b, where Z is defined as
*>
*>        Z = [ kron(In, A)  -kron(B**H, Im) ]             (2)
*>            [ kron(In, D)  -kron(E**H, Im) ],
*>
*> Ik is the identity matrix of size k and X**H is the conjuguate transpose of X.
*> kron(X, Y) is the Kronecker product between the matrices X and Y.
*>
*> If TRANS = 'C', y in the conjugate transposed system Z**H*y = scale*b
*> is solved for, which is equivalent to solve for R and L in
*>
*>             A**H * R  + D**H * L   = scale * C           (3)
*>             R  * B**H + L  * E**H  = scale * -F
*>
*> This case is used to compute an estimate of Dif[(A, D), (B, E)] =
*> = sigma_min(Z) using reverse communicaton with ZLACON.
*>
*> ZTGSY2 also (IJOB >= 1) contributes to the computation in ZTGSYL
*> of an upper bound on the separation between to matrix pairs. Then
*> the input (A, D), (B, E) are sub-pencils of two matrix pairs in
*> ZTGSYL.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          = 'N', solve the generalized Sylvester equation (1).
*>          = 'T': solve the 'transposed' system (3).
*> \endverbatim
*>
*> \param[in] IJOB
*> \verbatim
*>          IJOB is INTEGER
*>          Specifies what kind of functionality to be performed.
*>          =0: solve (1) only.
*>          =1: A contribution from this subsystem to a Frobenius
*>              norm-based estimate of the separation between two matrix
*>              pairs is computed. (look ahead strategy is used).
*>          =2: A contribution from this subsystem to a Frobenius
*>              norm-based estimate of the separation between two matrix
*>              pairs is computed. (DGECON on sub-systems is used.)
*>          Not referenced if TRANS = 'T'.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          On entry, M specifies the order of A and D, and the row
*>          dimension of C, F, R and L.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          On entry, N specifies the order of B and E, and the column
*>          dimension of C, F, R and L.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA, M)
*>          On entry, A contains an upper triangular matrix.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the matrix A. LDA >= max(1, M).
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is COMPLEX*16 array, dimension (LDB, N)
*>          On entry, B contains an upper triangular matrix.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the matrix B. LDB >= max(1, N).
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is COMPLEX*16 array, dimension (LDC, N)
*>          On entry, C contains the right-hand-side of the first matrix
*>          equation in (1).
*>          On exit, if IJOB = 0, C has been overwritten by the solution
*>          R.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the matrix C. LDC >= max(1, M).
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is COMPLEX*16 array, dimension (LDD, M)
*>          On entry, D contains an upper triangular matrix.
*> \endverbatim
*>
*> \param[in] LDD
*> \verbatim
*>          LDD is INTEGER
*>          The leading dimension of the matrix D. LDD >= max(1, M).
*> \endverbatim
*>
*> \param[in] E
*> \verbatim
*>          E is COMPLEX*16 array, dimension (LDE, N)
*>          On entry, E contains an upper triangular matrix.
*> \endverbatim
*>
*> \param[in] LDE
*> \verbatim
*>          LDE is INTEGER
*>          The leading dimension of the matrix E. LDE >= max(1, N).
*> \endverbatim
*>
*> \param[in,out] F
*> \verbatim
*>          F is COMPLEX*16 array, dimension (LDF, N)
*>          On entry, F contains the right-hand-side of the second matrix
*>          equation in (1).
*>          On exit, if IJOB = 0, F has been overwritten by the solution
*>          L.
*> \endverbatim
*>
*> \param[in] LDF
*> \verbatim
*>          LDF is INTEGER
*>          The leading dimension of the matrix F. LDF >= max(1, M).
*> \endverbatim
*>
*> \param[out] SCALE
*> \verbatim
*>          SCALE is DOUBLE PRECISION
*>          On exit, 0 <= SCALE <= 1. If 0 < SCALE < 1, the solutions
*>          R and L (C and F on entry) will hold the solutions to a
*>          slightly perturbed system but the input matrices A, B, D and
*>          E have not been changed. If SCALE = 0, R and L will hold the
*>          solutions to the homogeneous system with C = F = 0.
*>          Normally, SCALE = 1.
*> \endverbatim
*>
*> \param[in,out] RDSUM
*> \verbatim
*>          RDSUM is DOUBLE PRECISION
*>          On entry, the sum of squares of computed contributions to
*>          the Dif-estimate under computation by ZTGSYL, where the
*>          scaling factor RDSCAL (see below) has been factored out.
*>          On exit, the corresponding sum of squares updated with the
*>          contributions from the current sub-system.
*>          If TRANS = 'T' RDSUM is not touched.
*>          NOTE: RDSUM only makes sense when ZTGSY2 is called by
*>          ZTGSYL.
*> \endverbatim
*>
*> \param[in,out] RDSCAL
*> \verbatim
*>          RDSCAL is DOUBLE PRECISION
*>          On entry, scaling factor used to prevent overflow in RDSUM.
*>          On exit, RDSCAL is updated w.r.t. the current contributions
*>          in RDSUM.
*>          If TRANS = 'T', RDSCAL is not touched.
*>          NOTE: RDSCAL only makes sense when ZTGSY2 is called by
*>          ZTGSYL.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          On exit, if INFO is set to
*>            =0: Successful exit
*>            <0: If INFO = -i, input argument number i is illegal.
*>            >0: The matrix pairs (A, D) and (B, E) have common or very
*>                close eigenvalues.
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
*> \ingroup complex16SYauxiliary
*
*> \par Contributors:
*  ==================
*>
*>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
*>     Umea University, S-901 87 Umea, Sweden.
*
*  =====================================================================
      SUBROUTINE ZTGSY2( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D,
     $                   LDD, E, LDE, F, LDF, SCALE, RDSUM, RDSCAL,
     $                   INFO )
*
*  -- LAPACK auxiliary routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            IJOB, INFO, LDA, LDB, LDC, LDD, LDE, LDF, M, N
      DOUBLE PRECISION   RDSCAL, RDSUM, SCALE
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * ),
     $                   D( LDD, * ), E( LDE, * ), F( LDF, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      INTEGER            LDZ
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, LDZ = 2 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN
      INTEGER            I, IERR, J, K
      DOUBLE PRECISION   SCALOC
      COMPLEX*16         ALPHA
*     ..
*     .. Local Arrays ..
      INTEGER            IPIV( LDZ ), JPIV( LDZ )
      COMPLEX*16         RHS( LDZ ), Z( LDZ, LDZ )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZAXPY, ZGESC2, ZGETC2, ZLATDF, ZSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCMPLX, DCONJG, MAX
*     ..
*     .. Executable Statements ..
*
*     Decode and test input parameters
*
      INFO = 0
      IERR = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( NOTRAN ) THEN
         IF( ( IJOB.LT.0 ) .OR. ( IJOB.GT.2 ) ) THEN
            INFO = -2
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( M.LE.0 ) THEN
            INFO = -3
         ELSE IF( N.LE.0 ) THEN
            INFO = -4
         ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
            INFO = -6
         ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
            INFO = -8
         ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
            INFO = -10
         ELSE IF( LDD.LT.MAX( 1, M ) ) THEN
            INFO = -12
         ELSE IF( LDE.LT.MAX( 1, N ) ) THEN
            INFO = -14
         ELSE IF( LDF.LT.MAX( 1, M ) ) THEN
            INFO = -16
         END IF
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTGSY2', -INFO )
         RETURN
      END IF
*
      IF( NOTRAN ) THEN
*
*        Solve (I, J) - system
*           A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J)
*           D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J)
*        for I = M, M - 1, ..., 1; J = 1, 2, ..., N
*
         SCALE = ONE
         SCALOC = ONE
         DO 30 J = 1, N
            DO 20 I = M, 1, -1
*
*              Build 2 by 2 system
*
               Z( 1, 1 ) = A( I, I )
               Z( 2, 1 ) = D( I, I )
               Z( 1, 2 ) = -B( J, J )
               Z( 2, 2 ) = -E( J, J )
*
*              Set up right hand side(s)
*
               RHS( 1 ) = C( I, J )
               RHS( 2 ) = F( I, J )
*
*              Solve Z * x = RHS
*
               CALL ZGETC2( LDZ, Z, LDZ, IPIV, JPIV, IERR )
               IF( IERR.GT.0 )
     $            INFO = IERR
               IF( IJOB.EQ.0 ) THEN
                  CALL ZGESC2( LDZ, Z, LDZ, RHS, IPIV, JPIV, SCALOC )
                  IF( SCALOC.NE.ONE ) THEN
                     DO 10 K = 1, N
                        CALL ZSCAL( M, DCMPLX( SCALOC, ZERO ),
     $                              C( 1, K ), 1 )
                        CALL ZSCAL( M, DCMPLX( SCALOC, ZERO ),
     $                              F( 1, K ), 1 )
   10                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
               ELSE
                  CALL ZLATDF( IJOB, LDZ, Z, LDZ, RHS, RDSUM, RDSCAL,
     $                         IPIV, JPIV )
               END IF
*
*              Unpack solution vector(s)
*
               C( I, J ) = RHS( 1 )
               F( I, J ) = RHS( 2 )
*
*              Substitute R(I, J) and L(I, J) into remaining equation.
*
               IF( I.GT.1 ) THEN
                  ALPHA = -RHS( 1 )
                  CALL ZAXPY( I-1, ALPHA, A( 1, I ), 1, C( 1, J ), 1 )
                  CALL ZAXPY( I-1, ALPHA, D( 1, I ), 1, F( 1, J ), 1 )
               END IF
               IF( J.LT.N ) THEN
                  CALL ZAXPY( N-J, RHS( 2 ), B( J, J+1 ), LDB,
     $                        C( I, J+1 ), LDC )
                  CALL ZAXPY( N-J, RHS( 2 ), E( J, J+1 ), LDE,
     $                        F( I, J+1 ), LDF )
               END IF
*
   20       CONTINUE
   30    CONTINUE
      ELSE
*
*        Solve transposed (I, J) - system:
*           A(I, I)**H * R(I, J) + D(I, I)**H * L(J, J) = C(I, J)
*           R(I, I) * B(J, J) + L(I, J) * E(J, J)   = -F(I, J)
*        for I = 1, 2, ..., M, J = N, N - 1, ..., 1
*
         SCALE = ONE
         SCALOC = ONE
         DO 80 I = 1, M
            DO 70 J = N, 1, -1
*
*              Build 2 by 2 system Z**H
*
               Z( 1, 1 ) = DCONJG( A( I, I ) )
               Z( 2, 1 ) = -DCONJG( B( J, J ) )
               Z( 1, 2 ) = DCONJG( D( I, I ) )
               Z( 2, 2 ) = -DCONJG( E( J, J ) )
*
*
*              Set up right hand side(s)
*
               RHS( 1 ) = C( I, J )
               RHS( 2 ) = F( I, J )
*
*              Solve Z**H * x = RHS
*
               CALL ZGETC2( LDZ, Z, LDZ, IPIV, JPIV, IERR )
               IF( IERR.GT.0 )
     $            INFO = IERR
               CALL ZGESC2( LDZ, Z, LDZ, RHS, IPIV, JPIV, SCALOC )
               IF( SCALOC.NE.ONE ) THEN
                  DO 40 K = 1, N
                     CALL ZSCAL( M, DCMPLX( SCALOC, ZERO ), C( 1, K ),
     $                           1 )
                     CALL ZSCAL( M, DCMPLX( SCALOC, ZERO ), F( 1, K ),
     $                           1 )
   40             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
*
*              Unpack solution vector(s)
*
               C( I, J ) = RHS( 1 )
               F( I, J ) = RHS( 2 )
*
*              Substitute R(I, J) and L(I, J) into remaining equation.
*
               DO 50 K = 1, J - 1
                  F( I, K ) = F( I, K ) + RHS( 1 )*DCONJG( B( K, J ) ) +
     $                        RHS( 2 )*DCONJG( E( K, J ) )
   50          CONTINUE
               DO 60 K = I + 1, M
                  C( K, J ) = C( K, J ) - DCONJG( A( I, K ) )*RHS( 1 ) -
     $                        DCONJG( D( I, K ) )*RHS( 2 )
   60          CONTINUE
*
   70       CONTINUE
   80    CONTINUE
      END IF
      RETURN
*
*     End of ZTGSY2
*
      END
*> \brief \b ZTGSYL
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*> \htmlonly
*> Download ZTGSYL + dependencies 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztgsyl.f"> 
*> [TGZ]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztgsyl.f"> 
*> [ZIP]</a> 
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztgsyl.f"> 
*> [TXT]</a>
*> \endhtmlonly 
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZTGSYL( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D,
*                          LDD, E, LDE, F, LDF, SCALE, DIF, WORK, LWORK,
*                          IWORK, INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          TRANS
*       INTEGER            IJOB, INFO, LDA, LDB, LDC, LDD, LDE, LDF,
*      $                   LWORK, M, N
*       DOUBLE PRECISION   DIF, SCALE
*       ..
*       .. Array Arguments ..
*       INTEGER            IWORK( * )
*       COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * ),
*      $                   D( LDD, * ), E( LDE, * ), F( LDF, * ),
*      $                   WORK( * )
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZTGSYL solves the generalized Sylvester equation:
*>
*>             A * R - L * B = scale * C            (1)
*>             D * R - L * E = scale * F
*>
*> where R and L are unknown m-by-n matrices, (A, D), (B, E) and
*> (C, F) are given matrix pairs of size m-by-m, n-by-n and m-by-n,
*> respectively, with complex entries. A, B, D and E are upper
*> triangular (i.e., (A,D) and (B,E) in generalized Schur form).
*>
*> The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1
*> is an output scaling factor chosen to avoid overflow.
*>
*> In matrix notation (1) is equivalent to solve Zx = scale*b, where Z
*> is defined as
*>
*>        Z = [ kron(In, A)  -kron(B**H, Im) ]        (2)
*>            [ kron(In, D)  -kron(E**H, Im) ],
*>
*> Here Ix is the identity matrix of size x and X**H is the conjugate
*> transpose of X. Kron(X, Y) is the Kronecker product between the
*> matrices X and Y.
*>
*> If TRANS = 'C', y in the conjugate transposed system Z**H *y = scale*b
*> is solved for, which is equivalent to solve for R and L in
*>
*>             A**H * R + D**H * L = scale * C           (3)
*>             R * B**H + L * E**H = scale * -F
*>
*> This case (TRANS = 'C') is used to compute an one-norm-based estimate
*> of Dif[(A,D), (B,E)], the separation between the matrix pairs (A,D)
*> and (B,E), using ZLACON.
*>
*> If IJOB >= 1, ZTGSYL computes a Frobenius norm-based estimate of
*> Dif[(A,D),(B,E)]. That is, the reciprocal of a lower bound on the
*> reciprocal of the smallest singular value of Z.
*>
*> This is a level-3 BLAS algorithm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          = 'N': solve the generalized sylvester equation (1).
*>          = 'C': solve the "conjugate transposed" system (3).
*> \endverbatim
*>
*> \param[in] IJOB
*> \verbatim
*>          IJOB is INTEGER
*>          Specifies what kind of functionality to be performed.
*>          =0: solve (1) only.
*>          =1: The functionality of 0 and 3.
*>          =2: The functionality of 0 and 4.
*>          =3: Only an estimate of Dif[(A,D), (B,E)] is computed.
*>              (look ahead strategy is used).
*>          =4: Only an estimate of Dif[(A,D), (B,E)] is computed.
*>              (ZGECON on sub-systems is used).
*>          Not referenced if TRANS = 'C'.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The order of the matrices A and D, and the row dimension of
*>          the matrices C, F, R and L.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices B and E, and the column dimension
*>          of the matrices C, F, R and L.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA, M)
*>          The upper triangular matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A. LDA >= max(1, M).
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is COMPLEX*16 array, dimension (LDB, N)
*>          The upper triangular matrix B.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B. LDB >= max(1, N).
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is COMPLEX*16 array, dimension (LDC, N)
*>          On entry, C contains the right-hand-side of the first matrix
*>          equation in (1) or (3).
*>          On exit, if IJOB = 0, 1 or 2, C has been overwritten by
*>          the solution R. If IJOB = 3 or 4 and TRANS = 'N', C holds R,
*>          the solution achieved during the computation of the
*>          Dif-estimate.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C. LDC >= max(1, M).
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is COMPLEX*16 array, dimension (LDD, M)
*>          The upper triangular matrix D.
*> \endverbatim
*>
*> \param[in] LDD
*> \verbatim
*>          LDD is INTEGER
*>          The leading dimension of the array D. LDD >= max(1, M).
*> \endverbatim
*>
*> \param[in] E
*> \verbatim
*>          E is COMPLEX*16 array, dimension (LDE, N)
*>          The upper triangular matrix E.
*> \endverbatim
*>
*> \param[in] LDE
*> \verbatim
*>          LDE is INTEGER
*>          The leading dimension of the array E. LDE >= max(1, N).
*> \endverbatim
*>
*> \param[in,out] F
*> \verbatim
*>          F is COMPLEX*16 array, dimension (LDF, N)
*>          On entry, F contains the right-hand-side of the second matrix
*>          equation in (1) or (3).
*>          On exit, if IJOB = 0, 1 or 2, F has been overwritten by
*>          the solution L. If IJOB = 3 or 4 and TRANS = 'N', F holds L,
*>          the solution achieved during the computation of the
*>          Dif-estimate.
*> \endverbatim
*>
*> \param[in] LDF
*> \verbatim
*>          LDF is INTEGER
*>          The leading dimension of the array F. LDF >= max(1, M).
*> \endverbatim
*>
*> \param[out] DIF
*> \verbatim
*>          DIF is DOUBLE PRECISION
*>          On exit DIF is the reciprocal of a lower bound of the
*>          reciprocal of the Dif-function, i.e. DIF is an upper bound of
*>          Dif[(A,D), (B,E)] = sigma-min(Z), where Z as in (2).
*>          IF IJOB = 0 or TRANS = 'C', DIF is not referenced.
*> \endverbatim
*>
*> \param[out] SCALE
*> \verbatim
*>          SCALE is DOUBLE PRECISION
*>          On exit SCALE is the scaling factor in (1) or (3).
*>          If 0 < SCALE < 1, C and F hold the solutions R and L, resp.,
*>          to a slightly perturbed system but the input matrices A, B,
*>          D and E have not been changed. If SCALE = 0, R and L will
*>          hold the solutions to the homogenious system with C = F = 0.
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
*>          The dimension of the array WORK. LWORK > = 1.
*>          If IJOB = 1 or 2 and TRANS = 'N', LWORK >= max(1,2*M*N).
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (M+N+2)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>            =0: successful exit
*>            <0: If INFO = -i, the i-th argument had an illegal value.
*>            >0: (A, D) and (B, E) have common or very close
*>                eigenvalues.
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
*> \ingroup complex16SYcomputational
*
*> \par Contributors:
*  ==================
*>
*>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
*>     Umea University, S-901 87 Umea, Sweden.
*
*> \par References:
*  ================
*>
*>  [1] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software
*>      for Solving the Generalized Sylvester Equation and Estimating the
*>      Separation between Regular Matrix Pairs, Report UMINF - 93.23,
*>      Department of Computing Science, Umea University, S-901 87 Umea,
*>      Sweden, December 1993, Revised April 1994, Also as LAPACK Working
*>      Note 75.  To appear in ACM Trans. on Math. Software, Vol 22,
*>      No 1, 1996.
*> \n
*>  [2] B. Kagstrom, A Perturbation Analysis of the Generalized Sylvester
*>      Equation (AR - LB, DR - LE ) = (C, F), SIAM J. Matrix Anal.
*>      Appl., 15(4):1045-1060, 1994.
*> \n
*>  [3] B. Kagstrom and L. Westin, Generalized Schur Methods with
*>      Condition Estimators for Solving the Generalized Sylvester
*>      Equation, IEEE Transactions on Automatic Control, Vol. 34, No. 7,
*>      July 1989, pp 745-751.
*>
*  =====================================================================
      SUBROUTINE ZTGSYL( TRANS, IJOB, M, N, A, LDA, B, LDB, C, LDC, D,
     $                   LDD, E, LDE, F, LDF, SCALE, DIF, WORK, LWORK,
     $                   IWORK, INFO )
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            IJOB, INFO, LDA, LDB, LDC, LDD, LDE, LDF,
     $                   LWORK, M, N
      DOUBLE PRECISION   DIF, SCALE
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * ),
     $                   D( LDD, * ), E( LDE, * ), F( LDF, * ),
     $                   WORK( * )
*     ..
*
*  =====================================================================
*  Replaced various illegal calls to CCOPY by calls to CLASET.
*  Sven Hammarling, 1/5/02.
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO
      PARAMETER          ( CZERO = (0.0D+0, 0.0D+0) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, NOTRAN
      INTEGER            I, IE, IFUNC, IROUND, IS, ISOLVE, J, JE, JS, K,
     $                   LINFO, LWMIN, MB, NB, P, PQ, Q
      DOUBLE PRECISION   DSCALE, DSUM, SCALE2, SCALOC
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEMM, ZLACPY, ZLASET, ZSCAL, ZTGSY2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Decode and test input parameters
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
*
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( NOTRAN ) THEN
         IF( ( IJOB.LT.0 ) .OR. ( IJOB.GT.4 ) ) THEN
            INFO = -2
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( M.LE.0 ) THEN
            INFO = -3
         ELSE IF( N.LE.0 ) THEN
            INFO = -4
         ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
            INFO = -6
         ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
            INFO = -8
         ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
            INFO = -10
         ELSE IF( LDD.LT.MAX( 1, M ) ) THEN
            INFO = -12
         ELSE IF( LDE.LT.MAX( 1, N ) ) THEN
            INFO = -14
         ELSE IF( LDF.LT.MAX( 1, M ) ) THEN
            INFO = -16
         END IF
      END IF
*
      IF( INFO.EQ.0 ) THEN
         IF( NOTRAN ) THEN
            IF( IJOB.EQ.1 .OR. IJOB.EQ.2 ) THEN
               LWMIN = MAX( 1, 2*M*N )
            ELSE
               LWMIN = 1
            END IF
         ELSE
            LWMIN = 1
         END IF
         WORK( 1 ) = LWMIN
*
         IF( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) THEN
            INFO = -20
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTGSYL', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
         SCALE = 1
         IF( NOTRAN ) THEN
            IF( IJOB.NE.0 ) THEN
               DIF = 0
            END IF
         END IF
         RETURN
      END IF
*
*     Determine  optimal block sizes MB and NB
*
      MB = ILAENV( 2, 'ZTGSYL', TRANS, M, N, -1, -1 )
      NB = ILAENV( 5, 'ZTGSYL', TRANS, M, N, -1, -1 )
*
      ISOLVE = 1
      IFUNC = 0
      IF( NOTRAN ) THEN
         IF( IJOB.GE.3 ) THEN
            IFUNC = IJOB - 2
            CALL ZLASET( 'F', M, N, CZERO, CZERO, C, LDC )
            CALL ZLASET( 'F', M, N, CZERO, CZERO, F, LDF )
         ELSE IF( IJOB.GE.1 .AND. NOTRAN ) THEN
            ISOLVE = 2
         END IF
      END IF
*
      IF( ( MB.LE.1 .AND. NB.LE.1 ) .OR. ( MB.GE.M .AND. NB.GE.N ) )
     $     THEN
*
*        Use unblocked Level 2 solver
*
         DO 30 IROUND = 1, ISOLVE
*
            SCALE = ONE
            DSCALE = ZERO
            DSUM = ONE
            PQ = M*N
            CALL ZTGSY2( TRANS, IFUNC, M, N, A, LDA, B, LDB, C, LDC, D,
     $                   LDD, E, LDE, F, LDF, SCALE, DSUM, DSCALE,
     $                   INFO )
            IF( DSCALE.NE.ZERO ) THEN
               IF( IJOB.EQ.1 .OR. IJOB.EQ.3 ) THEN
                  DIF = SQRT( DBLE( 2*M*N ) ) / ( DSCALE*SQRT( DSUM ) )
               ELSE
                  DIF = SQRT( DBLE( PQ ) ) / ( DSCALE*SQRT( DSUM ) )
               END IF
            END IF
            IF( ISOLVE.EQ.2 .AND. IROUND.EQ.1 ) THEN
               IF( NOTRAN ) THEN
                  IFUNC = IJOB
               END IF
               SCALE2 = SCALE
               CALL ZLACPY( 'F', M, N, C, LDC, WORK, M )
               CALL ZLACPY( 'F', M, N, F, LDF, WORK( M*N+1 ), M )
               CALL ZLASET( 'F', M, N, CZERO, CZERO, C, LDC )
               CALL ZLASET( 'F', M, N, CZERO, CZERO, F, LDF )
            ELSE IF( ISOLVE.EQ.2 .AND. IROUND.EQ.2 ) THEN
               CALL ZLACPY( 'F', M, N, WORK, M, C, LDC )
               CALL ZLACPY( 'F', M, N, WORK( M*N+1 ), M, F, LDF )
               SCALE = SCALE2
            END IF
   30    CONTINUE
*
         RETURN
*
      END IF
*
*     Determine block structure of A
*
      P = 0
      I = 1
   40 CONTINUE
      IF( I.GT.M )
     $   GO TO 50
      P = P + 1
      IWORK( P ) = I
      I = I + MB
      IF( I.GE.M )
     $   GO TO 50
      GO TO 40
   50 CONTINUE
      IWORK( P+1 ) = M + 1
      IF( IWORK( P ).EQ.IWORK( P+1 ) )
     $   P = P - 1
*
*     Determine block structure of B
*
      Q = P + 1
      J = 1
   60 CONTINUE
      IF( J.GT.N )
     $   GO TO 70
*
      Q = Q + 1
      IWORK( Q ) = J
      J = J + NB
      IF( J.GE.N )
     $   GO TO 70
      GO TO 60
*
   70 CONTINUE
      IWORK( Q+1 ) = N + 1
      IF( IWORK( Q ).EQ.IWORK( Q+1 ) )
     $   Q = Q - 1
*
      IF( NOTRAN ) THEN
         DO 150 IROUND = 1, ISOLVE
*
*           Solve (I, J) - subsystem
*               A(I, I) * R(I, J) - L(I, J) * B(J, J) = C(I, J)
*               D(I, I) * R(I, J) - L(I, J) * E(J, J) = F(I, J)
*           for I = P, P - 1, ..., 1; J = 1, 2, ..., Q
*
            PQ = 0
            SCALE = ONE
            DSCALE = ZERO
            DSUM = ONE
            DO 130 J = P + 2, Q
               JS = IWORK( J )
               JE = IWORK( J+1 ) - 1
               NB = JE - JS + 1
               DO 120 I = P, 1, -1
                  IS = IWORK( I )
                  IE = IWORK( I+1 ) - 1
                  MB = IE - IS + 1
                  CALL ZTGSY2( TRANS, IFUNC, MB, NB, A( IS, IS ), LDA,
     $                         B( JS, JS ), LDB, C( IS, JS ), LDC,
     $                         D( IS, IS ), LDD, E( JS, JS ), LDE,
     $                         F( IS, JS ), LDF, SCALOC, DSUM, DSCALE,
     $                         LINFO )
                  IF( LINFO.GT.0 )
     $               INFO = LINFO
                  PQ = PQ + MB*NB
                  IF( SCALOC.NE.ONE ) THEN
                     DO 80 K = 1, JS - 1
                        CALL ZSCAL( M, DCMPLX( SCALOC, ZERO ),
     $                              C( 1, K ), 1 )
                        CALL ZSCAL( M, DCMPLX( SCALOC, ZERO ),
     $                              F( 1, K ), 1 )
   80                CONTINUE
                     DO 90 K = JS, JE
                        CALL ZSCAL( IS-1, DCMPLX( SCALOC, ZERO ),
     $                              C( 1, K ), 1 )
                        CALL ZSCAL( IS-1, DCMPLX( SCALOC, ZERO ),
     $                              F( 1, K ), 1 )
   90                CONTINUE
                     DO 100 K = JS, JE
                        CALL ZSCAL( M-IE, DCMPLX( SCALOC, ZERO ),
     $                              C( IE+1, K ), 1 )
                        CALL ZSCAL( M-IE, DCMPLX( SCALOC, ZERO ),
     $                              F( IE+1, K ), 1 )
  100                CONTINUE
                     DO 110 K = JE + 1, N
                        CALL ZSCAL( M, DCMPLX( SCALOC, ZERO ),
     $                              C( 1, K ), 1 )
                        CALL ZSCAL( M, DCMPLX( SCALOC, ZERO ),
     $                              F( 1, K ), 1 )
  110                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
*
*                 Substitute R(I,J) and L(I,J) into remaining equation.
*
                  IF( I.GT.1 ) THEN
                     CALL ZGEMM( 'N', 'N', IS-1, NB, MB,
     $                           DCMPLX( -ONE, ZERO ), A( 1, IS ), LDA,
     $                           C( IS, JS ), LDC, DCMPLX( ONE, ZERO ),
     $                           C( 1, JS ), LDC )
                     CALL ZGEMM( 'N', 'N', IS-1, NB, MB,
     $                           DCMPLX( -ONE, ZERO ), D( 1, IS ), LDD,
     $                           C( IS, JS ), LDC, DCMPLX( ONE, ZERO ),
     $                           F( 1, JS ), LDF )
                  END IF
                  IF( J.LT.Q ) THEN
                     CALL ZGEMM( 'N', 'N', MB, N-JE, NB,
     $                           DCMPLX( ONE, ZERO ), F( IS, JS ), LDF,
     $                           B( JS, JE+1 ), LDB,
     $                           DCMPLX( ONE, ZERO ), C( IS, JE+1 ),
     $                           LDC )
                     CALL ZGEMM( 'N', 'N', MB, N-JE, NB,
     $                           DCMPLX( ONE, ZERO ), F( IS, JS ), LDF,
     $                           E( JS, JE+1 ), LDE,
     $                           DCMPLX( ONE, ZERO ), F( IS, JE+1 ),
     $                           LDF )
                  END IF
  120          CONTINUE
  130       CONTINUE
            IF( DSCALE.NE.ZERO ) THEN
               IF( IJOB.EQ.1 .OR. IJOB.EQ.3 ) THEN
                  DIF = SQRT( DBLE( 2*M*N ) ) / ( DSCALE*SQRT( DSUM ) )
               ELSE
                  DIF = SQRT( DBLE( PQ ) ) / ( DSCALE*SQRT( DSUM ) )
               END IF
            END IF
            IF( ISOLVE.EQ.2 .AND. IROUND.EQ.1 ) THEN
               IF( NOTRAN ) THEN
                  IFUNC = IJOB
               END IF
               SCALE2 = SCALE
               CALL ZLACPY( 'F', M, N, C, LDC, WORK, M )
               CALL ZLACPY( 'F', M, N, F, LDF, WORK( M*N+1 ), M )
               CALL ZLASET( 'F', M, N, CZERO, CZERO, C, LDC )
               CALL ZLASET( 'F', M, N, CZERO, CZERO, F, LDF )
            ELSE IF( ISOLVE.EQ.2 .AND. IROUND.EQ.2 ) THEN
               CALL ZLACPY( 'F', M, N, WORK, M, C, LDC )
               CALL ZLACPY( 'F', M, N, WORK( M*N+1 ), M, F, LDF )
               SCALE = SCALE2
            END IF
  150    CONTINUE
      ELSE
*
*        Solve transposed (I, J)-subsystem
*            A(I, I)**H * R(I, J) + D(I, I)**H * L(I, J) = C(I, J)
*            R(I, J) * B(J, J)  + L(I, J) * E(J, J) = -F(I, J)
*        for I = 1,2,..., P; J = Q, Q-1,..., 1
*
         SCALE = ONE
         DO 210 I = 1, P
            IS = IWORK( I )
            IE = IWORK( I+1 ) - 1
            MB = IE - IS + 1
            DO 200 J = Q, P + 2, -1
               JS = IWORK( J )
               JE = IWORK( J+1 ) - 1
               NB = JE - JS + 1
               CALL ZTGSY2( TRANS, IFUNC, MB, NB, A( IS, IS ), LDA,
     $                      B( JS, JS ), LDB, C( IS, JS ), LDC,
     $                      D( IS, IS ), LDD, E( JS, JS ), LDE,
     $                      F( IS, JS ), LDF, SCALOC, DSUM, DSCALE,
     $                      LINFO )
               IF( LINFO.GT.0 )
     $            INFO = LINFO
               IF( SCALOC.NE.ONE ) THEN
                  DO 160 K = 1, JS - 1
                     CALL ZSCAL( M, DCMPLX( SCALOC, ZERO ), C( 1, K ),
     $                           1 )
                     CALL ZSCAL( M, DCMPLX( SCALOC, ZERO ), F( 1, K ),
     $                           1 )
  160             CONTINUE
                  DO 170 K = JS, JE
                     CALL ZSCAL( IS-1, DCMPLX( SCALOC, ZERO ),
     $                           C( 1, K ), 1 )
                     CALL ZSCAL( IS-1, DCMPLX( SCALOC, ZERO ),
     $                           F( 1, K ), 1 )
  170             CONTINUE
                  DO 180 K = JS, JE
                     CALL ZSCAL( M-IE, DCMPLX( SCALOC, ZERO ),
     $                           C( IE+1, K ), 1 )
                     CALL ZSCAL( M-IE, DCMPLX( SCALOC, ZERO ),
     $                           F( IE+1, K ), 1 )
  180             CONTINUE
                  DO 190 K = JE + 1, N
                     CALL ZSCAL( M, DCMPLX( SCALOC, ZERO ), C( 1, K ),
     $                           1 )
                     CALL ZSCAL( M, DCMPLX( SCALOC, ZERO ), F( 1, K ),
     $                           1 )
  190             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
*
*              Substitute R(I,J) and L(I,J) into remaining equation.
*
               IF( J.GT.P+2 ) THEN
                  CALL ZGEMM( 'N', 'C', MB, JS-1, NB,
     $                        DCMPLX( ONE, ZERO ), C( IS, JS ), LDC,
     $                        B( 1, JS ), LDB, DCMPLX( ONE, ZERO ),
     $                        F( IS, 1 ), LDF )
                  CALL ZGEMM( 'N', 'C', MB, JS-1, NB,
     $                        DCMPLX( ONE, ZERO ), F( IS, JS ), LDF,
     $                        E( 1, JS ), LDE, DCMPLX( ONE, ZERO ),
     $                        F( IS, 1 ), LDF )
               END IF
               IF( I.LT.P ) THEN
                  CALL ZGEMM( 'C', 'N', M-IE, NB, MB,
     $                        DCMPLX( -ONE, ZERO ), A( IS, IE+1 ), LDA,
     $                        C( IS, JS ), LDC, DCMPLX( ONE, ZERO ),
     $                        C( IE+1, JS ), LDC )
                  CALL ZGEMM( 'C', 'N', M-IE, NB, MB,
     $                        DCMPLX( -ONE, ZERO ), D( IS, IE+1 ), LDD,
     $                        F( IS, JS ), LDF, DCMPLX( ONE, ZERO ),
     $                        C( IE+1, JS ), LDC )
               END IF
  200       CONTINUE
  210    CONTINUE
      END IF
*
      WORK( 1 ) = LWMIN
*
      RETURN
*
*     End of ZTGSYL
*
      END
