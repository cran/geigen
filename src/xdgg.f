c     dummy interface to force availability of dggev and dsygv

      SUBROUTINE xDGGEV(JOBVL,JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI,
     *                  BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
c
c  -- LAPACK driver routine (version 3.4.1) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     April 2012
c
c     .. Scalar Arguments ..
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
     *                   B( LDB, * ), BETA( * ), VL( LDVL, * ),
     *                   VR( LDVR, * ), WORK( * )


      call dggev(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai,
     *                  beta, vl, ldvl, vr, ldvr, work, lwork, info)

      return
      end

      SUBROUTINE xDSYGV(ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,
     *                  LWORK, INFO )
c
c  -- LAPACK driver routine (version 3.4.0) --
c  -- LAPACK is a software package provided by Univ. of Tennessee,    --
c  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
c     November 2011
c
c     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, ITYPE, LDA, LDB, LWORK, N
c     ..
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), W( * ), WORK( * )

      call dsygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work,
     *                  lwork, info)
      return
      end
