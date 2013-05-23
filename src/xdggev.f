c interface to force availability of dggev and dsygv
c and to accept character arguments with declared length 255

      subroutine xdggev(kjobvl,kjobvr, n, a, lda, b, ldb,alphar,alphai,
     *                  beta, vl, ldvl, vr, ldvr, work, lwork, info )
c
c     .. Scalar Arguments ..
      integer            kjobvl, kjobvr
      integer            info, lda, ldb, ldvl, ldvr, lwork, n
c     ..
c     .. Array Arguments ..
      double precision   a(lda, *), alphai(*), alphar(*),
     *                   b(ldb, *), beta(*), vl( ldvl, *),
     *                   vr( ldvr, *), work(*)

      character          jobvl, jobvr

      jobvl = 'NV'(kjobvl:kjobvl)
      jobvr = 'NV'(kjobvr:kjobvr)

      call dggev(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai,
     *                  beta, vl, ldvl, vr, ldvr, work, lwork, info)

      return
      end

      subroutine xdsygv(itype, kjobz, kuplo, n, a, lda, b, ldb,w,work,
     *                  lwork, info )
c
c     .. Scalar Arguments ..
      integer            kjobz, kuplo
      integer            info, itype, lda, ldb, lwork, n
c     ..
c     .. Array Arguments ..
      double precision   a(lda, *), b( ldb, *), w(*), work(*)

      character          jobz, uplo

      jobz = 'NV'(kjobz:kjobz)
      uplo = 'UL'(kuplo:kuplo)

      call dsygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work,
     *                  lwork, info)
      return
      end
