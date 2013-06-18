c interface to force availability of dggev and dsygv

      subroutine xdggev(kjobvl,kjobvr, n, a, b, alphar, alphai,
     *                  beta, vl, ldvl, vr, ldvr, work, lwork, info )
c
c     .. Scalar Arguments ..
      integer            kjobvl, kjobvr
      integer            info, ldvl, ldvr, lwork, n
c     ..
c     .. Array Arguments ..
      double precision   a(n, *), alphai(*), alphar(*),
     *                   b(n, *), beta(*), vl(ldvl, *),
     *                   vr(ldvr, *), work(*)

      character          jobvl, jobvr

      jobvl = 'NV'(kjobvl:kjobvl)
      jobvr = 'NV'(kjobvr:kjobvr)

      call dggev(jobvl, jobvr, n, a, n, b, n, alphar, alphai,
     *                  beta, vl, ldvl, vr, ldvr, work, lwork, info)

      return
      end

      subroutine xdsygv(itype, kjobz, kuplo, n, a, b,w,work,lwork,info)
c
c     .. Scalar Arguments ..
      integer            kjobz, kuplo
      integer            info, itype, lwork, n
c     ..
c     .. Array Arguments ..
      double precision   a(n, *), b(n, *), w(*), work(*)

      character          jobz, uplo

      jobz = 'NV'(kjobz:kjobz)
      uplo = 'UL'(kuplo:kuplo)

      call dsygv(itype, jobz, uplo, n, a, n, b, n, w, work,
     *                  lwork, info)
      return
      end
