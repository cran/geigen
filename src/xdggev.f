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

      character*2        cjobv
      parameter(cjobv='NV')
      character          jobvl, jobvr

c     if you change the cxxx values don't forget to adjust the R functions

      jobvl = cjobv(kjobvl:kjobvl)
      jobvr = cjobv(kjobvr:kjobvr)

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

      character*2        cjobz,cjobul
      parameter(cjobz='NV', cjobul='UL')
      character          jobz, uplo

c     if you change the cxxx values don't forget to adjust the R functions

      jobz = cjobz(kjobz:kjobz)
      uplo = cjobul(kuplo:kuplo)

      call dsygv(itype, jobz, uplo, n, a, n, b, n, w, work,
     *                  lwork, info)
      return
      end
