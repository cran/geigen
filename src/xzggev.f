 
c interface to Lapack's zggev and zhegv
 
      subroutine xzggev(kjobvl, kjobvr, n, a, b, alpha, beta,
     *                  vl, ldvl, vr, ldvr, work, lwork, rwork, info )
c
c     .. Scalar Arguments ..
      integer            kjobvl, kjobvr
      integer            info, ldvl, ldvr, lwork, n
c     ..
c     .. Array Arguments ..
      double precision   rwork(*) 
      complex*16         a(n, *) , alpha(*), b(n, *),
     *                   beta(*), vl(ldvl, *), vr(ldvr, *),
     *                   work(*)

      character*2        cjobv
      parameter(cjobv='NV')
      character          jobvl, jobvr

c     if you change the cxxx values don't forget to adjust the R functions

      jobvl = cjobv(kjobvl:kjobvl)
      jobvr = cjobv(kjobvr:kjobvr)

      call zggev(jobvl, jobvr, n, a, n, b, n, alpha, beta,
     *           vl, ldvl, vr, ldvr, work, lwork, rwork, info )
     
      return
      end
      
      subroutine xzhegv(itype, kjobz, kuplo, n, a, b, w, work,
     *                  lwork, rwork, info )

c     .. Scalar Arguments ..
      integer            kjobz, kuplo
      integer            info, itype, lwork, n
c     ..
c     .. Array Arguments ..
      double precision   rwork(*), w(*)
      complex*16         a(n, *), b(n, *), work(*)
      
      character*2        cjobz, cjobul
      parameter(cjobz='NV', cjobul='UL')
      character          jobz, uplo

c     if you change the cxxx values don't forget to adjust the R functions

      jobz = cjobz(kjobz:kjobz)
      uplo = cjobul(kuplo:kuplo)

      call zhegv(itype, jobz, uplo, n, a, n, b, n, w, work,
     *           lwork, rwork, info )

      return
      end
