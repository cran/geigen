 
c interface to Lapack's zggev and zhegv
c input character arguments with length 255 converted to character*1.
 
      subroutine xzggev(kjobvl, kjobvr, n, a, lda, b, ldb, alpha, beta,
     *                  vl, ldvl, vr, ldvr, work, lwork, rwork, info )
c
c     .. Scalar Arguments ..
      integer            kjobvl, kjobvr
      integer            info, lda, ldb, ldvl, ldvr, lwork, n
c     ..
c     .. Array Arguments ..
      double precision   rwork(*) 
      complex*16         a(lda, *) , alpha(*), b(ldb, *),
     *                   beta(*), vl(ldvl, *), vr(ldvr, *),
     *                   work(*)

      character          jobvl, jobvr

      jobvl = 'NV'(kjobvl:kjobvl)
      jobvr = 'NV'(kjobvr:kjobvr)

      call zggev(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta,
     *           vl, ldvl, vr, ldvr, work, lwork, rwork, info )
     
      return
      end
      
      subroutine xzhegv(itype, kjobz, kuplo, n, a, lda, b, ldb, w, work,
     *                  lwork, rwork, info )

c     .. Scalar Arguments ..
      integer            kjobz, kuplo
      integer            info, itype, lda, ldb, lwork, n
c     ..
c     .. Array Arguments ..
      double precision   rwork(*), w(*)
      complex*16         a(lda, *), b(ldb, *), work(*)
      
      character          jobz, uplo

      jobz = 'NV'(kjobz:kjobz)
      uplo = 'UL'(kuplo:kuplo)

      call zhegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work,
     *           lwork, rwork, info )

      return
      end
