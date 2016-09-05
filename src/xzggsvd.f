c Wrapper for Lapack's zggsvd

      subroutine xzggsvd(kjobu, kjobv, kjobq, m, n, p, k, l, a, lda, b,
     *                   ldb, alpha, beta, u, ldu, v, ldv, q, ldq,
     *                   work, lwork, rwork, iwork, info)
c
c     .. scalar arguments ..
      integer            kjobu, kjobv, kjobq
      integer            info, k, l, lda, ldb, ldq, ldu, ldv, m, n, p
      integer lwork
c     ..
c     .. array arguments ..
      integer            iwork(*)
      double precision   alpha(*), beta(*), rwork(*)
      double complex     a(lda, *), b(ldb, *),
     *                   q(ldq, *), u(ldu, *),
     *                   v(ldv, *), work( *)

      character(2)  cjobu, cjobv,cjobq
      parameter(cjobu='UN',cjobv='VN', cjobq='QN')
      character     jobq, jobu, jobv

c     select job parameters
      jobu = cjobu(kjobu:kjobu)
      jobv = cjobv(kjobv:kjobv)
      jobq = cjobq(kjobq:kjobq)
c     jobu = 'U'
c     jobv = 'V'
c     jobq = 'Q'
      call zggsvd3(jobu, jobv, jobq, m, n, p, k, l, A, lda, B, ldb,
     *            alpha, beta, u, ldu, v, ldv, q, ldq,
     *            work, lwork, rwork, iwork, info)

      return
      end
