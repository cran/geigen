c Wrapper for Lapack's dggsvd

      subroutine xdggsvd(kjobu, kjobv, kjobq, m, n, p, k, l, a, lda, b,
     *                   ldb, alpha, beta, u, ldu, v, ldv, q, ldq,
     *                   work, iwork, info)
c
c  -- lapack driver routine (version 3.4.0) --
c  -- lapack is a software package provided by univ. of tennessee,    --
c  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
c     november 2011
c
c     .. scalar arguments ..
      integer            kjobu, kjobv, kjobq
      integer            info, k, l, lda, ldb, ldq, ldu, ldv, m, n, p
c     ..
c     .. array arguments ..
      integer            iwork(*)
      double precision   a(lda, *), alpha(*), b(ldb, *),
     *                   beta(*), q(ldq, *), u(ldu, *),
     *                   v(ldv, *), work( *)

      character*2  cjobu, cjobv,cjobq
      parameter(cjobu='UN',cjobv='VN', cjobq='QN')
      character*1  jobq, jobu, jobv

c     select job parameters
      jobu = cjobu(kjobu:kjobu)
      jobv = cjobv(kjobv:kjobv)
      jobq = cjobq(kjobq:kjobq)
c     jobu = 'U'
c     jobv = 'V'
c     jobq = 'Q'
      call dggsvd(jobu, jobv, jobq, m, n, p, k, l, A, lda, B, ldb,
     *            alpha, beta, u, ldu, v, ldv, q, ldq,
     *            work, iwork, info)

      return
      end
