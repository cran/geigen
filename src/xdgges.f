
c interface to Lapack's dgges
c accept character arguments with declared length 255

      subroutine xdgges(kjobvsl, kjobvsr, kevsort, n, a, lda, b, ldb,
     *                  sdim, alphar, alphai, beta, vsl, ldvsl, vsr,
     *                  ldvsr, work, lwork, bwork, info )

c     copied from dgges with argument selctg removed
c
c     .. Scalar Arguments ..
      integer            kjobvsl, kjobvsr, kevsort
      integer            info, lda, ldb, ldvsl, ldvsr, lwork, n, sdim
c     ..
c     .. Array Arguments ..
      logical            bwork(*)
      double precision   a(lda,*), alphai(*), alphar(*),
     *                   b(ldb,*), beta(*), vsl(ldvsl,*),
     *                   vsr(ldvsr,*), work(*)
c     ..
c     .. Function Arguments ..
      logical            selctg,evzero,revneg,revpos,evudi,evudo
      external           selctg,evzero,revneg,revpos,evudi,evudo
      
      character         jobvsl, jobvsr, evsort
      
      jobvsl = 'NV'(kjobvsl:kjobvsl)
      jobvsr = 'NV'(kjobvsr:kjobvsr)
      evsort = 'N-+SBR'(kevsort:kevsort)
      
      select case (evsort)
      case ('N')
          call dgges(jobvsl, jobvsr, 'N', selctg, n, a, lda, b, ldb,
     *               sdim, alphar, alphai, beta, vsl, ldvsl, vsr,
     *               ldvsr, work, lwork, bwork, info)
      case ('-')
          call dgges(jobvsl, jobvsr, 'S', revneg, n, a, lda, b, ldb,
     *               sdim, alphar, alphai, beta, vsl, ldvsl, vsr,
     *               ldvsr, work, lwork, bwork, info)
      case ('+')
          call dgges(jobvsl, jobvsr, 'S', revpos, n, a, lda, b, ldb,
     *               sdim, alphar, alphai, beta, vsl, ldvsl, vsr,
     *               ldvsr, work, lwork, bwork, info)
      case ('S')
          call dgges(jobvsl, jobvsr, 'S', evudi,  n, a, lda, b, ldb,
     *               sdim, alphar, alphai, beta, vsl, ldvsl, vsr,
     *               ldvsr, work, lwork, bwork, info)
      case ('B')
          call dgges(jobvsl, jobvsr, 'S', evudo,  n, a, lda, b, ldb,
     *               sdim, alphar, alphai, beta, vsl, ldvsl, vsr,
     *               ldvsr, work, lwork, bwork, info)
      case ('R')
          call dgges(jobvsl, jobvsr, 'S', evzero, n, a, lda, b, ldb,
     *               sdim, alphar, alphai, beta, vsl, ldvsl, vsr,
     *               ldvsr, work, lwork, bwork, info)
      
      end select
      return
      end

c for unordered result
      logical function selctg(alphar,alphai,beta)
      double precision alphar,alphai,beta

      selctg = .false.
      return
      end

c real eigenvalue
      logical function evzero(alphar,alphai,beta)
      double precision alphar,alphai,beta
      double precision Rzero
      parameter(Rzero=0.0d0)
      
      evzero = alphai .eq. Rzero
      return
      end

c real(ev) < 0
      logical function revneg(alphar,alphai,beta)
      double precision alphar,alphai,beta
      double precision Rzero
      parameter(Rzero=0.0d0)

      revneg = alphar * beta .lt. Rzero
      return
      end

c real(ev) > 0
      logical function revpos(alphar,alphai,beta)
      double precision alphar,alphai,beta
      double precision Rzero
      parameter(Rzero=0.0d0)

      revpos = alphar * beta .gt. Rzero
      return
      end

c abs(ev) < 1
      logical function evudi(alphar,alphai,beta)
      double precision alphar,alphai,beta

      evudi = abs(cmplx(alphar,alphai)) .lt. abs(beta)
      return
      end

c abs(ev) > 1
      logical function evudo(alphar,alphai,beta)
      double precision alphar,alphai,beta

      evudo = abs(cmplx(alphar,alphai)) .gt. abs(beta)
      return
      end
