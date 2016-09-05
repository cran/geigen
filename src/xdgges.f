
c interface to Lapack's dgges

      subroutine xdgges(kjobvsl, kjobvsr, kevsort, n, a, b,
     *                  sdim, alphar, alphai, beta, vsl, vsr,
     *                  work, lwork, bwork, info)

c     .. Scalar Arguments ..
      integer            kjobvsl, kjobvsr, kevsort
      integer            info, lwork, n, sdim
c     ..
c     .. Array Arguments ..
      logical            bwork(*)
      double precision   a(n,*), alphai(*), alphar(*),
     *                   b(n,*), beta(*), vsl(n,*),
     *                   vsr(n,*), work(*)
c     ..
c     .. Function Arguments ..
      logical            selctg,evzero,revneg,revpos,evudi,evudo
      external           selctg,evzero,revneg,revpos,evudi,evudo

      character(2)       cjobv
      character(6)       cevsort
      parameter(cjobv='NV', cevsort='N-+SBR')
      character          jobvsl, jobvsr, evsort

c     if you change the cxxx values don't forget to adjust the R functions

      jobvsl = cjobv(kjobvsl:kjobvsl)
      jobvsr = cjobv(kjobvsr:kjobvsr)
      evsort = cevsort(kevsort:kevsort)

      select case (evsort)
      case ('N')
          call dgges(jobvsl, jobvsr, 'N', selctg, n, a, n, b,n,
     *               sdim, alphar, alphai, beta, vsl, n, vsr,
     *               n, work, lwork, bwork, info)
      case ('-')
          call dgges(jobvsl, jobvsr, 'S', revneg, n, a, n, b, n,
     *               sdim, alphar, alphai, beta, vsl, n, vsr,
     *               n, work, lwork, bwork, info)
      case ('+')
          call dgges(jobvsl, jobvsr, 'S', revpos, n, a, n, b, n,
     *               sdim, alphar, alphai, beta, vsl, n, vsr,
     *               n, work, lwork, bwork, info)
      case ('S')
          call dgges(jobvsl, jobvsr, 'S', evudi,  n, a, n, b, n,
     *               sdim, alphar, alphai, beta, vsl, n, vsr,
     *               n, work, lwork, bwork, info)
      case ('B')
          call dgges(jobvsl, jobvsr, 'S', evudo,  n, a, n, b, n,
     *               sdim, alphar, alphai, beta, vsl, n, vsr,
     *               n, work, lwork, bwork, info)
      case ('R')
          call dgges(jobvsl, jobvsr, 'S', evzero, n, a, n, b, n,
     *               sdim, alphar, alphai, beta, vsl, n, vsr,
     *               n, work, lwork, bwork, info)

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
      double precision Rzero, Rtol
      parameter(Rzero=0.0d0, Rtol=100d0)
      double precision dlamch

      if( beta .eq. Rzero ) then
          evzero = .false.
      else
          evzero = abs(alphai/beta) .le. Rtol*dlamch('E')
      endif
      return
      end

c real(ev) < 0
      logical function revneg(alphar,alphai,beta)
      double precision alphar,alphai,beta
      double precision Rzero
      parameter(Rzero=0.0d0)

      if( beta .eq. Rzero ) then
          revneg = .false.
      else
          revneg = alphar * beta .lt. Rzero
      endif
      return
      end

c real(ev) > 0
      logical function revpos(alphar,alphai,beta)
      double precision alphar,alphai,beta
      double precision Rzero
      parameter(Rzero=0.0d0)

      if( beta .eq. Rzero ) then
          revpos = .false.
      else
          revpos = alphar * beta .gt. Rzero
      endif
      return
      end

c abs(ev) < 1
      logical function evudi(alphar,alphai,beta)
      double precision alphar,alphai,beta
      double precision Rzero
      parameter(Rzero=0.0d0)

      if( beta .eq. Rzero ) then
          evudi = .false.
      else
          evudi = abs(dcmplx(alphar,alphai)) .lt. abs(beta)
      endif
      return
      end

c abs(ev) > 1
      logical function evudo(alphar,alphai,beta)
      double precision alphar,alphai,beta
      double precision Rzero
      parameter(Rzero=0.0d0)

      if( beta .eq. Rzero ) then
          evudo = .false.
      else
          evudo = abs(dcmplx(alphar,alphai)) .gt. abs(beta)
      endif
      return
      end
