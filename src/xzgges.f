
c interface to Lapack's zgges

      subroutine xzgges(kjobvsl, kjobvsr, kevsort, n, a, b,
     *                  sdim, alpha, beta, vsl, vsr,
     *                  work, lwork, rwork, bwork, info)

c     .. Scalar Arguments ..
      integer            kjobvsl, kjobvsr, kevsort
      integer            info, lwork, n, sdim
c     ..
c     .. Array Arguments ..
      logical            bwork(*)
      double precision   rwork(*)
      double complex     a(n,*), alpha(*),
     *                   b(n,*), beta(*), vsl(n,*),
     *                   vsr(n,*), work(*)
c     ..
c     .. Function Arguments ..
      logical            zelctg,zevzero,zrevneg,zrevpos,zevudi,zevudo
      external           zelctg,zevzero,zrevneg,zrevpos,zevudi,zevudo

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
          call zgges(jobvsl, jobvsr, 'N', zelctg, n, a, n, b, n,
     *               sdim, alpha, beta, vsl, n, vsr,
     *               n, work, lwork, rwork, bwork, info)
      case ('-')
          call zgges(jobvsl, jobvsr, 'S', zrevneg, n, a, n, b, n,
     *               sdim, alpha, beta, vsl, n, vsr,
     *               n, work, lwork, rwork, bwork, info)
      case ('+')
          call zgges(jobvsl, jobvsr, 'S', zrevpos, n, a, n, b, n,
     *               sdim, alpha, beta, vsl, n, vsr,
     *               n, work, lwork, rwork, bwork, info)
      case ('S')
          call zgges(jobvsl, jobvsr, 'S', zevudi,  n, a, n, b, n,
     *               sdim, alpha, beta, vsl, n, vsr,
     *               n, work, lwork, rwork, bwork, info)
      case ('B')
          call zgges(jobvsl, jobvsr, 'S', zevudo,  n, a, n, b, n,
     *               sdim, alpha, beta, vsl, n, vsr,
     *               n, work, lwork, rwork, bwork, info)
      case ('R')
          call zgges(jobvsl, jobvsr, 'S', zevzero, n, a, n, b, n,
     *               sdim, alpha, beta, vsl, n, vsr,
     *               n, work, lwork, rwork, bwork, info)

      end select
      return
      end

c The beta array is non-negative real according to the documentation
c sign(alpha/beta) == sign(real(alpha)*real(beta))

c for unordered result
      logical function zelctg(alpha,beta)
      double complex alpha,beta

      zelctg = .false.
      return
      end

c real eigenvalue
      logical function zevzero(alpha,beta)
      double complex alpha,beta
      double precision Rzero, Rtol
      parameter(Rzero=0.0d0, Rtol=100.0d0)
      double precision dlamch

      if(dble(beta) .eq. Rzero .and. dimag(beta) .eq. Rzero ) then
          zevzero = .false.
      else
          zevzero = abs(dimag(alpha/beta)) .le. Rtol*dlamch('E')
      endif
      return
      end

c real(ev) < 0
      logical function zrevneg(alpha,beta)
      double complex alpha,beta
      double precision Rzero
      parameter(Rzero=0.0d0)

      zrevneg = dble(alpha) * dble(beta) .lt. Rzero
      return
      end

c real(ev) > 0
      logical function zrevpos(alpha,beta)
      double complex alpha,beta
      double precision Rzero
      parameter(Rzero=0.0d0)

      zrevpos = dble(alpha) * dble(beta) .gt. Rzero
      return
      end

c abs(ev) < 1
      logical function zevudi(alpha,beta)
      double complex alpha,beta
      double precision Rzero
      parameter(Rzero=0.0d0)

      if(dble(beta) .eq. Rzero .and. dimag(beta) .eq. Rzero ) then
          zevudi = .false.
      else
          zevudi = abs(alpha) .lt. abs(beta)
      endif
      return
      end

c abs(ev) > 1
      logical function zevudo(alpha,beta)
      double complex alpha,beta
      double precision Rzero
      parameter(Rzero=0.0d0)

      if(dble(beta) .eq. Rzero .and. dimag(beta) .eq. Rzero ) then
          zevudo = .false.
      else
          zevudo = abs(alpha) .gt. abs(beta)
      endif
      return
      end
