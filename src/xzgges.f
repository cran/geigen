
c interface to Lapack's zgges

      subroutine xzgges(jobvsl, jobvsr, evsort, n, a, lda, b, ldb,
     *                  sdim, alpha, beta, vsl, ldvsl, vsr,
     *                  ldvsr, work, lwork, rwork, bwork, info )

c     copied from zgges with argument selctg removed
c
c     .. Scalar Arguments ..
      character*1        jobvsl, jobvsr, evsort
      integer            info, lda, ldb, ldvsl, ldvsr, lwork, n, sdim
c     ..
c     .. Array Arguments ..
      logical            bwork(*)
      double precision   rwork(*)
      double complex     a(lda,*), alpha(*),
     *                   b(ldb,*), beta(*), vsl(ldvsl,*),
     *                   vsr(ldvsr,*), work(*)
c     ..
c     .. Function Arguments ..
      logical            zelctg,zevzero,zrevneg,zrevpos,zevudi,zevudo
      external           zelctg,zevzero,zrevneg,zrevpos,zevudi,zevudo
      
      select case (evsort)
      case ("N","n")
          call zgges(jobvsl, jobvsr, "N", zelctg, n, a, lda, b, ldb,
     *               sdim, alpha, beta, vsl, ldvsl, vsr,
     *               ldvsr, work, lwork, rwork, bwork, info)
      case ("-")
          call zgges(jobvsl, jobvsr, "S", zrevneg, n, a, lda, b, ldb,
     *               sdim, alpha, beta, vsl, ldvsl, vsr,
     *               ldvsr, work, lwork, rwork, bwork, info)
      case ("+")
          call zgges(jobvsl, jobvsr, "S", zrevpos, n, a, lda, b, ldb,
     *               sdim, alpha, beta, vsl, ldvsl, vsr,
     *               ldvsr, work, lwork, rwork, bwork, info)
      case ("S","s")
          call zgges(jobvsl, jobvsr, "S", zevudi,  n, a, lda, b, ldb,
     *               sdim, alpha, beta, vsl, ldvsl, vsr,
     *               ldvsr, work, lwork, rwork, bwork, info)
      case ("B","b")
          call zgges(jobvsl, jobvsr, "S", zevudo,  n, a, lda, b, ldb,
     *               sdim, alpha, beta, vsl, ldvsl, vsr,
     *               ldvsr, work, lwork, rwork, bwork, info)
      case ("R","r")
          call zgges(jobvsl, jobvsr, "S", zevzero, n, a, lda, b, ldb,
     *               sdim, alpha, beta, vsl, ldvsl, vsr,
     *               ldvsr, work, lwork, rwork, bwork, info)
      
      end select
      return
      end

c for unordered result
      logical function zelctg(alpha,beta)
      double complex alpha,beta

      zelctg = .false.
      return
      end

c real eigenvalue
      logical function zevzero(alpha,beta)
      double complex alpha,beta
      double precision Rzero
      parameter(Rzero=0.0d0)
      
      zevzero = aimag(alpha) .eq. Rzero
      return
      end

c real(ev) < 0
      logical function zrevneg(alpha,beta)
      double complex alpha,beta
      double precision Rzero
      parameter(Rzero=0.0d0)

      zrevneg = real(alpha) * real(beta) .lt. Rzero
      return
      end

c real(ev) > 0
      logical function zrevpos(alpha,beta)
      double complex alpha,beta
      double precision Rzero
      parameter(Rzero=0.0d0)

      zrevpos = real(alpha) * real(beta) .gt. Rzero
      return
      end

c abs(ev) < 1
      logical function zevudi(alpha,beta)
      double complex alpha,beta

      zevudi = abs(alpha) .lt. abs(beta)
      return
      end

c abs(ev) > 1
      logical function zevudo(alpha,beta)
      double complex alpha,beta

      zevudo = abs(alpha) .gt. abs(beta)
      return
      end
