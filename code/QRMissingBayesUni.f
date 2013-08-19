c===========================================================
c$$$
C$$$  Time-stamp: <liuminzhao 08/19/2013 00:59:08>
c$$$  2013/08/16 Bayesian MCMC for QRMissing Univariate single normal
c$$$
c===========================================================

CCCCCCCCCCCCCCCCCCCC
C     TARGET DELTA EQUATION 1
CCCCCCCCCCCCCCCCCCCC

      real*8 function TargetEqnf(delta,gamma,beta,sigma,p,tau,x,xdim)
      integer xdim, i
      real*8 delta, gamma(xdim), beta(xdim), sigma(2), p, tau, x(xdim)
      real*8 targeteqnf
C     OTHER C FUNCTION
      real*8 pnrm

C     TEMP
      real*8 quan, lp

      quan = 0
      lp = 0

      do i = 1, xdim
         quan = quan + gamma(i) * x(i)
         lp = lp + beta(i) * x(i)
      end do
      targeteqnf=tau-p*pnrm(quan-delta - lp,0.d0,sigma(1),1,0)-
     &     (1-p)*pnrm(quan-delta+lp,0.d0,sigma(2),1,0)

      return
      end

CCCCCCCCCCCCCCCCCCCC
C     BRACKET 1
CCCCCCCCCCCCCCCCCCCC
      SUBROUTINE zbrac1(func,x1,x2,success,gamma,
     &     beta,sigma,p, tau,x,xdim)
      implicit none
      integer ntry, xdim
      real*8 x1, x2, func
      real*8 delta, gamma(xdim), beta(xdim), sigma(2), p, tau, x(xdim)
      external func
      real*8 factor
      integer j
      real*8 f1, f2
      logical success

      factor = 1.6
      ntry = 50
      f1 = func(x1,gamma,beta,sigma ,p, tau, x, xdim)
      f2 = func(x2, gamma,beta,sigma ,p, tau, x, xdim)
      success = .true.
      do j = 1, ntry
         if (f1*f2 .lt. 0.) return
         if (abs(f1).lt.abs(f2)) then
            x1 = x1 + factor*(x1-x2)
            f1 = func(x1, gamma,beta,sigma, p,tau, x, xdim)
         else
            x2 = x2 + factor*(x2-x1)
            f2 = func(x2, gamma,beta,sigma, p,tau, x, xdim)
         endif
      enddo
      success = .false.
      print*, dot_product(gamma, x), dot_product(beta, x)
      print*, sigma, p, tau, x
      stop 0
      return
      end

CCCCCCCCCCCCCCCCCCCC
C     SOLVE DELTA 1
CCCCCCCCCCCCCCCCCCCC
      real*8 function root1(gamma,beta, sigma, p,tau, x, xdim)

      integer imax
      integer xdim, i
      real*8 delta, gamma(xdim), beta(xdim), sigma(2), p, tau, x(xdim)
      real*8 root1
      real*8 targeteqnf
C     TEMP

      real*8 a, b, fa, fb, m, fm, tol
      logical success
      real*8 dx

C     INITIAL

      a = -100
      b = 100
      fa = 0
      fb = 0
      m = (a + b)/2
      fm = 0
      tol = 0.00001
      success = .true.
      imax = 40

c$$$      call zbrac1(targeteqnf, a, b, success,
c$$$     &     gamma, beta, sigma, p,tau, x, xdim)

      if (.not. success) print*, 'fail to bracket delta1'

      dx = b - a

      fm = TargetEqnf(m, gamma,beta,sigma, p,tau, x,xdim)
      fb = TargetEqnf(b, gamma,beta,sigma, p,tau, x,xdim)
      fa = TargetEqnf(a, gamma,beta,sigma, p,tau, x,xdim)

C      print*, fa, fb, fm

      if (fa*fb .ge. 0) print*, 'root must be bracketed'

      if (abs(fm) .lt. tol) then
         root1 = m
         return
      end if

      do i = 1, imax
         if (fm * fb < 0) then
            a = m
         else
            b = m
            fb = fm
         end if
         m = (a + b)/2
         root1 = m
         fm = TargetEqnf(m, gamma,beta,sigma, p, tau, x, xdim)
         if (abs(fm) .lt. tol) return
      end do
      print*, 'too many bisections in rootfind'
      print*, dot_product(gamma, x), dot_product(beta, x)
      print*,  a, b, m, x, beta, gamma
      stop 0
      end


C------------------------------
C     likelihood
C------------------------------

      real*8 function likelihood(gamma, beta, sigma, p, tau,
     &     x, xdim, y, R, n)

      integer xdim, n
      real*8 y(n), R(n), x(n, xdim)
      real*8 gamma(xdim), beta(xdim), sigma(2), p, tau

      real*8 delta(n), mu

      integer i, j, k
      real*8 likelihood, root1, d1, tmp, lp

      real*8 dnrm

      likelihood = 0.d0

      do i = 1, n
c$$$         print*, gamma, beta, sigma, p, tau, x(i,:), xdim

         d1 = root1(gamma,beta, sigma, p, tau, x(i,:), xdim)
         lp = 0.d0
         do j  = 1, xdim
            lp = lp + x(i, j) * beta(j)
         end do
         if (R(i) .eq. 1) then
            mu = d1 + lp
            likelihood = likelihood + dnrm(y(i), mu, sigma(1), 1)
         else
            mu = d1 - lp
            likelihood = likelihood + dnrm(y(i), mu, sigma(2), 1)
         end if
      end do
      end



C------------------------------
C     Main function
C------------------------------

      subroutine qrmissingbayesuni(n, xdim, x, y, R,tau,
     &     gammapm, gammapv, betapm, betapv,
     &     a, b, c, d,
     &     mcmc, nsave,
     &     gammasave, betasave, sigmasave, psave)

      implicit none

C     DATA
      integer n, xdim
      real*8 x(n, xdim), y(n), tau
      integer R(n)

C     prior
      real*8 gammapm(xdim), gammapv(xdim)
      real*8 betapm(xdim), betapv(xdim)
      real*8 a, b, c, d

C     mcmc
      integer mcmc(4), nburn, nskip, nsave, ndisp

C     store
      real*8 gammasave(nsave, xdim), betasave(nsave, xdim)
      real*8 sigmasave(nsave, 2), psave(nsave)

C     parameters:
      real*8 gamma(xdim), beta(xdim), sigma(2), p

C     function
      real*8 likelihood, myrnorm, dnrm, myrunif, dgamma2, mydbeta

C     tmp
      integer i, j, k, isave, skipcount, dispcount, nscan,iscan
      real*8 loglikeo, loglikec, logprioro, logpriorc
      real*8 logcgkc, logcgko
      real*8 gammac(xdim), betac(xdim), sigmac(2), pc
      real*8 tmp
      real*8 thetac(2), theta(2)
      real*8 ratio

C     tuning
      real*8 sdgamma(xdim), sdbeta(xdim), arate
      integer attgamma(xdim), attbeta(xdim), accgamma(xdim),
     &     accbeta(xdim)

C     time
      real*8 sec00, sec0, sec1, sec


C------------------------------
C     initial
C------------------------------

      nsave=mcmc(1)
      nskip=mcmc(2)
      nburn=mcmc(3)
      ndisp=mcmc(4)

      do i = 1, xdim
         gamma(i) = 0.d0
         beta(i) = 0.d0
         sdgamma(i) = 1.d0
         sdbeta(i) = 1.d0
         accgamma(i) = 0
         accbeta(i) = 0
         attgamma(i) = 0
         attbeta(i) = 0
      end do

      sigma(1) = 1
      sigma(2) = 1

      p = 0.5

      arate = 0.25

c$$$      print*, gamma
c$$$      print*, beta, sigma, p, tau, xdim, n

C------------------------------
C     start mcmc
C------------------------------

      isave = 0
      skipcount = 0
      dispcount = 0
      nscan = nburn + nskip * nsave

      call cpu_time(sec0)
      sec00 = 0.d0

C     first
      loglikeo = 0.d0
      loglikeo = likelihood(gamma, beta, sigma, p, tau,
     &     x, xdim, y, R, n)

C      print*, loglikeo

C------------------------------
C     rolling
C------------------------------

      do iscan = 1, nscan

C------------------------------
C     gamma
C------------------------------

         do i = 1, xdim
            attgamma(i) = attgamma(i) + 1
            gammac = gamma
            gammac(i) = myrnorm(gamma(i), sdgamma(i))
            logpriorc = dnrm(gammac(i), gammapm(i),
     &           gammapv(i), 1)
            logprioro = dnrm(gamma(i), gammapm(i),
     &           gammapv(i), 1)

c$$$         print*, gammac

            loglikec = likelihood(gammac, beta, sigma, p, tau,
     &           x, xdim, y, R, n)

            ratio = loglikec + logpriorc - loglikeo - logprioro

c$$$         print*, loglikec, loglikeo, logpriorc, logprioro
c$$$         print*, exp(ratio)

            if (log(dble(myrunif(0.d0, 1.d0))).lt.ratio) then
               loglikeo = loglikec
               accgamma(i) = accgamma(i) + 1
               gamma = gammac
            end if

         end do
C------------------------------
C     beta
C------------------------------

         do i = 1, xdim
            attbeta(i) = attbeta(i) + 1
            betac = beta
            betac(i) = myrnorm(beta(i), sdbeta(i))
            logpriorc = dnrm(betac(i), betapm(i),
     &           betapv(i), 1)
            logprioro = dnrm(beta(i), betapm(i),
     &           betapv(i), 1)

            loglikec = likelihood(gamma, betac, sigma, p, tau,
     &           x, xdim, y, R, n)

c$$$         print*, loglikec, loglikeo, logpriorc, logprioro

            ratio = loglikec + logpriorc - loglikeo - logprioro

            if (log(dble(myrunif(0.d0, 1.d0))).lt.ratio) then
               loglikeo = loglikec
               accbeta(i) = accbeta(i) + 1
               beta = betac
            end if
         end do

         print*, beta

C------------------------------
C     sigma
C------------------------------

         theta(1) = log(sigma(1))
         theta(2) = log(sigma(2))
         thetac(1) = myrnorm(theta(1), 1.d0)
         thetac(2) = myrnorm(theta(2), 1.d0)
         logcgkc = -sum(theta)
         logcgko = -sum(thetac)
         sigmac(1) = exp(thetac(1))
         sigmac(2) = exp(thetac(2))

         sigmac(1) = 1
         sigmac(2) = 1

C         print*, sigmac

C     likelihood

         loglikec = likelihood(gamma, beta, sigmac, p, tau,
     &        x, xdim, y, R, n)

C     log prior

         logpriorc = dgamma2(sigmac(1), a/2, b/2, 1)
     &        + dgamma2(sigmac(2), a/2, b/2, 1)
         logprioro = dgamma2(sigma(1), a/2, b/2, 1)
     &        + dgamma2(sigma(2), a/2, b/2, 1)

c$$$         print*, loglikec, loglikeo, logpriorc, logprioro
c$$$         print*, exp(ratio)

         ratio=loglikec + logpriorc -loglikeo -logprioro+
     &        logcgkc - logcgko

         if (log(dble(myrunif(0.d0, 1.d0))).lt.ratio) then
            loglikeo = loglikec
            do i = 1, 2
               sigma(i) = sigmac(i)
            end do
         end if

C------------------------------
C     p
C------------------------------

C         pc = myrnorm(p, 0.1d0)
         pc = 0.46
         pc = max(pc, 0.01)
         pc = min(pc, 0.99)

C         print*, pc

         loglikec = likelihood(gamma, beta, sigma, pc, tau,
     &        x, xdim, y, R, n)

C     log prior
         logpriorc = mydbeta(pc, c/2, d/2, 1)
         logprioro = mydbeta(p, c/2, d/2, 1)

c$$$         print*, loglikec, loglikeo, logpriorc, logprioro
c$$$         print*, exp(ratio)

         ratio=loglikec + logpriorc -loglikeo -logprioro

         if (log(dble(myrunif(0.d0, 1.d0))).lt.ratio) then
            loglikeo = loglikec
            p = pc
         end if

C------------------------------
C     tuning
C------------------------------

         if ((attgamma(1) .ge. 100) .and. (iscan .le. nburn)) then
            do i = 1, xdim
               if (dble(accgamma(i))/dble(attgamma(i)) .gt. arate) then
                  sdgamma(i) = sdgamma(i) * 2
               else
                  sdgamma(i) = sdgamma(i) / 2.d0
               end if
               if (dble(accbeta(i))/dble(attbeta(i)) .gt. arate) then
                  sdbeta(i) = sdbeta(i) * 2
               else
                  sdbeta(i) = sdbeta(i) / 2.d0
               end if
c$$$               tmp = dble(accgamma)/dble(attgamma)
c$$$               print*, tmp
c$$$               tmp = dble(accbeta)/dble(attbeta)
c$$$               print*, tmp
c$$$               print*, sdgamma, sdbeta
               attgamma(i) = 0
               accgamma(i) = 0
               attbeta(i) = 0
               accbeta(i) = 0
            end do
         end if

C------------------------------
C     save
C------------------------------

         if (iscan .gt. nburn) then
            skipcount = skipcount + 1
            if (skipcount .ge. nskip) then
               isave = isave + 1
               dispcount = dispcount + 1
               gammasave(isave, :) = gamma
               betasave(isave, :) = beta
               sigmasave(isave, :) = sigma
               psave(isave) = p
               skipcount = 0
               if (dispcount .ge. ndisp) then
                  call cpu_time(sec1)
                  sec00=sec00 + (sec1-sec0)
                  sec=sec00
                  sec0 = sec1
                  print*,  isave, " out of " , nsave,
     &                 " for time ", floor(sec)
                  dispcount = 0
               end if
            end if
         end if
      end do

      return
      end
