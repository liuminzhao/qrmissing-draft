c===========================================================
c$$$
C$$$  Time-stamp: <liuminzhao 08/22/2013 17:46:29>
c$$$  2013/08/16 Bayesian MCMC for QRMissing Univariate single normal
c$$$
c===========================================================

CCCCCCCCCCCCCCCCCCCC
C     TARGET DELTA EQUATION 1
CCCCCCCCCCCCCCCCCCCC

      real*8 function TargetEqnf(delta,gamma,beta,sigma,p,tau,x,xdim)
      integer xdim, i
      real*8 delta, gamma(xdim), beta(xdim), sigma(2), p, tau, x(xdim)
C     OTHER C FUNCTION
      real*8 pnrm

C     TEMP
      real*8 quan, lp

      quan = dot_product(gamma, x)
      lp = dot_product(beta, x)

      targeteqnf=tau-p*pnrm(quan-delta - lp,0.d0,sigma(1),1,0)-
     &     (1-p)*pnrm(quan-delta+lp,0.d0,sigma(2),1,0)

      return
      end


C------------------------------
C     myzeroin
C------------------------------

      real*8 function myzeroin(gamma,beta, sigma, p,tau, x, xdim)
      integer xdim, i
      real*8 delta, gamma(xdim), beta(xdim), sigma(2), p, tau, x(xdim)
      real*8 targeteqnf
      real*8 myzeroin

      real*8 a, b, fa, fb, c, fc, tol, prevstep, newstep
      logical success
      real*8 dx, t1, cb, t2, pp, q
      integer maxit

      tol = 0.00001
      maxit = 40

      a = -100
      b = 100

      fa = targeteqnf(a, gamma, beta, sigma, p, tau, x, xdim)
      fb = targeteqnf(b, gamma, beta, sigma, p, tau, x, xdim)
      c = a
      fc = fa

C     First test if root is an endpoint
      if (fa .eq. 0.d0) then
         myzeroin = a
         return
      end if

      if (fb .eq. 0.d0) then
         myzeroin = b
         return
      end if

      do while (maxit .gt. 0)
         prevstep = b - a
         if (abs(fc) .lt. abs(fb)) then
            a = b; b = c; c = a;
            fa = fb; fb = fc; fc = fa;
         end if
         newstep = (c - b)/2

         if ((abs(newstep) .le. tol) .or. (fb .eq. 0.d0)) then
            myzeroin = b
            return
         end if

         if ((abs(prevstep) .ge. tol) .and. (abs(fa) .gt. abs(fb))) then
            cb = c - b
            if (a .eq. c) then
               t1 = fb/fa
               pp = cb*t1
               q = 1.d0 - t1
            else
               q = fa/fc; t1 = fb/fc; t2 = fb/fa;
               pp = t2 *(cb*q*(q -t1) - (b-a)*(t1-1.d0))
               q = (q - 1.d0) *(t1-1.d0) *(t2-1.d0)
            end if
            if (pp .gt. 0.d0) then
               q = -q
            else
               pp = -pp
            end if
            if ((pp .lt. (0.75*cb*q - abs(tol*q)/2)) .and.
     &           pp .le. abs(prevstep*q/2)) newstep = pp/q
         end if

         if (abs(newstep) .le. tol) then
            if (newstep .gt. 0.d0) then
               newstep = tol
            else
               newstep = -tol
            end if
         end if
         a = b
         fa = fb
         b = b + newstep
         fb = targeteqnf(b, gamma, beta, sigma, p, tau, x, xdim)
         if (((fb > 0.d0) .and. (fc > 0.d0)) .or.
     &        ((fb < 0.d0) .and. (fc < 0.d0))) then
            c = a
            fc = fa
         end if

         maxit = maxit - 1
      end do

      myzeroin = b
      return
      end


C------------------------------
C     mydelta
C------------------------------

      SUBROUTINE mydelta(x, gamma, beta, sigma, p, tau,n, xdim, delta)
      implicit none
      integer xdim, n
      real*8 x(n, xdim), gamma(xdim), beta(xdim),sigma(2)
      real*8 p, tau, delta(n)
      real*8 myzeroin
      integer i
      do i = 1, n
         delta(i) = myzeroin(gamma,beta, sigma, p,tau, x(i,:), xdim)
      end do
      return
      end
