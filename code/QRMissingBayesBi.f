c===========================================================
c$$$
C$$$  Time-stamp: <liuminzhao 08/22/2013 16:34:04>
c$$$  2013/08/22 Bayesian MCMC for QRMissing Bivariate single normal
c$$$
c===========================================================

CCCCCCCCCCCCCCCCCCCC
C     TARGET DELTA1 EQUATION
CCCCCCCCCCCCCCCCCCCC

      real*8 function TargetEqn1f(delta1,gamma1,beta1,sigma1,
     &     p,tau,x,xdim)
      integer xdim, i
      real*8 delta1, gamma1(xdim), beta1(xdim), sigma1(2),
     &     p, tau, x(xdim)
      real*8 pnrm

      real*8 quan, lp

      quan = dot_product(gamma1, x)
      lp = dot_product(beta1, x)

      targeteqn1f=tau-p*pnrm(quan-delta1 - lp,0.d0,sigma1(1),1,0)-
     &     (1-p)*pnrm(quan-delta1+lp,0.d0,sigma1(2),1,0)

      return
      end


C------------------------------
C     myzeroin1 : solve for delta1
C------------------------------

      real*8 function myzeroin1(gamma1,beta1, sigma1, p,tau, x, xdim)
      integer xdim, i
      real*8 delta1, gamma1(xdim), beta1(xdim), sigma1(2), p,
     &     tau, x(xdim)
      real*8 targeteqn1f
      real*8 myzeroin1

      real*8 a, b, fa, fb, c, fc, tol, prevstep, newstep
      logical success
      real*8 dx, t1, cb, t2, pp, q
      integer maxit

c$$$      print*, gamma1, beta1, sigma1
c$$$      print*, p, tau, xdim

      tol = 0.00001
      maxit = 40

      a = -100
      b = 100

      fa = targeteqn1f(a, gamma1, beta1, sigma1, p, tau, x, xdim)
      fb = targeteqn1f(b, gamma1, beta1, sigma1, p, tau, x, xdim)
      c = a
      fc = fa

C     First test if root is an endpoint
      if (fa .eq. 0.d0) then
         myzeroin1 = a
         return
      end if

      if (fb .eq. 0.d0) then
         myzeroin1 = b
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
            myzeroin1 = b
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
         fb = targeteqn1f(b, gamma1, beta1, sigma1, p, tau, x, xdim)
         if (((fb > 0.d0) .and. (fc > 0.d0)) .or.
     &        ((fb < 0.d0) .and. (fc < 0.d0))) then
            c = a
            fc = fa
         end if

         maxit = maxit - 1
      end do

      myzeroin1 = b
      return
      end


C------------------------------
C     Target function for Delta2 given other parameters
C------------------------------

      real*8 function targeteqn2f(d2, d1, gamma1, beta1, sigma1,
     &     gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
     &     p, tau, x, xdim)
      integer xdim, i
      real*8 d2, gamma1(xdim), beta1(xdim), sigma1(2)
      real*8 gamma2(xdim), beta2sp(xdim), sigma21, sigma20, sigma21sp
      real*8 betay, betay0, betaysp, p, tau, x(xdim), d1
      real*8 myzeroin1
      real*8 p1, p2, pnrm

c$$$      print*, gamma1, beta1, sigma1, gamma2, beta2sp, sigma21
c$$$      print*, sigma21sp, betay, betaysp, p, tau, xdim
c$$$      print*, d1, 'x', x

      sigma20 = sigma21 * exp(sigma21sp)
      betay0 = betay + betaysp

      lp1 = dot_product(beta1, x)
      lp2 = dot_product(beta2sp, x)
      quan2 = dot_product(gamma2, x)

      if (betay .ne. 0) then
         p1=pnrm(((-d2+ quan2
     &        )/betay-(d1+lp1))/
     &        sigma1(1)/sqrt(
     &        sigma21**2/
     &        sigma1(1)**2/betay**2+1),
     &        0.d0,1.d0,1,0)
      else
         p1=pnrm((quan2-d2
     &        )/sigma21, 0.d0, 1.d0, 1, 0)
      end if

      if (betay0 .ne. 0) then
         p2=pnrm(((-d2+quan2 - lp2
     &        )/betay0-(d1-lp1))/
     &        sigma1(2)/sqrt(
     &        sigma20**2/
     &        sigma1(2)**2/betay0**2+1),
     &        0.d0,1.d0,1,0)
      else
         p2=pnrm((quan2-d2-lp2
     &        )/sigma20, 0.d0, 1.d0, 1, 0)
      end if

      if (betay < 0) then
         p1 = 1 - p1
      end if

      if (betay0 < 0) then
         p2 = 1 - p2
      end if

      targeteqn2f = tau - p * p1 - (1 - p) * p2

c$$$      print*, p1, p2, p, targeteqn2f

      return
      end

C------------------------------
C     Solver for Delta2
C------------------------------

      real*8 function myzeroin2(gamma1, beta1, sigma1,
     &     gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
     &     p, tau, x, xdim, d1)
      integer xdim, i
      real*8 d2, gamma1(xdim), beta1(xdim), sigma1(2)
      real*8 gamma2(xdim), beta2sp(xdim), sigma21, sigma20, sigma21sp
      real*8 betay, betay0, betaysp, p, tau, x(xdim), d1
      real*8 p1, p2, pnrm
      real*8 a, b, fa, fb, c, fc, tol, prevstep, newstep
      real*8 targeteqn2f
      logical success
      real*8 dx, t1, cb, t2, pp, q
      integer maxit

c$$$      print*, gamma1, beta1, sigma1, gamma2, beta2sp, sigma21
c$$$      print*, sigma21sp, betay, betaysp, p, tau, xdim
c$$$
c$$$      print*, d1
c$$$      myzeroin2 = d1

      tol = 0.00001
      maxit = 40

      a = -100
      b = 100

      fa = targeteqn2f(a, d1, gamma1, beta1, sigma1,
     &     gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
     &     p, tau, x, xdim)
      fb = targeteqn2f(b, d1, gamma1, beta1, sigma1,
     &     gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
     &     p, tau, x, xdim)
      c = a
      fc = fa

c$$$      print*, fa, fb

C     First test if root is an endpoint
      if (fa .eq. 0.d0) then
         myzeroin2 = a
c$$$         print*, myzeroin2
c$$$         stop 0
         return
      end if

      if (fb .eq. 0.d0) then
         myzeroin2 = b
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
            myzeroin2 = b
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
         fb = targeteqn2f(b, d1, gamma1, beta1, sigma1,
     &     gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
     &     p, tau, x, xdim)
         if (((fb > 0.d0) .and. (fc > 0.d0)) .or.
     &        ((fb < 0.d0) .and. (fc < 0.d0))) then
            c = a
            fc = fa
         end if

         maxit = maxit - 1
      end do

      myzeroin2 = b
      return
      end

C------------------------------
C     mydelta : solve delta1, delta2 for all X
C------------------------------

      SUBROUTINE mydelta2(x, gamma1, beta1, sigma1,
     &     gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
     &     p, tau, n, xdim, delta)
      implicit none
      integer xdim, n
      real*8 x(n, xdim), gamma1(xdim), beta1(xdim),sigma1(2)
      real*8 gamma2(xdim), beta2sp(xdim), sigma21, sigma21sp
      real*8 betay, betaysp
      real*8 p, tau, delta(n, 2)
      real*8 myzeroin2, myzeroin1
      integer i

c$$$      print*, gamma1, beta1, sigma1, gamma2, beta2sp, sigma21
c$$$      print*, sigma21sp, betay, betaysp, p, tau, n, xdim

      do i = 1, n
c$$$         print*, x(i, :), delta(i, :)
         delta(i,1) = myzeroin1(gamma1,beta1,sigma1, p,
     &        tau, x(i,:), xdim)
         delta(i,2) = myzeroin2(gamma1, beta1, sigma1,
     &        gamma2, beta2sp, sigma21, sigma21sp, betay, betaysp,
     &        p, tau, x(i,:), xdim, delta(i, 1))
      end do
      return
      end
