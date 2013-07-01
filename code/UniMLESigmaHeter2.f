c===========================================================
c$$$
C$$$  Time-stamp: <liuminzhao 07/01/2013 16:45:03>
c$$$  Univariate MLE using sigma
c$$$  Heter by Exponential exp(alpha0 + alpha1*x)
c===========================================================

CCCCCCCCCCCCCCCCCCCCC
C     TARGET DELTA EQUATION
CCCCCCCCCCCCCCCCCCCCC
      real*8 function targeteqnf(delta, param, x, tau)
      integer JMAX
      real*8 param(9)
      real*8 delta, gamma(2), beta(2), sigma(2), tau, p, x
      real*8 targeteqnf, alpha1(2), alpha2(2)
C     OTHER C FUNCTION
      real*8 pnrm
C     TEMP
      real*8 quan, lp

      gamma(1) = param(1)
      gamma(2) = param(2)
      beta(1) = param(3)
      beta(2) = param(4)
      alpha1(1) = param(5)
      alpha1(2) = param(6)
      alpha2(1) = param(7)
      alpha2(2) = param(8)
      p = param(9)

      sigma(1) = exp(alpha1(1) + alpha1(2) * x)
      sigma(2) = exp(alpha2(1) + alpha2(2) * x)

      quan = gamma(1) + gamma(2) * x
      lp = beta(1) + beta(2) * x

      targeteqnf = tau - p * pnrm((quan - delta - lp)/sigma(1),
     &     0.d0,1.d0,1,0)-
     &     (1 - p) * pnrm((quan - delta + lp)/sigma(2),0.d0,1.d0,1,0)

      return
      end

CCCCCCCCCCCCCCCCCCCC
C     BRACKET
CCCCCCCCCCCCCCCCCCCC
      SUBROUTINE zbrac(func, x1, x2, success, param, tau, x)
      implicit none
      integer ntry
      real*8 x1, x2, func, tau, x, param(9)
      external func
      real*8 factor
      integer j
      real*8 f1, f2
      logical success
      factor = 1.6
      f1 = func(x1, param, x, tau)
      f2 = func(x2, param, x, tau)
      success = .true.
      do j = 1, ntry
         if (f1*f2 .lt. 0.) return
         if (abs(f1).lt.abs(f2)) then
            x1 = x1 + factor*(x1-x2)
            f1 = func(x1, param, x, tau)
         else
            x2 = x2 + factor*(x2-x1)
            f2 = func(x2, param, x, tau)
         endif
      enddo
      success = .false.
      return
      end

CCCCCCCCCCCCCCCCCCCC
C     SOLVE DELTA
C     2013/07/01 change to function
CCCCCCCCCCCCCCCCCCCC
C      SUBROUTINE SolveDeltaH2f(param, tau, x, root)
      real*8 function root(param, tau, x)
      integer jmax
      real*8 root, param(9)
      real*8 tau, x
      real*8 tol
C     TEMP
      real*8 targeteqnf
      integer j
      real*8 a, b, fa, fb, m, fm
      logical success
      real*8 dx
C     INITIAL

      a = -0.1
      b = 0.1
      fa = 0
      fb = 0
      m = (a + b)/2
      fm = 0
      tol = 0.00001
      success = .true.
      jmax = 40

      call zbrac(targeteqnf, a, b, success, param, tau, x)

      dx = b - a

      fm = TargetEqnf(m, param, x, tau)
      fb = TargetEqnf(b, param, x, tau)
      fa = TargetEqnf(a, param, x, tau)

      if (fa*fb .ge. 0) print*, 'root must be bracketed'

      do j = 1, jmax
         dx = dx*0.5
         if (fm * fb < 0) then
            a = m
         else
            b = m
            fb = fm
         end if
         m = (a + b)/2
         root = m
         fm = TargetEqnf(m, param, x, tau)
         if (abs(dx) .lt. tol .or. abs(fm) .lt. tol) return
      end do
      print*, 'too many bisections in rootfind'
      end

CCCCCCCCCCCCCCCCCCCC
C     NEGLOGLIKELIHOOD
CCCCCCCCCCCCCCCCCCCC

      SUBROUTINE NegLogLikelihoodH2f(y, S, x, delta, param, n, nll)

      implicit none

      integer n
      integer S(n)
      real*8 y(n), x(n), delta(n), beta(2), nll, p, param(9)

      real*8 dnrm

C     TEMP
      integer i

      nll = 0
      beta(1) = param(3)
      beta(2) = param(4)
      p = param(9)

      do i = 1, n
         if (S(i) .eq. 1) then
            nll=nll+dnrm(y(i),delta(i)+beta(1)+beta(2)*x(i),
     &           exp(param(5) + param(6)*x(i)),1) + log(p)
         else
            nll=nll+dnrm(y(i),delta(i)-beta(1)-beta(2)*x(i),
     &           exp(param(7) + param(8)*x(i)),1)+ log(1 - p)
         end if
      end do

      nll = -nll

      return
      end

CCCCCCCCCCCCCCCCCCCC
C     PARTIAL ALL
CCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PartialH2f(param, tau, x, y, S, n, pp)

      implicit none

      integer n
      real*8 param(9)
      integer S(n)
      real*8 tau, x(n), y(n), pp(9)
      real*8 root

C     TEMP
      real*8 epsilon, delta1(n), delta2(n), param1(9), param2(9)
      real*8 ll1, ll2
      integer i, j

      epsilon = 0.0003

      do i = 1, 9
         do j = 1, 9
            param1(j) = param(j)
            param2(j) = param(j)
         end do
         param1(i) = param(i) + epsilon
         param2(i) = param(i) - epsilon

         do j = 1, n
            delta1(j) = root(param1, tau, x(j))
            delta2(j) = root(param2, tau, x(j))
         end do

         call NegLogLikelihoodH2f(y, S, x, delta1, param1, n, ll1)
         call NegLogLikelihoodH2f(y, S, x, delta2, param2, n, ll2)

         pp(i) = (ll1 - ll2)/2/epsilon

      end do

      return
      end

CCCCCCCCCCCCCCCCCCCC
C      MAIN FUNCTION
CCCCCCCCCCCCCCCCCCCC

      SUBROUTINE QRGradientH2f(y, S, x, tau, n, niter, param, paramsave)

      implicit none

      integer n, niter
      integer S(n)
      real*8 param(9), y(n), x(n), tau, delta(n), paramsave(niter, 10)
      real*8 root
C     TEMP
      integer i, iter, j
      real*8 dif, nll0, nll, pp(9), ppp(9)
      real*8 alpha(9), alphamax, alphamin, etap, etam

C     INITIAL
      param(1) = 0
      param(2) = 0
      param(3) = 0
      param(4) = 0
      param(5) = 0
      param(6) = 0
      param(7) = 0
      param(8) = 0
      param(9) = 0.5

      do i = 1, 9
         alpha(i) = 0.1
         ppp(i) = 1
      end do
      alphamax = 1
      alphamin = 0.000001
      etap = 1.2
      etam = 0.5

      dif = 1
      iter = 1
      nll0 = 0

      do i = 1, niter
         do j = 1, 10
            paramsave(i, j) = 0
         end do
      end do

      do while (dif > 0.00001 .and. iter .le. niter)
         call PartialH2f(param, tau, x, y, S, n, pp)
         do i = 1, 9
            if (pp(i) * ppp(i) > 0) then
               alpha(i) = min(alpha(i)*etap, alphamax)
               param(i) = param(i) - alpha(i)*pp(i)/abs(pp(i))
            else if (pp(i)*ppp(i) < 0) then
               alpha(i) = max(alpha(i)*etam, alphamin)
               param(i) = param(i) - alpha(i)*pp(i)/abs(pp(i))
               pp(i) = 0
            else if (ppp(i)*pp(i) .eq. 0) then
               if (pp(i) .eq. 0.d0) then
                  param(i) = param(i)
               else
                  param(i) = param(i) - alpha(i)*pp(i)/abs(pp(i))
               end if
            end if
            ppp(i) = pp(i)
         end do
         param(9) = max(min(param(9), 0.99), 0.001)
         do i = 1, n
            delta(i) =  root(param, tau, x(i))
         end do
         call NegLogLikelihoodH2f(y, S, x, delta, param, n, nll)
         dif = abs(nll - nll0)
         nll0 = nll
         do i = 1, 9
            paramsave(iter, i) = param(i)
         end do
         paramsave(iter, 10) = nll
         call progress(iter, niter)
         iter = iter + 1
      end do

      print*, '\n'

      return
      end

CCCCCCCCCCCCCCCCCCCC
C   PROGRESS BAR
CCCCCCCCCCCCCCCCCCCC

      subroutine progress(j, n)

      implicit none
      integer(kind=4) :: j,k,n
      character(len=18) :: bar="\r???% |          |"

      write(unit=bar(2:4),fmt="(i3)") 100*j/n
      do k = 1, j*10/n
         bar(7+k:7+k)="*"
      enddo

      write(*,'(a)',advance='no') bar

      return
      end
