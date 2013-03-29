c===========================================================
c$$$  
C$$$  Time-stamp: <liuminzhao 03/29/2013 09:44:20>
c$$$  Bivariate MLE using sigma
c$$$  
c===========================================================


CCCCCCCCCCCCCCCCCCCC
C        TARGET DELTA EQUATION 1
CCCCCCCCCCCCCCCCCCCC

      SUBROUTINE TargetEqn1f(delta, gamma, beta, sigma, tau, p, x, f)
      
      implicit none

      real*8 delta, gamma(2), beta(2), sigma(2), tau, p, x, f

      real*8 ans

C     OTHER C FUNCTION
      real*8 pnrm
      
C     TEMP
      real*8 quan, lp

      quan = gamma(1) + gamma(2) * x
      lp = beta(1) + beta(2) * x

      f = tau - p * pnrm((quan - delta - lp)/sigma(1),0.d0,1.d0,1,0)-
     &     (1 - p) * pnrm((quan - delta + lp)/sigma(2),0.d0,1.d0,1,0)

      return 
      end 

CCCCCCCCCCCCCCCCCCCC
C     SOLVE DELTA 1
CCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SolveDelta1f(param, tau, x, root)

      implicit none

      real*8 param(14)
      real*8 gamma(2), beta(2), sigma(2), tau, p, x, root

C     TEMP

      real*8 a, b, fa, fb, m, fm, tol

C     INITIAL 
      
      a = -30 
      b = 30 
      fa = 0
      fb = 0
      m = (a + b)/2
      fm = 0
      tol = 0.00001
      
      gamma(1) = param(1)
      gamma(2) = param(2)
      beta(1) = param(3)
      beta(2) = param(4)
      sigma(1) = param(5)
      sigma(2) = param(6)
      p = param(14)

      call TargetEqn1f(m, gamma, beta, sigma, tau, p, x, fm)
      call TargetEqn1f(b, gamma, beta, sigma, tau, p, x, fb)
      call TargetEqn1f(a, gamma, beta, sigma, tau, p, x, fa)

      do while(abs(fm) > tol)
         if (fm * fb < 0) then
            a = m
         else
            b = m
            fb = fm
         end if
         m = (a + b)/2
         call TargetEqn1f(m, gamma, beta, sigma, tau, p, x, fm)
      end do

      root = m
      return 
      end 


CCCCCCCCCCCCCCCCCCCC
C TARGET EQN FOR DELTA2     
CCCCCCCCCCCCCCCCCCCC

      subroutine targeteqn2f(delta2, gamma2, beta1, beta2, sigma1,
     &     sigma2, lambda, tau, p, x, delta1, f)
      implicit none
      real*8 delta2, gamma2(2), beta1(2), beta2(3)
      real*8 sigma1(2), sigma2, lambda, tau, p, x, delta1
      real*8 f

      real*8 pnrm
      real*8 p1, p2
      
      if (beta2(3) .ne. 0) then
         p1=1-pnrm(((delta2-gamma2(1)-x*gamma2(2)-beta2(1)-x*beta2(2))
     &        /beta2(3)-(delta1+beta1(1)+x*beta1(2)))/sigma1(1)/sqrt(
     &        sigma2**2/sigma1(1)**2/beta2(3)**2+1),
     &        0.d0,1.d0,1,0)
         p2=pnrm(((-delta2+gamma2(1)+x*gamma2(2)-beta2(1)-x*beta2(2))
     &        /beta2(3)-(delta1-beta1(1)-x*beta1(2)))/sigma1(2)/sqrt(
     &        lambda**2*sigma2**2/sigma1(2)**2/beta2(3)**2+1),
     &        0.d0,1.d0,1,0)
         if (beta2(3) > 0) then
            f = tau - p*p1 - (1 - p) * p2
         else
            f = tau - p * (1 - p1) - (1 - p) * (1 - p2)
         end if
      else
         p1=pnrm((gamma2(1)+x*gamma2(2)-delta2+beta2(1)+x*beta2(2))
     &        /sigma2, 0.d0, 1.d0, 1, 0)
         p2=pnrm((gamma2(1)+x*gamma2(2)-delta2-beta2(1)-x*beta2(2))
     &        /sigma2/lambda, 0.d0, 1.d0, 1, 0)
         f = tau - p * p1 - (1 - p) * p2
      end if

      return 
      end


CCCCCCCCCCCCCCCCCCCC
C SOLVE DELTA2     
CCCCCCCCCCCCCCCCCCCC

      subroutine solvedelta2f(param, tau, x, delta1, root)
      implicit none
      real*8 param(14)
      real*8 gamma1(2), beta1(2), sigma1(2)
      real*8 gamma2(2), beta2(3), sigma2, lambda, p
      real*8 tau, x, delta1, root

      real*8 a, b, fa, fb, m, fm, tol

C     INITIAL 
      
      a = -30 
      b = 30 
      fa = 0
      fb = 0
      m = (a + b)/2
      fm = 0
      tol = 0.00001      

      gamma1(1) = param(1)
      gamma1(2) = param(2)
      beta1(1) = param(3)
      beta1(2) = param(4)
      sigma1(1) = param(5)
      sigma1(2) = param(6)
      gamma2(1) = param(7)
      gamma2(2) = param(8)
      beta2(1) = param(9)
      beta2(2) = param(10)
      beta2(3) = param(11)
      sigma2 = param(12)
      lambda = param(13)
      p = param(14)

      call targeteqn2f(m, gamma2, beta1, beta2, sigma1,
     &     sigma2, lambda, tau, p, x, delta1, fm)
      call targeteqn2f(b, gamma2, beta1, beta2, sigma1,
     &     sigma2, lambda, tau, p, x, delta1, fb)
      call targeteqn2f(a, gamma2, beta1, beta2, sigma1,
     &     sigma2, lambda, tau, p, x, delta1, fa)

      do while(abs(fm) > tol)
         if (fm * fb < 0) then
            a = m
         else
            b = m
            fb = fm
         end if
         m = (a + b)/2
         call targeteqn2f(m, gamma2, beta1, beta2, sigma1,
     &        sigma2, lambda, tau, p, x, delta1, fm)
      end do

      root = m
      return 
      end 


CCCCCCCCCCCCCCCCCCCC
C NEGLOGLIKELIHOOD     
CCCCCCCCCCCCCCCCCCCC

      subroutine negloglikelihood(y, R, x, delta1, delta2,param,n,nll)
      implicit none
      integer n
      integer R(n)
      real*8 y(n, 2),x(n), delta1(n), delta2(n), param(14), nll
      real*8 dnrm
      real*8 gamma1(2), beta1(2), sigma1(2)
      real*8 gamma2(2), beta2(3), sigma2, lambda, p
      integer  i

      gamma1(1) = param(1)
      gamma1(2) = param(2)
      beta1(1) = param(3)
      beta1(2) = param(4)
      sigma1(1) = param(5)
      sigma1(2) = param(6)
      gamma2(1) = param(7)
      gamma2(2) = param(8)
      beta2(1) = param(9)
      beta2(2) = param(10)
      beta2(3) = param(11)
      sigma2 = param(12)
      lambda = param(13)
      p = param(14)

      nll = 0
      do i = 1, n
         if (R(i) .eq. 1) then
            nll=nll+dnrm(y(i,1),delta1(i)+beta1(1)+x(i)*beta1(2),
     &           sigma1(1),1)+dnrm(y(i,2),delta2(i)-beta2(1)-
     &           x(i)*beta2(2)-y(i,1)*beta2(3),sigma2,1)+log(p)
         else
            nll=nll+dnrm(y(i,1),delta1(i)-beta1(1)-x(i)*beta1(2),
     &           sigma1(2),1) + log(1 - p)
         end if
      end do
      nll = -nll
      return 
      end


CCCCCCCCCCCCCCCCCCCC
C       PARTIAL ALL
CCCCCCCCCCCCCCCCCCCC

      subroutine Partialf(param, tau, x, y, R, n, pp)
      implicit none
      integer n
      real*8 param(14)
      integer R(n)
      real*8 tau, x(n), y(n, 2), pp(14)

      real*8 epsilon, delta1p(n), delta2p(n), param1(14), param2(14)
      real*8 delta1m(n), delta2m(n)
      real*8 ll1, ll2
      integer i, j

      epsilon = 0.0003

      do i = 1, 14
         if (i .ne. 9 .and. i.ne.10.and.i.ne.11.and.i.ne.13) then
            do j = 1, 14
               param1(j) = param(j)
               param2(j) = param(j)
            end do
            param1(i) = param(i) + epsilon
            param2(i) = param(i) - epsilon

            do j = 1,n 
               call SolveDelta1f(param1, tau, x(j), delta1p(j))
               call SolveDelta1f(param2, tau, x(j), delta1m(j))
               call solvedelta2f(param1,tau,x(j),delta1p(j),delta2p(j))
               call solvedelta2f(param2,tau,x(j),delta1m(j),delta2m(j))
            end do

         call negloglikelihood(y, R, x, delta1p, delta2p,param1,n,ll1)
         call negloglikelihood(y, R, x, delta1m, delta2m,param2,n,ll2)

         pp(i) = (ll1 - ll2)/2/epsilon
         end if 
      end do
      return
      end

    

CCCCCCCCCCCCCCCCCCCC
C MAIN FUNCTION     
CCCCCCCCCCCCCCCCCCCC

      subroutine BiQRGradientf(y,R,x,tau,n,niter,param,paramsave)
      implicit none
      integer n, niter
      integer R(n)
      real*8 param(14), y(n,2), x(n),tau, paramsave(niter, 15)
      real*8 delta1(n), delta2(n)

      integer i, iter, j
      real*8 dif, nll0, nll, pp(14), alpha(14)

      
      do i = 1, 14
         pp(i) = 0
         alpha(i) = 0.0003
      end do
      alpha(14) = 0.0001

      dif = 1
      iter = 1
      nll0 = 0

      do i = 1, niter
         do j = 1, 15
            paramsave(i, j) = 0
         end do
      end do

      do while (dif > 0.00001 .and. iter .le. niter)
         call Partialf(param, tau, x, y, R, n, pp)
         do i = 1, 14
            param(i) = param(i) - alpha(i) * pp(i)
         end do
         param(5) = max(param(5), 0.01)
         param(6) = max(param(6), 0.01)
         param(12) = max(param(12), 0.01)
         param(14) = max(min(param(14), 0.99), 0.01)
         do j = 1, n
            call SolveDelta1f(param, tau, x(j), delta1(j))
            call solvedelta2f(param,tau,x(j),delta1(j),delta2(j))
         end do
         call negloglikelihood(y, R, x, delta1, delta2,param,n,nll)
         dif = abs(nll - nll0)
         nll0 = nll
         do i = 1, 14
            paramsave(iter, i) = param(i)
         end do
         paramsave(iter, 15) = nll
         iter = iter + 1
      end do
      return
      end

