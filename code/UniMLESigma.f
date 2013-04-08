c===========================================================
c$$$  
C$$$  Time-stamp: <liuminzhao 04/08/2013 11:26:30>
c$$$  Univariate MLE using sigma
c$$$  
c===========================================================

CCCCCCCCCCCCCCCCCCCCC      
C     TARGET DELTA EQUATION
CCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE TargetEqnf(delta,gamma,beta,sigma,heter, tau, p, x, f)
      
      implicit none

      real*8 delta, gamma(2), beta(2), sigma(2), tau, p, x, f
      real*8 heter(2)

      real*8 ans

C     OTHER C FUNCTION
      real*8 pnrm
      
C     TEMP
      real*8 quan, lp

      quan = gamma(1) + gamma(2) * x
      lp = beta(1) + beta(2) * x

      f = tau - p * pnrm((quan - delta - lp)/(1+heter(1)*x)/sigma(1)
     &     ,0.d0,1.d0,1,0)-(1 - p) * pnrm((quan - delta + lp)/
     &     (1+heter(2)*x)/sigma(2),0.d0,1.d0,1,0)

      return 
      end 


CCCCCCCCCCCCCCCCCCCC
C     SOLVE DELTA
CCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SolveDeltaf(param, tau, x, root)

      implicit none

      real*8 param(9)
      real*8 gamma(2), beta(2), sigma(2), tau, p, x, root
      real*8 heter(2)
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
      p = param(7)
      heter(1) = param(8)
      heter(2) = param(9)

      call TargetEqnf(m,gamma,beta,sigma,heter,tau, p, x, fm)
      call TargetEqnf(b,gamma,beta,sigma,heter,tau, p, x, fb)
      call TargetEqnf(a,gamma,beta,sigma,heter,tau, p, x, fa)

      do while(abs(fm) > tol)
         if (fm * fb < 0) then
            a = m
         else
            b = m
            fb = fm
         end if
         m = (a + b)/2
         call TargetEqnf(m,gamma,beta,sigma,heter,tau, p, x, fm)
      end do

      root = m
      return 
      end 


CCCCCCCCCCCCCCCCCCCC
C     NEGLOGLIKELIHOOD
CCCCCCCCCCCCCCCCCCCC

      SUBROUTINE NegLogLikelihoodf(y, S, x, delta, param, n, nll)
      
      implicit none
      
      integer n
      integer S(n)
      real*8 y(n), x(n), delta(n), beta(2), sigma(2), nll, p, param(9)
      real*8 heter(2)
      real*8 dnrm
      
C     TEMP
      integer i

      nll = 0
      beta(1) = param(3)
      beta(2) = param(4)
      sigma(1) = param(5)
      sigma(2) = param(6)
      p = param(7)
      heter(1) = param(8)
      heter(2) = param(9)

      do i = 1, n
         if (S(i) .eq. 1) then
            nll=nll+dnrm(y(i),delta(i)+beta(1)+beta(2)*x(i),
     &           (1+heter(1)*x(i))*sigma(1),1)+ log(p)
         else
            nll=nll+dnrm(y(i),delta(i)-beta(1)-beta(2)*x(i),
     &           (1+heter(2)*x(i))*sigma(2),1)+ log(1 - p)
         end if
      end do

      nll = -nll
      
      return 
      end


CCCCCCCCCCCCCCCCCCCC
C     PARTIAL ALL
CCCCCCCCCCCCCCCCCCCC
      SUBROUTINE Partialf(param, tau, x, y, S, n, pp)

      implicit none

      integer n
      real*8 param(9)
      integer S(n)
      real*8 tau, x(n), y(n), pp(9)

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
            call SolveDeltaf(param1, tau, x(j), delta1(j))
            call SolveDeltaf(param2, tau, x(j), delta2(j))
         end do
         
         call NegLogLikelihoodf(y, S, x, delta1, param1, n, ll1)
         call NegLogLikelihoodf(y, S, x, delta2, param2, n, ll2)
         
         pp(i) = (ll1 - ll2)/2/epsilon

      end do
      
      return 
      end


CCCCCCCCCCCCCCCCCCCC
C      MAIN FUNCTION
CCCCCCCCCCCCCCCCCCCC

      SUBROUTINE QRGradientf(y, S, x, tau, n, niter, param, paramsave)

      implicit none

      integer n, niter
      integer S(n)
      real*8 param(9), y(n), x(n), tau, delta(n), paramsave(niter, 10)

C     TEMP
      integer i, iter, j
      real*8 dif, nll0, nll, pp(9), ppp(9)
      real*8 alpha(9), alphamax, alphamin, etap, etam
      real*8 maxx

C     INITIAL
      param(1) = 0
      param(2) = 0
      param(3) = 0
      param(4) = 0
      param(5) = 1
      param(6) = 1
      param(7) = 0.5
      param(8) = 0
      param(9) = 0
      
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

      maxx = 0
      do i = 1, n
         if (x(i) > maxx) then
            maxx = x(i)
         end if
      end do

      do while (dif > 0.00001 .and. iter .le. niter)
         call Partialf(param, tau, x, y, S, n, pp)
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
         param(5) = max(param(5), 0.01)
         param(6) = max(param(6), 0.01)
         param(7) = max(min(param(7), 0.99), 0.001)
         param(8) = max(param(8), -1/maxx)
         param(9) = max(param(9), -1/maxx)

         do i = 1, n
            call SolveDeltaf(param, tau, x(i), delta(i))        
         end do
         call NegLogLikelihoodf(y, S, x, delta, param, n, nll)
         dif = abs(nll - nll0)
         nll0 = nll
         do i = 1, 7
            paramsave(iter, i) = param(i)
         end do
         paramsave(iter, 8) = nll
         paramsave(iter, 9) = param(8)
         paramsave(iter, 10) = param(9)

         call progress(iter, niter)
         iter = iter + 1
      end do

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
