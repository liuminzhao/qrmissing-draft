c===========================================================
c$$$  
C$$$  Time-stamp: <liuminzhao 04/21/2013 15:31:47>
c$$$  Bivariate MLE using sigma
c$$$  exp(a0 + a1*x) as sigma
c===========================================================


CCCCCCCCCCCCCCCCCCCC
C        TARGET DELTA EQUATION 1
CCCCCCCCCCCCCCCCCCCC

      SUBROUTINE TargetEqn1H2f(delta,gamma,beta,sigma,tau,p,x,f,xdim)
      
      implicit none
      integer xdim, i
      real*8 delta, gamma(xdim), beta(xdim), sigma(2), tau, p, x(xdim)
      real*8 ans, f

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
      f=tau-p*pnrm(quan-delta - lp,0.d0,sigma(1),1,0)-
     &     (1-p)*pnrm(quan-delta+lp,0.d0,sigma(2),1,0)

      return 
      end 

CCCCCCCCCCCCCCCCCCCC
C     SOLVE DELTA 1
CCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SolveDelta1H2f(param, tau, x, root,xdim)

      implicit none
      integer xdim, i
      real*8 param(8*xdim + 3)
      real*8 gamma(xdim), beta(xdim), sigma(2), tau, p, x(xdim), root

C     TEMP

      real*8 a, b, fa, fb, m, fm, tol

C     INITIAL 
      
      a = -100 
      b = 100 
      fa = 0
      fb = 0
      m = (a + b)/2
      fm = 0
      tol = 0.00001
      
      sigma(1) = 0
      sigma(2) = 0

      do i = 1, xdim
         gamma(i) = param(i)
         beta(i) = param(xdim + i)
         sigma(1) = sigma(1) +  param(2*xdim + i) * x(i)
         sigma(2) = sigma(2) +  param(3*xdim + i) * x(i)
      end do
      p = param(8*xdim + 3)

      sigma(1) = exp(sigma(1))
      sigma(2) = exp(sigma(2))


      call TargetEqn1H2f(m,gamma,beta,sigma,tau, p, x, fm,xdim)
      call TargetEqn1H2f(b,gamma,beta,sigma,tau, p, x, fb,xdim)
      call TargetEqn1H2f(a,gamma,beta,sigma,tau, p, x, fa,xdim)

      do while(abs(fm) > tol)
         if (fm * fb < 0) then
            a = m
         else
            b = m
            fb = fm
         end if
         m = (a + b)/2
         call TargetEqn1H2f(m,gamma,beta,sigma,tau, p, x, fm,xdim)
      end do

      root = m
      return 
      end 


CCCCCCCCCCCCCCCCCCCC
C TARGET EQN FOR DELTA2     
CCCCCCCCCCCCCCCCCCCC

      subroutine targeteqn2h2f(delta2, gamma2, beta1, beta2, h, sigma1,
     &     sigma2, tau, p, x, delta1, f, xdim)
      implicit none
      integer xdim, i
      real*8 delta2, gamma2(xdim), beta1(xdim), beta2(xdim+1), h
      real*8 sigma1(2), sigma2(2), tau, p, x(xdim), delta1
      real*8 f

      real*8 pnrm
      real*8 p1, p2
      real*8 beta22

      real*8 beta1lp, beta2lp, gamma2lp

      beta22 = beta2(xdim+1) + h

      beta1lp = 0
      beta2lp = 0
      gamma2lp = 0

      do i  = 1, xdim
         beta1lp = beta1lp + beta1(i)*x(i)
         beta2lp = beta2lp + beta2(i)*x(i)
         gamma2lp = gamma2lp + gamma2(i)*x(i)
      end do

      if (beta22 .ne. 0) then
         p1=pnrm(((-delta2+ gamma2lp+beta2lp
     &        )/beta22-(delta1+beta1lp ))/
     &        sigma1(1)/sqrt(
     &        sigma2(1)**2/
     &        sigma1(1)**2/beta22**2+1),
     &        0.d0,1.d0,1,0)
      else 
         p1=pnrm((gamma2lp-delta2+ beta2lp 
     &        )/sigma2(1), 0.d0, 1.d0, 1, 0)
      end if
      
      if (beta2(3) .ne. 0) then
         p2=pnrm(((-delta2+gamma2lp - beta2lp
     &        )/beta2(3)-(delta1-beta1lp ))/
     &        sigma1(2)/sqrt(
     &        sigma2(2)**2/
     &        sigma1(2)**2/beta2(3)**2+1),
     &        0.d0,1.d0,1,0)
      else
         p2=pnrm((gamma2lp -delta2-beta2lp
     &        )/sigma2(2), 0.d0, 1.d0, 1, 0)
      end if 

      if (beta22 < 0) then
         p1 = 1 - p1
      end if

      if (beta2(3) < 0) then
         p2 = 1 - p2
      end if

      f = tau - p * p1 - (1 - p) * p2

      return 
      end


CCCCCCCCCCCCCCCCCCCC
C SOLVE DELTA2     
CCCCCCCCCCCCCCCCCCCC

      subroutine solvedelta2h2f(param, tau, x, delta1, root,xdim)
      implicit none
      integer xdim
      real*8 param(8*xdim + 3)
      real*8 gamma1(xdim), beta1(xdim), sigma1(2)
      real*8 gamma2(xdim), beta2(xdim + 1), sigma2(2), p, h
      real*8 tau, x(xdim), delta1, root
      real*8 a, b, fa, fb, m, fm, tol
      integer i, j
C     INITIAL 
      
      a = -100 
      b = 100 
      fa = 0
      fb = 0
      m = (a + b)/2
      fm = 0
      tol = 0.00001      

      sigma1(1) = 0
      sigma1(2) = 0
      sigma2(1) = 0
      sigma2(2) = 0

      do i = 1, xdim
         gamma1(i) = param(i)
         beta1(i) = param(xdim + i)
         sigma1(1) = sigma1(1) + param(2*xdim + i)*x(i)
         sigma1(2) = sigma1(2) + param(3*xdim + i)*x(i)
         gamma2(i) = param(4*xdim + i)
         beta2(i) = param(5*xdim + i)
         sigma2(1) = sigma2(1) + param(6*xdim + i)*x(i)
        sigma2(2) = sigma2(2) + (param(6*xdim + i)+param(7*xdim+i))*x(i)
      end do
      beta2(xdim + 1) = param(8*xdim + 1)
      h = param(xdim*8 + 2)
      p = param(xdim*8 + 3)

      sigma1(1) = exp(sigma1(1))
      sigma1(2) = exp(sigma1(2))
      sigma2(1) = exp(sigma2(1))
      sigma2(2) = exp(sigma2(2))

      call targeteqn2h2f(m, gamma2, beta1, beta2,h, sigma1,
     &     sigma2, tau, p, x, delta1, fm,xdim)
      call targeteqn2h2f(b, gamma2, beta1, beta2,h, sigma1,
     &     sigma2, tau, p, x, delta1, fb,xdim)
      call targeteqn2h2f(a, gamma2, beta1, beta2,h, sigma1,
     &     sigma2, tau, p, x, delta1, fa,xdim)

      do while(abs(fm) > tol)
         if (fm * fb < 0) then
            a = m
         else
            b = m
            fb = fm
         end if
         m = (a + b)/2
         call targeteqn2h2f(m, gamma2, beta1, beta2,h, sigma1,
     &        sigma2, tau, p, x, delta1, fm,xdim)
      end do

      root = m
      return 
      end 


CCCCCCCCCCCCCCCCCCCC
C NEGLOGLIKELIHOOD     
CCCCCCCCCCCCCCCCCCCC

      subroutine negloglikelihoodh2(y, R, x, delta1, delta2,
     &     param,n,nll,xdim)
      implicit none
      integer n, xdim
      integer R(n)
      real*8 y(n, 2),x(n,xdim), delta1(n), delta2(n), param(8*xdim+3)
      real*8 dnrm, nll
      real*8 gamma1(xdim), beta1(xdim), sigma1(2)
      real*8 gamma2(xdim), beta2(xdim + 1), sigma2(2), p
      integer  i,j,k
      real*8 beta22, h

      real*8 beta1lp, beta2lp
      real*8 alpha11(xdim), alpha01(xdim)
      real*8 alpha12(xdim)

      do i = 1, xdim
         gamma1(i) = param(i)
         beta1(i) = param(xdim+i)
         gamma2(i) = param(4*xdim + i)
         beta2(i) = param(5*xdim + i)
         alpha11(i) = param(2*xdim + i)
         alpha01(i) = param(3*xdim + i)
         alpha12(i) = param(6*xdim + i)         
      end do
      beta2(xdim + 1) = param(8*xdim + 1)
      h = param(8*xdim + 2)
      p = param(8* xdim + 3)

      beta22 = beta2(xdim + 1) + h

      nll = 0
      do i = 1, n
         beta1lp = 0
         beta2lp = 0
         sigma1(1) = 0
         sigma1(2) = 0
         sigma2(1) = 0
         do j = 1, xdim
            beta1lp = beta1lp + beta1(j)*x(i,j)
            beta2lp = beta2lp + beta2(j)*x(i,j)
            sigma1(1) = sigma1(1) + alpha11(j)*x(i,j)
            sigma1(2) = sigma1(2) + alpha01(j)*x(i,j)
            sigma2(1) = sigma2(1) + alpha12(j)*x(i,j)
         end do
         if (R(i) .eq. 1) then
            nll=nll+dnrm(y(i,1),delta1(i)+beta1lp
     &           ,exp(sigma1(1)),1)
     &           +dnrm(y(i,2),delta2(i)-beta2lp + 
     &           y(i,1)*beta22,
     &           exp(sigma2(1)),1)+log(p)
         else
            nll=nll+dnrm(y(i,1),delta1(i)-beta1lp
     &           ,exp(sigma1(2)),1) 
     &           + log(1 - p)
         end if
      end do
      nll = -nll
      return 
      end


CCCCCCCCCCCCCCCCCCCC
C       PARTIAL ALL
CCCCCCCCCCCCCCCCCCCC

      subroutine PartialH2f(param, tau, x, y, R, n, pp,xdim)
      implicit none
      integer n, xdim
      real*8 param(8*xdim+3)
      integer R(n)
      real*8 tau, x(n,xdim), y(n, 2), pp(8*xdim + 3)

      real*8 epsilon, delta1p(n), delta2p(n), param1(8*xdim+3)
      real*8 delta1m(n), delta2m(n), param2(8*xdim+3)
      real*8 ll1, ll2
      integer i, j, k

      real*8 tmpx(xdim)

      epsilon = 0.0003

      do i = 1, 8*xdim + 3
         if (.not. ((i>5*xdim.and.i.le.6*xdim).or.
     &        (i>7*xdim.and.i.le.8*xdim).or.i.eq.8*xdim+2)) 
     &        then
            do j = 1, 8*xdim+3
               param1(j) = param(j)
               param2(j) = param(j)
            end do
            param1(i) = param(i) + epsilon
            param2(i) = param(i) - epsilon

            do j = 1,n 
               do k = 1, xdim
                  tmpx(k) = x(j, k)
               end do
               call SolveDelta1H2f(param1, tau, tmpx, delta1p(j), xdim)
               call SolveDelta1H2f(param2, tau, tmpx, delta1m(j),xdim)
               call solvedelta2h2f(param1,tau,tmpx,delta1p(j),delta2p(j)
     &              ,xdim)
               call solvedelta2h2f(param2,tau,tmpx,delta1m(j),delta2m(j)
     &              ,xdim)
            end do

            call negloglikelihoodh2(y, R, x, delta1p, delta2p
     &           ,param1,n,ll1,xdim)
            call negloglikelihoodh2(y, R, x, delta1m, delta2m
     &           ,param2,n,ll2,xdim)

            pp(i) = (ll1 - ll2)/2/epsilon
         end if 
      end do
      return
      end

    

CCCCCCCCCCCCCCCCCCCC
C MAIN FUNCTION     
CCCCCCCCCCCCCCCCCCCC

      subroutine BiQRGradientH2f(y,R,x,tau,n,niter,param,paramsave,xdim)
      implicit none
      integer n, niter, xdim
      integer R(n)
      real*8 param(8*xdim+3), y(n,2), x(n,xdim),tau
      real*8 delta1(n), delta2(n),paramsave(niter,8*xdim+4)

      integer i, iter, j, k
      real*8 dif, nll0, nll,pp(8*xdim+3), alpha(8*xdim+3), ppp(8*xdim+3)
      real*8 alphamax, alphamin, etap, etam

      real*8 tmpx(xdim)

      do i = 1, 8*xdim+3
         pp(i) = 0
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
         do j = 1, 8*xdim + 4
            paramsave(i, j) = 0
         end do
      end do

      do while (dif > 0.00001 .and. iter .le. niter)
         call PartialH2f(param, tau, x, y, R, n, pp,xdim)
         do i = 1, 8*xdim + 3
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

         param(8*xdim+3) = max(min(param(8*xdim+3), 0.99), 0.01)

         do j = 1, n
            do k = 1, xdim
               tmpx(k) = x(j, k)
            end do
            call SolveDelta1H2f(param, tau, tmpx, delta1(j),xdim)
            call solvedelta2h2f(param,tau,tmpx,delta1(j),delta2(j),xdim)
         end do
         call negloglikelihoodh2(y,R,x,delta1, delta2,param,n,nll,xdim)
         dif = abs(nll - nll0)
         nll0 = nll
         do i = 1, 8*xdim+3
            paramsave(iter, i) = param(i)
         end do
         paramsave(iter, 8*xdim+4) = nll
C         print*, param, '\n' ,pp, '\n', iter

         call progress(iter, niter)
         iter = iter + 1
      end do
      return
      end

CCCCCCCCCCCCCCCCCCCC
C     PROGRESS BAR
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

