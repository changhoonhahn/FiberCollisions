      implicit none
      real mu
      real LegendreP2 
      external LegendreP2
      ! read power multipoles: could be from theory or from data measurements
      mu=1.d0
      write(*,*) LegendreP2(mu)

      ! interpolation 
      end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      REAL function f_l_lprime(x,y,l,lp) ! calculate f_l,l' values 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      implicit none 
      integer l,lp
      real x,y

      !if (l.eq.0) then 
      !    if (lp.eq.0) then 



      !    endif 

      !elseif (l.eq.1) then 

      !elseif 

      !endif 

      f_l_lprime = (2.d0)
      end function f_l_lprime

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!legendre polynomials 
      function LegendreP2(mu)        
      implicit none
      real*8 mu,LegendreP2
      LegendreP2=(3._dp*mu**2 - 1._dp)/2._dp
      end function LegendreP2

      function LegendreP4(mu)        
      implicit none
      real*8 mu,LegendreP4
      LegendreP4=(35._dp*mu**4-30._dp*mu**2+3._dp)/8._dp
      return
      end

!window functions 

      function wnext(x) !window for (krc)^2 term
      implicit none
      real*8 x,bessj2
      bessj2= 2.d0*bessj1(x)/x - bessj0(x) !from J_nu's recursion relation
      wnext = bessj2/x**2 - 0.5d0*bessj1(x)/x
      return
      end 

      function w2d(x) !window for leading order term
      implicit none
      double x
      w2d = 2.d0*bessj1(x)/x
      end function w2d

!bessel functions 
      function bessj2(x)
      implicit none
      double x
      bessj2 = (2.d0/x**3-0.5d0/x)*bessj1(x)-bessj0(x)/x**2
      end function bessj2

      function bessj0(x)
      implicit none
      double :: bessj0,x,ax,xx,z
      double :: p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6, &
                s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,s5,s6
      DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4, &
      -.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-1, &
      .1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      DATA r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0, &
      651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/,s1,s2, &
      s3,s4,s5,s6/57568490411.d0,1029532985.d0,9494680.718d0, &
      59272.64853d0,267.8532712d0,1.d0/
      if(dabs(x).lt.8.d0)then
        y=x**2
        bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y* &
               (s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8.d0/ax
        y=z**2
        xx=ax-0.785398164d0
        bessj0=dsqrt(0.636619772d0/ax)*(dcos(xx)*(p1+y*(p2+y*(p3+y*(p4+y* &
               p5))))-z*dsin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      end function bessj0
      
      FUNCTION bessj1(x)
      implicit none
      double :: bessj1,x,ax,xx,z
      double ::  p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6, &
                 s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,s5,s6
      DATA r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0, &
      242396853.1d0,-2972611.439d0,15704.48260d0,-30.16036606d0/,s1,s2, & 
      s3,s4,s5,s6/144725228442.d0,2300535178.d0,18583304.74d0, & 
      99447.43394d0,376.9991397d0,1.d0/ &
      DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4, &
      .2457520174d-5,-.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0, &
      -.2002690873d-3,.8449199096d-5,-.88228987d-6,.105787412d-6/ &
      if(dabs(x).lt.8.d0)then
        y=x**2
        bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+ &
               y*(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8.d0/ax
        y=z**2
        xx=ax-2.356194491d0
        bessj1=dsqrt(0.636619772d0/ax)*(dcos(xx)*(p1+y*(p2+y*(p3+y*(p4+y* &
          p5))))-z*dsin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*dsign(1.d0,x)
      endif
      return
      end function bessj1
