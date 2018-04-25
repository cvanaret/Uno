C     hristen this file wdotd.f

      subroutine wdotd (n, d, ws, lws, v)

c     ==========================================================
c     Computes W.d where W is Hessian and d is a vector for AMPL
c     Assumes v=0 on entry (OK, if called from gdotx, see QPsolve*.f)
c     ==========================================================

      implicit none

c     ... declaration of passed parameters
      integer          n, lws(0:*)
      double precision d(n), v(n), ws(*)

c     ... declaration of internal variables
      integer i, j, pjp, row

      double precision alpha
      common /kktalphac/ alpha

c     ========================  procedure body  =========================

c     ... form v = W.d from sparse, upper triangular Hessian
      pjp = lws(0)
      do i=1,n
         do j=lws(pjp+i-1),lws(pjp+i)-1
            row  = lws(j)
            v(i) = v(i) + ws(j)*d(row)
            if (row.ne.i) then
               v(row) = v(row) + ws(j)*d(i)
            else
               v(i)   = v(i) + alpha*d(row)
            endif
         enddo
      enddo

c     print * '  d = ',(d(i),i=1,n)
c     print * 'W.d = ',(v(i),i=1,n)

      return
      end

c     ******************************************************************

      subroutine gdotx (n, x, ws, lws, v)

      implicit none

c     ... declaration of passed parameters
      integer n, lws(*)
      double precision    x(n), v(n), ws(*)

c     ... declaration of internal variables
      integer i

c     ... storage map for hessian and scale_mode
      integer         scale_mode, phe
      common /scalec/ scale_mode, phe

c     ========================  procedure body  =========================

c     ... set v = 0
      do i=1,n
         v(i) = 0.D0
      enddo

c     ... allow for scaling of variables 
      if ((scale_mode.eq.1).or.(scale_mode.eq.3)) then
         do i=1,n
            x(i) = x(i) * ws(i)
         enddo
      endif

c     ... form v = W.d from sparse, upper triangular Hessian
      call Wdotd (n, x, ws(phe+1), lws, v)

c     ... allow for scaling of variables 
      if ((scale_mode.eq.1).or.(scale_mode.eq.3)) then
         do i=1,n
            v(i) = v(i) * ws(i)
            x(i) = x(i) / ws(i)
         enddo
      endif

      return
      end
c     ******************************************************************
      subroutine saipy2(s,a,la,i,y,n)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),y(*)
c     ========================  procedure body  =========================
c     saxpy with column i of A
      if(s.eq.0.D0)return
      jp=la(0)+i
      do j=la(jp),la(jp+1)-1
         ir=la(j)
         y(ir)=y(ir)+s*a(j)
      enddo
      return
      end

c     **************************** E N D *********************************
      function daiscpr2(n,a,la,i,x,b)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),x(*)
      DOUBLE PRECISION daiscpr2
      daiscpr2=dble(b)
      jp=la(0)+i
      do j=la(jp),la(jp+1)-1
         ir=la(j)
         daiscpr2=daiscpr2+dble(x(ir))*dble(a(j))
      enddo
      return
      end
