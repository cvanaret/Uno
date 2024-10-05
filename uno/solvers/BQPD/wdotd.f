C Copyright (c) 2018-2024 Sven Leyffer
C Licensed under the MIT license. See LICENSE file in the project directory for details.

C     hristen this file wdotd.f

      subroutine wdotd (n, x, ws, lws, result)

c     ==========================================================
c     Computes result = W.x where W is Hessian and x is a vector for AMPL
c     Assumes v=0 on entry (OK, if called from gdotx, see QPsolve*.f)
c     ==========================================================

      implicit none

c     ... declaration of passed parameters
      integer          n, lws(0:*)
      double precision x(n), result(n), ws(*)

c     ... declaration of internal variables
      integer i, j, k, footer_start

c     inertia control for diagonal terms
      double precision alpha
      common /kktalphac/ alpha

c     ========================  procedure body  =========================

c     ... form result = W.x from sparse, upper triangular Hessian
      footer_start = lws(0)
      do i=1,n
         do k=lws(footer_start + i - 1), lws(footer_start + i)-1
            j = lws(k)
            result(i) = result(i) + ws(k)*x(j)
            if (j.ne.i) then
c               off-diagonal term
               result(j) = result(j) + ws(k)*x(i)
            else
c               diagonal term
               result(i) = result(i) + alpha*x(j)
            endif
         enddo
      enddo
      return
      end

c     ******************************************************************

      subroutine gdotx (n, x, ws, lws, result)

      implicit none

c     ... declaration of passed parameters
      integer n, lws(*)
      double precision    x(n), result(n), ws(*)

c     ... declaration of internal variables
      integer i

c     ... storage map for hessian and scale_mode
      integer         scale_mode, phe
      common /scalec/ scale_mode, phe

c     ========================  procedure body  =========================

c     ... set result = 0
      do i=1,n
         result(i) = 0.D0
      enddo

c     ... allow for scaling of variables 
      if ((scale_mode.eq.1).or.(scale_mode.eq.3)) then
         do i=1,n
            x(i) = x(i) * ws(i)
         enddo
      endif

c     ... form v = W.d from sparse, upper triangular Hessian
      call Wdotd (n, x, ws(phe+1), lws, result)

c     ... allow for scaling of variables 
      if ((scale_mode.eq.1).or.(scale_mode.eq.3)) then
         do i=1,n
            result(i) = result(i) * ws(i)
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
c     saxpy with column i of A: y + s*A_{i, :}
      if(s.eq.0.D0) return
      j_column_start = la(0) + i
      do j = la(j_column_start), la(j_column_start+1)-1
         i_variable = la(j)
         y(i_variable) = y(i_variable) + s*a(j)
      enddo
      return
      end

c     **************************** E N D *********************************
      function daiscpr2(n,a,la,i,x,b)
      implicit double precision (a-h,o-z)
      dimension a(*),la(0:*),x(*)
      DOUBLE PRECISION daiscpr2
c     dot product of x and row i of A
      daiscpr2 = dble(b)
      j_column_start = la(0) + i
      do j = la(j_column_start), la(j_column_start+1)-1
         i_variable = la(j)
         daiscpr2 = daiscpr2 + dble(x(i_variable))*dble(a(j))
      enddo
      return
      end
