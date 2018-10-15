c                                                                                      
c  L-BFGS-B is released under the “New BSD License” (aka “Modified BSD License”        
c  or “3-clause license”)                                                              
c  Please read attached file License.txt                                               
c                                        
c                             DRIVER 2 in Fortran 77
c     --------------------------------------------------------------
c              CUSTOMIZED DRIVER FOR L-BFGS-B (version 3.0)
c     --------------------------------------------------------------
c
c        L-BFGS-B is a code for solving large nonlinear optimization
c             problems with simple bounds on the variables.
c
c        The code can also be used for unconstrained problems and is
c        as efficient for these problems as the earlier limited memory
c                          code L-BFGS.
c
c        This driver illustrates how to control the termination of the
c        run and how to design customized output.
c
c     References:
c
c        [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
c        memory algorithm for bound constrained optimization'',
c        SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
c
c        [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
c        Subroutines for Large Scale Bound Constrained Optimization''
c        Tech. Report, NAM-11, EECS Department, Northwestern University,
c        1994.
c
c
c          (Postscript files of these papers are available via anonymous
c           ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
c
c                              *  *  *
c
c         February 2011   (latest revision)
c         Optimization Center at Northwestern University
c         Instituto Tecnologico Autonomo de Mexico
c
c         Jorge Nocedal and Jose Luis Morales
c         Jorge Nocedal and Jose Luis Morales, Remark on "Algorithm 778: 
c         L-BFGS-B: Fortran Subroutines for Large-Scale Bound Constrained 
c         Optimization"  (2011). To appear in  ACM Transactions on 
c         Mathematical Software,
c
c     **************

      program driver
 
c     This driver shows how to replace the default stopping test
c       by other termination criteria. It also illustrates how to
c       print the values of several parameters during the course of
c       the iteration. The sample problem used here is the same as in 
c       DRIVER1 (the extended Rosenbrock function with bounds on the 
c       variables).
 
      integer          nmax, mmax
      parameter        (nmax=1024, mmax=17)
c        nmax is the dimension of the largest problem to be solved.
c        mmax is the maximum number of limited memory corrections.
 
c     Declare the variables needed by the code.
c       A description of all these variables is given at the end of 
c       driver1.
 
      character*60     task, csave
      logical          lsave(4)
      integer          n, m, iprint, 
     +                 nbd(nmax), iwa(3*nmax), isave(44)
      double precision f, factr, pgtol, 
     +                 x(nmax), l(nmax), u(nmax), g(nmax), dsave(29), 
     +                 wa(2*mmax*nmax+5*nmax+11*mmax*mmax+8*mmax)

c     Declare a few additional variables for the sample problem.

      double precision t1, t2
      integer          i
 
c     We suppress the default output.

      iprint = -1

c     We suppress both code-supplied stopping tests because the
c        user is providing his own stopping criteria.

      factr=0.0d0
      pgtol=0.0d0

c     We specify the dimension n of the sample problem and the number
c        m of limited memory corrections stored.  (n and m should not
c        exceed the limits nmax and mmax respectively.)
 
      n=25
      m=5
 
c     We now specify nbd which defines the bounds on the variables:
c                    l   specifies the lower bounds,
c                    u   specifies the upper bounds. 
 
c     First set bounds on the odd numbered variables.

      do 10 i=1,n,2
         nbd(i)=2
         l(i)=1.0d0
         u(i)=1.0d2
  10  continue

c     Next set bounds on the even numbered variables.

      do 12 i=2,n,2
         nbd(i)=2
         l(i)=-1.0d2
         u(i)=1.0d2
  12   continue

c     We now define the starting point.

      do 14 i=1,n
         x(i)=3.0d0
  14  continue
 
c     We now write the heading of the output.

      write (6,16)
  16  format(/,5x, 'Solving sample problem.',
     +       /,5x, ' (f = 0.0 at the optimal solution.)',/) 

c     We start the iteration by initializing task.
c 
      task = 'START'

c        ------- the beginning of the loop ----------
 
 111  continue
      
c     This is the call to the L-BFGS-B code.
 
      call setulb(n,m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint,
     +            csave,lsave,isave,dsave)
 
      if (task(1:2) .eq. 'FG') then
c        the minimization routine has returned to request the
c        function f and gradient g values at the current x.

c        Compute function value f for the sample problem.

         f=.25d0*(x(1)-1.d0)**2
         do 20 i=2,n
            f=f+(x(i)-x(i-1)**2)**2
  20     continue
         f=4.d0*f

c        Compute gradient g for the sample problem.

         t1=x(2)-x(1)**2
         g(1)=2.d0*(x(1)-1.d0)-1.6d1*x(1)*t1
         do 22 i=2,n-1
            t2=t1
            t1=x(i+1)-x(i)**2
            g(i)=8.d0*t2-1.6d1*x(i)*t1
  22     continue
         g(n)=8.d0*t1

c          go back to the minimization routine.
         goto 111
      endif
c
      if (task(1:5) .eq. 'NEW_X') then   
c     
c        the minimization routine has returned with a new iterate.
c        At this point have the opportunity of stopping the iteration 
c        or observing the values of certain parameters
c
c        First are two examples of stopping tests.

c        Note: task(1:4) must be assigned the value 'STOP' to terminate  
c          the iteration and ensure that the final results are
c          printed in the default format. The rest of the character
c          string TASK may be used to store other information.

c        1) Terminate if the total number of f and g evaluations
c             exceeds 99.

         if (isave(34) .ge. 99)
     +      task='STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIMIT'

c        2) Terminate if  |proj g|/(1+|f|) < 1.0d-10, where 
c           "proj g" denoted the projected gradient

         if (dsave(13) .le. 1.d-10*(1.0d0 + abs(f)))
     +      task='STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL'

c        We now wish to print the following information at each
c        iteration:
c        
c          1) the current iteration number, isave(30),
c          2) the total number of f and g evaluations, isave(34),
c          3) the value of the objective function f,
c          4) the norm of the projected gradient,  dsve(13)
c
c        See the comments at the end of driver1 for a description
c        of the variables isave and dsave.
         
         write (6,'(2(a,i5,4x),a,1p,d12.5,4x,a,1p,d12.5)') 'Iterate'
     +      ,isave(30),'nfg =',isave(34),'f =',f,'|proj g| =',dsave(13)

c        If the run is to be terminated, we print also the information
c        contained in task as well as the final value of x.

         if (task(1:4) .eq. 'STOP') then
            write (6,*) task  
            write (6,*) 'Final X='
            write (6,'((1x,1p, 6(1x,d11.4)))') (x(i),i = 1,n)
         endif

c          go back to the minimization routine.
         goto 111

      endif

c           ---------- the end of the loop -------------
 
c     If task is neither FG nor NEW_X we terminate execution.

      stop
 
      end

c======================= The end of driver2 ============================

