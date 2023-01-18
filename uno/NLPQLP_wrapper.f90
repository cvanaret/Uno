! This is pyOptSparse wrapper to NLPQLP. It currently only supports
! the serial version of NLPQLP. This file is used to catch the
! "reverse communication" call backs (that are not *actually* reverse
! mode call backs) and to call the *actual* user supplied
! function/gradient routines like all other *normal* optimizers. This
! is F90 wrapper and as such some of the work arrays etc are done as
! allocatable arrays here instead of in python.

subroutine wrapper(np, m, me, mmax, n, nmax, mnn2, x, f, g, df, dg, u, xl, xu, &
     c, d, acc, accqp, stpmin, maxfun, maxit, maxnm, rho, iprint, mode, iout, &
     ifile, ifail, lql, nlfunc, nlgrad)

  implicit none

  ! Input Variables
  integer :: np, m, me, mmax, n, nmax, mnn2
  double precision :: X(nmax, np), F(np), G(mmax, np), df(nmax), dg(mmax, nmax)
  double precision :: u(mnn2), xl(n), xu(n), c(nmax, nmax), d(nmax)
  double precision :: acc, accqp, stpmin
  integer :: maxit, maxnm, maxfun
  double precision :: rho
  character*(*) ifile
  integer :: mode, iprint, iout, ifail
  logical :: lql
  external nlfunc, nlgrad

  ! Local Variables
  logical, allocatable, dimension(:) :: active
  integer, allocatable :: kwa(:)
  double precision, allocatable, dimension(:) :: wa
  integer :: lwa, lkwa, lactiv
  logical :: fail
  external ql

  ! The structure of this routine is based on NLP_DEMOA.f90. It really
  ! pains me to use goto statements, but they kind of make sense here. 

  ifail  = 0
  lwa = 23*N+4*M+3*MMAX+NP*(N+M+1)+150 + 3*NMAX*NMAX/2+10*NMAX+MMAX+M+1 
  lkwa   = n + 30
  lactiv = 2*m + 10      

  ! Allocate some workspace
  allocate (wa(lwa), kwa(lkwa), active(lactiv))

  ! Open the output file
  open(unit=iout, file=trim(ifile), status='unknown')

  !  ------------ Call the user supplied function ----------------
1 continue

  call nlfunc(m, me, mmax, n, f, g, x, active, fail)
  call flush(iout)
  if (fail) then 
     ! NLPQLP says to set ifail to -10 and it will back-off on the
     ! step and try again
     ifail = -10
     print *, "+===========================================+"
     print *, "|  Failed user supplied function in NLPQLP  |" 
     print *, "+===========================================+"
  end if

  ! Now go back to NLPQLP IF this ISN'T the first pass when ifail == 0
  if (ifail == -1 .or. ifail == -10) goto 4
   
  !  ------------ Call the user supplied gradient ----------------
2 continue
  call nlgrad(m, me, mmax, n, f, g, df, dg, x, active, wa)
  call flush(iout)
  ! Now go back to NLPQLP
  if (ifail == -2) goto 4
    
  ! --------------------------------------------------------------
4 continue

  ! The actual NLPQLP run command. 
  call nlpqlp (np, m, me, mmax, n, nmax, mnn2, x, f, g, &
       df, dg, u, xl, xu, c, d, acc, accqp, stpmin, maxfun, maxit, &
       maxnm, rho, iprint, mode, iout, ifail, wa, lwa, kwa, lkwa, &
       active, lactiv, lql, ql)

   ! Now check to see what NLPQL has told us with the ifail call
   if (ifail == -1) goto 1
   if (ifail == -2) goto 2

   ! If we're here than NLPQL has finished. Clean up this routine and
   ! return to python.
   deallocate(wa, kwa, active)
   close(iout)

end subroutine wrapper
