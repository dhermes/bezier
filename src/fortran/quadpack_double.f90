subroutine dgtsl ( n, c, d, e, b, info )

!*****************************************************************************80
!
!! DGTSL solves a general tridiagonal linear system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 May 2005
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the tridiagonal matrix.
!
!    Input/output, real ( kind = 8 ) C(N), contains the subdiagonal of the
!    tridiagonal matrix in entries C(2:N).  On output, C is destroyed.
!
!    Input/output, real ( kind = 8 ) D(N).  On input, the diagonal of the
!    matrix.  On output, D is destroyed.
!
!    Input/output, real ( kind = 8 ) E(N), contains the superdiagonal of the
!    tridiagonal matrix in entries E(1:N-1).  On output E is destroyed.
!
!    Input/output, real ( kind = 8 ) B(N).  On input, the right hand side.
!    On output, the solution.
!
!    Output, integer ( kind = 4 ) INFO, error flag.
!    0, normal value.
!    K, the K-th element of the diagonal becomes exactly zero.  The
!    routine returns if this error condition is detected.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) b(n)
  real ( kind = 8 ) c(n)
  real ( kind = 8 ) d(n)
  real ( kind = 8 ) e(n)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) k
  real ( kind = 8 ) t

  info = 0
  c(1) = d(1)

  if ( 2 <= n ) then

    d(1) = e(1)
    e(1) = 0.0D+00
    e(n) = 0.0D+00

    do k = 1, n - 1
!
!  Find the larger of the two rows.
!
  if ( abs ( c(k) ) <= abs ( c(k+1) ) ) then
!
!  Interchange rows.
!
    t = c(k+1)
    c(k+1) = c(k)
    c(k) = t

    t = d(k+1)
    d(k+1) = d(k)
    d(k) = t

    t = e(k+1)
    e(k+1) = e(k)
    e(k) = t

    t = b(k+1)
    b(k+1) = b(k)
    b(k) = t

  end if
!
!  Zero elements.
!
  if ( c(k) == 0.0D+00 ) then
    info = k
    return
  end if

  t = -c(k+1) / c(k)
  c(k+1) = d(k+1) + t * d(k)
  d(k+1) = e(k+1) + t * e(k)
  e(k+1) = 0.0D+00
  b(k+1) = b(k+1) + t * b(k)

    end do

  end if

  if ( c(n) == 0.0D+00 ) then
    info = n
    return
  end if
!
!  Back solve.
!
  b(n) = b(n) / c(n)

  if ( 1 < n ) then

    b(n-1) = ( b(n-1) - d(n-1) * b(n) ) / c(n-1)

    do k = n-2, 1, -1
  b(k) = ( b(k) - d(k) * b(k+1) - e(k) * b(k+2) ) / c(k)
    end do

  end if

  return
end
subroutine dqage ( f, a, b, epsabs, epsrel, key, limit, result, abserr, &
  neval, ier, alist, blist, rlist, elist, iord, last )

!*****************************************************************************80
!
!! DQAGE estimates a definite integral.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  the routine calculates an approximation result to a given
!      definite integral   i = integral of f over (a,b),
!      hopefully satisfying following claim for accuracy
!      abs(i-reslt).le.max(epsabs,epsrel*abs(i)).
!
!  Parameters:
!
!   on entry
!      f      - real ( kind = 8 )
!               function subprogram defining the integrand
!               function f(x). the actual name for f needs to be
!               declared e x t e r n a l in the driver program.
!
!      a      - real ( kind = 8 )
!               lower limit of integration
!
!      b      - real ( kind = 8 )
!               upper limit of integration
!
!      epsabs - real ( kind = 8 )
!               absolute accuracy requested
!      epsrel - real ( kind = 8 )
!               relative accuracy requested
!               if  epsabs.le.0
!               and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!               the routine will end with ier = 6.
!
!      key    - integer ( kind = 4 )
!               key for choice of local integration rule
!               a gauss-kronrod pair is used with
!                    7 - 15 points if key.lt.2,
!                   10 - 21 points if key = 2,
!                   15 - 31 points if key = 3,
!                   20 - 41 points if key = 4,
!                   25 - 51 points if key = 5,
!                   30 - 61 points if key.gt.5.
!
!      limit  - integer ( kind = 4 )
!               gives an upperbound on the number of subintervals
!               in the partition of (a,b), limit.ge.1.
!
!   on return
!      result - real ( kind = 8 )
!               approximation to the integral
!
!      abserr - real ( kind = 8 )
!               estimate of the modulus of the absolute error,
!               which should equal or exceed abs(i-result)
!
!      neval  - integer ( kind = 4 )
!               number of integrand evaluations
!
!      ier    - integer ( kind = 4 )
!               ier = 0 normal and reliable termination of the
!                       routine. it is assumed that the requested
!                       accuracy has been achieved.
!               ier.gt.0 abnormal termination of the routine
!                       the estimates for result and error are
!                       less reliable. it is assumed that the
!                       requested accuracy has not been achieved.
!      error messages
!               ier = 1 maximum number of subdivisions allowed
!                       has been achieved. one can allow more
!                       subdivisions by increasing the value
!                       of limit.
!                       however, if this yields no improvement it
!                       is rather advised to analyze the integrand
!                       in order to determine the integration
!                       difficulties. if the position of a local
!                       difficulty can be determined(e.g.
!                       singularity, discontinuity within the
!                       interval) one will probably gain from
!                       splitting up the interval at this point
!                       and calling the integrator on the
!                       subranges. if possible, an appropriate
!                       special-purpose integrator should be used
!                       which is designed for handling the type of
!                       difficulty involved.
!                   = 2 the occurrence of roundoff error is
!                       detected, which prevents the requested
!                       tolerance from being achieved.
!                   = 3 extremely bad integrand behaviour occurs
!                       at some points of the integration
!                       interval.
!                   = 6 the input is invalid, because
!                       (epsabs.le.0 and
!                        epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!                       result, abserr, neval, last, rlist(1) ,
!                       elist(1) and iord(1) are set to zero.
!                       alist(1) and blist(1) are set to a and b
!                       respectively.
!
!      alist   - real ( kind = 8 )
!                vector of dimension at least limit, the first
!                 last  elements of which are the left
!                end points of the subintervals in the partition
!                of the given integration range (a,b)
!
!      blist   - real ( kind = 8 )
!                vector of dimension at least limit, the first
!                 last  elements of which are the right
!                end points of the subintervals in the partition
!                of the given integration range (a,b)
!
!      rlist   - real ( kind = 8 )
!                vector of dimension at least limit, the first
!                 last  elements of which are the
!                integral approximations on the subintervals
!
!      elist   - real ( kind = 8 )
!                vector of dimension at least limit, the first
!                 last  elements of which are the moduli of the
!                absolute error estimates on the subintervals
!
!      iord    - integer ( kind = 4 )
!                vector of dimension at least limit, the first k
!                elements of which are pointers to the
!                error estimates over the subintervals,
!                such that elist(iord(1)), ...,
!                elist(iord(k)) form a decreasing sequence,
!                with k = last if last.le.(limit/2+2), and
!                k = limit+1-last otherwise
!
!      last    - integer ( kind = 4 )
!                number of subintervals actually produced in the
!                subdivision process
!
!  Local Parameters:
!
!     alist     - list of left end points of all subintervals
!                 considered up to now
!     blist     - list of right end points of all subintervals
!                 considered up to now
!     rlist(i)  - approximation to the integral over
!                (alist(i),blist(i))
!     elist(i)  - error estimate applying to rlist(i)
!     maxerr    - pointer to the interval with largest
!                 error estimate
!     errmax    - elist(maxerr)
!     area      - sum of the integrals over the subintervals
!     errsum    - sum of the errors over the subintervals
!     errbnd    - requested accuracy max(epsabs,epsrel*
!                 abs(result))
!     *****1    - variable for the left subinterval
!     *****2    - variable for the right subinterval
!     last      - index for subdivision
!
!
!     machine dependent constants
!
!     epmach  is the largest relative spacing.
!     uflow  is the smallest positive magnitude.
!
  implicit none

  real ( kind = 8 ) a,abserr,alist,area,area1,area12,area2,a1,a2,b, &
    blist,b1,b2,defabs,defab1,defab2,elist,epmach, &
    epsabs,epsrel,errbnd,errmax,error1,error2,erro12,errsum,f, &
    resabs,result,rlist,uflow
  integer ( kind = 4 ) ier,iord,iroff1,iroff2,k,key,keyf,last,limit, &
    maxerr, nrmax, neval

  dimension alist(limit),blist(limit),elist(limit),iord(limit), &
    rlist(limit)

  external f

  epmach = epsilon ( epmach )
  uflow = tiny ( uflow )
!
!  test on validity of parameters
!
  ier = 0
  neval = 0
  last = 0
  result = 0.0D+00
  abserr = 0.0D+00
  alist(1) = a
  blist(1) = b
  rlist(1) = 0.0D+00
  elist(1) = 0.0D+00
  iord(1) = 0

  if(epsabs.le.0.0D+00.and. &
    epsrel.lt. max ( 0.5D+02*epmach,0.5d-28)) then
    ier = 6
    return
  end if
!
!  first approximation to the integral
!
  keyf = key
  if(key.le.0) keyf = 1
  if(key.ge.7) keyf = 6
  neval = 0
  if(keyf.eq.1) call dqk15(f,a,b,result,abserr,defabs,resabs)
  if(keyf.eq.2) call dqk21(f,a,b,result,abserr,defabs,resabs)
  if(keyf.eq.3) call dqk31(f,a,b,result,abserr,defabs,resabs)
  if(keyf.eq.4) call dqk41(f,a,b,result,abserr,defabs,resabs)
  if(keyf.eq.5) call dqk51(f,a,b,result,abserr,defabs,resabs)
  if(keyf.eq.6) call dqk61(f,a,b,result,abserr,defabs,resabs)
  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1
!
!  test on accuracy.
!
  errbnd =  max ( epsabs, epsrel* abs ( result ) )

  if(abserr.le.0.5D+02* epmach * defabs .and. &
    abserr.gt.errbnd) then
    ier = 2
  end if

  if(limit.eq.1) then
    ier = 1
  end if

  if ( ier .ne. 0 .or. &
    (abserr .le. errbnd .and. abserr .ne. resabs ) .or. &
    abserr .eq. 0.0D+00 ) then

    if(keyf.ne.1) then
      neval = (10*keyf+1)*(2*neval+1)
    else
      neval = 30*neval+15
    end if

    return

  end if
!
!  initialization
!
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  nrmax = 1
  iroff1 = 0
  iroff2 = 0
!
!  main do-loop
!
  do last = 2, limit
!
!  bisect the subinterval with the largest error estimate.
!
    a1 = alist(maxerr)
    b1 = 0.5D+00*(alist(maxerr)+blist(maxerr))
    a2 = b1
    b2 = blist(maxerr)

    if(keyf.eq.1) call dqk15(f,a1,b1,area1,error1,resabs,defab1)
    if(keyf.eq.2) call dqk21(f,a1,b1,area1,error1,resabs,defab1)
    if(keyf.eq.3) call dqk31(f,a1,b1,area1,error1,resabs,defab1)
    if(keyf.eq.4) call dqk41(f,a1,b1,area1,error1,resabs,defab1)
    if(keyf.eq.5) call dqk51(f,a1,b1,area1,error1,resabs,defab1)
    if(keyf.eq.6) call dqk61(f,a1,b1,area1,error1,resabs,defab1)

    if(keyf.eq.1) call dqk15(f,a2,b2,area2,error2,resabs,defab2)
    if(keyf.eq.2) call dqk21(f,a2,b2,area2,error2,resabs,defab2)
    if(keyf.eq.3) call dqk31(f,a2,b2,area2,error2,resabs,defab2)
    if(keyf.eq.4) call dqk41(f,a2,b2,area2,error2,resabs,defab2)
    if(keyf.eq.5) call dqk51(f,a2,b2,area2,error2,resabs,defab2)
    if(keyf.eq.6) call dqk61(f,a2,b2,area2,error2,resabs,defab2)
!
!  improve previous approximations to integral
!  and error and test for accuracy.
!
    neval = neval+1
    area12 = area1+area2
    erro12 = error1+error2
    errsum = errsum+erro12-errmax
    area = area+area12-rlist(maxerr)

    if ( defab1 .ne. error1 .and. defab2 .ne. error2 ) then

      if( abs ( rlist(maxerr)-area12).le.0.1D-04* abs ( area12) &
        .and. erro12.ge.0.99D+00*errmax) then
        iroff1 = iroff1+1
      end if

      if(last.gt.10.and.erro12.gt.errmax) then
        iroff2 = iroff2+1
      end if

    end if

    rlist(maxerr) = area1
    rlist(last) = area2
    errbnd =  max ( epsabs,epsrel* abs ( area))

    if ( errbnd .lt. errsum ) then
!
!  test for roundoff error and eventually set error flag.
!
      if(iroff1.ge.6.or.iroff2.ge.20) then
        ier = 2
      end if
!
!  set error flag in the case that the number of subintervals
!  equals limit.
!
      if(last.eq.limit) then
        ier = 1
      end if
!
!  set error flag in the case of bad integrand behaviour
!  at a point of the integration range.
!
      if( max (  abs ( a1), abs ( b2)).le.(0.1D+01+0.1D+03* &
        epmach)*( abs ( a2)+0.1D+04*uflow)) then
        ier = 3
      end if

    end if
!
!  append the newly-created intervals to the list.
!
    if(error2.le.error1) then
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
    else
      alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
    end if
!
!  call dqpsrt to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with the largest error estimate (to be bisected next).
!
    call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)

    if(ier.ne.0.or.errsum.le.errbnd) then
      exit
    end if

  end do
!
!  compute final result.
!
  result = 0.0D+00
  do k=1,last
    result = result+rlist(k)
  end do
  abserr = errsum

  if(keyf.ne.1) then
    neval = (10*keyf+1)*(2*neval+1)
  else
    neval = 30*neval+15
  end if

  return
end
subroutine dqag ( f, a, b, epsabs, epsrel, key, result, abserr, neval, ier, &
  limit, lenw, last, iwork, work )

!*****************************************************************************80
!
!! DQAG approximates an integral over a finite interval.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  the routine calculates an approximation result to a given
!      definite integral i = integral of f over (a,b),
!      hopefully satisfying following claim for accuracy
!      abs(i-result)le.max(epsabs,epsrel*abs(i)).
!
!  Parameters:
!
!      f      - real ( kind = 8 )
!               function subprogam defining the integrand
!               function f(x). the actual name for f needs to be
!               declared e x t e r n a l in the driver program.
!
!      a      - real ( kind = 8 )
!               lower limit of integration
!
!      b      - real ( kind = 8 )
!               upper limit of integration
!
!      epsabs - real ( kind = 8 )
!               absolute accoracy requested
!      epsrel - real ( kind = 8 )
!               relative accuracy requested
!               if  epsabs.le.0
!               and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!               the routine will end with ier = 6.
!
!      key    - integer ( kind = 4 )
!               key for choice of local integration rule
!               a gauss-kronrod pair is used with
!                 7 - 15 points if key.lt.2,
!                10 - 21 points if key = 2,
!                15 - 31 points if key = 3,
!                20 - 41 points if key = 4,
!                25 - 51 points if key = 5,
!                30 - 61 points if key.gt.5.
!
!   on return
!      result - real ( kind = 8 )
!               approximation to the integral
!
!      abserr - real ( kind = 8 )
!               estimate of the modulus of the absolute error,
!               which should equal or exceed abs(i-result)
!
!      neval  - integer ( kind = 4 )
!               number of integrand evaluations
!
!      ier    - integer ( kind = 4 )
!               ier = 0 normal and reliable termination of the
!                       routine. it is assumed that the requested
!                       accuracy has been achieved.
!               ier.gt.0 abnormal termination of the routine
!                       the estimates for result and error are
!                       less reliable. it is assumed that the
!                       requested accuracy has not been achieved.
!                error messages
!               ier = 1 maximum number of subdivisions allowed
!                       has been achieved. one can allow more
!                       subdivisions by increasing the value of
!                       limit (and taking the according dimension
!                       adjustments into account). however, if
!                       this yield no improvement it is advised
!                       to analyze the integrand in order to
!                       determine the integration difficulaties.
!                       if the position of a local difficulty can
!                       be determined (i.e.singularity,
!                       discontinuity within the interval) one
!                       will probably gain from splitting up the
!                       interval at this point and calling the
!                       integrator on the subranges. if possible,
!                       an appropriate special-purpose integrator
!                       should be used which is designed for
!                       handling the type of difficulty involved.
!                   = 2 the occurrence of roundoff error is
!                       detected, which prevents the requested
!                       tolerance from being achieved.
!                   = 3 extremely bad integrand behaviour occurs
!                       at some points of the integration
!                       interval.
!                   = 6 the input is invalid, because
!                       (epsabs.le.0 and
!                        epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
!                       or limit.lt.1 or lenw.lt.limit*4.
!                       result, abserr, neval, last are set
!                       to zero.
!                       except when lenw is invalid, iwork(1),
!                       work(limit*2+1) and work(limit*3+1) are
!                       set to zero, work(1) is set to a and
!                       work(limit+1) to b.
!
!   dimensioning parameters
!      limit - integer ( kind = 4 )
!              dimensioning parameter for iwork
!              limit determines the maximum number of subintervals
!              in the partition of the given integration interval
!              (a,b), limit.ge.1.
!              if limit.lt.1, the routine will end with ier = 6.
!
!      lenw  - integer ( kind = 4 )
!              dimensioning parameter for work
!              lenw must be at least limit*4.
!              if lenw.lt.limit*4, the routine will end with
!              ier = 6.
!
!      last  - integer ( kind = 4 )
!              on return, last equals the number of subintervals
!              produced in the subdivision process, which
!              determines the number of significant elements
!              actually in the work arrays.
!
!   work arrays
!      iwork - integer ( kind = 4 )
!              vector of dimension at least limit, the first k
!              elements of which contain pointers to the error
!              estimates over the subintervals, such that
!              work(limit*3+iwork(1)),... , work(limit*3+iwork(k))
!              form a decreasing sequence with k = last if
!              last.le.(limit/2+2), and k = limit+1-last otherwise
!
!      work  - real ( kind = 8 )
!              vector of dimension at least lenw
!              on return
!              work(1), ..., work(last) contain the left end
!              points of the subintervals in the partition of
!               (a,b),
!              work(limit+1), ..., work(limit+last) contain the
!               right end points,
!              work(limit*2+1), ..., work(limit*2+last) contain
!               the integral approximations over the subintervals,
!              work(limit*3+1), ..., work(limit*3+last) contain
!               the error estimates.
!
  implicit none

  integer ( kind = 4 ) lenw
  integer ( kind = 4 ) limit

  real ( kind = 8 ) a
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ) epsabs
  real ( kind = 8 ) epsrel
  real ( kind = 8 ), external :: f
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) iwork(limit)
  integer ( kind = 4 ) key
  integer ( kind = 4 ) last
  integer ( kind = 4 ) lvl
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) l3
  integer ( kind = 4 ) neval
  real ( kind = 8 ) result
  real ( kind = 8 ) work(lenw)
!
!  check validity of lenw.
!
  ier = 6
  neval = 0
  last = 0
  result = 0.0D+00
  abserr = 0.0D+00
  if(limit.lt.1.or.lenw.lt.limit*4) go to 10
!
!  prepare call for dqage.
!
  l1 = limit+1
  l2 = limit+l1
  l3 = limit+l2

  call dqage(f,a,b,epsabs,epsrel,key,limit,result,abserr,neval, &
    ier,work(1),work(l1),work(l2),work(l3),iwork,last)
!
!  call error handler if necessary.
!
  lvl = 0
10    continue

  if(ier.eq.6) lvl = 1
  if(ier.ne.0) call xerror('abnormal return from dqag ',26,ier,lvl)

  return
end
subroutine dqagie ( f, bound, inf, epsabs, epsrel, limit, result, abserr, &
  neval, ier, alist, blist, rlist, elist, iord, last )

!*****************************************************************************80
!
!! DQAGIE estimates an integral over a semi-infinite or infinite interval.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  the routine calculates an approximation result to a given
!      integral   i = integral of f over (bound,+infinity)
!      or i = integral of f over (-infinity,bound)
!      or i = integral of f over (-infinity,+infinity),
!      hopefully satisfying following claim for accuracy
!      abs(i-result).le.max(epsabs,epsrel*abs(i))
!
!  Parameters:
!
!      f      - real ( kind = 8 )
!               function subprogram defining the integrand
!               function f(x). the actual name for f needs to be
!               declared e x t e r n a l in the driver program.
!
!      bound  - real ( kind = 8 )
!               finite bound of integration range
!               (has no meaning if interval is doubly-infinite)
!
!      inf    - real ( kind = 8 )
!               indicating the kind of integration range involved
!               inf = 1 corresponds to  (bound,+infinity),
!               inf = -1            to  (-infinity,bound),
!               inf = 2             to (-infinity,+infinity).
!
!      epsabs - real ( kind = 8 )
!               absolute accuracy requested
!      epsrel - real ( kind = 8 )
!               relative accuracy requested
!               if  epsabs.le.0
!               and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!               the routine will end with ier = 6.
!
!      limit  - integer ( kind = 4 )
!               gives an upper bound on the number of subintervals
!               in the partition of (a,b), limit.ge.1
!
!   on return
!      result - real ( kind = 8 )
!               approximation to the integral
!
!      abserr - real ( kind = 8 )
!               estimate of the modulus of the absolute error,
!               which should equal or exceed abs(i-result)
!
!      neval  - integer ( kind = 4 )
!               number of integrand evaluations
!
!      ier    - integer ( kind = 4 )
!               ier = 0 normal and reliable termination of the
!                       routine. it is assumed that the requested
!                       accuracy has been achieved.
!             - ier.gt.0 abnormal termination of the routine. the
!                       estimates for result and error are less
!                       reliable. it is assumed that the requested
!                       accuracy has not been achieved.
!      error messages
!               ier = 1 maximum number of subdivisions allowed
!                       has been achieved. one can allow more
!                       subdivisions by increasing the value of
!                       limit (and taking the according dimension
!                       adjustments into account). however,if
!                       this yields no improvement it is advised
!                       to analyze the integrand in order to
!                       determine the integration difficulties.
!                       if the position of a local difficulty can
!                       be determined (e.g. singularity,
!                       discontinuity within the interval) one
!                       will probably gain from splitting up the
!                       interval at this point and calling the
!                       integrator on the subranges. if possible,
!                       an appropriate special-purpose integrator
!                       should be used, which is designed for
!                       handling the type of difficulty involved.
!                   = 2 the occurrence of roundoff error is
!                       detected, which prevents the requested
!                       tolerance from being achieved.
!                       the error may be under-estimated.
!                   = 3 extremely bad integrand behaviour occurs
!                       at some points of the integration
!                       interval.
!                   = 4 the algorithm does not converge.
!                       roundoff error is detected in the
!                       extrapolation table.
!                       it is assumed that the requested tolerance
!                       cannot be achieved, and that the returned
!                       result is the best which can be obtained.
!                   = 5 the integral is probably divergent, or
!                       slowly convergent. it must be noted that
!                       divergence can occur with any other value
!                       of ier.
!                   = 6 the input is invalid, because
!                       (epsabs.le.0 and
!                        epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!                       result, abserr, neval, last, rlist(1),
!                       elist(1) and iord(1) are set to zero.
!                       alist(1) and blist(1) are set to 0
!                       and 1 respectively.
!
!      alist  - real ( kind = 8 )
!               vector of dimension at least limit, the first
!                last  elements of which are the left
!               end points of the subintervals in the partition
!               of the transformed integration range (0,1).
!
!      blist  - real ( kind = 8 )
!               vector of dimension at least limit, the first
!                last  elements of which are the right
!               end points of the subintervals in the partition
!               of the transformed integration range (0,1).
!
!      rlist  - real ( kind = 8 )
!               vector of dimension at least limit, the first
!                last  elements of which are the integral
!               approximations on the subintervals
!
!      elist  - real ( kind = 8 )
!               vector of dimension at least limit,  the first
!               last elements of which are the moduli of the
!               absolute error estimates on the subintervals
!
!      iord   - integer ( kind = 4 )
!               vector of dimension limit, the first k
!               elements of which are pointers to the
!               error estimates over the subintervals,
!               such that elist(iord(1)), ..., elist(iord(k))
!               form a decreasing sequence, with k = last
!               if last.le.(limit/2+2), and k = limit+1-last
!               otherwise
!
!      last   - integer ( kind = 4 )
!               number of subintervals actually produced
!               in the subdivision process
!
!  Local Parameters:
!
!      the dimension of rlist2 is determined by the value of
!      limexp in routine dqelg.
!
!     alist     - list of left end points of all subintervals
!                 considered up to now
!     blist     - list of right end points of all subintervals
!                 considered up to now
!     rlist(i)  - approximation to the integral over
!                 (alist(i),blist(i))
!     rlist2    - array of dimension at least (limexp+2),
!                 containing the part of the epsilon table
!                 wich is still needed for further computations
!     elist(i)  - error estimate applying to rlist(i)
!     maxerr    - pointer to the interval with largest error
!                 estimate
!     errmax    - elist(maxerr)
!     erlast    - error on the interval currently subdivided
!                 (before that subdivision has taken place)
!     area      - sum of the integrals over the subintervals
!     errsum    - sum of the errors over the subintervals
!     errbnd    - requested accuracy max(epsabs,epsrel*
!                 abs(result))
!     *****1    - variable for the left subinterval
!     *****2    - variable for the right subinterval
!     last      - index for subdivision
!     nres      - number of calls to the extrapolation routine
!     numrl2    - number of elements currently in rlist2. if an
!                 appropriate approximation to the compounded
!                 integral has been obtained, it is put in
!                 rlist2(numrl2) after numrl2 has been increased
!                 by one.
!     small     - length of the smallest interval considered up
!                 to now, multiplied by 1.5
!     erlarg    - sum of the errors over the intervals larger
!                 than the smallest interval considered up to now
!     extrap    - logical variable denoting that the routine
!                 is attempting to perform extrapolation. i.e.
!                 before subdividing the smallest interval we
!                 try to decrease the value of erlarg.
!     noext     - logical variable denoting that extrapolation
!                 is no longer allowed (true-value)
!
!      machine dependent constants
!
!     epmach is the largest relative spacing.
!     uflow is the smallest positive magnitude.
!     oflow is the largest positive magnitude.
!
  implicit none

  real ( kind = 8 ) abseps,abserr,alist,area,area1,area12,area2,a1, &
    a2,blist,boun,bound,b1,b2,correc,defabs,defab1,defab2, &
    dres,elist,epmach,epsabs,epsrel,erlarg,erlast, &
    errbnd,errmax,error1,error2,erro12,errsum,ertest,f,oflow,resabs, &
    reseps,result,res3la,rlist,rlist2,small,uflow
  integer ( kind = 4 ) id,ier,ierro,inf,iord,iroff1,iroff2, &
    iroff3,jupbnd,k,ksgn, &
    ktmin,last,limit,maxerr,neval,nres,nrmax,numrl2
  logical extrap,noext
  dimension alist(limit),blist(limit),elist(limit),iord(limit), &
    res3la(3),rlist(limit),rlist2(52)

  external f

  epmach = epsilon ( epmach )
!
!  test on validity of parameters
!
  ier = 0
  neval = 0
  last = 0
  result = 0.0D+00
  abserr = 0.0D+00
  alist(1) = 0.0D+00
  blist(1) = 0.1D+01
  rlist(1) = 0.0D+00
  elist(1) = 0.0D+00
  iord(1) = 0

  if(epsabs.le.0.0D+00.and.epsrel.lt. max ( 0.5D+02*epmach,0.5D-28)) then
    ier = 6
  end if

  if(ier.eq.6) then
    return
  end if
!
!  first approximation to the integral
!
!  determine the interval to be mapped onto (0,1).
!  if inf = 2 the integral is computed as i = i1+i2, where
!  i1 = integral of f over (-infinity,0),
!  i2 = integral of f over (0,+infinity).
!
  boun = bound
  if(inf.eq.2) boun = 0.0D+00
  call dqk15i(f,boun,inf,0.0D+00,0.1D+01,result,abserr, &
    defabs,resabs)
!
!  test on accuracy
!
  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1
  dres =  abs ( result)
  errbnd =  max ( epsabs,epsrel*dres)
  if(abserr.le.1.0D+02*epmach*defabs.and.abserr.gt.errbnd) ier = 2
  if(limit.eq.1) ier = 1
  if(ier.ne.0.or.(abserr.le.errbnd.and.abserr.ne.resabs).or. &
    abserr.eq.0.0D+00) go to 130
!
!  initialization
!
  uflow = tiny ( uflow )
  oflow = huge ( oflow )
  rlist2(1) = result
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  abserr = oflow
  nrmax = 1
  nres = 0
  ktmin = 0
  numrl2 = 2
  extrap = .false.
  noext = .false.
  ierro = 0
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0
  ksgn = -1
  if(dres.ge.(0.1D+01-0.5D+02*epmach)*defabs) ksgn = 1
!
!  main do-loop
!
  do 90 last = 2,limit
!
!  bisect the subinterval with nrmax-th largest error estimate.
!
    a1 = alist(maxerr)
    b1 = 0.5D+00*(alist(maxerr)+blist(maxerr))
    a2 = b1
    b2 = blist(maxerr)
    erlast = errmax
    call dqk15i(f,boun,inf,a1,b1,area1,error1,resabs,defab1)
    call dqk15i(f,boun,inf,a2,b2,area2,error2,resabs,defab2)
!
!  improve previous approximations to integral
!  and error and test for accuracy.
!
    area12 = area1+area2
    erro12 = error1+error2
    errsum = errsum+erro12-errmax
    area = area+area12-rlist(maxerr)
    if(defab1.eq.error1.or.defab2.eq.error2)go to 15
    if( abs ( rlist(maxerr)-area12).gt.0.1D-04* abs ( area12) &
    .or.erro12.lt.0.99D+00*errmax) go to 10
    if(extrap) iroff2 = iroff2+1
    if(.not.extrap) iroff1 = iroff1+1
   10   if(last.gt.10.and.erro12.gt.errmax) iroff3 = iroff3+1
   15   rlist(maxerr) = area1
    rlist(last) = area2
    errbnd =  max ( epsabs,epsrel* abs ( area))
!
!  test for roundoff error and eventually set error flag.
!
    if(iroff1+iroff2.ge.10.or.iroff3.ge.20) ier = 2
    if(iroff2.ge.5) ierro = 3
!
!  set error flag in the case that the number of
!  subintervals equals limit.
!
    if(last.eq.limit) ier = 1
!
!  set error flag in the case of bad integrand behaviour
!  at some points of the integration range.
!
    if( max (  abs ( a1), abs ( b2)).le.(0.1D+01+0.1D+03*epmach)* &
    ( abs ( a2)+0.1D+04*uflow)) then
      ier = 4
    end if
!
!  append the newly-created intervals to the list.
!
    if(error2.gt.error1) go to 20
    alist(last) = a2
    blist(maxerr) = b1
    blist(last) = b2
    elist(maxerr) = error1
    elist(last) = error2
    go to 30
   20   continue

    alist(maxerr) = a2
    alist(last) = a1
    blist(last) = b1
    rlist(maxerr) = area2
    rlist(last) = area1
    elist(maxerr) = error2
    elist(last) = error1
!
!  call dqpsrt to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with nrmax-th largest error estimate (to be bisected next).
!
   30   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
    if(errsum.le.errbnd) go to 115
    if(ier.ne.0) go to 100
    if(last.eq.2) go to 80
    if(noext) go to 90
    erlarg = erlarg-erlast
    if( abs ( b1-a1).gt.small) erlarg = erlarg+erro12
    if(extrap) go to 40
!
!  test whether the interval to be bisected next is the
!  smallest interval.
!
    if( abs ( blist(maxerr)-alist(maxerr)).gt.small) go to 90
    extrap = .true.
    nrmax = 2
   40   if(ierro.eq.3.or.erlarg.le.ertest) go to 60
!
!  the smallest interval has the largest error.
!  before bisecting decrease the sum of the errors over the
!  larger intervals (erlarg) and perform extrapolation.
!
    id = nrmax
    jupbnd = last
    if(last.gt.(2+limit/2)) jupbnd = limit+3-last

    do k = id,jupbnd
      maxerr = iord(nrmax)
      errmax = elist(maxerr)
      if( abs ( blist(maxerr)-alist(maxerr)).gt.small) go to 90
      nrmax = nrmax+1
    end do
!
!  perform extrapolation.
!
   60   numrl2 = numrl2+1
    rlist2(numrl2) = area
    call dqelg(numrl2,rlist2,reseps,abseps,res3la,nres)
    ktmin = ktmin+1
    if(ktmin.gt.5.and.abserr.lt.0.1D-02*errsum) ier = 5
    if(abseps.ge.abserr) go to 70
    ktmin = 0
    abserr = abseps
    result = reseps
    correc = erlarg
    ertest =  max ( epsabs,epsrel* abs ( reseps))
    if(abserr.le.ertest) go to 100
!
!  prepare bisection of the smallest interval.
!
   70   if(numrl2.eq.1) noext = .true.
    if(ier.eq.5) go to 100
    maxerr = iord(1)
    errmax = elist(maxerr)
    nrmax = 1
    extrap = .false.
    small = small*0.5D+00
    erlarg = errsum
    go to 90
   80   small = 0.375D+00
    erlarg = errsum
    ertest = errbnd
    rlist2(2) = area
   90 continue
!
!  set final result and error estimate.
!
  100 if(abserr.eq.oflow) go to 115
  if((ier+ierro).eq.0) go to 110
  if(ierro.eq.3) abserr = abserr+correc
  if(ier.eq.0) ier = 3
  if(result.ne.0.0D+00.and.area.ne.0.0D+00)go to 105
  if(abserr.gt.errsum)go to 115
  if(area.eq.0.0D+00) go to 130
  go to 110
  105 if(abserr/ abs ( result).gt.errsum/ abs ( area))go to 115
!
!  test on divergence
!
  110 continue

  if ( ksgn .eq. (-1) .and. &
    max ( abs ( result), abs ( area)) .le. defabs*0.1D-01 ) then
    go to 130
  end if

  if ( 0.1D-01 .gt. (result/area) .or. &
    (result/area) .gt. 0.1D+03 .or. &
    errsum .gt. abs ( area) ) then
    ier = 6
  end if

  go to 130
!
!  compute global integral sum.
!
  115 result = 0.0D+00
  do k = 1,last
    result = result+rlist(k)
  end do
  abserr = errsum
  130 continue

  neval = 30*last-15
  if(inf.eq.2) neval = 2*neval
  if(ier.gt.2) ier=ier-1

  return
end
subroutine dqagi ( f, bound, inf, epsabs, epsrel, result, abserr, neval, &
  ier,limit,lenw,last,iwork,work)

!*****************************************************************************80
!
!! DQAGI estimates an integral over a semi-infinite or infinite interval.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  the routine calculates an approximation result to a given
!      integral   i = integral of f over (bound,+infinity)
!      or i = integral of f over (-infinity,bound)
!      or i = integral of f over (-infinity,+infinity)
!      hopefully satisfying following claim for accuracy
!      abs(i-result).le.max(epsabs,epsrel*abs(i)).
!
!  Parameters:
!
!   on entry
!      f      - real ( kind = 8 )
!               function subprogram defining the integrand
!               function f(x). the actual name for f needs to be
!               declared e x t e r n a l in the driver program.
!
!      bound  - real ( kind = 8 )
!               finite bound of integration range
!               (has no meaning if interval is doubly-infinite)
!
!      inf    - integer ( kind = 4 )
!               indicating the kind of integration range involved
!               inf = 1 corresponds to  (bound,+infinity),
!               inf = -1            to  (-infinity,bound),
!               inf = 2             to (-infinity,+infinity).
!
!      epsabs - real ( kind = 8 )
!               absolute accuracy requested
!      epsrel - real ( kind = 8 )
!               relative accuracy requested
!               if  epsabs.le.0
!               and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!               the routine will end with ier = 6.
!
!
!   on return
!      result - real ( kind = 8 )
!               approximation to the integral
!
!      abserr - real ( kind = 8 )
!               estimate of the modulus of the absolute error,
!               which should equal or exceed abs(i-result)
!
!      neval  - integer ( kind = 4 )
!               number of integrand evaluations
!
!      ier    - integer ( kind = 4 )
!               ier = 0 normal and reliable termination of the
!                       routine. it is assumed that the requested
!                       accuracy has been achieved.
!             - ier.gt.0 abnormal termination of the routine. the
!                       estimates for result and error are less
!                       reliable. it is assumed that the requested
!                       accuracy has not been achieved.
!      error messages
!               ier = 1 maximum number of subdivisions allowed
!                       has been achieved. one can allow more
!                       subdivisions by increasing the value of
!                       limit (and taking the according dimension
!                       adjustments into account). however, if
!                       this yields no improvement it is advised
!                       to analyze the integrand in order to
!                       determine the integration difficulties. if
!                       the position of a local difficulty can be
!                       determined (e.g. singularity,
!                       discontinuity within the interval) one
!                       will probably gain from splitting up the
!                       interval at this point and calling the
!                       integrator on the subranges. if possible,
!                       an appropriate special-purpose integrator
!                       should be used, which is designed for
!                       handling the type of difficulty involved.
!                   = 2 the occurrence of roundoff error is
!                       detected, which prevents the requested
!                       tolerance from being achieved.
!                       the error may be under-estimated.
!                   = 3 extremely bad integrand behaviour occurs
!                       at some points of the integration
!                       interval.
!                   = 4 the algorithm does not converge.
!                       roundoff error is detected in the
!                       extrapolation table.
!                       it is assumed that the requested tolerance
!                       cannot be achieved, and that the returned
!                       result is the best which can be obtained.
!                   = 5 the integral is probably divergent, or
!                       slowly convergent. it must be noted that
!                       divergence can occur with any other value
!                       of ier.
!                   = 6 the input is invalid, because
!                       (epsabs.le.0 and
!                        epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
!                        or limit.lt.1 or leniw.lt.limit*4.
!                       result, abserr, neval, last are set to
!                       zero. exept when limit or leniw is
!                       invalid, iwork(1), work(limit*2+1) and
!                       work(limit*3+1) are set to zero, work(1)
!                       is set to a and work(limit+1) to b.
!
!   dimensioning parameters
!      limit - integer ( kind = 4 )
!              dimensioning parameter for iwork
!              limit determines the maximum number of subintervals
!              in the partition of the given integration interval
!              (a,b), limit.ge.1.
!              if limit.lt.1, the routine will end with ier = 6.
!
!      lenw  - integer ( kind = 4 )
!              dimensioning parameter for work
!              lenw must be at least limit*4.
!              if lenw.lt.limit*4, the routine will end
!              with ier = 6.
!
!      last  - integer ( kind = 4 )
!              on return, last equals the number of subintervals
!              produced in the subdivision process, which
!              determines the number of significant elements
!              actually in the work arrays.
!
!   work arrays
!      iwork - integer ( kind = 4 )
!              vector of dimension at least limit, the first
!              k elements of which contain pointers
!              to the error estimates over the subintervals,
!              such that work(limit*3+iwork(1)),... ,
!              work(limit*3+iwork(k)) form a decreasing
!              sequence, with k = last if last.le.(limit/2+2), and
!              k = limit+1-last otherwise
!
!      work  - real ( kind = 8 )
!              vector of dimension at least lenw
!              on return
!              work(1), ..., work(last) contain the left
!               end points of the subintervals in the
!               partition of (a,b),
!              work(limit+1), ..., work(limit+last) contain
!               the right end points,
!              work(limit*2+1), ...,work(limit*2+last) contain the
!               integral approximations over the subintervals,
!              work(limit*3+1), ..., work(limit*3)
!               contain the error estimates.
!
  implicit none

  real ( kind = 8 ) abserr,bound,epsabs,epsrel,f,result,work
  integer ( kind = 4 ) ier,inf,iwork,last,lenw,limit,lvl,l1,l2,l3,neval

  dimension iwork(limit),work(lenw)

  external f
!
!  check validity of limit and lenw.
!
  ier = 6
  neval = 0
  last = 0
  result = 0.0D+00
  abserr = 0.0D+00
  if(limit.lt.1.or.lenw.lt.limit*4) go to 10
!
!  prepare call for dqagie.
!
  l1 = limit+1
  l2 = limit+l1
  l3 = limit+l2

  call dqagie(f,bound,inf,epsabs,epsrel,limit,result,abserr, &
    neval,ier,work(1),work(l1),work(l2),work(l3),iwork,last)
!
!  call error handler if necessary.
!
   lvl = 0
10    if(ier.eq.6) lvl = 1

  if(ier.ne.0) then
    call xerror('abnormal return from dqagi',26,ier,lvl)
  end if

  return
end
subroutine dqagpe(f,a,b,npts2,points,epsabs,epsrel,limit,result, &
  abserr,neval,ier,alist,blist,rlist,elist,pts,iord,level,ndin, &
  last)

!*****************************************************************************80
!
!! DQAGPE computes a definite integral.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  the routine calculates an approximation result to a given
!      definite integral i = integral of f over (a,b), hopefully
!      satisfying following claim for accuracy abs(i-result).le.
!      max(epsabs,epsrel*abs(i)). break points of the integration
!      interval, where local difficulties of the integrand may
!      occur(e.g. singularities,discontinuities),provided by user.
!
!  Parameters:
!
!   on entry
!      f      - real ( kind = 8 )
!               function subprogram defining the integrand
!               function f(x). the actual name for f needs to be
!               declared e x t e r n a l in the driver program.
!
!      a      - real ( kind = 8 )
!               lower limit of integration
!
!      b      - real ( kind = 8 )
!               upper limit of integration
!
!      npts2  - integer ( kind = 4 )
!               number equal to two more than the number of
!               user-supplied break points within the integration
!               range, npts2.ge.2.
!               if npts2.lt.2, the routine will end with ier = 6.
!
!      points - real ( kind = 8 )
!               vector of dimension npts2, the first (npts2-2)
!               elements of which are the user provided break
!               points. if these points do not constitute an
!               ascending sequence there will be an automatic
!               sorting.
!
!      epsabs - real ( kind = 8 )
!               absolute accuracy requested
!      epsrel - real ( kind = 8 )
!               relative accuracy requested
!               if  epsabs.le.0
!               and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!               the routine will end with ier = 6.
!
!      limit  - integer ( kind = 4 )
!               gives an upper bound on the number of subintervals
!               in the partition of (a,b), limit.ge.npts2
!               if limit.lt.npts2, the routine will end with
!               ier = 6.
!
!   on return
!      result - real ( kind = 8 )
!               approximation to the integral
!
!      abserr - real ( kind = 8 )
!               estimate of the modulus of the absolute error,
!               which should equal or exceed abs(i-result)
!
!      neval  - integer ( kind = 4 )
!               number of integrand evaluations
!
!      ier    - integer ( kind = 4 )
!               ier = 0 normal and reliable termination of the
!                       routine. it is assumed that the requested
!                       accuracy has been achieved.
!               ier.gt.0 abnormal termination of the routine.
!                       the estimates for integral and error are
!                       less reliable. it is assumed that the
!                       requested accuracy has not been achieved.
!      error messages
!               ier = 1 maximum number of subdivisions allowed
!                       has been achieved. one can allow more
!                       subdivisions by increasing the value of
!                       limit (and taking the according dimension
!                       adjustments into account). however, if
!                       this yields no improvement it is advised
!                       to analyze the integrand in order to
!                       determine the integration difficulties. if
!                       the position of a local difficulty can be
!                       determined (i.e. singularity,
!                       discontinuity within the interval), it
!                       should be supplied to the routine as an
!                       element of the vector points. if necessary
!                       an appropriate special-purpose integrator
!                       must be used, which is designed for
!                       handling the type of difficulty involved.
!                   = 2 the occurrence of roundoff error is
!                       detected, which prevents the requested
!                       tolerance from being achieved.
!                       the error may be under-estimated.
!                   = 3 extremely bad integrand behaviour occurs
!                       at some points of the integration
!                       interval.
!                   = 4 the algorithm does not converge.
!                       roundoff error is detected in the
!                       extrapolation table. it is presumed that
!                       the requested tolerance cannot be
!                       achieved, and that the returned result is
!                       the best which can be obtained.
!                   = 5 the integral is probably divergent, or
!                       slowly convergent. it must be noted that
!                       divergence can occur with any other value
!                       of ier.gt.0.
!                   = 6 the input is invalid because
!                       npts2.lt.2 or
!                       break points are specified outside
!                       the integration range or
!                       (epsabs.le.0 and
!                        epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
!                       or limit.lt.npts2.
!                       result, abserr, neval, last, rlist(1),
!                       and elist(1) are set to zero. alist(1) and
!                       blist(1) are set to a and b respectively.
!
!      alist  - real ( kind = 8 )
!               vector of dimension at least limit, the first
!                last  elements of which are the left end points
!               of the subintervals in the partition of the given
!               integration range (a,b)
!
!      blist  - real ( kind = 8 )
!               vector of dimension at least limit, the first
!                last  elements of which are the right end points
!               of the subintervals in the partition of the given
!               integration range (a,b)
!
!      rlist  - real ( kind = 8 )
!               vector of dimension at least limit, the first
!                last  elements of which are the integral
!               approximations on the subintervals
!
!      elist  - real ( kind = 8 )
!               vector of dimension at least limit, the first
!                last  elements of which are the moduli of the
!               absolute error estimates on the subintervals
!
!      pts    - real ( kind = 8 )
!               vector of dimension at least npts2, containing the
!               integration limits and the break points of the
!               interval in ascending sequence.
!
!      level  - integer ( kind = 4 )
!               vector of dimension at least limit, containing the
!               subdivision levels of the subinterval, i.e. if
!               (aa,bb) is a subinterval of (p1,p2) where p1 as
!               well as p2 is a user-provided break point or
!               integration limit, then (aa,bb) has level l if
!               abs(bb-aa) = abs(p2-p1)*2**(-l).
!
!      ndin   - integer ( kind = 4 )
!               vector of dimension at least npts2, after first
!               integration over the intervals (pts(i)),pts(i+1),
!               i = 0,1, ..., npts2-2, the error estimates over
!               some of the intervals may have been increased
!               artificially, in order to put their subdivision
!               forward. if this happens for the subinterval
!               numbered k, ndin(k) is put to 1, otherwise
!               ndin(k) = 0.
!
!      iord   - integer ( kind = 4 )
!               vector of dimension at least limit, the first k
!               elements of which are pointers to the
!               error estimates over the subintervals,
!               such that elist(iord(1)), ..., elist(iord(k))
!               form a decreasing sequence, with k = last
!               if last.le.(limit/2+2), and k = limit+1-last
!               otherwise
!
!      last   - integer ( kind = 4 )
!               number of subintervals actually produced in the
!               subdivisions process
!
!  Local Parameters:
!
!      the dimension of rlist2 is determined by the value of
!      limexp in routine epsalg (rlist2 should be of dimension
!      (limexp+2) at least).
!
!     alist     - list of left end points of all subintervals
!                 considered up to now
!     blist     - list of right end points of all subintervals
!                 considered up to now
!     rlist(i)  - approximation to the integral over
!                 (alist(i),blist(i))
!     rlist2    - array of dimension at least limexp+2
!                 containing the part of the epsilon table which
!                 is still needed for further computations
!     elist(i)  - error estimate applying to rlist(i)
!     maxerr    - pointer to the interval with largest error
!                 estimate
!     errmax    - elist(maxerr)
!     erlast    - error on the interval currently subdivided
!                 (before that subdivision has taken place)
!     area      - sum of the integrals over the subintervals
!     errsum    - sum of the errors over the subintervals
!     errbnd    - requested accuracy max(epsabs,epsrel*
!                 abs(result))
!     *****1    - variable for the left subinterval
!     *****2    - variable for the right subinterval
!     last      - index for subdivision
!     nres      - number of calls to the extrapolation routine
!     numrl2    - number of elements in rlist2. if an appropriate
!                 approximation to the compounded integral has
!                 been obtained, it is put in rlist2(numrl2) after
!                 numrl2 has been increased by one.
!     erlarg    - sum of the errors over the intervals larger
!                 than the smallest interval considered up to now
!     extrap    - logical variable denoting that the routine
!                 is attempting to perform extrapolation. i.e.
!                 before subdividing the smallest interval we
!                 try to decrease the value of erlarg.
!     noext     - logical variable denoting that extrapolation is
!                 no longer allowed (true-value)
!
!      machine dependent constants
!
!     epmach is the largest relative spacing.
!     uflow is the smallest positive magnitude.
!     oflow is the largest positive magnitude.
!
  implicit none

  real ( kind = 8 ) a,abseps,abserr,alist,area,area1,area12,area2,a1, &
    a2,b,blist,b1,b2,correc,defabs,defab1,defab2, &
    dres,elist,epmach,epsabs,epsrel,erlarg,erlast,errbnd, &
    errmax,error1,erro12,error2,errsum,ertest,f,oflow,points,pts, &
    resa,resabs,reseps,result,res3la,rlist,rlist2,sgn,temp,uflow
  integer ( kind = 4 ) i,id,ier,ierro,ind1,ind2,iord,ip1, &
    iroff1,iroff2,iroff3,j, &
    jlow,jupbnd,k,ksgn,ktmin,last,levcur,level,levmax,limit,maxerr, &
    ndin,neval,nint,nintp1,npts,npts2,nres,nrmax,numrl2
  logical extrap,noext

  dimension alist(limit),blist(limit),elist(limit),iord(limit), &
    level(limit),ndin(npts2),points(npts2),pts(npts2),res3la(3), &
    rlist(limit),rlist2(52)

  external f

  epmach = epsilon ( epmach )
!
!  test on validity of parameters
!
  ier = 0
  neval = 0
  last = 0
  result = 0.0D+00
  abserr = 0.0D+00
  alist(1) = a
  blist(1) = b
  rlist(1) = 0.0D+00
  elist(1) = 0.0D+00
  iord(1) = 0
  level(1) = 0
  npts = npts2-2
  if(npts2.lt.2.or.limit.le.npts.or.(epsabs.le.0.0D+00.and. &
    epsrel.lt. max ( 0.5D+02*epmach,0.5d-28))) ier = 6

  if(ier.eq.6) then
    return
  end if
!
!  if any break points are provided, sort them into an
!  ascending sequence.
!
  sgn = 1.0D+00
  if(a.gt.b) sgn = -1.0D+00
  pts(1) =  min (a,b)
  if(npts.eq.0) go to 15
  do i = 1,npts
    pts(i+1) = points(i)
  end do
   15 pts(npts+2) =  max ( a,b)
  nint = npts+1
  a1 = pts(1)
  if(npts.eq.0) go to 40
  nintp1 = nint+1
  do i = 1,nint
    ip1 = i+1
    do j = ip1,nintp1
      if(pts(i).gt.pts(j)) then
        temp = pts(i)
        pts(i) = pts(j)
        pts(j) = temp
      end if
    end do
  end do
  if(pts(1).ne. min (a,b).or.pts(nintp1).ne. max ( a,b)) ier = 6

  if(ier.eq.6) then
    return
  end if
!
!  compute first integral and error approximations.
!
   40 resabs = 0.0D+00

  do i = 1,nint
    b1 = pts(i+1)
    call dqk21(f,a1,b1,area1,error1,defabs,resa)
    abserr = abserr+error1
    result = result+area1
    ndin(i) = 0
    if(error1.eq.resa.and.error1.ne.0.0D+00) ndin(i) = 1
    resabs = resabs+defabs
    level(i) = 0
    elist(i) = error1
    alist(i) = a1
    blist(i) = b1
    rlist(i) = area1
    iord(i) = i
    a1 = b1
  end do

  errsum = 0.0D+00
  do i = 1,nint
    if(ndin(i).eq.1) elist(i) = abserr
    errsum = errsum+elist(i)
  end do
!
!  test on accuracy.
!
  last = nint
  neval = 21*nint
  dres =  abs ( result)
  errbnd =  max ( epsabs,epsrel*dres)
  if(abserr.le.0.1D+03*epmach*resabs.and.abserr.gt.errbnd) ier = 2
  if(nint.eq.1) go to 80

  do i = 1,npts
    jlow = i+1
    ind1 = iord(i)
    do j = jlow,nint
      ind2 = iord(j)
      if(elist(ind1).le.elist(ind2)) then
        ind1 = ind2
        k = j
      end if
    end do
    if(ind1.ne.iord(i)) then
      iord(k) = iord(i)
      iord(i) = ind1
    end if
  end do

  if(limit.lt.npts2) ier = 1
   80 if(ier.ne.0.or.abserr.le.errbnd) go to 210
!
!  initialization
!
  rlist2(1) = result
  maxerr = iord(1)
  errmax = elist(maxerr)
  area = result
  nrmax = 1
  nres = 0
  numrl2 = 1
  ktmin = 0
  extrap = .false.
  noext = .false.
  erlarg = errsum
  ertest = errbnd
  levmax = 1
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0
  ierro = 0
  uflow = tiny ( uflow )
  oflow = huge ( oflow )
  abserr = oflow
  ksgn = -1
  if(dres.ge.(0.1D+01-0.5D+02*epmach)*resabs) ksgn = 1
!
!  main do-loop
!
  do 160 last = npts2,limit
!
!  bisect the subinterval with the nrmax-th largest error estimate.
!
    levcur = level(maxerr)+1
    a1 = alist(maxerr)
    b1 = 0.5D+00*(alist(maxerr)+blist(maxerr))
    a2 = b1
    b2 = blist(maxerr)
    erlast = errmax
    call dqk21(f,a1,b1,area1,error1,resa,defab1)
    call dqk21(f,a2,b2,area2,error2,resa,defab2)
!
!  improve previous approximations to integral
!  and error and test for accuracy.
!
    neval = neval+42
    area12 = area1+area2
    erro12 = error1+error2
    errsum = errsum+erro12-errmax
    area = area+area12-rlist(maxerr)
    if(defab1.eq.error1.or.defab2.eq.error2) go to 95
    if( abs ( rlist(maxerr)-area12).gt.0.1D-04* abs ( area12) &
    .or.erro12.lt.0.99D+00*errmax) go to 90
    if(extrap) iroff2 = iroff2+1
    if(.not.extrap) iroff1 = iroff1+1
   90   if(last.gt.10.and.erro12.gt.errmax) iroff3 = iroff3+1
   95   level(maxerr) = levcur
    level(last) = levcur
    rlist(maxerr) = area1
    rlist(last) = area2
    errbnd =  max ( epsabs,epsrel* abs ( area))
!
!  test for roundoff error and eventually set error flag.
!
    if(iroff1+iroff2.ge.10.or.iroff3.ge.20) ier = 2
    if(iroff2.ge.5) ierro = 3
!
!  set error flag in the case that the number of
!  subintervals equals limit.
!
    if(last.eq.limit) ier = 1
!
!  set error flag in the case of bad integrand behaviour
!  at a point of the integration range
!
    if( max (  abs ( a1), abs ( b2)).le.(0.1D+01+0.1D+03*epmach)* &
    ( abs ( a2)+0.1D+04*uflow)) ier = 4
!
!  append the newly-created intervals to the list.
!
    if(error2.gt.error1) go to 100
    alist(last) = a2
    blist(maxerr) = b1
    blist(last) = b2
    elist(maxerr) = error1
    elist(last) = error2
    go to 110
  100   alist(maxerr) = a2
    alist(last) = a1
    blist(last) = b1
    rlist(maxerr) = area2
    rlist(last) = area1
    elist(maxerr) = error2
    elist(last) = error1
!
!  call dqpsrt to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with nrmax-th largest error estimate (to be bisected next).
!
  110   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
    if(errsum.le.errbnd) go to 190
    if(ier.ne.0) go to 170
    if(noext) go to 160
    erlarg = erlarg-erlast
    if(levcur+1.le.levmax) erlarg = erlarg+erro12
    if(extrap) go to 120
!
!     test whether the interval to be bisected next is the
!     smallest interval.
!
    if(level(maxerr)+1.le.levmax) go to 160
    extrap = .true.
    nrmax = 2
  120   if(ierro.eq.3.or.erlarg.le.ertest) go to 140
!
!  the smallest interval has the largest error.
!  before bisecting decrease the sum of the errors over
!  the larger intervals (erlarg) and perform extrapolation.
!
    id = nrmax
    jupbnd = last
    if(last.gt.(2+limit/2)) jupbnd = limit+3-last

    do k = id,jupbnd
      maxerr = iord(nrmax)
      errmax = elist(maxerr)
      if(level(maxerr)+1.le.levmax) go to 160
      nrmax = nrmax+1
    end do
!
!  perform extrapolation.
!
  140   numrl2 = numrl2+1
    rlist2(numrl2) = area
    if(numrl2.le.2) go to 155
    call dqelg(numrl2,rlist2,reseps,abseps,res3la,nres)
    ktmin = ktmin+1
    if(ktmin.gt.5.and.abserr.lt.0.1D-02*errsum) ier = 5
    if(abseps.ge.abserr) go to 150
    ktmin = 0
    abserr = abseps
    result = reseps
    correc = erlarg
    ertest =  max ( epsabs,epsrel* abs ( reseps))
    if(abserr.lt.ertest) go to 170
!
!  prepare bisection of the smallest interval.
!
  150   if(numrl2.eq.1) noext = .true.
    if(ier.ge.5) go to 170
  155   maxerr = iord(1)
    errmax = elist(maxerr)
    nrmax = 1
    extrap = .false.
    levmax = levmax+1
    erlarg = errsum
  160 continue
!
!  set the final result.
!
  170 continue

  if(abserr.eq.oflow) go to 190
  if((ier+ierro).eq.0) go to 180
  if(ierro.eq.3) abserr = abserr+correc
  if(ier.eq.0) ier = 3
  if(result.ne.0.0D+00.and.area.ne.0.0D+00)go to 175
  if(abserr.gt.errsum)go to 190
  if(area.eq.0.0D+00) go to 210
  go to 180
  175 if(abserr/ abs ( result).gt.errsum/ abs ( area))go to 190
!
!  test on divergence.
!
  180 if(ksgn.eq.(-1).and. max (  abs ( result), abs ( area)).le. &
    resabs*0.1D-01) go to 210
  if(0.1D-01.gt.(result/area).or.(result/area).gt.0.1D+03.or. &
    errsum.gt. abs ( area)) ier = 6
  go to 210
!
!  compute global integral sum.
!
  190 result = 0.0D+00
  do k = 1,last
    result = result+rlist(k)
  end do

  abserr = errsum
  210 if(ier.gt.2) ier = ier-1
  result = result*sgn

  return
end
subroutine dqagp ( f, a, b, npts2, points, epsabs, epsrel, result, abserr, &
  neval, ier, leniw, lenw, last, iwork, work )

!*****************************************************************************80
!
!! DQAGP computes a definite integral.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  the routine calculates an approximation result to a given
!      definite integral i = integral of f over (a,b),
!      hopefully satisfying following claim for accuracy
!      break points of the integration interval, where local
!      difficulties of the integrand may occur (e.g.
!      singularities, discontinuities), are provided by the user.
!
!  Parameters:
!
!   on entry
!      f      - real ( kind = 8 )
!               function subprogram defining the integrand
!               function f(x). the actual name for f needs to be
!               declared e x t e r n a l in the driver program.
!
!      a      - real ( kind = 8 )
!               lower limit of integration
!
!      b      - real ( kind = 8 )
!               upper limit of integration
!
!      npts2  - integer ( kind = 4 )
!               number equal to two more than the number of
!               user-supplied break points within the integration
!               range, npts.ge.2.
!               if npts2.lt.2, the routine will end with ier = 6.
!
!      points - real ( kind = 8 )
!               vector of dimension npts2, the first (npts2-2)
!               elements of which are the user provided break
!               points. if these points do not constitute an
!               ascending sequence there will be an automatic
!               sorting.
!
!      epsabs - real ( kind = 8 )
!               absolute accuracy requested
!      epsrel - real ( kind = 8 )
!               relative accuracy requested
!               if  epsabs.le.0
!               and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!               the routine will end with ier = 6.
!
!   on return
!      result - real ( kind = 8 )
!               approximation to the integral
!
!      abserr - real ( kind = 8 )
!               estimate of the modulus of the absolute error,
!               which should equal or exceed abs(i-result)
!
!      neval  - integer ( kind = 4 )
!               number of integrand evaluations
!
!      ier    - integer ( kind = 4 )
!               ier = 0 normal and reliable termination of the
!                       routine. it is assumed that the requested
!                       accuracy has been achieved.
!               ier.gt.0 abnormal termination of the routine.
!                       the estimates for integral and error are
!                       less reliable. it is assumed that the
!                       requested accuracy has not been achieved.
!      error messages
!               ier = 1 maximum number of subdivisions allowed
!                       has been achieved. one can allow more
!                       subdivisions by increasing the value of
!                       limit (and taking the according dimension
!                       adjustments into account). however, if
!                       this yields no improvement it is advised
!                       to analyze the integrand in order to
!                       determine the integration difficulties. if
!                       the position of a local difficulty can be
!                       determined (i.e. singularity,
!                       discontinuity within the interval), it
!                       should be supplied to the routine as an
!                       element of the vector points. if necessary
!                       an appropriate special-purpose integrator
!                       must be used, which is designed for
!                       handling the type of difficulty involved.
!                   = 2 the occurrence of roundoff error is
!                       detected, which prevents the requested
!                       tolerance from being achieved.
!                       the error may be under-estimated.
!                   = 3 extremely bad integrand behaviour occurs
!                       at some points of the integration
!                       interval.
!                   = 4 the algorithm does not converge.
!                       roundoff error is detected in the
!                       extrapolation table.
!                       it is presumed that the requested
!                       tolerance cannot be achieved, and that
!                       the returned result is the best which
!                       can be obtained.
!                   = 5 the integral is probably divergent, or
!                       slowly convergent. it must be noted that
!                       divergence can occur with any other value
!                       of ier.gt.0.
!                   = 6 the input is invalid because
!                       npts2.lt.2 or
!                       break points are specified outside
!                       the integration range or
!                       (epsabs.le.0 and
!                        epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
!                       result, abserr, neval, last are set to
!                       zero. exept when leniw or lenw or npts2 is
!                       invalid, iwork(1), iwork(limit+1),
!                       work(limit*2+1) and work(limit*3+1)
!                       are set to zero.
!                       work(1) is set to a and work(limit+1)
!                       to b (where limit = (leniw-npts2)/2).
!
!   dimensioning parameters
!      leniw - integer ( kind = 4 )
!              dimensioning parameter for iwork
!              leniw determines limit = (leniw-npts2)/2,
!              which is the maximum number of subintervals in the
!              partition of the given integration interval (a,b),
!              leniw.ge.(3*npts2-2).
!              if leniw.lt.(3*npts2-2), the routine will end with
!              ier = 6.
!
!      lenw  - integer ( kind = 4 )
!              dimensioning parameter for work
!              lenw must be at least leniw*2-npts2.
!              if lenw.lt.leniw*2-npts2, the routine will end
!              with ier = 6.
!
!      last  - integer ( kind = 4 )
!              on return, last equals the number of subintervals
!              produced in the subdivision process, which
!              determines the number of significant elements
!              actually in the work arrays.
!
!   work arrays
!      iwork - integer ( kind = 4 )
!              vector of dimension at least leniw. on return,
!              the first k elements of which contain
!              pointers to the error estimates over the
!              subintervals, such that work(limit*3+iwork(1)),...,
!              work(limit*3+iwork(k)) form a decreasing
!              sequence, with k = last if last.le.(limit/2+2), and
!              k = limit+1-last otherwise
!              iwork(limit+1), ...,iwork(limit+last) contain the
!               subdivision levels of the subintervals, i.e.
!               if (aa,bb) is a subinterval of (p1,p2)
!               where p1 as well as p2 is a user-provided
!               break point or integration limit, then (aa,bb) has
!               level l if abs(bb-aa) = abs(p2-p1)*2**(-l),
!              iwork(limit*2+1), ..., iwork(limit*2+npts2) have
!               no significance for the user,
!              note that limit = (leniw-npts2)/2.
!
!      work  - real ( kind = 8 )
!              vector of dimension at least lenw
!              on return
!              work(1), ..., work(last) contain the left
!               end points of the subintervals in the
!               partition of (a,b),
!              work(limit+1), ..., work(limit+last) contain
!               the right end points,
!              work(limit*2+1), ..., work(limit*2+last) contain
!               the integral approximations over the subintervals,
!              work(limit*3+1), ..., work(limit*3+last)
!               contain the corresponding error estimates,
!              work(limit*4+1), ..., work(limit*4+npts2)
!               contain the integration limits and the
!               break points sorted in an ascending sequence.
!              note that limit = (leniw-npts2)/2.
!
  implicit none

  real ( kind = 8 ) a,abserr,b,epsabs,epsrel,f,points,result,work
  integer ( kind = 4 ) ier,iwork,last,leniw,lenw,limit,lvl,l1,l2,l3, &
    l4,neval,npts2

  dimension iwork(leniw),points(npts2),work(lenw)

  external f
!
!  check validity of limit and lenw.
!
  ier = 6
  neval = 0
  last = 0
  result = 0.0D+00
  abserr = 0.0D+00
  if(leniw.lt.(3*npts2-2).or.lenw.lt.(leniw*2-npts2).or.npts2.lt.2) &
    go to 10
!
!  prepare call for dqagpe.
!
  limit = (leniw-npts2)/2
  l1 = limit+1
  l2 = limit+l1
  l3 = limit+l2
  l4 = limit+l3

  call dqagpe(f,a,b,npts2,points,epsabs,epsrel,limit,result,abserr, &
    neval,ier,work(1),work(l1),work(l2),work(l3),work(l4), &
    iwork(1),iwork(l1),iwork(l2),last)
!
!  call error handler if necessary.
!
  lvl = 0
10    if(ier.eq.6) lvl = 1

  if(ier.ne.0) then
    call xerror('abnormal return from dqagp',26,ier,lvl)
  end if

  return
end
subroutine dqagse(f,a,b,epsabs,epsrel,limit,result,abserr,neval, &
     ier,alist,blist,rlist,elist,iord,last)

!*****************************************************************************80
!
!! DQAGSE estimates the integral of a function.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  the routine calculates an approximation result to a given
!      definite integral i = integral of f over (a,b),
!      hopefully satisfying following claim for accuracy
!      abs(i-result).le.max(epsabs,epsrel*abs(i)).
!
!  Parameters:
!
!   on entry
!      f      - real ( kind = 8 )
!               function subprogram defining the integrand
!               function f(x). the actual name for f needs to be
!               declared e x t e r n a l in the driver program.
!
!      a      - real ( kind = 8 )
!               lower limit of integration
!
!      b      - real ( kind = 8 )
!               upper limit of integration
!
!      epsabs - real ( kind = 8 )
!               absolute accuracy requested
!      epsrel - real ( kind = 8 )
!               relative accuracy requested
!               if  epsabs.le.0
!               and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!               the routine will end with ier = 6.
!
!      limit  - integer ( kind = 4 )
!               gives an upperbound on the number of subintervals
!               in the partition of (a,b)
!
!   on return
!      result - real ( kind = 8 )
!               approximation to the integral
!
!      abserr - real ( kind = 8 )
!               estimate of the modulus of the absolute error,
!               which should equal or exceed abs(i-result)
!
!      neval  - integer ( kind = 4 )
!               number of integrand evaluations
!
!      ier    - integer ( kind = 4 )
!               ier = 0 normal and reliable termination of the
!                       routine. it is assumed that the requested
!                       accuracy has been achieved.
!               ier.gt.0 abnormal termination of the routine
!                       the estimates for integral and error are
!                       less reliable. it is assumed that the
!                       requested accuracy has not been achieved.
!      error messages
!                   = 1 maximum number of subdivisions allowed
!                       has been achieved. one can allow more sub-
!                       divisions by increasing the value of limit
!                       (and taking the according dimension
!                       adjustments into account). however, if
!                       this yields no improvement it is advised
!                       to analyze the integrand in order to
!                       determine the integration difficulties. if
!                       the position of a local difficulty can be
!                       determined (e.g. singularity,
!                       discontinuity within the interval) one
!                       will probably gain from splitting up the
!                       interval at this point and calling the
!                       integrator on the subranges. if possible,
!                       an appropriate special-purpose integrator
!                       should be used, which is designed for
!                       handling the type of difficulty involved.
!                   = 2 the occurrence of roundoff error is detec-
!                       ted, which prevents the requested
!                       tolerance from being achieved.
!                       the error may be under-estimated.
!                   = 3 extremely bad integrand behaviour
!                       occurs at some points of the integration
!                       interval.
!                   = 4 the algorithm does not converge.
!                       roundoff error is detected in the
!                       extrapolation table.
!                       it is presumed that the requested
!                       tolerance cannot be achieved, and that the
!                       returned result is the best which can be
!                       obtained.
!                   = 5 the integral is probably divergent, or
!                       slowly convergent. it must be noted that
!                       divergence can occur with any other value
!                       of ier.
!                   = 6 the input is invalid, because
!                       epsabs.le.0 and
!                       epsrel.lt.max(50*rel.mach.acc.,0.5d-28).
!                       result, abserr, neval, last, rlist(1),
!                       iord(1) and elist(1) are set to zero.
!                       alist(1) and blist(1) are set to a and b
!                       respectively.
!
!      alist  - real ( kind = 8 )
!               vector of dimension at least limit, the first
!                last  elements of which are the left end points
!               of the subintervals in the partition of the
!               given integration range (a,b)
!
!      blist  - real ( kind = 8 )
!               vector of dimension at least limit, the first
!                last  elements of which are the right end points
!               of the subintervals in the partition of the given
!               integration range (a,b)
!
!      rlist  - real ( kind = 8 )
!               vector of dimension at least limit, the first
!                last  elements of which are the integral
!               approximations on the subintervals
!
!      elist  - real ( kind = 8 )
!               vector of dimension at least limit, the first
!                last  elements of which are the moduli of the
!               absolute error estimates on the subintervals
!
!      iord   - integer ( kind = 4 )
!               vector of dimension at least limit, the first k
!               elements of which are pointers to the
!               error estimates over the subintervals,
!               such that elist(iord(1)), ..., elist(iord(k))
!               form a decreasing sequence, with k = last
!               if last.le.(limit/2+2), and k = limit+1-last
!               otherwise
!
!      last   - integer ( kind = 4 )
!               number of subintervals actually produced in the
!               subdivision process
!
!  Local parameters:
!
!      the dimension of rlist2 is determined by the value of
!      limexp in routine dqelg (rlist2 should be of dimension
!      (limexp+2) at least).
!
!      list of major variables
!
!     alist     - list of left end points of all subintervals
!                 considered up to now
!     blist     - list of right end points of all subintervals
!                 considered up to now
!     rlist(i)  - approximation to the integral over
!                 (alist(i),blist(i))
!     rlist2    - array of dimension at least limexp+2 containing
!                 the part of the epsilon table which is still
!                 needed for further computations
!     elist(i)  - error estimate applying to rlist(i)
!     maxerr    - pointer to the interval with largest error
!                 estimate
!     errmax    - elist(maxerr)
!     erlast    - error on the interval currently subdivided
!                 (before that subdivision has taken place)
!     area      - sum of the integrals over the subintervals
!     errsum    - sum of the errors over the subintervals
!     errbnd    - requested accuracy max(epsabs,epsrel*
!                 abs(result))
!     *****1    - variable for the left interval
!     *****2    - variable for the right interval
!     last      - index for subdivision
!     nres      - number of calls to the extrapolation routine
!     numrl2    - number of elements currently in rlist2. if an
!                 appropriate approximation to the compounded
!                 integral has been obtained it is put in
!                 rlist2(numrl2) after numrl2 has been increased
!                 by one.
!     small     - length of the smallest interval considered up
!                 to now, multiplied by 1.5
!     erlarg    - sum of the errors over the intervals larger
!                 than the smallest interval considered up to now
!     extrap    - logical variable denoting that the routine is
!                 attempting to perform extrapolation i.e. before
!                 subdividing the smallest interval we try to
!                 decrease the value of erlarg.
!     noext     - logical variable denoting that extrapolation
!                 is no longer allowed (true value)
!
!      machine dependent constants
!
!     epmach is the largest relative spacing.
!     uflow is the smallest positive magnitude.
!     oflow is the largest positive magnitude.
!
  implicit none

  real ( kind = 8 ) a,abseps,abserr,alist,area,area1,area12,area2,a1, &
    a2,b,blist,b1,b2,correc,defabs,defab1,defab2, &
    dres,elist,epmach,epsabs,epsrel,erlarg,erlast,errbnd,errmax, &
    error1,error2,erro12,errsum,ertest,f,oflow,resabs,reseps,result, &
    res3la,rlist,rlist2,small,uflow
  integer ( kind = 4 ) id,ier,ierro,iord,iroff1,iroff2,iroff3,jupbnd, &
    k,ksgn,ktmin,last,limit,maxerr,neval,nres,nrmax,numrl2
  logical extrap,noext
  dimension alist(limit),blist(limit),elist(limit),iord(limit), &
   res3la(3),rlist(limit),rlist2(52)

  external f

  epmach = epsilon ( epmach )
!
!  test on validity of parameters
!
  ier = 0
  neval = 0
  last = 0
  result = 0.0D+00
  abserr = 0.0D+00
  alist(1) = a
  blist(1) = b
  rlist(1) = 0.0D+00
  elist(1) = 0.0D+00

  if(epsabs.le.0.0D+00.and.epsrel.lt. max ( 0.5D+02*epmach,0.5d-28)) then
    ier = 6
    return
  end if
!
!  first approximation to the integral
!
  uflow = tiny ( uflow )
  oflow = huge ( oflow )
  ierro = 0
  call dqk21(f,a,b,result,abserr,defabs,resabs)
!
!  test on accuracy.
!
  dres =  abs ( result)
  errbnd =  max ( epsabs,epsrel*dres)
  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1
  if(abserr.le.1.0D+02*epmach*defabs.and.abserr.gt.errbnd) ier = 2
  if(limit.eq.1) ier = 1
  if(ier.ne.0.or.(abserr.le.errbnd.and.abserr.ne.resabs).or. &
    abserr.eq.0.0D+00) go to 140
!
!  initialization
!
  rlist2(1) = result
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  abserr = oflow
  nrmax = 1
  nres = 0
  numrl2 = 2
  ktmin = 0
  extrap = .false.
  noext = .false.
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0
  ksgn = -1
  if(dres.ge.(0.1D+01-0.5D+02*epmach)*defabs) ksgn = 1
!
!  main do-loop
!
  do 90 last = 2,limit
!
!  bisect the subinterval with the nrmax-th largest error estimate.
!
    a1 = alist(maxerr)
    b1 = 0.5D+00*(alist(maxerr)+blist(maxerr))
    a2 = b1
    b2 = blist(maxerr)
    erlast = errmax
    call dqk21(f,a1,b1,area1,error1,resabs,defab1)
    call dqk21(f,a2,b2,area2,error2,resabs,defab2)
!
!  improve previous approximations to integral
!  and error and test for accuracy.
!
    area12 = area1+area2
    erro12 = error1+error2
    errsum = errsum+erro12-errmax
    area = area+area12-rlist(maxerr)
    if(defab1.eq.error1.or.defab2.eq.error2) go to 15
    if( abs ( rlist(maxerr)-area12).gt.0.1D-04* abs ( area12) &
    .or.erro12.lt.0.99D+00*errmax) go to 10
    if(extrap) iroff2 = iroff2+1
    if(.not.extrap) iroff1 = iroff1+1
   10   if(last.gt.10.and.erro12.gt.errmax) iroff3 = iroff3+1
   15   rlist(maxerr) = area1
    rlist(last) = area2
    errbnd =  max ( epsabs,epsrel* abs ( area))
!
!  test for roundoff error and eventually set error flag.
!
    if(iroff1+iroff2.ge.10.or.iroff3.ge.20) ier = 2
    if(iroff2.ge.5) ierro = 3
!
!  set error flag in the case that the number of subintervals
!  equals limit.
!
    if(last.eq.limit) ier = 1
!
!  set error flag in the case of bad integrand behaviour
!  at a point of the integration range.
!
    if( max (  abs ( a1), abs ( b2)).le.(0.1D+01+0.1D+03*epmach)* &
    ( abs ( a2)+0.1D+04*uflow)) ier = 4
!
!  append the newly-created intervals to the list.
!
    if(error2.gt.error1) go to 20
    alist(last) = a2
    blist(maxerr) = b1
    blist(last) = b2
    elist(maxerr) = error1
    elist(last) = error2
    go to 30
   20   alist(maxerr) = a2
    alist(last) = a1
    blist(last) = b1
    rlist(maxerr) = area2
    rlist(last) = area1
    elist(maxerr) = error2
    elist(last) = error1
!
!  call dqpsrt to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with nrmax-th largest error estimate (to be bisected next).
!
   30   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
    if(errsum.le.errbnd) go to 115
    if(ier.ne.0) go to 100
    if(last.eq.2) go to 80
    if(noext) go to 90
    erlarg = erlarg-erlast
    if( abs ( b1-a1).gt.small) erlarg = erlarg+erro12
    if(extrap) go to 40
!
!  test whether the interval to be bisected next is the
!  smallest interval.
!
    if( abs ( blist(maxerr)-alist(maxerr)).gt.small) go to 90
    extrap = .true.
    nrmax = 2
   40   if(ierro.eq.3.or.erlarg.le.ertest) go to 60
!
!  the smallest interval has the largest error.
!  before bisecting decrease the sum of the errors over the
!  larger intervals (erlarg) and perform extrapolation.
!
    id = nrmax
    jupbnd = last
    if(last.gt.(2+limit/2)) jupbnd = limit+3-last
    do k = id,jupbnd
      maxerr = iord(nrmax)
      errmax = elist(maxerr)
      if( abs ( blist(maxerr)-alist(maxerr)).gt.small) go to 90
      nrmax = nrmax+1
    end do
!
!  perform extrapolation.
!
   60   numrl2 = numrl2+1
    rlist2(numrl2) = area
    call dqelg(numrl2,rlist2,reseps,abseps,res3la,nres)
    ktmin = ktmin+1
    if(ktmin.gt.5.and.abserr.lt.0.1D-02*errsum) ier = 5
    if(abseps.ge.abserr) go to 70
    ktmin = 0
    abserr = abseps
    result = reseps
    correc = erlarg
    ertest =  max ( epsabs,epsrel* abs ( reseps))
    if(abserr.le.ertest) go to 100
!
!  prepare bisection of the smallest interval.
!
   70   if(numrl2.eq.1) noext = .true.
    if(ier.eq.5) go to 100
    maxerr = iord(1)
    errmax = elist(maxerr)
    nrmax = 1
    extrap = .false.
    small = small*0.5D+00
    erlarg = errsum
    go to 90
   80   small =  abs ( b-a)*0.375D+00
    erlarg = errsum
    ertest = errbnd
    rlist2(2) = area
   90 continue
!
!  set final result and error estimate.
!
  100 if(abserr.eq.oflow) go to 115
  if(ier+ierro.eq.0) go to 110
  if(ierro.eq.3) abserr = abserr+correc
  if(ier.eq.0) ier = 3
  if(result.ne.0.0D+00.and.area.ne.0.0D+00) go to 105
  if(abserr.gt.errsum) go to 115
  if(area.eq.0.0D+00) go to 130
  go to 110
  105 if(abserr/ abs ( result).gt.errsum/ abs ( area)) go to 115
!
!  test on divergence.
!
  110 if(ksgn.eq.(-1).and. max (  abs ( result), abs ( area)).le. &
   defabs*0.1D-01) go to 130
  if(0.1D-01.gt.(result/area).or.(result/area).gt.0.1D+03 &
   .or.errsum.gt. abs ( area)) ier = 6
  go to 130
!
!  compute global integral sum.
!
  115 result = 0.0D+00
  do k = 1,last
     result = result+rlist(k)
  end do
  abserr = errsum
  130 if(ier.gt.2) ier = ier-1
  140 neval = 42*last-21

  return
end
subroutine dqags ( f, a, b, epsabs, epsrel, result, abserr, neval, ier, &
  limit, lenw, last, iwork, work )

!*****************************************************************************80
!
!! DQAGS estimates the integral of a function.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  the routine calculates an approximation result to a given
!      definite integral  i = integral of f over (a,b),
!      hopefully satisfying following claim for accuracy
!      abs(i-result).le.max(epsabs,epsrel*abs(i)).
!
!  Parameters:
!
!   on entry
!      f      - real ( kind = 8 )
!               function subprogram defining the integrand
!               function f(x). the actual name for f needs to be
!               declared e x t e r n a l in the driver program.
!
!      a      - real ( kind = 8 )
!               lower limit of integration
!
!      b      - real ( kind = 8 )
!               upper limit of integration
!
!      epsabs - real ( kind = 8 )
!               absolute accuracy requested
!      epsrel - real ( kind = 8 )
!               relative accuracy requested
!               if  epsabs.le.0
!               and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!               the routine will end with ier = 6.
!
!   on return
!      result - real ( kind = 8 )
!               approximation to the integral
!
!      abserr - real ( kind = 8 )
!               estimate of the modulus of the absolute error,
!               which should equal or exceed abs(i-result)
!
!      neval  - integer ( kind = 4 )
!               number of integrand evaluations
!
!      ier    - integer ( kind = 4 )
!               ier = 0 normal and reliable termination of the
!                       routine. it is assumed that the requested
!                       accuracy has been achieved.
!               ier.gt.0 abnormal termination of the routine
!                       the estimates for integral and error are
!                       less reliable. it is assumed that the
!                       requested accuracy has not been achieved.
!      error messages
!               ier = 1 maximum number of subdivisions allowed
!                       has been achieved. one can allow more sub-
!                       divisions by increasing the value of limit
!                       (and taking the according dimension
!                       adjustments into account. however, if
!                       this yields no improvement it is advised
!                       to analyze the integrand in order to
!                       determine the integration difficulties. if
!                       the position of a local difficulty can be
!                       determined (e.g. singularity,
!                       discontinuity within the interval) one
!                       will probably gain from splitting up the
!                       interval at this point and calling the
!                       integrator on the subranges. if possible,
!                       an appropriate special-purpose integrator
!                       should be used, which is designed for
!                       handling the type of difficulty involved.
!                   = 2 the occurrence of roundoff error is detec-
!                       ted, which prevents the requested
!                       tolerance from being achieved.
!                       the error may be under-estimated.
!                   = 3 extremely bad integrand behaviour
!                       occurs at some points of the integration
!                       interval.
!                   = 4 the algorithm does not converge.
!                       roundoff error is detected in the
!                       extrapolation table. it is presumed that
!                       the requested tolerance cannot be
!                       achieved, and that the returned result is
!                       the best which can be obtained.
!                   = 5 the integral is probably divergent, or
!                       slowly convergent. it must be noted that
!                       divergence can occur with any other value
!                       of ier.
!                   = 6 the input is invalid, because
!                       (epsabs.le.0 and
!                        epsrel.lt.max(50*rel.mach.acc.,0.5d-28)
!                       or limit.lt.1 or lenw.lt.limit*4.
!                       result, abserr, neval, last are set to
!                       zero.except when limit or lenw is invalid,
!                       iwork(1), work(limit*2+1) and
!                       work(limit*3+1) are set to zero, work(1)
!                       is set to a and work(limit+1) to b.
!
!   dimensioning parameters
!      limit - integer ( kind = 4 )
!              dimensioning parameter for iwork
!              limit determines the maximum number of subintervals
!              in the partition of the given integration interval
!              (a,b), limit.ge.1.
!              if limit.lt.1, the routine will end with ier = 6.
!
!      lenw  - integer ( kind = 4 )
!              dimensioning parameter for work
!              lenw must be at least limit*4.
!              if lenw.lt.limit*4, the routine will end
!              with ier = 6.
!
!      last  - integer ( kind = 4 )
!              on return, last equals the number of subintervals
!              produced in the subdivision process, detemines the
!              number of significant elements actually in the work
!              arrays.
!
!   work arrays
!      iwork - integer ( kind = 4 )
!              vector of dimension at least limit, the first k
!              elements of which contain pointers
!              to the error estimates over the subintervals
!              such that work(limit*3+iwork(1)),... ,
!              work(limit*3+iwork(k)) form a decreasing
!              sequence, with k = last if last.le.(limit/2+2),
!              and k = limit+1-last otherwise
!
!      work  - real ( kind = 8 )
!              vector of dimension at least lenw
!              on return
!              work(1), ..., work(last) contain the left
!               end-points of the subintervals in the
!               partition of (a,b),
!              work(limit+1), ..., work(limit+last) contain
!               the right end-points,
!              work(limit*2+1), ..., work(limit*2+last) contain
!               the integral approximations over the subintervals,
!              work(limit*3+1), ..., work(limit*3+last)
!               contain the error estimates.
!
  implicit none

  real ( kind = 8 ) a,abserr,b,epsabs,epsrel,f,result,work
  integer ( kind = 4 ) ier,iwork,last,lenw,limit,lvl,l1,l2,l3,neval
  dimension iwork(limit),work(lenw)

  external f
!
!  check validity of limit and lenw.
!
  ier = 6
  neval = 0
  last = 0
  result = 0.0D+00
  abserr = 0.0D+00
  if(limit.lt.1.or.lenw.lt.limit*4) go to 10
!
!  prepare call for dqagse.
!
  l1 = limit+1
  l2 = limit+l1
  l3 = limit+l2

  call dqagse(f,a,b,epsabs,epsrel,limit,result,abserr,neval, &
    ier,work(1),work(l1),work(l2),work(l3),iwork,last)
!
!  call error handler if necessary.
!
  lvl = 0
10    if(ier.eq.6) lvl = 1
  if(ier.ne.0) call xerror('abnormal return from dqags',26,ier,lvl)

  return
end
subroutine dqawce(f,a,b,c,epsabs,epsrel,limit,result,abserr,neval, &
  ier,alist,blist,rlist,elist,iord,last)

!*****************************************************************************80
!
!! DQAWCE computes a Cauchy principal value.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***  purpose  the routine calculates an approximation result to a
!        cauchy principal value i = integral of f*w over (a,b)
!        (w(x) = 1/(x-c), (c.ne.a, c.ne.b), hopefully satisfying
!        following claim for accuracy
!        abs(i-result).le.max(epsabs,epsrel*abs(i))
!
!  Parameters:
!
!   on entry
!      f      - real ( kind = 8 )
!               function subprogram defining the integrand
!               function f(x). the actual name for f needs to be
!               declared e x t e r n a l in the driver program.
!
!      a      - real ( kind = 8 )
!               lower limit of integration
!
!      b      - real ( kind = 8 )
!               upper limit of integration
!
!      c      - real ( kind = 8 )
!               parameter in the weight function, c.ne.a, c.ne.b
!               if c = a or c = b, the routine will end with
!               ier = 6.
!
!      epsabs - real ( kind = 8 )
!               absolute accuracy requested
!      epsrel - real ( kind = 8 )
!               relative accuracy requested
!               if  epsabs.le.0
!               and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!               the routine will end with ier = 6.
!
!      limit  - integer ( kind = 4 )
!               gives an upper bound on the number of subintervals
!               in the partition of (a,b), limit.ge.1
!
!   on return
!      result - real ( kind = 8 )
!               approximation to the integral
!
!      abserr - real ( kind = 8 )
!               estimate of the modulus of the absolute error,
!               which should equal or exceed abs(i-result)
!
!      neval  - integer ( kind = 4 )
!               number of integrand evaluations
!
!      ier    - integer ( kind = 4 )
!               ier = 0 normal and reliable termination of the
!                       routine. it is assumed that the requested
!                       accuracy has been achieved.
!               ier.gt.0 abnormal termination of the routine
!                       the estimates for integral and error are
!                       less reliable. it is assumed that the
!                       requested accuracy has not been achieved.
!      error messages
!               ier = 1 maximum number of subdivisions allowed
!                       has been achieved. one can allow more sub-
!                       divisions by increasing the value of
!                       limit. however, if this yields no
!                       improvement it is advised to analyze the
!                       the integrand, in order to determine the
!                       the integration difficulties. if the
!                       position of a local difficulty can be
!                       determined (e.g. singularity,
!                       discontinuity within the interval) one
!                       will probably gain from splitting up the
!                       interval at this point and calling
!                       appropriate integrators on the subranges.
!                   = 2 the occurrence of roundoff error is detec-
!                       ted, which prevents the requested
!                       tolerance from being achieved.
!                   = 3 extremely bad integrand behaviour
!                       occurs at some interior points of
!                       the integration interval.
!                   = 6 the input is invalid, because
!                       c = a or c = b or
!                       (epsabs.le.0 and
!                        epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
!                       or limit.lt.1.
!                       result, abserr, neval, rlist(1), elist(1),
!                       iord(1) and last are set to zero. alist(1)
!                       and blist(1) are set to a and b
!                       respectively.
!
!      alist   - real ( kind = 8 )
!                vector of dimension at least limit, the first
!                 last  elements of which are the left
!                end points of the subintervals in the partition
!                of the given integration range (a,b)
!
!      blist   - real ( kind = 8 )
!                vector of dimension at least limit, the first
!                 last  elements of which are the right
!                end points of the subintervals in the partition
!                of the given integration range (a,b)
!
!      rlist   - real ( kind = 8 )
!                vector of dimension at least limit, the first
!                 last  elements of which are the integral
!                approximations on the subintervals
!
!      elist   - real ( kind = 8 )
!                vector of dimension limit, the first  last
!                elements of which are the moduli of the absolute
!                error estimates on the subintervals
!
!      iord    - integer ( kind = 4 )
!                vector of dimension at least limit, the first k
!                elements of which are pointers to the error
!                estimates over the subintervals, so that
!                elist(iord(1)), ..., elist(iord(k)) with k = last
!                if last.le.(limit/2+2), and k = limit+1-last
!                otherwise, form a decreasing sequence
!
!      last    - integer ( kind = 4 )
!                number of subintervals actually produced in
!                the subdivision process
!
!  Local Parameters:
!
!     alist     - list of left end points of all subintervals
!                 considered up to now
!     blist     - list of right end points of all subintervals
!                 considered up to now
!     rlist(i)  - approximation to the integral over
!                 (alist(i),blist(i))
!     elist(i)  - error estimate applying to rlist(i)
!     maxerr    - pointer to the interval with largest
!                 error estimate
!     errmax    - elist(maxerr)
!     area      - sum of the integrals over the subintervals
!     errsum    - sum of the errors over the subintervals
!     errbnd    - requested accuracy max(epsabs,epsrel*
!                 abs(result))
!     *****1    - variable for the left subinterval
!     *****2    - variable for the right subinterval
!     last      - index for subdivision
!
!
!      machine dependent constants
!
!     epmach is the largest relative spacing.
!     uflow is the smallest positive magnitude.
!
  implicit none

  real ( kind = 8 ) a,aa,abserr,alist,area,area1,area12,area2,a1,a2, &
    b,bb,blist,b1,b2,c,elist,epmach,epsabs,epsrel, &
    errbnd,errmax,error1,erro12,error2,errsum,f,result,rlist,uflow
  integer ( kind = 4 ) ier,iord,iroff1,iroff2,k,krule,last,limit,&
    maxerr, nev, &
    neval,nrmax
  dimension alist(limit),blist(limit),rlist(limit),elist(limit), &
    iord(limit)

  external f

  epmach = epsilon ( epmach )
  uflow = tiny ( uflow )
!
!  test on validity of parameters
!
  ier = 6
  neval = 0
  last = 0
  alist(1) = a
  blist(1) = b
  rlist(1) = 0.0D+00
  elist(1) = 0.0D+00
  iord(1) = 0
  result = 0.0D+00
  abserr = 0.0D+00

  if ( c.eq.a .or. &
    c.eq.b .or. &
    (epsabs.le.0.0D+00 .and. epsrel.lt. max ( 0.5D+02*epmach,0.5d-28)) ) then
    ier = 6
    return
  end if
!
!  first approximation to the integral
!
  if ( a <= b ) then
    aa=a
    bb=b
  else
    aa=b
    bb=a
  end if

  ier=0
  krule = 1
  call dqc25c(f,aa,bb,c,result,abserr,krule,neval)
  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1
  alist(1) = a
  blist(1) = b
!
!  test on accuracy
!
  errbnd =  max ( epsabs,epsrel* abs ( result))
  if(limit.eq.1) ier = 1

  if(abserr.lt. min (0.1D-01* abs ( result),errbnd) &
    .or.ier.eq.1) go to 70
!
!  initialization
!
  alist(1) = aa
  blist(1) = bb
  rlist(1) = result
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  nrmax = 1
  iroff1 = 0
  iroff2 = 0
!
!  main do-loop
!
  do 40 last = 2,limit
!
!  bisect the subinterval with nrmax-th largest error estimate.
!
    a1 = alist(maxerr)
    b1 = 0.5D+00*(alist(maxerr)+blist(maxerr))
    b2 = blist(maxerr)
    if(c.le.b1.and.c.gt.a1) b1 = 0.5D+00*(c+b2)
    if(c.gt.b1.and.c.lt.b2) b1 = 0.5D+00*(a1+c)
    a2 = b1
    krule = 2
    call dqc25c(f,a1,b1,c,area1,error1,krule,nev)
    neval = neval+nev
    call dqc25c(f,a2,b2,c,area2,error2,krule,nev)
    neval = neval+nev
!
!  improve previous approximations to integral
!  and error and test for accuracy.
!
    area12 = area1+area2
    erro12 = error1+error2
    errsum = errsum+erro12-errmax
    area = area+area12-rlist(maxerr)
    if( abs ( rlist(maxerr)-area12).lt.0.1D-04* abs ( area12) &
      .and.erro12.ge.0.99D+00*errmax.and.krule.eq.0) &
      iroff1 = iroff1+1
    if(last.gt.10.and.erro12.gt.errmax.and.krule.eq.0) &
      iroff2 = iroff2+1
    rlist(maxerr) = area1
    rlist(last) = area2
    errbnd =  max ( epsabs,epsrel* abs ( area))
    if(errsum.le.errbnd) go to 15
!
!  test for roundoff error and eventually set error flag.
!
    if(iroff1.ge.6.and.iroff2.gt.20) ier = 2
!
!  set error flag in the case that number of interval bisections exceeds limit.
!
    if(last.eq.limit) ier = 1
!
!  set error flag in the case of bad integrand behaviour
!  at a point of the integration range.
!
    if( max (  abs ( a1), abs ( b2)).le.(0.1D+01+0.1D+03*epmach) &
      *( abs ( a2)+0.1D+04*uflow)) ier = 3
!
!  append the newly-created intervals to the list.
!
   15   continue

    if ( error2 .le. error1 ) then
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
    else
      alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
    end if
!
!  call dqpsrt to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with nrmax-th largest error estimate (to be bisected next).
!
    call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)

    if(ier.ne.0.or.errsum.le.errbnd) go to 50

   40 continue
!
!  compute final result.
!
   50 continue

  result = 0.0D+00
  do k=1,last
    result = result+rlist(k)
  end do

  abserr = errsum
   70 if (aa.eq.b) result=-result

  return
end
subroutine dqawc ( f, a, b, c, epsabs, epsrel, result, abserr, neval, ier, &
  limit, lenw, last, iwork, work )

!*****************************************************************************80
!
!! DQAWC computes a Cauchy principal value.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  the routine calculates an approximation result to a
!      cauchy principal value i = integral of f*w over (a,b)
!      (w(x) = 1/((x-c), c.ne.a, c.ne.b), hopefully satisfying
!      following claim for accuracy
!      abs(i-result).le.max(epsabe,epsrel*abs(i)).
!
!  Parameters:
!
!   on entry
!      f      - real ( kind = 8 )
!               function subprogram defining the integrand
!               function f(x). the actual name for f needs to be
!               declared e x t e r n a l in the driver program.
!
!      a      - real ( kind = 8 )
!               under limit of integration
!
!      b      - real ( kind = 8 )
!               upper limit of integration
!
!      c      - parameter in the weight function, c.ne.a, c.ne.b.
!               if c = a or c = b, the routine will end with
!               ier = 6 .
!
!      epsabs - real ( kind = 8 )
!               absolute accuracy requested
!      epsrel - real ( kind = 8 )
!               relative accuracy requested
!               if  epsabs.le.0
!               and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!               the routine will end with ier = 6.
!
!   on return
!      result - real ( kind = 8 )
!               approximation to the integral
!
!      abserr - real ( kind = 8 )
!               estimate or the modulus of the absolute error,
!               which should equal or exceed abs(i-result)
!
!      neval  - integer ( kind = 4 )
!               number of integrand evaluations
!
!      ier    - integer ( kind = 4 )
!               ier = 0 normal and reliable termination of the
!                       routine. it is assumed that the requested
!                       accuracy has been achieved.
!               ier.gt.0 abnormal termination of the routine
!                       the estimates for integral and error are
!                       less reliable. it is assumed that the
!                       requested accuracy has not been achieved.
!      error messages
!               ier = 1 maximum number of subdivisions allowed
!                       has been achieved. one can allow more sub-
!                       divisions by increasing the value of limit
!                       (and taking the according dimension
!                       adjustments into account). however, if
!                       this yields no improvement it is advised
!                       to analyze the integrand in order to
!                       determine the integration difficulties.
!                       if the position of a local difficulty
!                       can be determined (e.g. singularity,
!                       discontinuity within the interval) one
!                       will probably gain from splitting up the
!                       interval at this point and calling
!                       appropriate integrators on the subranges.
!                   = 2 the occurrence of roundoff error is detec-
!                       ted, which prevents the requested
!                       tolerance from being achieved.
!                   = 3 extremely bad integrand behaviour occurs
!                       at some points of the integration
!                       interval.
!                   = 6 the input is invalid, because
!                       c = a or c = b or
!                       (epsabs.le.0 and
!                        epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
!                       or limit.lt.1 or lenw.lt.limit*4.
!                       result, abserr, neval, last are set to
!                       zero. exept when lenw or limit is invalid,
!                       iwork(1), work(limit*2+1) and
!                       work(limit*3+1) are set to zero, work(1)
!                       is set to a and work(limit+1) to b.
!
!   dimensioning parameters
!      limit - integer ( kind = 4 )
!              dimensioning parameter for iwork
!              limit determines the maximum number of subintervals
!              in the partition of the given integration interval
!              (a,b), limit.ge.1.
!              if limit.lt.1, the routine will end with ier = 6.
!
!     lenw   - integer ( kind = 4 )
!              dimensioning parameter for work
!              lenw must be at least limit*4.
!              if lenw.lt.limit*4, the routine will end with
!              ier = 6.
!
!      last  - integer ( kind = 4 )
!              on return, last equals the number of subintervals
!              produced in the subdivision process, which
!              determines the number of significant elements
!              actually in the work arrays.
!
!   work arrays
!      iwork - integer ( kind = 4 )
!              vector of dimension at least limit, the first k
!              elements of which contain pointers
!              to the error estimates over the subintervals,
!              such that work(limit*3+iwork(1)), ... ,
!              work(limit*3+iwork(k)) form a decreasing
!              sequence, with k = last if last.le.(limit/2+2),
!              and k = limit+1-last otherwise
!
!      work  - real ( kind = 8 )
!              vector of dimension at least lenw
!              on return
!              work(1), ..., work(last) contain the left
!               end points of the subintervals in the
!               partition of (a,b),
!              work(limit+1), ..., work(limit+last) contain
!               the right end points,
!              work(limit*2+1), ..., work(limit*2+last) contain
!               the integral approximations over the subintervals,
!              work(limit*3+1), ..., work(limit*3+last)
!               contain the error estimates.
!
  implicit none

  real ( kind = 8 ) a,abserr,b,c,epsabs,epsrel,f,result,work
  integer ( kind = 4 ) ier,iwork,last,lenw,limit,lvl,l1,l2,l3,neval
  dimension iwork(limit),work(lenw)

  external f
!
!  check validity of limit and lenw.
!
  ier = 6
  neval = 0
  last = 0
  result = 0.0D+00
  abserr = 0.0D+00
  if(limit.lt.1.or.lenw.lt.limit*4) go to 10
!
!  prepare call for dqawce.
!
  l1 = limit+1
  l2 = limit+l1
  l3 = limit+l2
  call dqawce(f,a,b,c,epsabs,epsrel,limit,result,abserr,neval,ier, &
    work(1),work(l1),work(l2),work(l3),iwork,last)
!
!  call error handler if necessary.
!
  lvl = 0
10    if(ier.eq.6) lvl = 1

  if(ier.ne.0) then
    call xerror('abnormal return from dqawc',26,ier,lvl)
  end if

  return
end
subroutine dqawfe(f,a,omega,integr,epsabs,limlst,limit,maxp1, &
     result,abserr,neval,ier,rslst,erlst,ierlst,lst,alist,blist, &
     rlist,elist,iord,nnlog,chebmo)

!*****************************************************************************80
!
!! DQAWFE computes Fourier integrals.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  the routine calculates an approximation result to a
!      given fourier integal
!      i = integral of f(x)*w(x) over (a,infinity)
!      where w(x)=cos(omega*x) or w(x)=sin(omega*x),
!      hopefully satisfying following claim for accuracy
!      abs(i-result).le.epsabs.
!
!  Parameters:
!
!   on entry
!      f      - real ( kind = 8 )
!               function subprogram defining the integrand
!               function f(x). the actual name for f needs to
!               be declared e x t e r n a l in the driver program.
!
!      a      - real ( kind = 8 )
!               lower limit of integration
!
!      omega  - real ( kind = 8 )
!               parameter in the weight function
!
!      integr - integer ( kind = 4 )
!               indicates which weight function is used
!               integr = 1      w(x) = cos(omega*x)
!               integr = 2      w(x) = sin(omega*x)
!               if integr.ne.1.and.integr.ne.2, the routine will
!               end with ier = 6.
!
!      epsabs - real ( kind = 8 )
!               absolute accuracy requested, epsabs.gt.0
!               if epsabs.le.0, the routine will end with ier = 6.
!
!      limlst - integer ( kind = 4 )
!               limlst gives an upper bound on the number of
!               cycles, limlst.ge.1.
!               if limlst.lt.3, the routine will end with ier = 6.
!
!      limit  - integer ( kind = 4 )
!               gives an upper bound on the number of subintervals
!               allowed in the partition of each cycle, limit.ge.1
!               each cycle, limit.ge.1.
!
!      maxp1  - integer ( kind = 4 )
!               gives an upper bound on the number of
!               chebyshev moments which can be stored, i.e.
!               for the intervals of lengths abs(b-a)*2**(-l),
!               l=0,1, ..., maxp1-2, maxp1.ge.1
!
!   on return
!      result - real ( kind = 8 )
!               approximation to the integral x
!
!      abserr - real ( kind = 8 )
!               estimate of the modulus of the absolute error,
!               which should equal or exceed abs(i-result)
!
!      neval  - integer ( kind = 4 )
!               number of integrand evaluations
!
!      ier    - ier = 0 normal and reliable termination of
!                       the routine. it is assumed that the
!                       requested accuracy has been achieved.
!               ier.gt.0 abnormal termination of the routine. the
!                       estimates for integral and error are less
!                       reliable. it is assumed that the requested
!                       accuracy has not been achieved.
!      error messages
!              if omega.ne.0
!               ier = 1 maximum number of  cycles  allowed
!                       has been achieved., i.e. of subintervals
!                       (a+(k-1)c,a+kc) where
!                       c = (2*int(abs(omega))+1)*pi/abs(omega),
!                       for k = 1, 2, ..., lst.
!                       one can allow more cycles by increasing
!                       the value of limlst (and taking the
!                       according dimension adjustments into
!                       account).
!                       examine the array iwork which contains
!                       the error flags on the cycles, in order to
!                       look for eventual local integration
!                       difficulties. if the position of a local
!                       difficulty can be determined (e.g.
!                       singularity, discontinuity within the
!                       interval) one will probably gain from
!                       splitting up the interval at this point
!                       and calling appropriate integrators on
!                       the subranges.
!                   = 4 the extrapolation table constructed for
!                       convergence acceleration of the series
!                       formed by the integral contributions over
!                       the cycles, does not converge to within
!                       the requested accuracy. as in the case of
!                       ier = 1, it is advised to examine the
!                       array iwork which contains the error
!                       flags on the cycles.
!                   = 6 the input is invalid because
!                       (integr.ne.1 and integr.ne.2) or
!                        epsabs.le.0 or limlst.lt.3.
!                        result, abserr, neval, lst are set
!                        to zero.
!                   = 7 bad integrand behaviour occurs within one
!                       or more of the cycles. location and type
!                       of the difficulty involved can be
!                       determined from the vector ierlst. here
!                       lst is the number of cycles actually
!                       needed (see below).
!                       ierlst(k) = 1 the maximum number of
!                                     subdivisions (= limit) has
!                                     been achieved on the k th
!                                     cycle.
!                                 = 2 occurrence of roundoff error
!                                     is detected and prevents the
!                                     tolerance imposed on the
!                                     k th cycle, from being
!                                     achieved.
!                                 = 3 extremely bad integrand
!                                     behaviour occurs at some
!                                     points of the k th cycle.
!                                 = 4 the integration procedure
!                                     over the k th cycle does
!                                     not converge (to within the
!                                     required accuracy) due to
!                                     roundoff in the
!                                     extrapolation procedure
!                                     invoked on this cycle. it
!                                     is assumed that the result
!                                     on this interval is the
!                                     best which can be obtained.
!                                 = 5 the integral over the k th
!                                     cycle is probably divergent
!                                     or slowly convergent. it
!                                     must be noted that
!                                     divergence can occur with
!                                     any other value of
!                                     ierlst(k).
!              if omega = 0 and integr = 1,
!              the integral is calculated by means of dqagie
!              and ier = ierlst(1) (with meaning as described
!              for ierlst(k), k = 1).
!
!      rslst  - real ( kind = 8 )
!               vector of dimension at least limlst
!               rslst(k) contains the integral contribution
!               over the interval (a+(k-1)c,a+kc) where
!               c = (2*int(abs(omega))+1)*pi/abs(omega),
!               k = 1, 2, ..., lst.
!               note that, if omega = 0, rslst(1) contains
!               the value of the integral over (a,infinity).
!
!      erlst  - real ( kind = 8 )
!               vector of dimension at least limlst
!               erlst(k) contains the error estimate corresponding
!               with rslst(k).
!
!      ierlst - integer ( kind = 4 )
!               vector of dimension at least limlst
!               ierlst(k) contains the error flag corresponding
!               with rslst(k). for the meaning of the local error
!               flags see description of output parameter ier.
!
!      lst    - integer ( kind = 4 )
!               number of subintervals needed for the integration
!               if omega = 0 then lst is set to 1.
!
!      alist, blist, rlist, elist - real ( kind = 8 )
!               vector of dimension at least limit,
!
!      iord, nnlog - integer ( kind = 4 )
!               vector of dimension at least limit, providing
!               space for the quantities needed in the subdivision
!               process of each cycle
!
!      chebmo - real ( kind = 8 )
!               array of dimension at least (maxp1,25), providing
!               space for the chebyshev moments needed within the
!               cycles
!
!  Local Parameters:
!
!      the dimension of  psum  is determined by the value of
!      limexp in routine dqelg (psum must be of dimension
!      (limexp+2) at least).
!
!     c1, c2    - end points of subinterval (of length cycle)
!     cycle     - (2*int(abs(omega))+1)*pi/abs(omega)
!     psum      - vector of dimension at least (limexp+2)
!                 (see routine dqelg)
!                 psum contains the part of the epsilon table
!                 which is still needed for further computations.
!                 each element of psum is a partial sum of the
!                 series which should sum to the value of the
!                 integral.
!     errsum    - sum of error estimates over the subintervals,
!                 calculated cumulatively
!     epsa      - absolute tolerance requested over current
!                 subinterval
!     chebmo    - array containing the modified chebyshev
!                 moments (see also routine dqc25f)
!
  implicit none

  real ( kind = 8 ) a,abseps,abserr,alist,blist,chebmo,correc,cycle, &
    c1,c2,dl,dla,drl,elist,erlst,ep,eps,epsa, &
    epsabs,errsum,f,fact,omega,p,pi,p1,psum,reseps,result,res3la, &
    rlist,rslst,uflow
  integer ( kind = 4 ) ier,ierlst,integr,iord,ktmin,l,last,lst,limit
  integer ( kind = 4 ) limlst
  integer ( kind = 4 ) ll
  integer ( kind = 4 ) maxp1,momcom,nev,neval,nnlog,nres,numrl2

  dimension alist(limit),blist(limit),chebmo(maxp1,25),elist(limit), &
    erlst(limlst),ierlst(limlst),iord(limit),nnlog(limit),psum(52), &
    res3la(3),rlist(limit),rslst(limlst)

  external f

  data p / 0.9D+00 /
  data pi / 3.14159265358979323846264338327950D+00 /
!
!  test on validity of parameters
!
  result = 0.0D+00
  abserr = 0.0D+00
  neval = 0
  lst = 0
  ier = 0

  if((integr.ne.1.and.integr.ne.2).or.epsabs.le.0.0D+00.or. &
    limlst.lt.3) then
    ier = 6
    return
  end if

  if(omega.ne.0.0D+00) go to 10
!
!  integration by dqagie if omega is zero
!
  if(integr.eq.1) then
    call dqagie(f,0.0D+00,1,epsabs,0.0D+00,limit, &
      result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)
  end if

  rslst(1) = result
  erlst(1) = abserr
  ierlst(1) = ier
  lst = 1
  go to 999
!
!  initializations
!
   10 l =  abs ( omega)
  dl = 2*l+1
  cycle = dl*pi/ abs ( omega)
  ier = 0
  ktmin = 0
  neval = 0
  numrl2 = 0
  nres = 0
  c1 = a
  c2 = cycle+a
  p1 = 0.1D+01-p
  uflow = tiny ( uflow )
  eps = epsabs
  if(epsabs.gt.uflow/p1) eps = epsabs*p1
  ep = eps
  fact = 0.1D+01
  correc = 0.0D+00
  abserr = 0.0D+00
  errsum = 0.0D+00
!
!  main do-loop
!
  do lst = 1,limlst
!
!  integrate over current subinterval.
!
    dla = lst
    epsa = eps*fact
    call dqawoe(f,c1,c2,omega,integr,epsa,0.0D+00,limit,lst,maxp1, &
    rslst(lst),erlst(lst),nev,ierlst(lst),last,alist,blist,rlist, &
    elist,iord,nnlog,momcom,chebmo)
    neval = neval+nev
    fact = fact*p
    errsum = errsum+erlst(lst)
    drl = 0.5D+02* abs ( rslst(lst))
!
!  test on accuracy with partial sum
!
    if((errsum+drl).le.epsabs.and.lst.ge.6) go to 80
    correc =  max ( correc,erlst(lst))
    if(ierlst(lst).ne.0) eps =  max ( ep,correc*p1)
    if(ierlst(lst).ne.0) ier = 7
    if(ier.eq.7.and.(errsum+drl).le.correc*0.1D+02.and. &
    lst.gt.5) go to 80
    numrl2 = numrl2+1
    if(lst.gt.1) go to 20
    psum(1) = rslst(1)
    go to 40
   20   psum(numrl2) = psum(ll)+rslst(lst)
    if(lst.eq.2) go to 40
!
!  test on maximum number of subintervals
!
    if(lst.eq.limlst) ier = 1
!
!  perform new extrapolation
!
    call dqelg(numrl2,psum,reseps,abseps,res3la,nres)
!
!  test whether extrapolated result is influenced by roundoff
!
    ktmin = ktmin+1
    if(ktmin.ge.15.and.abserr.le.0.1D-02*(errsum+drl)) ier = 4
    if(abseps.gt.abserr.and.lst.ne.3) go to 30
    abserr = abseps
    result = reseps
    ktmin = 0
!
!  if ier is not 0, check whether direct result (partial sum)
!  or extrapolated result yields the best integral
!  approximation
!
    if((abserr+0.1D+02*correc).le.epsabs.or. &
    (abserr.le.epsabs.and.0.1D+02*correc.ge.epsabs)) go to 60
   30   if(ier.ne.0.and.ier.ne.7) go to 60
   40   ll = numrl2
    c1 = c2
    c2 = c2+cycle
  end do
!
!  set final result and error estimate
!
   60 abserr = abserr+0.1D+02*correc
  if(ier.eq.0) go to 999
  if(result.ne.0.0D+00.and.psum(numrl2).ne.0.0D+00) go to 70
  if(abserr.gt.errsum) go to 80
  if(psum(numrl2).eq.0.0D+00) go to 999
   70 if(abserr/ abs ( result).gt.(errsum+drl)/ abs ( psum(numrl2))) &
    go to 80
  if(ier.ge.1.and.ier.ne.7) abserr = abserr+drl
  go to 999
   80 result = psum(numrl2)
  abserr = errsum+drl
  999 continue

  return
end
subroutine dqawf ( f, a, omega, integr, epsabs, result, abserr, neval, ier, &
  limlst, lst, leniw, maxp1, lenw, iwork, work )

!*****************************************************************************80
!
!! DQAWF computes Fourier integrals over the interval [ A, +Infinity ).
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  the routine calculates an approximation result to a given
!      fourier integral i=integral of f(x)*w(x) over (a,infinity)
!      where w(x) = cos(omega*x) or w(x) = sin(omega*x).
!      hopefully satisfying following claim for accuracy
!      abs(i-result).le.epsabs.
!
!  Parameters:
!
!   on entry
!      f      - real ( kind = 8 )
!               function subprogram defining the integrand
!               function f(x). the actual name for f needs to be
!               declared e x t e r n a l in the driver program.
!
!      a      - real ( kind = 8 )
!               lower limit of integration
!
!      omega  - real ( kind = 8 )
!               parameter in the integrand weight function
!
!      integr - integer ( kind = 4 )
!               indicates which of the weight functions is used
!               integr = 1      w(x) = cos(omega*x)
!               integr = 2      w(x) = sin(omega*x)
!               if integr.ne.1.and.integr.ne.2, the routine
!               will end with ier = 6.
!
!      epsabs - real ( kind = 8 )
!               absolute accuracy requested, epsabs.gt.0.
!               if epsabs.le.0, the routine will end with ier = 6.
!
!   on return
!      result - real ( kind = 8 )
!               approximation to the integral
!
!      abserr - real ( kind = 8 )
!               estimate of the modulus of the absolute error,
!               which should equal or exceed abs(i-result)
!
!      neval  - integer ( kind = 4 )
!               number of integrand evaluations
!
!      ier    - integer ( kind = 4 )
!               ier = 0 normal and reliable termination of the
!                       routine. it is assumed that the requested
!                       accuracy has been achieved.
!               ier.gt.0 abnormal termination of the routine.
!                       the estimates for integral and error are
!                       less reliable. it is assumed that the
!                       requested accuracy has not been achieved.
!      error messages
!              if omega.ne.0
!               ier = 1 maximum number of cycles allowed
!                       has been achieved, i.e. of subintervals
!                       (a+(k-1)c,a+kc) where
!                       c = (2*int(abs(omega))+1)*pi/abs(omega),
!                       for k = 1, 2, ..., lst.
!                       one can allow more cycles by increasing
!                       the value of limlst (and taking the
!                       according dimension adjustments into
!                       account). examine the array iwork which
!                       contains the error flags on the cycles, in
!                       order to look for eventual local
!                       integration difficulties.
!                       if the position of a local difficulty
!                       can be determined (e.g. singularity,
!                       discontinuity within the interval) one
!                       will probably gain from splitting up the
!                       interval at this point and calling
!                       appropriate integrators on the subranges.
!                   = 4 the extrapolation table constructed for
!                       convergence accelaration of the series
!                       formed by the integral contributions over
!                       the cycles, does not converge to within
!                       the requested accuracy.
!                       as in the case of ier = 1, it is advised
!                       to examine the array iwork which contains
!                       the error flags on the cycles.
!                   = 6 the input is invalid because
!                       (integr.ne.1 and integr.ne.2) or
!                        epsabs.le.0 or limlst.lt.1 or
!                        leniw.lt.(limlst+2) or maxp1.lt.1 or
!                        lenw.lt.(leniw*2+maxp1*25).
!                        result, abserr, neval, lst are set to
!                        zero.
!                   = 7 bad integrand behaviour occurs within
!                       one or more of the cycles. location and
!                       type of the difficulty involved can be
!                       determined from the first lst elements of
!                       vector iwork.  here lst is the number of
!                       cycles actually needed (see below).
!                       iwork(k) = 1 the maximum number of
!                                    subdivisions (=(leniw-limlst)
!                                    /2) has been achieved on the
!                                    k th cycle.
!                                = 2 occurrence of roundoff error
!                                    is detected and prevents the
!                                    tolerance imposed on the k th
!                                    cycle, from being achieved
!                                    on this cycle.
!                                = 3 extremely bad integrand
!                                    behaviour occurs at some
!                                    points of the k th cycle.
!                                = 4 the integration procedure
!                                    over the k th cycle does
!                                    not converge (to within the
!                                    required accuracy) due to
!                                    roundoff in the extrapolation
!                                    procedure invoked on this
!                                    cycle. it is assumed that the
!                                    result on this interval is
!                                    the best which can be
!                                    obtained.
!                                = 5 the integral over the k th
!                                    cycle is probably divergent
!                                    or slowly convergent. it must
!                                    be noted that divergence can
!                                    occur with any other value of
!                                    iwork(k).
!              if omega = 0 and integr = 1,
!              the integral is calculated by means of dqagie,
!              and ier = iwork(1) (with meaning as described
!              for iwork(k),k = 1).
!
!   dimensioning parameters
!      limlst - integer ( kind = 4 )
!               limlst gives an upper bound on the number of
!               cycles, limlst.ge.3.
!               if limlst.lt.3, the routine will end with ier = 6.
!
!      lst    - integer ( kind = 4 )
!               on return, lst indicates the number of cycles
!               actually needed for the integration.
!               if omega = 0, then lst is set to 1.
!
!      leniw  - integer ( kind = 4 )
!               dimensioning parameter for iwork. on entry,
!               (leniw-limlst)/2 equals the maximum number of
!               subintervals allowed in the partition of each
!               cycle, leniw.ge.(limlst+2).
!               if leniw.lt.(limlst+2), the routine will end with
!               ier = 6.
!
!      maxp1  - integer ( kind = 4 )
!               maxp1 gives an upper bound on the number of
!               chebyshev moments which can be stored, i.e. for
!               the intervals of lengths abs(b-a)*2**(-l),
!               l = 0,1, ..., maxp1-2, maxp1.ge.1.
!               if maxp1.lt.1, the routine will end with ier = 6.
!      lenw   - integer ( kind = 4 )
!               dimensioning parameter for work
!               lenw must be at least leniw*2+maxp1*25.
!               if lenw.lt.(leniw*2+maxp1*25), the routine will
!               end with ier = 6.
!
!   work arrays
!      iwork  - integer ( kind = 4 )
!               vector of dimension at least leniw
!               on return, iwork(k) for k = 1, 2, ..., lst
!               contain the error flags on the cycles.
!
!      work   - real ( kind = 8 )
!               vector of dimension at least lenw
!               on return,
!               work(1), ..., work(lst) contain the integral
!                approximations over the cycles,
!               work(limlst+1), ..., work(limlst+lst) contain
!                the error extimates over the cycles.
!               further elements of work have no specific
!               meaning for the user.
!
  implicit none

  integer ( kind = 4 ) leniw
  integer ( kind = 4 ) lenw

  real ( kind = 8 ) a
  real ( kind = 8 ) abserr
  real ( kind = 8 ) epsabs
  real ( kind = 8 ), external :: f
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) integr
  integer ( kind = 4 ) iwork(leniw)
  integer ( kind = 4 ) last
  integer ( kind = 4 ) limit
  integer ( kind = 4 ) limlst
  integer ( kind = 4 ) ll2
  integer ( kind = 4 ) lst
  integer ( kind = 4 ) lvl
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) l3
  integer ( kind = 4 ) l4
  integer ( kind = 4 ) l5
  integer ( kind = 4 ) l6
  integer ( kind = 4 ) maxp1
  integer ( kind = 4 ) neval
  real ( kind = 8 ) omega
  real ( kind = 8 ) result
  real ( kind = 8 ) work(lenw)
!
!  check validity of limlst, leniw, maxp1 and lenw.
!
  ier = 6
  neval = 0
  last = 0
  result = 0.0D+00
  abserr = 0.0D+00
  if(limlst.lt.3.or.leniw.lt.(limlst+2).or.maxp1.lt.1.or.lenw.lt. &
     (leniw*2+maxp1*25)) go to 10
!
!  prepare call for dqawfe
!
  limit = (leniw-limlst)/2
  l1 = limlst+1
  l2 = limlst+l1
  l3 = limit+l2
  l4 = limit+l3
  l5 = limit+l4
  l6 = limit+l5
  ll2 = limit+l1
  call dqawfe(f,a,omega,integr,epsabs,limlst,limit,maxp1,result, &
    abserr,neval,ier,work(1),work(l1),iwork(1),lst,work(l2), &
    work(l3),work(l4),work(l5),iwork(l1),iwork(ll2),work(l6))
!
!  call error handler if necessary
!
  lvl = 0
10  continue

  if(ier.eq.6) lvl = 1
  if(ier.ne.0) call xerror('abnormal return from dqawf',26,ier,lvl)

  return
end
subroutine dqawoe ( f, a, b, omega, integr, epsabs, epsrel, limit, icall, &
  maxp1, result, abserr, neval, ier, last, alist, blist, rlist, elist, iord, &
  nnlog, momcom, chebmo )

!*****************************************************************************80
!
!! DQAWOE computes the integrals of oscillatory integrands.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  the routine calculates an approximation result to a given
!      definite integral
!      i = integral of f(x)*w(x) over (a,b)
!      where w(x) = cos(omega*x) or w(x)=sin(omega*x),
!      hopefully satisfying following claim for accuracy
!      abs(i-result).le.max(epsabs,epsrel*abs(i)).
!
!  Parameters:
!
!   on entry
!      f      - real ( kind = 8 )
!               function subprogram defining the integrand
!               function f(x). the actual name for f needs to be
!               declared e x t e r n a l in the driver program.
!
!      a      - real ( kind = 8 )
!               lower limit of integration
!
!      b      - real ( kind = 8 )
!               upper limit of integration
!
!      omega  - real ( kind = 8 )
!               parameter in the integrand weight function
!
!      integr - integer ( kind = 4 )
!               indicates which of the weight functions is to be
!               used
!               integr = 1      w(x) = cos(omega*x)
!               integr = 2      w(x) = sin(omega*x)
!               if integr.ne.1 and integr.ne.2, the routine
!               will end with ier = 6.
!
!      epsabs - real ( kind = 8 )
!               absolute accuracy requested
!      epsrel - real ( kind = 8 )
!               relative accuracy requested
!               if  epsabs.le.0
!               and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!               the routine will end with ier = 6.
!
!      limit  - integer ( kind = 4 )
!               gives an upper bound on the number of subdivisions
!               in the partition of (a,b), limit.ge.1.
!
!      icall  - integer ( kind = 4 )
!               if dqawoe is to be used only once, icall must
!               be set to 1.  assume that during this call, the
!               chebyshev moments (for clenshaw-curtis integration
!               of degree 24) have been computed for intervals of
!               lenghts (abs(b-a))*2**(-l), l=0,1,2,...momcom-1.
!               if icall.gt.1 this means that dqawoe has been
!               called twice or more on intervals of the same
!               length abs(b-a). the chebyshev moments already
!               computed are then re-used in subsequent calls.
!               if icall.lt.1, the routine will end with ier = 6.
!
!      maxp1  - integer ( kind = 4 )
!               gives an upper bound on the number of chebyshev
!               moments which can be stored, i.e. for the
!               intervals of lenghts abs(b-a)*2**(-l),
!               l=0,1, ..., maxp1-2, maxp1.ge.1.
!               if maxp1.lt.1, the routine will end with ier = 6.
!
!   on return
!      result - real ( kind = 8 )
!               approximation to the integral
!
!      abserr - real ( kind = 8 )
!               estimate of the modulus of the absolute error,
!               which should equal or exceed abs(i-result)
!
!      neval  - integer ( kind = 4 )
!               number of integrand evaluations
!
!      ier    - integer ( kind = 4 )
!               ier = 0 normal and reliable termination of the
!                       routine. it is assumed that the
!                       requested accuracy has been achieved.
!             - ier.gt.0 abnormal termination of the routine.
!                       the estimates for integral and error are
!                       less reliable. it is assumed that the
!                       requested accuracy has not been achieved.
!      error messages
!               ier = 1 maximum number of subdivisions allowed
!                       has been achieved. one can allow more
!                       subdivisions by increasing the value of
!                       limit (and taking according dimension
!                       adjustments into account). however, if
!                       this yields no improvement it is advised
!                       to analyze the integrand, in order to
!                       determine the integration difficulties.
!                       if the position of a local difficulty can
!                       be determined (e.g. singularity,
!                       discontinuity within the interval) one
!                       will probably gain from splitting up the
!                       interval at this point and calling the
!                       integrator on the subranges. if possible,
!                       an appropriate special-purpose integrator
!                       should be used which is designed for
!                       handling the type of difficulty involved.
!                   = 2 the occurrence of roundoff error is
!                       detected, which prevents the requested
!                       tolerance from being achieved.
!                       the error may be under-estimated.
!                   = 3 extremely bad integrand behaviour occurs
!                       at some points of the integration
!                       interval.
!                   = 4 the algorithm does not converge.
!                       roundoff error is detected in the
!                       extrapolation table.
!                       it is presumed that the requested
!                       tolerance cannot be achieved due to
!                       roundoff in the extrapolation table,
!                       and that the returned result is the
!                       best which can be obtained.
!                   = 5 the integral is probably divergent, or
!                       slowly convergent. it must be noted that
!                       divergence can occur with any other value
!                       of ier.gt.0.
!                   = 6 the input is invalid, because
!                       (epsabs.le.0 and
!                        epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
!                       or (integr.ne.1 and integr.ne.2) or
!                       icall.lt.1 or maxp1.lt.1.
!                       result, abserr, neval, last, rlist(1),
!                       elist(1), iord(1) and nnlog(1) are set
!                       to zero. alist(1) and blist(1) are set
!                       to a and b respectively.
!
!      last  -  integer ( kind = 4 )
!               on return, last equals the number of
!               subintervals produces in the subdivision
!               process, which determines the number of
!               significant elements actually in the
!               work arrays.
!      alist  - real ( kind = 8 )
!               vector of dimension at least limit, the first
!                last  elements of which are the left
!               end points of the subintervals in the partition
!               of the given integration range (a,b)
!
!      blist  - real ( kind = 8 )
!               vector of dimension at least limit, the first
!                last  elements of which are the right
!               end points of the subintervals in the partition
!               of the given integration range (a,b)
!
!      rlist  - real ( kind = 8 )
!               vector of dimension at least limit, the first
!                last  elements of which are the integral
!               approximations on the subintervals
!
!      elist  - real ( kind = 8 )
!               vector of dimension at least limit, the first
!                last  elements of which are the moduli of the
!               absolute error estimates on the subintervals
!
!      iord   - integer ( kind = 4 )
!               vector of dimension at least limit, the first k
!               elements of which are pointers to the error
!               estimates over the subintervals,
!               such that elist(iord(1)), ...,
!               elist(iord(k)) form a decreasing sequence, with
!               k = last if last.le.(limit/2+2), and
!               k = limit+1-last otherwise.
!
!      nnlog  - integer ( kind = 4 )
!               vector of dimension at least limit, containing the
!               subdivision levels of the subintervals, i.e.
!               iwork(i) = l means that the subinterval
!               numbered i is of length abs(b-a)*2**(1-l)
!
!   on entry and return
!      momcom - integer ( kind = 4 )
!               indicating that the chebyshev moments
!               have been computed for intervals of lengths
!               (abs(b-a))*2**(-l), l=0,1,2, ..., momcom-1,
!               momcom.lt.maxp1
!
!      chebmo - real ( kind = 8 )
!               array of dimension (maxp1,25) containing the
!               chebyshev moments
!
!  Local Parameters:
!
!      the dimension of rlist2 is determined by  the value of
!      limexp in routine dqelg (rlist2 should be of
!      dimension (limexp+2) at least).
!
!      list of major variables
!
!     alist     - list of left end points of all subintervals
!                 considered up to now
!     blist     - list of right end points of all subintervals
!                 considered up to now
!     rlist(i)  - approximation to the integral over
!                 (alist(i),blist(i))
!     rlist2    - array of dimension at least limexp+2
!                 containing the part of the epsilon table
!                 which is still needed for further computations
!     elist(i)  - error estimate applying to rlist(i)
!     maxerr    - pointer to the interval with largest
!                 error estimate
!     errmax    - elist(maxerr)
!     erlast    - error on the interval currently subdivided
!     area      - sum of the integrals over the subintervals
!     errsum    - sum of the errors over the subintervals
!     errbnd    - requested accuracy max(epsabs,epsrel*
!                 abs(result))
!     *****1    - variable for the left subinterval
!     *****2    - variable for the right subinterval
!     last      - index for subdivision
!     nres      - number of calls to the extrapolation routine
!     numrl2    - number of elements in rlist2. if an appropriate
!                 approximation to the compounded integral has
!                 been obtained it is put in rlist2(numrl2) after
!                 numrl2 has been increased by one
!     small     - length of the smallest interval considered
!                 up to now, multiplied by 1.5
!     erlarg    - sum of the errors over the intervals larger
!                 than the smallest interval considered up to now
!     extrap    - logical variable denoting that the routine is
!                 attempting to perform extrapolation, i.e. before
!                 subdividing the smallest interval we try to
!                 decrease the value of erlarg
!     noext     - logical variable denoting that extrapolation
!                 is no longer allowed (true  value)
!
!      machine dependent constants
!
!     epmach is the largest relative spacing.
!     uflow is the smallest positive magnitude.
!     oflow is the largest positive magnitude.
!
  implicit none

  real ( kind = 8 ) a,abseps,abserr,alist,area,area1,area12,area2,a1, &
    a2,b,blist,b1,b2,chebmo,correc,defab1,defab2,defabs, &
    domega,dres,elist,epmach,epsabs,epsrel,erlarg,erlast, &
    errbnd,errmax,error1,erro12,error2,errsum,ertest,f,oflow, &
    omega,resabs,reseps,result,res3la,rlist,rlist2,small,uflow,width
  integer ( kind = 4 ) icall,id,ier,ierro,integr,iord,iroff1,iroff2
  integer ( kind = 4 ) iroff3
  integer ( kind = 4 ) jupbnd
  integer ( kind = 4 ) k,ksgn,ktmin,last,limit,maxerr,maxp1,momcom,nev,neval, &
    nnlog,nres,nrmax,nrmom,numrl2
  logical extrap,noext,extall

  dimension alist(limit),blist(limit),rlist(limit),elist(limit), &
    iord(limit),rlist2(52),res3la(3),chebmo(maxp1,25),nnlog(limit)
  external f

  epmach = epsilon ( epmach )
!
!  test on validity of parameters
!
  ier = 0
  neval = 0
  last = 0
  result = 0.0D+00
  abserr = 0.0D+00
  alist(1) = a
  blist(1) = b
  rlist(1) = 0.0D+00
  elist(1) = 0.0D+00
  iord(1) = 0
  nnlog(1) = 0
  if((integr.ne.1.and.integr.ne.2).or.(epsabs.le.0.0D+00.and. &
    epsrel.lt. max ( 0.5D+02*epmach,0.5D-28)).or.icall.lt.1.or. &
    maxp1.lt.1) ier = 6
  if(ier.eq.6) go to 999
!
!  first approximation to the integral
!
  domega =  abs ( omega)
  nrmom = 0
  if (icall.gt.1) go to 5
  momcom = 0
    5 call dqc25f(f,a,b,domega,integr,nrmom,maxp1,0,result,abserr, &
    neval,defabs,resabs,momcom,chebmo)
!
!  test on accuracy.
!
  dres =  abs ( result)
  errbnd =  max ( epsabs,epsrel*dres)
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1
  if(abserr.le.0.1D+03*epmach*defabs.and.abserr.gt.errbnd) ier = 2
  if(limit.eq.1) ier = 1
  if(ier.ne.0.or.abserr.le.errbnd) go to 200
!
!  initializations
!
  uflow = tiny ( uflow )
  oflow = huge ( oflow )
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  abserr = oflow
  nrmax = 1
  extrap = .false.
  noext = .false.
  ierro = 0
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0
  ktmin = 0
  small =  abs ( b-a)*0.75D+00
  nres = 0
  numrl2 = 0
  extall = .false.
  if(0.5D+00* abs ( b-a)*domega.gt.0.2D+01) go to 10
  numrl2 = 1
  extall = .true.
  rlist2(1) = result
   10 if(0.25D+00* abs ( b-a)*domega.le.0.2D+01) extall = .true.
  ksgn = -1
  if(dres.ge.(0.1D+01-0.5D+02*epmach)*defabs) ksgn = 1
!
!  main do-loop
!
  do 140 last = 2,limit
!
!  bisect the subinterval with the nrmax-th largest error estimate.
!
    nrmom = nnlog(maxerr)+1
    a1 = alist(maxerr)
    b1 = 0.5D+00*(alist(maxerr)+blist(maxerr))
    a2 = b1
    b2 = blist(maxerr)
    erlast = errmax
    call dqc25f(f,a1,b1,domega,integr,nrmom,maxp1,0, &
    area1,error1,nev,resabs,defab1,momcom,chebmo)
    neval = neval+nev
    call dqc25f(f,a2,b2,domega,integr,nrmom,maxp1,1, &
    area2,error2,nev,resabs,defab2,momcom,chebmo)
    neval = neval+nev
!
!  improve previous approximations to integral
!  and error and test for accuracy.
!
    area12 = area1+area2
    erro12 = error1+error2
    errsum = errsum+erro12-errmax
    area = area+area12-rlist(maxerr)
    if(defab1.eq.error1.or.defab2.eq.error2) go to 25
    if( abs ( rlist(maxerr)-area12).gt.0.1D-04* abs ( area12) &
    .or.erro12.lt.0.99D+00*errmax) go to 20
    if(extrap) iroff2 = iroff2+1
    if(.not.extrap) iroff1 = iroff1+1
   20   if(last.gt.10.and.erro12.gt.errmax) iroff3 = iroff3+1
   25   rlist(maxerr) = area1
    rlist(last) = area2
    nnlog(maxerr) = nrmom
    nnlog(last) = nrmom
    errbnd =  max ( epsabs,epsrel* abs ( area))
!
!  test for roundoff error and eventually set error flag.
!
    if(iroff1+iroff2.ge.10.or.iroff3.ge.20) ier = 2
    if(iroff2.ge.5) ierro = 3
!
!  set error flag in the case that the number of
!  subintervals equals limit.
!
    if(last.eq.limit) ier = 1
!
!  set error flag in the case of bad integrand behaviour
!  at a point of the integration range.
!
    if( max (  abs ( a1), abs ( b2)).le.(0.1D+01+0.1D+03*epmach) &
    *( abs ( a2)+0.1D+04*uflow)) ier = 4
!
!  append the newly-created intervals to the list.
!
    if(error2.gt.error1) go to 30
    alist(last) = a2
    blist(maxerr) = b1
    blist(last) = b2
    elist(maxerr) = error1
    elist(last) = error2
    go to 40
   30   alist(maxerr) = a2
    alist(last) = a1
    blist(last) = b1
    rlist(maxerr) = area2
    rlist(last) = area1
    elist(maxerr) = error2
    elist(last) = error1
!
!  call dqpsrt to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with nrmax-th largest error estimate (to bisected next).
!
   40   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
  if(errsum.le.errbnd) go to 170
  if(ier.ne.0) go to 150
    if(last.eq.2.and.extall) go to 120
    if(noext) go to 140
    if(.not.extall) go to 50
    erlarg = erlarg-erlast
    if( abs ( b1-a1).gt.small) erlarg = erlarg+erro12
    if(extrap) go to 70
!
!  test whether the interval to be bisected next is the
!  smallest interval.
!
   50   width =  abs ( blist(maxerr)-alist(maxerr))
    if(width.gt.small) go to 140
    if(extall) go to 60
!
!  test whether we can start with the extrapolation procedure
!  (we do this if we integrate over the next interval with
!  use of a gauss-kronrod rule - see routine dqc25f).
!
    small = small*0.5D+00
    if(0.25D+00*width*domega.gt.0.2D+01) go to 140
    extall = .true.
    go to 130
   60   extrap = .true.
    nrmax = 2
   70   if(ierro.eq.3.or.erlarg.le.ertest) go to 90
!
!  the smallest interval has the largest error.
!  before bisecting decrease the sum of the errors over
!  the larger intervals (erlarg) and perform extrapolation.
!
    jupbnd = last
    if (last.gt.(limit/2+2)) jupbnd = limit+3-last
    id = nrmax
    do k = id,jupbnd
      maxerr = iord(nrmax)
      errmax = elist(maxerr)
      if( abs ( blist(maxerr)-alist(maxerr)).gt.small) go to 140
      nrmax = nrmax+1
    end do
!
!  perform extrapolation.
!
   90   numrl2 = numrl2+1
    rlist2(numrl2) = area
    if(numrl2.lt.3) go to 110
    call dqelg(numrl2,rlist2,reseps,abseps,res3la,nres)
    ktmin = ktmin+1
    if(ktmin.gt.5.and.abserr.lt.0.1D-02*errsum) ier = 5
    if(abseps.ge.abserr) go to 100
    ktmin = 0
    abserr = abseps
    result = reseps
    correc = erlarg
    ertest =  max ( epsabs,epsrel* abs ( reseps))
    if(abserr.le.ertest) go to 150
!
!  prepare bisection of the smallest interval.
!
  100   if(numrl2.eq.1) noext = .true.
    if(ier.eq.5) go to 150
  110   maxerr = iord(1)
    errmax = elist(maxerr)
    nrmax = 1
    extrap = .false.
    small = small*0.5D+00
    erlarg = errsum
    go to 140
  120   small = small*0.5D+00
    numrl2 = numrl2+1
    rlist2(numrl2) = area
  130   ertest = errbnd
    erlarg = errsum
  140 continue
!
!  set the final result.-
!
  150 if(abserr.eq.oflow.or.nres.eq.0) go to 170
  if(ier+ierro.eq.0) go to 165
  if(ierro.eq.3) abserr = abserr+correc
  if(ier.eq.0) ier = 3
  if(result.ne.0.0D+00.and.area.ne.0.0D+00) go to 160
  if(abserr.gt.errsum) go to 170
  if(area.eq.0.0D+00) go to 190
  go to 165
  160 if(abserr/ abs ( result).gt.errsum/ abs ( area)) go to 170
!
!  test on divergence.
!
  165 if(ksgn.eq.(-1).and. max (  abs ( result), abs ( area)).le. &
   defabs*0.1D-01) go to 190
  if(0.1D-01.gt.(result/area).or.(result/area).gt.0.1D+03 &
   .or.errsum.ge. abs ( area)) ier = 6
  go to 190
!
!  compute global integral sum.
!
  170 result = 0.0D+00
  do k=1,last
    result = result+rlist(k)
  end do
  abserr = errsum
  190 if (ier.gt.2) ier=ier-1
  200 if (integr.eq.2.and.omega.lt.0.0D+00) result=-result
  999 continue

  return
end
subroutine dqawo ( f, a, b, omega, integr, epsabs, epsrel, result, abserr, &
  neval, ier, leniw, maxp1, lenw, last, iwork, work )

!*****************************************************************************80
!
!! DQAWO computes the integrals of oscillatory integrands.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  the routine calculates an approximation result to a given
!      definite integral i=integral of f(x)*w(x) over (a,b)
!      where w(x) = cos(omega*x)
!      or w(x) = sin(omega*x),
!      hopefully satisfying following claim for accuracy
!      abs(i-result).le.max(epsabs,epsrel*abs(i)).
!
!  Parameters:
!
!   on entry
!      f      - real ( kind = 8 )
!               function subprogram defining the function
!               f(x).  the actual name for f needs to be
!               declared e x t e r n a l in the driver program.
!
!      a      - real ( kind = 8 )
!               lower limit of integration
!
!      b      - real ( kind = 8 )
!               upper limit of integration
!
!      omega  - real ( kind = 8 )
!               parameter in the integrand weight function
!
!      integr - integer ( kind = 4 )
!               indicates which of the weight functions is used
!               integr = 1      w(x) = cos(omega*x)
!               integr = 2      w(x) = sin(omega*x)
!               if integr.ne.1.and.integr.ne.2, the routine will
!               end with ier = 6.
!
!      epsabs - real ( kind = 8 )
!               absolute accuracy requested
!      epsrel - real ( kind = 8 )
!               relative accuracy requested
!               if epsabs.le.0 and
!               epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!               the routine will end with ier = 6.
!
!   on return
!      result - real ( kind = 8 )
!               approximation to the integral
!
!      abserr - real ( kind = 8 )
!               estimate of the modulus of the absolute error,
!               which should equal or exceed abs(i-result)
!
!      neval  - integer ( kind = 4 )
!               number of  integrand evaluations
!
!      ier    - integer ( kind = 4 )
!               ier = 0 normal and reliable termination of the
!                       routine. it is assumed that the requested
!                       accuracy has been achieved.
!             - ier.gt.0 abnormal termination of the routine.
!                       the estimates for integral and error are
!                       less reliable. it is assumed that the
!                       requested accuracy has not been achieved.
!      error messages
!               ier = 1 maximum number of subdivisions allowed
!                       (= leniw/2) has been achieved. one can
!                       allow more subdivisions by increasing the
!                       value of leniw (and taking the according
!                       dimension adjustments into account).
!                       however, if this yields no improvement it
!                       is advised to analyze the integrand in
!                       order to determine the integration
!                       difficulties. if the position of a local
!                       difficulty can be determined (e.g.
!                       singularity, discontinuity within the
!                       interval) one will probably gain from
!                       splitting up the interval at this point
!                       and calling the integrator on the
!                       subranges. if possible, an appropriate
!                       special-purpose integrator should be used
!                       which is designed for handling the type of
!                       difficulty involved.
!                   = 2 the occurrence of roundoff error is
!                       detected, which prevents the requested
!                       tolerance from being achieved.
!                       the error may be under-estimated.
!                   = 3 extremely bad integrand behaviour occurs
!                       at some interior points of the
!                       integration interval.
!                   = 4 the algorithm does not converge.
!                       roundoff error is detected in the
!                       extrapolation table. it is presumed that
!                       the requested tolerance cannot be achieved
!                       due to roundoff in the extrapolation
!                       table, and that the returned result is
!                       the best which can be obtained.
!                   = 5 the integral is probably divergent, or
!                       slowly convergent. it must be noted that
!                       divergence can occur with any other value
!                       of ier.
!                   = 6 the input is invalid, because
!                       (epsabs.le.0 and
!                        epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
!                       or (integr.ne.1 and integr.ne.2),
!                       or leniw.lt.2 or maxp1.lt.1 or
!                       lenw.lt.leniw*2+maxp1*25.
!                       result, abserr, neval, last are set to
!                       zero. except when leniw, maxp1 or lenw are
!                       invalid, work(limit*2+1), work(limit*3+1),
!                       iwork(1), iwork(limit+1) are set to zero,
!                       work(1) is set to a and work(limit+1) to
!                       b.
!
!   dimensioning parameters
!      leniw  - integer ( kind = 4 )
!               dimensioning parameter for iwork.
!               leniw/2 equals the maximum number of subintervals
!               allowed in the partition of the given integration
!               interval (a,b), leniw.ge.2.
!               if leniw.lt.2, the routine will end with ier = 6.
!
!      maxp1  - integer ( kind = 4 )
!               gives an upper bound on the number of chebyshev
!               moments which can be stored, i.e. for the
!               intervals of lengths abs(b-a)*2**(-l),
!               l=0,1, ..., maxp1-2, maxp1.ge.1
!               if maxp1.lt.1, the routine will end with ier = 6.
!
!      lenw   - integer ( kind = 4 )
!               dimensioning parameter for work
!               lenw must be at least leniw*2+maxp1*25.
!               if lenw.lt.(leniw*2+maxp1*25), the routine will
!               end with ier = 6.
!
!      last   - integer ( kind = 4 )
!               on return, last equals the number of subintervals
!               produced in the subdivision process, which
!               determines the number of significant elements
!               actually in the work arrays.
!
!   work arrays
!      iwork  - integer ( kind = 4 )
!               vector of dimension at least leniw
!               on return, the first k elements of which contain
!               pointers to the error estimates over the
!               subintervals, such that work(limit*3+iwork(1)), ..
!               work(limit*3+iwork(k)) form a decreasing
!               sequence, with limit = lenw/2 , and k = last
!               if last.le.(limit/2+2), and k = limit+1-last
!               otherwise.
!               furthermore, iwork(limit+1), ..., iwork(limit+
!               last) indicate the subdivision levels of the
!               subintervals, such that iwork(limit+i) = l means
!               that the subinterval numbered i is of length
!               abs(b-a)*2**(1-l).
!
!      work   - real ( kind = 8 )
!               vector of dimension at least lenw
!               on return
!               work(1), ..., work(last) contain the left
!                end points of the subintervals in the
!                partition of (a,b),
!               work(limit+1), ..., work(limit+last) contain
!                the right end points,
!               work(limit*2+1), ..., work(limit*2+last) contain
!                the integral approximations over the
!                subintervals,
!               work(limit*3+1), ..., work(limit*3+last)
!                contain the error estimates.
!               work(limit*4+1), ..., work(limit*4+maxp1*25)
!                provide space for storing the chebyshev moments.
!               note that limit = lenw/2.
!
  implicit none

  real ( kind = 8 ) a,abserr,b,epsabs,epsrel,f,omega,result,work
  integer ( kind = 4 ) ier,integr,iwork,last,limit,lenw,leniw,lvl,l
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) l3
  integer ( kind = 4 ) l4
  integer ( kind = 4 ) maxp1,momcom,neval
  dimension iwork(leniw),work(lenw)

  external f
!
!  check validity of leniw, maxp1 and lenw.
!
  ier = 6
  neval = 0
  last = 0
  result = 0.0D+00
  abserr = 0.0D+00
  if(leniw.lt.2.or.maxp1.lt.1.or.lenw.lt.(leniw*2+maxp1*25)) &
    go to 10
!
!  prepare call for dqawoe
!
  limit = leniw/2
  l1 = limit+1
  l2 = limit+l1
  l3 = limit+l2
  l4 = limit+l3
  call dqawoe(f,a,b,omega,integr,epsabs,epsrel,limit,1,maxp1,result, &
     abserr,neval,ier,last,work(1),work(l1),work(l2),work(l3), &
     iwork(1),iwork(l1),momcom,work(l4))
!
!  call error handler if necessary
!
  lvl = 0
10    if(ier.eq.6) lvl = 0
  if(ier.ne.0) call xerror('abnormal return from dqawo',26,ier,lvl)

  return
end
subroutine dqawse(f,a,b,alfa,beta,integr,epsabs,epsrel,limit, &
     result,abserr,neval,ier,alist,blist,rlist,elist,iord,last)

!*****************************************************************************80
!
!! DQAWSE estimates integrals with algebraico-logarithmic end singularities.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  the routine calculates an approximation result to a given
!      definite integral i = integral of f*w over (a,b),
!      (where w shows a singular behaviour at the end points,
!      see parameter integr).
!      hopefully satisfying following claim for accuracy
!      abs(i-result).le.max(epsabs,epsrel*abs(i)).
!
!  Parameters:
!
!   on entry
!      f      - real ( kind = 8 )
!               function subprogram defining the integrand
!               function f(x). the actual name for f needs to be
!               declared e x t e r n a l in the driver program.
!
!      a      - real ( kind = 8 )
!               lower limit of integration
!
!      b      - real ( kind = 8 )
!               upper limit of integration, b.gt.a
!               if b.le.a, the routine will end with ier = 6.
!
!      alfa   - real ( kind = 8 )
!               parameter in the weight function, alfa.gt.(-1)
!               if alfa.le.(-1), the routine will end with
!               ier = 6.
!
!      beta   - real ( kind = 8 )
!               parameter in the weight function, beta.gt.(-1)
!               if beta.le.(-1), the routine will end with
!               ier = 6.
!
!      integr - integer ( kind = 4 )
!               indicates which weight function is to be used
!               = 1  (x-a)**alfa*(b-x)**beta
!               = 2  (x-a)**alfa*(b-x)**beta*log(x-a)
!               = 3  (x-a)**alfa*(b-x)**beta*log(b-x)
!               = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
!               if integr.lt.1 or integr.gt.4, the routine
!               will end with ier = 6.
!
!      epsabs - real ( kind = 8 )
!               absolute accuracy requested
!      epsrel - real ( kind = 8 )
!               relative accuracy requested
!               if  epsabs.le.0
!               and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!               the routine will end with ier = 6.
!
!      limit  - integer ( kind = 4 )
!               gives an upper bound on the number of subintervals
!               in the partition of (a,b), limit.ge.2
!               if limit.lt.2, the routine will end with ier = 6.
!
!   on return
!      result - real ( kind = 8 )
!               approximation to the integral
!
!      abserr - real ( kind = 8 )
!               estimate of the modulus of the absolute error,
!               which should equal or exceed abs(i-result)
!
!      neval  - integer ( kind = 4 )
!               number of integrand evaluations
!
!      ier    - integer ( kind = 4 )
!               ier = 0 normal and reliable termination of the
!                       routine. it is assumed that the requested
!                       accuracy has been achieved.
!               ier.gt.0 abnormal termination of the routine
!                       the estimates for the integral and error
!                       are less reliable. it is assumed that the
!                       requested accuracy has not been achieved.
!      error messages
!                   = 1 maximum number of subdivisions allowed
!                       has been achieved. one can allow more
!                       subdivisions by increasing the value of
!                       limit. however, if this yields no
!                       improvement, it is advised to analyze the
!                       integrand in order to determine the
!                       integration difficulties which prevent the
!                       requested tolerance from being achieved.
!                       in case of a jump discontinuity or a local
!                       singularity of algebraico-logarithmic type
!                       at one or more interior points of the
!                       integration range, one should proceed by
!                       splitting up the interval at these
!                       points and calling the integrator on the
!                       subranges.
!                   = 2 the occurrence of roundoff error is
!                       detected, which prevents the requested
!                       tolerance from being achieved.
!                   = 3 extremely bad integrand behaviour occurs
!                       at some points of the integration
!                       interval.
!                   = 6 the input is invalid, because
!                       b.le.a or alfa.le.(-1) or beta.le.(-1), or
!                       integr.lt.1 or integr.gt.4, or
!                       (epsabs.le.0 and
!                        epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!                       or limit.lt.2.
!                       result, abserr, neval, rlist(1), elist(1),
!                       iord(1) and last are set to zero. alist(1)
!                       and blist(1) are set to a and b
!                       respectively.
!
!      alist  - real ( kind = 8 )
!               vector of dimension at least limit, the first
!                last  elements of which are the left
!               end points of the subintervals in the partition
!               of the given integration range (a,b)
!
!      blist  - real ( kind = 8 )
!               vector of dimension at least limit, the first
!                last  elements of which are the right
!               end points of the subintervals in the partition
!               of the given integration range (a,b)
!
!      rlist  - real ( kind = 8 )
!               vector of dimension at least limit,the first
!                last  elements of which are the integral
!               approximations on the subintervals
!
!      elist  - real ( kind = 8 )
!               vector of dimension at least limit, the first
!                last  elements of which are the moduli of the
!               absolute error estimates on the subintervals
!
!      iord   - integer ( kind = 4 )
!               vector of dimension at least limit, the first k
!               of which are pointers to the error
!               estimates over the subintervals, so that
!               elist(iord(1)), ..., elist(iord(k)) with k = last
!               if last.le.(limit/2+2), and k = limit+1-last
!               otherwise form a decreasing sequence
!
!      last   - integer ( kind = 4 )
!               number of subintervals actually produced in
!               the subdivision process
!
!  Local parameters:
!
!     alist     - list of left end points of all subintervals
!                 considered up to now
!     blist     - list of right end points of all subintervals
!                 considered up to now
!     rlist(i)  - approximation to the integral over
!                 (alist(i),blist(i))
!     elist(i)  - error estimate applying to rlist(i)
!     maxerr    - pointer to the interval with largest
!                 error estimate
!     errmax    - elist(maxerr)
!     area      - sum of the integrals over the subintervals
!     errsum    - sum of the errors over the subintervals
!     errbnd    - requested accuracy max(epsabs,epsrel*
!                 abs(result))
!     *****1    - variable for the left subinterval
!     *****2    - variable for the right subinterval
!     last      - index for subdivision
!
!      machine dependent constants
!
!     epmach is the largest relative spacing.
!     uflow is the smallest positive magnitude.
!
  implicit none

  real ( kind = 8 ) a,abserr,alfa,alist,area,area1,area12,area2,a1, &
    a2,b,beta,blist,b1,b2,centre,elist,epmach, &
    epsabs,epsrel,errbnd,errmax,error1,erro12,error2,errsum,f, &
    resas1,resas2,result,rg,rh,ri,rj,rlist,uflow
  integer ( kind = 4 ) ier,integr,iord,iroff1,iroff2,k,last,limit
  integer ( kind = 4 )maxerr
  integer ( kind = 4 ) nev
  integer ( kind = 4 ) neval,nrmax

  external f

  dimension alist(limit),blist(limit),rlist(limit),elist(limit), &
    iord(limit),ri(25),rj(25),rh(25),rg(25)

  epmach = epsilon ( epmach )
  uflow = tiny ( uflow )
!
!  test on validity of parameters
!
  neval = 0
  last = 0
  rlist(1) = 0.0D+00
  elist(1) = 0.0D+00
  iord(1) = 0
  result = 0.0D+00
  abserr = 0.0D+00

  if ( b.le.a .or. &
    (epsabs.eq.0.0D+00 .and. epsrel .lt. max ( 0.5D+02*epmach,0.5D-28) ) .or. &
    alfa .le. (-0.1D+01) .or. &
    beta .le. (-0.1D+01) .or. &
    integr.lt.1 .or. &
    integr.gt.4 .or. &
    limit.lt.2 ) then
    ier = 6
    return
  end if

  ier = 0
!
!  compute the modified chebyshev moments.
!
  call dqmomo(alfa,beta,ri,rj,rg,rh,integr)
!
!  integrate over the intervals (a,(a+b)/2) and ((a+b)/2,b).
!
  centre = 0.5D+00*(b+a)
  call dqc25s(f,a,b,a,centre,alfa,beta,ri,rj,rg,rh,area1, &
    error1,resas1,integr,nev)
  neval = nev
  call dqc25s(f,a,b,centre,b,alfa,beta,ri,rj,rg,rh,area2, &
    error2,resas2,integr,nev)
  last = 2
  neval = neval+nev
  result = area1+area2
  abserr = error1+error2
!
!  test on accuracy.
!
  errbnd = max ( epsabs,epsrel* abs ( result))
!
!  initialization
!
  if ( error2 .le. error1 ) then
    alist(1) = a
    alist(2) = centre
    blist(1) = centre
    blist(2) = b
    rlist(1) = area1
    rlist(2) = area2
    elist(1) = error1
    elist(2) = error2
  else
    alist(1) = centre
    alist(2) = a
    blist(1) = b
    blist(2) = centre
    rlist(1) = area2
    rlist(2) = area1
    elist(1) = error2
    elist(2) = error1
  end if

  iord(1) = 1
  iord(2) = 2
  if(limit.eq.2) ier = 1

  if(abserr.le.errbnd.or.ier.eq.1) then
    return
  end if

  errmax = elist(1)
  maxerr = 1
  nrmax = 1
  area = result
  errsum = abserr
  iroff1 = 0
  iroff2 = 0
!
!  main do-loop
!
  do 60 last = 3,limit
!
!  bisect the subinterval with largest error estimate.
!
    a1 = alist(maxerr)
    b1 = 0.5D+00*(alist(maxerr)+blist(maxerr))
    a2 = b1
    b2 = blist(maxerr)

    call dqc25s(f,a,b,a1,b1,alfa,beta,ri,rj,rg,rh,area1, &
    error1,resas1,integr,nev)
    neval = neval+nev
    call dqc25s(f,a,b,a2,b2,alfa,beta,ri,rj,rg,rh,area2, &
    error2,resas2,integr,nev)
    neval = neval+nev
!
!  improve previous approximations integral and error and test for accuracy.
!
    area12 = area1+area2
    erro12 = error1+error2
    errsum = errsum+erro12-errmax
    area = area+area12-rlist(maxerr)
    if(a.eq.a1.or.b.eq.b2) go to 30
    if(resas1.eq.error1.or.resas2.eq.error2) go to 30
!
!  test for roundoff error.
!
    if( abs ( rlist(maxerr)-area12).lt.0.1D-04* abs ( area12) &
    .and.erro12.ge.0.99D+00*errmax) iroff1 = iroff1+1
    if(last.gt.10.and.erro12.gt.errmax) iroff2 = iroff2+1
   30   rlist(maxerr) = area1
    rlist(last) = area2
!
!  test on accuracy.
!
    errbnd =  max ( epsabs,epsrel* abs ( area))
    if(errsum.le.errbnd) go to 35
!
!  set error flag in the case that the number of interval
!  bisections exceeds limit.
!
    if(last.eq.limit) ier = 1
!
!  set error flag in the case of roundoff error.
!
    if(iroff1.ge.6.or.iroff2.ge.20) ier = 2
!
!  set error flag in the case of bad integrand behaviour
!  at interior points of integration range.
!
    if( max (  abs ( a1), abs ( b2)).le.(0.1D+01+0.1D+03*epmach)* &
    ( abs ( a2)+0.1D+04*uflow)) ier = 3
!
!  append the newly-created intervals to the list.
!
   35   if(error2.gt.error1) go to 40
    alist(last) = a2
    blist(maxerr) = b1
    blist(last) = b2
    elist(maxerr) = error1
    elist(last) = error2
    go to 50

   40   alist(maxerr) = a2
    alist(last) = a1
    blist(last) = b1
    rlist(maxerr) = area2
    rlist(last) = area1
    elist(maxerr) = error2
    elist(last) = error1
!
!  call dqpsrt to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with largest error estimate (to be bisected next).
!
   50   call dqpsrt(limit,last,maxerr,errmax,elist,iord,nrmax)
    if (ier.ne.0.or.errsum.le.errbnd) go to 70
   60 continue
!
!  compute final result.
!
   70 continue

  result = 0.0D+00
  do k=1,last
    result = result+rlist(k)
  end do

  abserr = errsum
  999 continue

  return
end
subroutine dqaws ( f, a, b, alfa, beta, integr, epsabs, epsrel, result, &
  abserr, neval, ier, limit, lenw, last, iwork, work )

!*****************************************************************************80
!
!! DQAWS estimates integrals with algebraico-logarithmic endpoint singularities.
!
!  Modified:
!
!    12 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  the routine calculates an approximation result to a given
!      definite integral i = integral of f*w over (a,b),
!      (where w shows a singular behaviour at the end points
!      see parameter integr).
!      hopefully satisfying following claim for accuracy
!      abs(i-result).le.max(epsabs,epsrel*abs(i)).
!
!  Parameters:
!
!   on entry
!      f      - real ( kind = 8 )
!               function subprogram defining the integrand
!               function f(x). the actual name for f needs to be
!               declared e x t e r n a l in the driver program.
!
!      a      - real ( kind = 8 )
!               lower limit of integration
!
!      b      - real ( kind = 8 )
!               upper limit of integration, b.gt.a
!               if b.le.a, the routine will end with ier = 6.
!
!      alfa   - real ( kind = 8 )
!               parameter in the integrand function, alfa.gt.(-1)
!               if alfa.le.(-1), the routine will end with
!               ier = 6.
!
!      beta   - real ( kind = 8 )
!               parameter in the integrand function, beta.gt.(-1)
!               if beta.le.(-1), the routine will end with
!               ier = 6.
!
!      integr - integer ( kind = 4 )
!               indicates which weight function is to be used
!               = 1  (x-a)**alfa*(b-x)**beta
!               = 2  (x-a)**alfa*(b-x)**beta*log(x-a)
!               = 3  (x-a)**alfa*(b-x)**beta*log(b-x)
!               = 4  (x-a)**alfa*(b-x)**beta*log(x-a)*log(b-x)
!               if integr.lt.1 or integr.gt.4, the routine
!               will end with ier = 6.
!
!      epsabs - real ( kind = 8 )
!               absolute accuracy requested
!      epsrel - real ( kind = 8 )
!               relative accuracy requested
!               if  epsabs.le.0
!               and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!               the routine will end with ier = 6.
!
!   on return
!      result - real ( kind = 8 )
!               approximation to the integral
!
!      abserr - real ( kind = 8 )
!               estimate of the modulus of the absolute error,
!               which should equal or exceed abs(i-result)
!
!      neval  - integer ( kind = 4 )
!               number of integrand evaluations
!
!      ier    - integer ( kind = 4 )
!               ier = 0 normal and reliable termination of the
!                       routine. it is assumed that the requested
!                       accuracy has been achieved.
!               ier.gt.0 abnormal termination of the routine
!                       the estimates for the integral and error
!                       are less reliable. it is assumed that the
!                       requested accuracy has not been achieved.
!      error messages
!               ier = 1 maximum number of subdivisions allowed
!                       has been achieved. one can allow more
!                       subdivisions by increasing the value of
!                       limit (and taking the according dimension
!                       adjustments into account). however, if
!                       this yields no improvement it is advised
!                       to analyze the integrand, in order to
!                       determine the integration difficulties
!                       which prevent the requested tolerance from
!                       being achieved. in case of a jump
!                       discontinuity or a local singularity
!                       of algebraico-logarithmic type at one or
!                       more interior points of the integration
!                       range, one should proceed by splitting up
!                       the interval at these points and calling
!                       the integrator on the subranges.
!                   = 2 the occurrence of roundoff error is
!                       detected, which prevents the requested
!                       tolerance from being achieved.
!                   = 3 extremely bad integrand behaviour occurs
!                       at some points of the integration
!                       interval.
!                   = 6 the input is invalid, because
!                       b.le.a or alfa.le.(-1) or beta.le.(-1) or
!                       or integr.lt.1 or integr.gt.4 or
!                       (epsabs.le.0 and
!                        epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
!                       or limit.lt.2 or lenw.lt.limit*4.
!                       result, abserr, neval, last are set to
!                       zero. except when lenw or limit is invalid
!                       iwork(1), work(limit*2+1) and
!                       work(limit*3+1) are set to zero, work(1)
!                       is set to a and work(limit+1) to b.
!
!   dimensioning parameters
!      limit  - integer ( kind = 4 )
!               dimensioning parameter for iwork
!               limit determines the maximum number of
!               subintervals in the partition of the given
!               integration interval (a,b), limit.ge.2.
!               if limit.lt.2, the routine will end with ier = 6.
!
!      lenw   - integer ( kind = 4 )
!               dimensioning parameter for work
!               lenw must be at least limit*4.
!               if lenw.lt.limit*4, the routine will end
!               with ier = 6.
!
!      last   - integer ( kind = 4 )
!               on return, last equals the number of
!               subintervals produced in the subdivision process,
!               which determines the significant number of
!               elements actually in the work arrays.
!
!   work arrays
!      iwork  - integer ( kind = 4 )
!               vector of dimension limit, the first k
!               elements of which contain pointers
!               to the error estimates over the subintervals,
!               such that work(limit*3+iwork(1)), ...,
!               work(limit*3+iwork(k)) form a decreasing
!               sequence with k = last if last.le.(limit/2+2),
!               and k = limit+1-last otherwise
!
!      work   - real ( kind = 8 )
!               vector of dimension lenw
!               on return
!               work(1), ..., work(last) contain the left
!                end points of the subintervals in the
!                partition of (a,b),
!               work(limit+1), ..., work(limit+last) contain
!                the right end points,
!               work(limit*2+1), ..., work(limit*2+last)
!                contain the integral approximations over
!                the subintervals,
!               work(limit*3+1), ..., work(limit*3+last)
!                contain the error estimates.
!
  implicit none

  real ( kind = 8 ) a,abserr,alfa,b,beta,epsabs,epsrel,f,result,work
  integer ( kind = 4 ) ier,integr,iwork,last,lenw,limit,lvl,l1,l2,l3
  integer ( kind = 4 ) neval
  dimension iwork(limit),work(lenw)

  external f
!
!  check validity of limit and lenw.
!
  ier = 6
  neval = 0
  last = 0
  result = 0.0D+00
  abserr = 0.0D+00
  if(limit.lt.2.or.lenw.lt.limit*4) go to 10
!
!  prepare call for dqawse.
!
  l1 = limit+1
  l2 = limit+l1
  l3 = limit+l2

  call dqawse(f,a,b,alfa,beta,integr,epsabs,epsrel,limit,result, &
    abserr,neval,ier,work(1),work(l1),work(l2),work(l3),iwork,last)
!
!  call error handler if necessary.
!
  lvl = 0
10    if(ier.eq.6) lvl = 1
  if(ier.ne.0) call xerror('abnormal return from dqaws',26,ier,lvl)

  return
end
subroutine dqc25c(f,a,b,c,result,abserr,krul,neval)

!*****************************************************************************80
!
!! DQC25C returns integration rules for Cauchy Principal Value integrals.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  to compute i = integral of f*w over (a,b) with
!      error estimate, where w(x) = 1/(x-c)
!
!  Parameters:
!
!     f      - real ( kind = 8 )
!              function subprogram defining the integrand function
!              f(x). the actual name for f needs to be declared
!              e x t e r n a l  in the driver program.
!
!     a      - real ( kind = 8 )
!              left end point of the integration interval
!
!     b      - real ( kind = 8 )
!              right end point of the integration interval, b.gt.a
!
!     c      - real ( kind = 8 )
!              parameter in the weight function
!
!     result - real ( kind = 8 )
!              approximation to the integral
!              result is computed by using a generalized
!              clenshaw-curtis method if c lies within ten percent
!              of the integration interval. in the other case the
!              15-point kronrod rule obtained by optimal addition
!              of abscissae to the 7-point gauss rule, is applied.
!
!     abserr - real ( kind = 8 )
!              estimate of the modulus of the absolute error,
!              which should equal or exceed abs(i-result)
!
!     krul   - integer ( kind = 4 )
!              key which is decreased by 1 if the 15-point
!              gauss-kronrod scheme has been used
!
!     neval  - integer ( kind = 4 )
!              number of integrand evaluations
!
!  Local Parameters:
!
!     fval   - value of the function f at the points
!              cos(k*pi/24),  k = 0, ..., 24
!     cheb12 - chebyshev series expansion coefficients,
!              for the function f, of degree 12
!     cheb24 - chebyshev series expansion coefficients,
!              for the function f, of degree 24
!     res12  - approximation to the integral corresponding
!              to the use of cheb12
!     res24  - approximation to the integral corresponding
!              to the use of cheb24
!     dqwgtc - external function subprogram defining
!              the weight function
!     hlgth  - half-length of the interval
!     centr  - mid point of the interval
!
!     the vector x contains the values cos(k*pi/24),
!     k = 1, ..., 11, to be used for the chebyshev series
!     expansion of f
!
  implicit none

  real ( kind = 8 ) a,abserr,ak22,amom0,amom1,amom2,b,c,cc,centr, &
    cheb12,cheb24,dqwgtc,f,fval,hlgth,p2,p3,p4,resabs, &
    resasc,result,res12,res24,u,x
  integer ( kind = 4 ) i,isym,k,kp,krul,neval
  dimension x(11),fval(25),cheb12(13),cheb24(25)

  external f
  external dqwgtc

  data x(1) / 0.991444861373810411144557526928563d0 /
  data x(2) / 0.965925826289068286749743199728897d0 /
  data x(3) / 0.923879532511286756128183189396788d0 /
  data x(4) / 0.866025403784438646763723170752936d0 /
  data x(5) / 0.793353340291235164579776961501299d0 /
  data x(6) / 0.707106781186547524400844362104849d0 /
  data x(7) / 0.608761429008720639416097542898164d0 /
  data x(8) / 0.500000000000000000000000000000000d0 /
  data x(9) / 0.382683432365089771728459984030399d0 /
  data x(10) / 0.258819045102520762348898837624048d0 /
  data x(11) / 0.130526192220051591548406227895489d0 /
!
!  check the position of c.
!
  cc = (0.2D+01*c-b-a)/(b-a)
  if( abs ( cc).lt.0.11D+01) go to 10
!
!  apply the 15-point gauss-kronrod scheme.
!
  krul = krul-1
  call dqk15w(f,dqwgtc,c,p2,p3,p4,kp,a,b,result,abserr, &
    resabs,resasc)
  neval = 15
  if (resasc.eq.abserr) krul = krul+1
  go to 50
!
!  use the generalized clenshaw-curtis method.
!
   10 hlgth = 0.5D+00*(b-a)
  centr = 0.5D+00*(b+a)
  neval = 25
  fval(1) = 0.5D+00*f(hlgth+centr)
  fval(13) = f(centr)
  fval(25) = 0.5D+00*f(centr-hlgth)

  do i=2,12
    u = hlgth*x(i-1)
    isym = 26-i
    fval(i) = f(u+centr)
    fval(isym) = f(centr-u)
  end do
!
!  compute the chebyshev series expansion.
!
  call dqcheb(x,fval,cheb12,cheb24)
!
!  the modified chebyshev moments are computed by forward
!  recursion, using amom0 and amom1 as starting values.
!
  amom0 = log ( abs ( (0.1D+01-cc)/(0.1D+01+cc)))
  amom1 = 0.2D+01+cc*amom0
  res12 = cheb12(1)*amom0+cheb12(2)*amom1
  res24 = cheb24(1)*amom0+cheb24(2)*amom1

  do k=3,13
    amom2 = 0.2D+01*cc*amom1-amom0
    ak22 = (k-2)*(k-2)
    if((k/2)*2.eq.k) amom2 = amom2-0.4D+01/(ak22-0.1D+01)
    res12 = res12+cheb12(k)*amom2
    res24 = res24+cheb24(k)*amom2
    amom0 = amom1
    amom1 = amom2
  end do

  do k=14,25
    amom2 = 0.2D+01*cc*amom1-amom0
    ak22 = (k-2)*(k-2)
    if((k/2)*2.eq.k) amom2 = amom2-0.4D+01/(ak22-0.1D+01)
    res24 = res24+cheb24(k)*amom2
    amom0 = amom1
    amom1 = amom2
  end do

  result = res24
  abserr =  abs ( res24-res12)
   50 continue

  return
end
subroutine dqc25f(f,a,b,omega,integr,nrmom,maxp1,ksave,result, &
     abserr,neval,resabs,resasc,momcom,chebmo)

!*****************************************************************************80
!
!! DQC25F returns integration rules for functions with a COS or SIN factor.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  to compute the integral i=integral of f(x) over (a,b)
!      where w(x) = cos(omega*x) or w(x)=sin(omega*x) and to
!      compute j = integral of abs(f) over (a,b). for small value
!      of omega or small intervals (a,b) the 15-point gauss-kronro
!      rule is used. otherwise a generalized clenshaw-curtis
!      method is used.
!
!  Parameters:
!
!   on entry
!     f      - real ( kind = 8 )
!              function subprogram defining the integrand
!              function f(x). the actual name for f needs to
!              be declared e x t e r n a l in the calling program.
!
!     a      - real ( kind = 8 )
!              lower limit of integration
!
!     b      - real ( kind = 8 )
!              upper limit of integration
!
!     omega  - real ( kind = 8 )
!              parameter in the weight function
!
!     integr - integer ( kind = 4 )
!              indicates which weight function is to be used
!                 integr = 1   w(x) = cos(omega*x)
!                 integr = 2   w(x) = sin(omega*x)
!
!     nrmom  - integer ( kind = 4 )
!              the length of interval (a,b) is equal to the length
!              of the original integration interval divided by
!              2**nrmom (we suppose that the routine is used in an
!              adaptive integration process, otherwise set
!              nrmom = 0). nrmom must be zero at the first call.
!
!     maxp1  - integer ( kind = 4 )
!              gives an upper bound on the number of chebyshev
!              moments which can be stored, i.e. for the
!              intervals of lengths abs(bb-aa)*2**(-l),
!              l = 0,1,2, ..., maxp1-2.
!
!     ksave  - integer ( kind = 4 )
!              key which is one when the moments for the
!              current interval have been computed
!
!   on return
!     result - real ( kind = 8 )
!              approximation to the integral i
!
!     abserr - real ( kind = 8 )
!              estimate of the modulus of the absolute
!              error, which should equal or exceed abs(i-result)
!
!     neval  - integer ( kind = 4 )
!              number of integrand evaluations
!
!     resabs - real ( kind = 8 )
!              approximation to the integral j
!
!     resasc - real ( kind = 8 )
!              approximation to the integral of abs(f-i/(b-a))
!
!   on entry and return
!     momcom - integer ( kind = 4 )
!              for each interval length we need to compute the
!              chebyshev moments. momcom counts the number of
!              intervals for which these moments have already been
!              computed. if nrmom.lt.momcom or ksave = 1, the
!              chebyshev moments for the interval (a,b) have
!              already been computed and stored, otherwise we
!              compute them and we increase momcom.
!
!     chebmo - real ( kind = 8 )
!              array of dimension at least (maxp1,25) containing
!              the modified chebyshev moments for the first momcom
!              momcom interval lengths
!
!  Local Parameters:
!
!    the vector x contains the values cos(k*pi/24)
!    k = 1, ...,11, to be used for the chebyshev expansion of f
!
!     centr  - mid point of the integration interval
!     hlgth  - half-length of the integration interval
!     fval   - value of the function f at the points
!              (b-a)*0.5*cos(k*pi/12) + (b+a)*0.5, k = 0, ..., 24
!     cheb12 - coefficients of the chebyshev series expansion
!              of degree 12, for the function f, in the
!              interval (a,b)
!     cheb24 - coefficients of the chebyshev series expansion
!              of degree 24, for the function f, in the
!              interval (a,b)
!     resc12 - approximation to the integral of
!              cos(0.5*(b-a)*omega*x)*f(0.5*(b-a)*x+0.5*(b+a))
!              over (-1,+1), using the chebyshev series
!              expansion of degree 12
!     resc24 - approximation to the same integral, using the
!              chebyshev series expansion of degree 24
!     ress12 - the analogue of resc12 for the sine
!     ress24 - the analogue of resc24 for the sine
!
!
!     machine dependent constant
!
!     oflow is the largest positive magnitude.
!
  implicit none

  real ( kind = 8 ) a,abserr,ac,an,an2,as,asap,ass,b,centr,chebmo, &
    cheb12,cheb24,conc,cons,cospar,d,dqwgtf,d1, &
    d2,estc,ests,f,fval,hlgth,oflow,omega,parint,par2,par22, &
    p2,p3,p4,resabs,resasc,resc12,resc24,ress12,ress24,result, &
    sinpar,v,x
  integer ( kind = 4 ) i,iers,integr,isym,j,k,ksave,m,momcom,neval, maxp1,&
    noequ,noeq1,nrmom
  dimension chebmo(maxp1,25),cheb12(13),cheb24(25),d(25),d1(25), &
    d2(25),fval(25),v(28),x(11)

  external f,dqwgtf

  data x(1) / 0.991444861373810411144557526928563d0 /
  data x(2) / 0.965925826289068286749743199728897d0 /
  data x(3) / 0.923879532511286756128183189396788d0 /
  data x(4) / 0.866025403784438646763723170752936d0 /
  data x(5) / 0.793353340291235164579776961501299d0 /
  data x(6) / 0.707106781186547524400844362104849d0 /
  data x(7) / 0.608761429008720639416097542898164d0 /
  data x(8) / 0.500000000000000000000000000000000d0 /
  data x(9) / 0.382683432365089771728459984030399d0 /
  data x(10) / 0.258819045102520762348898837624048d0 /
  data x(11) / 0.130526192220051591548406227895489d0 /

  oflow = huge ( oflow )
  centr = 0.5D+00*(b+a)
  hlgth = 0.5D+00*(b-a)
  parint = omega*hlgth
!
!  compute the integral using the 15-point gauss-kronrod
!  formula if the value of the parameter in the integrand is small.
!
  if( abs ( parint).gt.0.2D+01) go to 10
  call dqk15w(f,dqwgtf,omega,p2,p3,p4,integr,a,b,result, &
    abserr,resabs,resasc)
  neval = 15
  go to 170
!
!  compute the integral using the generalized clenshaw-
!  curtis method.
!
   10 conc = hlgth*dcos(centr*omega)
  cons = hlgth*dsin(centr*omega)
  resasc = oflow
  neval = 25
!
!  check whether the chebyshev moments for this interval
!  have already been computed.
!
  if(nrmom.lt.momcom.or.ksave.eq.1) go to 120
!
!  compute a new set of chebyshev moments.
!
  m = momcom+1
  par2 = parint*parint
  par22 = par2+0.2D+01
  sinpar = dsin(parint)
  cospar = dcos(parint)
!
!  compute the chebyshev moments with respect to cosine.
!
  v(1) = 0.2D+01*sinpar/parint
  v(2) = (0.8D+01*cospar+(par2+par2-0.8D+01)*sinpar/parint)/par2
  v(3) = (0.32D+02*(par2-0.12D+02)*cospar+(0.2D+01* &
    ((par2-0.80D+02)*par2+0.192D+03)*sinpar)/parint)/(par2*par2)
  ac = 0.8D+01*cospar
  as = 0.24D+02*parint*sinpar
  if( abs ( parint).gt.0.24D+02) go to 30
!
!  compute the chebyshev moments as the solutions of a
!  boundary value problem with 1 initial value (v(3)) and 1
!  end value (computed using an asymptotic formula).
!
  noequ = 25
  noeq1 = noequ-1
  an = 0.6D+01

  do k = 1,noeq1
    an2 = an*an
    d(k) = -0.2D+01*(an2-0.4D+01)*(par22-an2-an2)
    d2(k) = (an-0.1D+01)*(an-0.2D+01)*par2
    d1(k+1) = (an+0.3D+01)*(an+0.4D+01)*par2
    v(k+3) = as-(an2-0.4D+01)*ac
    an = an+0.2D+01
  end do

  an2 = an*an
  d(noequ) = -0.2D+01*(an2-0.4D+01)*(par22-an2-an2)
  v(noequ+3) = as-(an2-0.4D+01)*ac
  v(4) = v(4)-0.56D+02*par2*v(3)
  ass = parint*sinpar
  asap = (((((0.210D+03*par2-0.1D+01)*cospar-(0.105D+03*par2 &
    -0.63D+02)*ass)/an2-(0.1D+01-0.15D+02*par2)*cospar &
    +0.15D+02*ass)/an2-cospar+0.3D+01*ass)/an2-cospar)/an2
  v(noequ+3) = v(noequ+3)-0.2D+01*asap*par2*(an-0.1D+01)* &
     (an-0.2D+01)
!
!  solve the tridiagonal system by means of gaussian
!  elimination with partial pivoting.
!
  call dgtsl(noequ,d1,d,d2,v(4),iers)
  go to 50
!
!  compute the chebyshev moments by means of forward recursion.
!
   30 an = 0.4D+01

  do i = 4,13
    an2 = an*an
    v(i) = ((an2-0.4D+01)*(0.2D+01*(par22-an2-an2)*v(i-1)-ac) &
    +as-par2*(an+0.1D+01)*(an+0.2D+01)*v(i-2))/ &
    (par2*(an-0.1D+01)*(an-0.2D+01))
    an = an+0.2D+01
  end do

   50 continue

  do j = 1,13
    chebmo(m,2*j-1) = v(j)
  end do
!
!  compute the chebyshev moments with respect to sine.
!
  v(1) = 0.2D+01*(sinpar-parint*cospar)/par2
  v(2) = (0.18D+02-0.48D+02/par2)*sinpar/par2 &
    +(-0.2D+01+0.48D+02/par2)*cospar/parint
  ac = -0.24D+02*parint*cospar
  as = -0.8D+01*sinpar
  if( abs ( parint).gt.0.24D+02) go to 80
!
!  compute the chebyshev moments as the solutions of a boundary
!  value problem with 1 initial value (v(2)) and 1 end value
!  (computed using an asymptotic formula).
!
  an = 0.5D+01

  do k = 1,noeq1
    an2 = an*an
    d(k) = -0.2D+01*(an2-0.4D+01)*(par22-an2-an2)
    d2(k) = (an-0.1D+01)*(an-0.2D+01)*par2
    d1(k+1) = (an+0.3D+01)*(an+0.4D+01)*par2
    v(k+2) = ac+(an2-0.4D+01)*as
    an = an+0.2D+01
  end do

  an2 = an*an
  d(noequ) = -0.2D+01*(an2-0.4D+01)*(par22-an2-an2)
  v(noequ+2) = ac+(an2-0.4D+01)*as
  v(3) = v(3)-0.42D+02*par2*v(2)
  ass = parint*cospar
  asap = (((((0.105D+03*par2-0.63D+02)*ass+(0.210D+03*par2 &
    -0.1D+01)*sinpar)/an2+(0.15D+02*par2-0.1D+01)*sinpar- &
    0.15D+02*ass)/an2-0.3D+01*ass-sinpar)/an2-sinpar)/an2
  v(noequ+2) = v(noequ+2)-0.2D+01*asap*par2*(an-0.1D+01) &
    *(an-0.2D+01)
!
!  solve the tridiagonal system by means of gaussian
!  elimination with partial pivoting.
!
  call dgtsl(noequ,d1,d,d2,v(3),iers)
  go to 100
!
!  compute the chebyshev moments by means of forward recursion.
!
   80 an = 0.3D+01

  do i = 3,12
    an2 = an*an
    v(i) = ((an2-0.4D+01)*(0.2D+01*(par22-an2-an2)*v(i-1)+as) &
    +ac-par2*(an+0.1D+01)*(an+0.2D+01)*v(i-2)) &
    /(par2*(an-0.1D+01)*(an-0.2D+01))
    an = an+0.2D+01
  end do

  100 continue

  do j = 1,12
    chebmo(m,2*j) = v(j)
  end do

  120 if (nrmom.lt.momcom) m = nrmom+1
   if (momcom.lt.(maxp1-1).and.nrmom.ge.momcom) momcom = momcom+1
!
!  compute the coefficients of the chebyshev expansions
!  of degrees 12 and 24 of the function f.
!
  fval(1) = 0.5D+00*f(centr+hlgth)
  fval(13) = f(centr)
  fval(25) = 0.5D+00*f(centr-hlgth)
  do i = 2,12
    isym = 26-i
    fval(i) = f(hlgth*x(i-1)+centr)
    fval(isym) = f(centr-hlgth*x(i-1))
  end do
  call dqcheb(x,fval,cheb12,cheb24)
!
!  compute the integral and error estimates.
!
  resc12 = cheb12(13)*chebmo(m,13)
  ress12 = 0.0D+00
  k = 11
  do j = 1,6
    resc12 = resc12+cheb12(k)*chebmo(m,k)
    ress12 = ress12+cheb12(k+1)*chebmo(m,k+1)
    k = k-2
  end do
  resc24 = cheb24(25)*chebmo(m,25)
  ress24 = 0.0D+00
  resabs =  abs ( cheb24(25))
  k = 23
  do j = 1,12
    resc24 = resc24+cheb24(k)*chebmo(m,k)
    ress24 = ress24+cheb24(k+1)*chebmo(m,k+1)
    resabs =  abs ( cheb24(k))+ abs ( cheb24(k+1))
    k = k-2
  end do
  estc =  abs ( resc24-resc12)
  ests =  abs ( ress24-ress12)
  resabs = resabs* abs ( hlgth)
  if(integr.eq.2) go to 160
  result = conc*resc24-cons*ress24
  abserr =  abs ( conc*estc)+ abs ( cons*ests)
  go to 170
  160 result = conc*ress24+cons*resc24
  abserr =  abs ( conc*ests)+ abs ( cons*estc)
  170 continue

  return
end
subroutine dqc25s(f,a,b,bl,br,alfa,beta,ri,rj,rg,rh,result, &
     abserr,resasc,integr,nev)

!*****************************************************************************80
!
!! DQC25S returns rules for algebraico-logarithmic end point singularities.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  to compute i = integral of f*w over (bl,br), with error
!      estimate, where the weight function w has a singular
!      behaviour of algebraico-logarithmic type at the points
!      a and/or b. (bl,br) is a part of (a,b).
!
!  Parameters:
!
!     f      - real ( kind = 8 )
!              function subprogram defining the integrand
!              f(x). the actual name for f needs to be declared
!              e x t e r n a l  in the driver program.
!
!     a      - real ( kind = 8 )
!              left end point of the original interval
!
!     b      - real ( kind = 8 )
!              right end point of the original interval, b.gt.a
!
!     bl     - real ( kind = 8 )
!              lower limit of integration, bl.ge.a
!
!     br     - real ( kind = 8 )
!              upper limit of integration, br.le.b
!
!     alfa   - real ( kind = 8 )
!              parameter in the weight function
!
!     beta   - real ( kind = 8 )
!              parameter in the weight function
!
!     ri,rj,rg,rh - real ( kind = 8 )
!              modified chebyshev moments for the application
!              of the generalized clenshaw-curtis
!              method (computed in routine dqmomo)
!
!     result - real ( kind = 8 )
!              approximation to the integral
!              result is computed by using a generalized
!              clenshaw-curtis method if b1 = a or br = b.
!              in all other cases the 15-point kronrod
!              rule is applied, obtained by optimal addition of
!              abscissae to the 7-point gauss rule.
!
!     abserr - real ( kind = 8 )
!              estimate of the modulus of the absolute error,
!              which should equal or exceed abs(i-result)
!
!     resasc - real ( kind = 8 )
!              approximation to the integral of abs(f*w-i/(b-a))
!
!     integr - integer ( kind = 4 )
!              which determines the weight function
!              = 1   w(x) = (x-a)**alfa*(b-x)**beta
!              = 2   w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)
!              = 3   w(x) = (x-a)**alfa*(b-x)**beta*log(b-x)
!              = 4   w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)*
!                           log(b-x)
!
!     nev    - integer ( kind = 4 )
!              number of integrand evaluations
!
!  Local Parameters:
!
!     the vector x contains the values cos(k*pi/24)
!     k = 1, ..., 11, to be used for the computation of the
!     chebyshev series expansion of f.
!
!     fval   - value of the function f at the points
!              (br-bl)*0.5*cos(k*pi/24)+(br+bl)*0.5
!              k = 0, ..., 24
!     cheb12 - coefficients of the chebyshev series expansion
!              of degree 12, for the function f, in the
!              interval (bl,br)
!     cheb24 - coefficients of the chebyshev series expansion
!              of degree 24, for the function f, in the
!              interval (bl,br)
!     res12  - approximation to the integral obtained from cheb12
!     res24  - approximation to the integral obtained from cheb24
!     dqwgts - external function subprogram defining
!              the four possible weight functions
!     hlgth  - half-length of the interval (bl,br)
!     centr  - mid point of the interval (bl,br)
!
  implicit none

  real ( kind = 8 ) a,abserr,alfa,b,beta,bl,br,centr,cheb12,cheb24, &
    dc,f,factor,fix,fval,hlgth,resabs,resasc,result,res12, &
    res24,rg,rh,ri,rj,u,dqwgts,x
  integer ( kind = 4 ) i,integr,isym,nev

  dimension cheb12(13),cheb24(25),fval(25),rg(25),rh(25),ri(25), &
    rj(25),x(11)

  external f,dqwgts

  data x(1) / 0.991444861373810411144557526928563d0 /
  data x(2) / 0.965925826289068286749743199728897d0 /
  data x(3) / 0.923879532511286756128183189396788d0 /
  data x(4) / 0.866025403784438646763723170752936d0 /
  data x(5) / 0.793353340291235164579776961501299d0 /
  data x(6) / 0.707106781186547524400844362104849d0 /
  data x(7) / 0.608761429008720639416097542898164d0 /
  data x(8) / 0.500000000000000000000000000000000d0 /
  data x(9) / 0.382683432365089771728459984030399d0 /
  data x(10) / 0.258819045102520762348898837624048d0 /
  data x(11) / 0.130526192220051591548406227895489d0 /

  nev = 25
  if(bl.eq.a.and.(alfa.ne.0.0D+00.or.integr.eq.2.or.integr.eq.4)) &
   go to 10
  if(br.eq.b.and.(beta.ne.0.0D+00.or.integr.eq.3.or.integr.eq.4)) &
   go to 140
!
!  if a.gt.bl and b.lt.br, apply the 15-point gauss-kronrod scheme.
!
!
  call dqk15w(f,dqwgts,a,b,alfa,beta,integr,bl,br, &
      result,abserr,resabs,resasc)
  nev = 15
  go to 270
!
!  this part of the program is executed only if a = bl.
!
!  compute the chebyshev series expansion of the
!  following function
!  f1 = (0.5*(b+b-br-a)-0.5*(br-a)*x)**beta
!         *f(0.5*(br-a)*x+0.5*(br+a))
!
   10 hlgth = 0.5D+00*(br-bl)
  centr = 0.5D+00*(br+bl)
  fix = b-centr
  fval(1) = 0.5D+00*f(hlgth+centr)*(fix-hlgth)**beta
  fval(13) = f(centr)*(fix**beta)
  fval(25) = 0.5D+00*f(centr-hlgth)*(fix+hlgth)**beta
  do i=2,12
    u = hlgth*x(i-1)
    isym = 26-i
    fval(i) = f(u+centr)*(fix-u)**beta
    fval(isym) = f(centr-u)*(fix+u)**beta
  end do

  factor = hlgth**(alfa+0.1D+01)
  result = 0.0D+00
  abserr = 0.0D+00
  res12 = 0.0D+00
  res24 = 0.0D+00
  if(integr.gt.2) go to 70
  call dqcheb(x,fval,cheb12,cheb24)
!
!  integr = 1  (or 2)
!
  do i=1,13
    res12 = res12+cheb12(i)*ri(i)
    res24 = res24+cheb24(i)*ri(i)
  end do

  do i=14,25
    res24 = res24+cheb24(i)*ri(i)
  end do

  if(integr.eq.1) go to 130
!
!  integr = 2
!
  dc = log (br-bl)
  result = res24*dc
  abserr =  abs ( (res24-res12)*dc)
  res12 = 0.0D+00
  res24 = 0.0D+00
  do i=1,13
    res12 = res12+cheb12(i)*rg(i)
    res24 = res12+cheb24(i)*rg(i)
  end do
  do i=14,25
    res24 = res24+cheb24(i)*rg(i)
  end do
  go to 130
!
!  compute the chebyshev series expansion of the
!  following function
!  f4 = f1*log(0.5*(b+b-br-a)-0.5*(br-a)*x)
!
   70 fval(1) = fval(1)* log (fix-hlgth)
  fval(13) = fval(13)* log (fix)
  fval(25) = fval(25)* log (fix+hlgth)
  do i=2,12
    u = hlgth*x(i-1)
    isym = 26-i
    fval(i) = fval(i)* log (fix-u)
    fval(isym) = fval(isym)* log (fix+u)
  end do
  call dqcheb(x,fval,cheb12,cheb24)
!
!  integr = 3  (or 4)
!
  do i=1,13
    res12 = res12+cheb12(i)*ri(i)
    res24 = res24+cheb24(i)*ri(i)
  end do

  do i=14,25
    res24 = res24+cheb24(i)*ri(i)
  end do
  if(integr.eq.3) go to 130
!
!  integr = 4
!
  dc = log (br-bl)
  result = res24*dc
  abserr =  abs ( (res24-res12)*dc)
  res12 = 0.0D+00
  res24 = 0.0D+00
  do i=1,13
    res12 = res12+cheb12(i)*rg(i)
    res24 = res24+cheb24(i)*rg(i)
  end do
  do i=14,25
    res24 = res24+cheb24(i)*rg(i)
  end do
  130 result = (result+res24)*factor
  abserr = (abserr+ abs ( res24-res12))*factor
  go to 270
!
!  this part of the program is executed only if b = br.
!
!  compute the chebyshev series expansion of the following function:
!
!    f2 = (0.5*(b+bl-a-a)+0.5*(b-bl)*x)**alfa*f(0.5*(b-bl)*x+0.5*(b+bl))
!
  140 hlgth = 0.5D+00*(br-bl)
  centr = 0.5D+00*(br+bl)
  fix = centr-a
  fval(1) = 0.5D+00*f(hlgth+centr)*(fix+hlgth)**alfa
  fval(13) = f(centr)*(fix**alfa)
  fval(25) = 0.5D+00*f(centr-hlgth)*(fix-hlgth)**alfa
  do i=2,12
    u = hlgth*x(i-1)
    isym = 26-i
    fval(i) = f(u+centr)*(fix+u)**alfa
    fval(isym) = f(centr-u)*(fix-u)**alfa
  end do
  factor = hlgth**(beta+0.1D+01)
  result = 0.0D+00
  abserr = 0.0D+00
  res12 = 0.0D+00
  res24 = 0.0D+00
  if(integr.eq.2.or.integr.eq.4) go to 200
!
!  integr = 1  (or 3)
!
  call dqcheb(x,fval,cheb12,cheb24)

  do i=1,13
    res12 = res12+cheb12(i)*rj(i)
    res24 = res24+cheb24(i)*rj(i)
  end do

  do i=14,25
    res24 = res24+cheb24(i)*rj(i)
  end do

  if(integr.eq.1) go to 260
!
! integr = 3
!
  dc = log (br-bl)
  result = res24*dc
  abserr =  abs ( (res24-res12)*dc)
  res12 = 0.0D+00
  res24 = 0.0D+00
  do i=1,13
    res12 = res12+cheb12(i)*rh(i)
    res24 = res24+cheb24(i)*rh(i)
  end do

  do i=14,25
    res24 = res24+cheb24(i)*rh(i)
  end do
  go to 260
!
!  compute the chebyshev series expansion of the
!  following function
!  f3 = f2*log(0.5*(b-bl)*x+0.5*(b+bl-a-a))
!
  200 fval(1) = fval(1)* log (hlgth+fix)
  fval(13) = fval(13)* log (fix)
  fval(25) = fval(25)* log (fix-hlgth)
  do i=2,12
    u = hlgth*x(i-1)
    isym = 26-i
    fval(i) = fval(i)* log (u+fix)
    fval(isym) = fval(isym)* log (fix-u)
  end do
  call dqcheb(x,fval,cheb12,cheb24)
!
!  integr = 2  (or 4)
!
  do i=1,13
    res12 = res12+cheb12(i)*rj(i)
    res24 = res24+cheb24(i)*rj(i)
  end do

  do i=14,25
    res24 = res24+cheb24(i)*rj(i)
  end do

  if(integr.eq.2) go to 260
  dc = log (br-bl)
  result = res24*dc
  abserr =  abs ( (res24-res12)*dc)
  res12 = 0.0D+00
  res24 = 0.0D+00
!
!  integr = 4
!
  do i=1,13
    res12 = res12+cheb12(i)*rh(i)
    res24 = res24+cheb24(i)*rh(i)
  end do

  do i=14,25
    res24 = res24+cheb24(i)*rh(i)
  end do

  260 result = (result+res24)*factor
  abserr = (abserr+ abs ( res24-res12))*factor
  270 return
end
subroutine dqcheb ( x, fval, cheb12, cheb24 )

!*****************************************************************************80
!
!! DQCHEB computes the Chebyshev series expansion.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  this routine computes the chebyshev series expansion
!      of degrees 12 and 24 of a function using a
!      fast fourier transform method
!      f(x) = sum(k=1,..,13) (cheb12(k)*t(k-1,x)),
!      f(x) = sum(k=1,..,25) (cheb24(k)*t(k-1,x)),
!      where t(k,x) is the chebyshev polynomial of degree k.
!
!  Parameters:
!
!    on entry
!     x      - real ( kind = 8 )
!              vector of dimension 11 containing the
!              values cos(k*pi/24), k = 1, ..., 11
!
!     fval   - real ( kind = 8 )
!              vector of dimension 25 containing the
!              function values at the points
!              (b+a+(b-a)*cos(k*pi/24))/2, k = 0, ...,24,
!              where (a,b) is the approximation interval.
!              fval(1) and fval(25) are divided by two
!              (these values are destroyed at output).
!
!    on return
!     cheb12 - real ( kind = 8 )
!              vector of dimension 13 containing the
!              chebyshev coefficients for degree 12
!
!     cheb24 - real ( kind = 8 )
!              vector of dimension 25 containing the
!              chebyshev coefficients for degree 24
!
  implicit none

  real ( kind = 8 ) alam,alam1,alam2,cheb12,cheb24,fval,part1,part2, &
    part3,v,x
  integer ( kind = 4 ) i,j

  dimension cheb12(13),cheb24(25),fval(25),v(12),x(11)

  do i=1,12
    j = 26-i
    v(i) = fval(i)-fval(j)
    fval(i) = fval(i)+fval(j)
  end do

  alam1 = v(1)-v(9)
  alam2 = x(6)*(v(3)-v(7)-v(11))
  cheb12(4) = alam1+alam2
  cheb12(10) = alam1-alam2
  alam1 = v(2)-v(8)-v(10)
  alam2 = v(4)-v(6)-v(12)
  alam = x(3)*alam1+x(9)*alam2
  cheb24(4) = cheb12(4)+alam
  cheb24(22) = cheb12(4)-alam
  alam = x(9)*alam1-x(3)*alam2
  cheb24(10) = cheb12(10)+alam
  cheb24(16) = cheb12(10)-alam
  part1 = x(4)*v(5)
  part2 = x(8)*v(9)
  part3 = x(6)*v(7)
  alam1 = v(1)+part1+part2
  alam2 = x(2)*v(3)+part3+x(10)*v(11)
  cheb12(2) = alam1+alam2
  cheb12(12) = alam1-alam2
  alam = x(1)*v(2)+x(3)*v(4)+x(5)*v(6)+x(7)*v(8) &
    +x(9)*v(10)+x(11)*v(12)
  cheb24(2) = cheb12(2)+alam
  cheb24(24) = cheb12(2)-alam
  alam = x(11)*v(2)-x(9)*v(4)+x(7)*v(6)-x(5)*v(8) &
    +x(3)*v(10)-x(1)*v(12)
  cheb24(12) = cheb12(12)+alam
  cheb24(14) = cheb12(12)-alam
  alam1 = v(1)-part1+part2
  alam2 = x(10)*v(3)-part3+x(2)*v(11)
  cheb12(6) = alam1+alam2
  cheb12(8) = alam1-alam2
  alam = x(5)*v(2)-x(9)*v(4)-x(1)*v(6) &
    -x(11)*v(8)+x(3)*v(10)+x(7)*v(12)
  cheb24(6) = cheb12(6)+alam
  cheb24(20) = cheb12(6)-alam
  alam = x(7)*v(2)-x(3)*v(4)-x(11)*v(6)+x(1)*v(8) &
    -x(9)*v(10)-x(5)*v(12)
  cheb24(8) = cheb12(8)+alam
  cheb24(18) = cheb12(8)-alam

  do i=1,6
    j = 14-i
    v(i) = fval(i)-fval(j)
    fval(i) = fval(i)+fval(j)
  end do

  alam1 = v(1)+x(8)*v(5)
  alam2 = x(4)*v(3)
  cheb12(3) = alam1+alam2
  cheb12(11) = alam1-alam2
  cheb12(7) = v(1)-v(5)
  alam = x(2)*v(2)+x(6)*v(4)+x(10)*v(6)
  cheb24(3) = cheb12(3)+alam
  cheb24(23) = cheb12(3)-alam
  alam = x(6)*(v(2)-v(4)-v(6))
  cheb24(7) = cheb12(7)+alam
  cheb24(19) = cheb12(7)-alam
  alam = x(10)*v(2)-x(6)*v(4)+x(2)*v(6)
  cheb24(11) = cheb12(11)+alam
  cheb24(15) = cheb12(11)-alam

  do i=1,3
    j = 8-i
    v(i) = fval(i)-fval(j)
    fval(i) = fval(i)+fval(j)
  end do

  cheb12(5) = v(1)+x(8)*v(3)
  cheb12(9) = fval(1)-x(8)*fval(3)
  alam = x(4)*v(2)
  cheb24(5) = cheb12(5)+alam
  cheb24(21) = cheb12(5)-alam
  alam = x(8)*fval(2)-fval(4)
  cheb24(9) = cheb12(9)+alam
  cheb24(17) = cheb12(9)-alam
  cheb12(1) = fval(1)+fval(3)
  alam = fval(2)+fval(4)
  cheb24(1) = cheb12(1)+alam
  cheb24(25) = cheb12(1)-alam
  cheb12(13) = v(1)-v(3)
  cheb24(13) = cheb12(13)
  alam = 0.1D+01/0.6D+01

  do i=2,12
    cheb12(i) = cheb12(i)*alam
  end do

  alam = 0.5D+00*alam
  cheb12(1) = cheb12(1)*alam
  cheb12(13) = cheb12(13)*alam

  do i=2,24
    cheb24(i) = cheb24(i)*alam
  end do

  cheb24(1) = 0.5D+00*alam*cheb24(1)
  cheb24(25) = 0.5D+00*alam*cheb24(25)

  return
end
subroutine dqelg ( n, epstab, result, abserr, res3la, nres )

!*****************************************************************************80
!
!! DQELG carries out the Epsilon extrapolation algorithm.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  the routine determines the limit of a given sequence of
!      approximations, by means of the epsilon algorithm of
!      p.wynn. an estimate of the absolute error is also given.
!      the condensed epsilon table is computed. only those
!      elements needed for the computation of the next diagonal
!      are preserved.
!
!  Parameters:
!
!        n      - integer ( kind = 4 )
!                 epstab(n) contains the new element in the
!                 first column of the epsilon table.
!
!        epstab - real ( kind = 8 )
!                 vector of dimension 52 containing the elements
!                 of the two lower diagonals of the triangular
!                 epsilon table. the elements are numbered
!                 starting at the right-hand corner of the
!                 triangle.
!
!        result - real ( kind = 8 )
!                 resulting approximation to the integral
!
!        abserr - real ( kind = 8 )
!                 estimate of the absolute error computed from
!                 result and the 3 previous results
!
!        res3la - real ( kind = 8 )
!                 vector of dimension 3 containing the last 3
!                 results
!
!        nres   - integer ( kind = 4 )
!                 number of calls to the routine
!                 (should be zero at first call)
!
!  Local Parameters:
!
!     e0     - the 4 elements on which the computation of a new
!     e1       element in the epsilon table is based
!     e2
!     e3                 e0
!                  e3    e1    new
!                        e2
!     newelm - number of elements to be computed in the new
!              diagonal
!     error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
!     result - the element in the new diagonal with least value
!              of error
!
!     machine dependent constants
!
!     epmach is the largest relative spacing.
!     oflow is the largest positive magnitude.
!     limexp is the maximum number of elements the epsilon
!     table can contain. if this number is reached, the upper
!     diagonal of the epsilon table is deleted.
!
  implicit none

  real ( kind = 8 ) abserr,delta1,delta2,delta3, &
    epmach,epsinf,epstab,error,err1,err2,err3,e0,e1,e1abs,e2,e3, &
    oflow,res,result,res3la,ss,tol1,tol2,tol3
  integer ( kind = 4 ) i,ib,ib2,ie,indx,k1,k2,k3,limexp,n,newelm
  integer ( kind = 4 ) nres
  integer ( kind = 4 ) num
  dimension epstab(52),res3la(3)

  epmach = epsilon ( epmach )
  oflow = huge ( oflow )
  nres = nres+1
  abserr = oflow
  result = epstab(n)
  if(n.lt.3) go to 100
  limexp = 50
  epstab(n+2) = epstab(n)
  newelm = (n-1)/2
  epstab(n) = oflow
  num = n
  k1 = n

  do 40 i = 1,newelm

    k2 = k1-1
    k3 = k1-2
    res = epstab(k1+2)
    e0 = epstab(k3)
    e1 = epstab(k2)
    e2 = res
    e1abs =  abs ( e1)
    delta2 = e2-e1
    err2 =  abs ( delta2)
    tol2 =  max (  abs ( e2),e1abs)*epmach
    delta3 = e1 - e0
    err3 =  abs ( delta3)
    tol3 =  max ( e1abs, abs ( e0))*epmach
    if(err2.gt.tol2.or.err3.gt.tol3) go to 10
!
!  if e0, e1 and e2 are equal to machine accuracy, convergence is assumed.
!
    result = res
    abserr = err2+err3
    go to 100
   10   e3 = epstab(k1)
    epstab(k1) = e1
    delta1 = e1-e3
    err1 =  abs ( delta1)
    tol1 =  max ( e1abs, abs ( e3))*epmach
!
!  if two elements are very close to each other, omit
!  a part of the table by adjusting the value of n
!
    if(err1.le.tol1.or.err2.le.tol2.or.err3.le.tol3) go to 20
    ss = 0.1D+01/delta1+0.1D+01/delta2-0.1D+01/delta3
    epsinf =  abs ( ss*e1)
!
!  test to detect irregular behaviour in the table, and
!  eventually omit a part of the table adjusting the value
!  of n.
!
    if(epsinf.gt.0.1D-03) go to 30
   20   n = i+i-1
    go to 50
!
!  compute a new element and eventually adjust
!  the value of result.
!
   30   res = e1+0.1D+01/ss
    epstab(k1) = res
    k1 = k1-2
    error = err2 + abs ( res-e2 ) + err3

    if ( error .le. abserr ) then
      abserr = error
      result = res
    end if

   40 continue
!
!  shift the table.
!
   50 if(n.eq.limexp) n = 2*(limexp/2)-1
  ib = 1
  if((num/2)*2.eq.num) ib = 2
  ie = newelm+1
  do i=1,ie
    ib2 = ib+2
    epstab(ib) = epstab(ib2)
    ib = ib2
  end do
  if(num.eq.n) go to 80
  indx = num-n+1
  do i = 1,n
    epstab(i)= epstab(indx)
    indx = indx+1
  end do
   80 if(nres.ge.4) go to 90
  res3la(nres) = result
  abserr = oflow
  go to 100
!
!  compute error estimate
!
   90 abserr =  abs ( result-res3la(3))+ abs ( result-res3la(2)) &
    + abs ( result-res3la(1))
  res3la(1) = res3la(2)
  res3la(2) = res3la(3)
  res3la(3) = result
  100 continue

  abserr =  max ( abserr, 0.5D+01*epmach* abs ( result))

  return
end
subroutine dqk15(f,a,b,result,abserr,resabs,resasc)

!*****************************************************************************80
!
!! DQK15 carries out a 15 point Gauss-Kronrod quadrature rule.
!
!     the abscissae and weights are given for the interval (-1,1).
!     because of symmetry only the positive abscissae and their
!     corresponding weights are given.
!
!     xgk    - abscissae of the 15-point kronrod rule
!              xgk(2), xgk(4), ...  abscissae of the 7-point
!              gauss rule
!              xgk(1), xgk(3), ...  abscissae which are optimally
!              added to the 7-point gauss rule
!
!     wgk    - weights of the 15-point kronrod rule
!
!     wg     - weights of the 7-point gauss rule
!
!
!   gauss quadrature weights and kronron quadrature abscissae and weights
!   as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
!   bell labs, nov. 1981.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  to compute i = integral of f over (a,b), with error
!                     estimate
!                 j = integral of abs(f) over (a,b)
!  Parameters:
!
!      on entry
!        f      - real ( kind = 8 )
!                 function subprogram defining the integrand
!                 function f(x). the actual name for f needs to be
!                 declared e x t e r n a l in the calling program.
!
!        a      - real ( kind = 8 )
!                 lower limit of integration
!
!        b      - real ( kind = 8 )
!                 upper limit of integration
!
!      on return
!        result - real ( kind = 8 )
!                 approximation to the integral i
!                 result is computed by applying the 15-point
!                 kronrod rule (resk) obtained by optimal addition
!                 of abscissae to the7-point gauss rule(resg).
!
!        abserr - real ( kind = 8 )
!                 estimate of the modulus of the absolute error,
!                 which should not exceed abs(i-result)
!
!        resabs - real ( kind = 8 )
!                 approximation to the integral j
!
!        resasc - real ( kind = 8 )
!                 approximation to the integral of abs(f-i/(b-a))
!                 over (a,b)
!
!  Local Parameters:
!
!     centr  - mid point of the interval
!     hlgth  - half-length of the interval
!     absc   - abscissa
!     fval*  - function value
!     resg   - result of the 7-point gauss formula
!     resk   - result of the 15-point kronrod formula
!     reskh  - approximation to the mean value of f over (a,b),
!              i.e. to i/(b-a)
!
!     machine dependent constants
!
!     epmach is the largest relative spacing.
!     uflow is the smallest positive magnitude.
!
  implicit none

  real ( kind = 8 ) a,absc,abserr,b,centr,dhlgth, &
    epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc, &
    resg,resk,reskh,result,uflow,wg,wgk,xgk
  integer ( kind = 4 ) j,jtw,jtwm1
  external f
  dimension fv1(7),fv2(7),wg(4),wgk(8),xgk(8)

  data wg  (  1) / 0.129484966168869693270611432679082d0 /
  data wg  (  2) / 0.279705391489276667901467771423780d0 /
  data wg  (  3) / 0.381830050505118944950369775488975d0 /
  data wg  (  4) / 0.417959183673469387755102040816327d0 /

  data xgk (  1) / 0.991455371120812639206854697526329d0 /
  data xgk (  2) / 0.949107912342758524526189684047851d0 /
  data xgk (  3) / 0.864864423359769072789712788640926d0 /
  data xgk (  4) / 0.741531185599394439863864773280788d0 /
  data xgk (  5) / 0.586087235467691130294144838258730d0 /
  data xgk (  6) / 0.405845151377397166906606412076961d0 /
  data xgk (  7) / 0.207784955007898467600689403773245d0 /
  data xgk (  8) / 0.000000000000000000000000000000000d0 /

  data wgk (  1) / 0.022935322010529224963732008058970d0 /
  data wgk (  2) / 0.063092092629978553290700663189204d0 /
  data wgk (  3) / 0.104790010322250183839876322541518d0 /
  data wgk (  4) / 0.140653259715525918745189590510238d0 /
  data wgk (  5) / 0.169004726639267902826583426598550d0 /
  data wgk (  6) / 0.190350578064785409913256402421014d0 /
  data wgk (  7) / 0.204432940075298892414161999234649d0 /
  data wgk (  8) / 0.209482141084727828012999174891714d0 /

  epmach = epsilon ( epmach )
  uflow = tiny ( uflow )
  centr = 0.5D+00*(a+b)
  hlgth = 0.5D+00*(b-a)
  dhlgth =  abs ( hlgth)
!
!  compute the 15-point kronrod approximation to
!  the integral, and estimate the absolute error.
!
  fc = f(centr)
  resg = fc*wg(4)
  resk = fc*wgk(8)
  resabs =  abs ( resk)

  do j=1,3
    jtw = j*2
    absc = hlgth*xgk(jtw)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*( abs ( fval1)+ abs ( fval2))
  end do

  do j = 1,4
    jtwm1 = j*2-1
    absc = hlgth*xgk(jtwm1)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*( abs ( fval1)+ abs ( fval2))
  end do

  reskh = resk*0.5D+00
  resasc = wgk(8)* abs ( fc-reskh)
  do j=1,7
    resasc = resasc+wgk(j)*( abs ( fv1(j)-reskh)+ abs ( fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr =  abs ( (resk-resg)*hlgth)
  if(resasc.ne.0.0D+00.and.abserr.ne.0.0D+00) &
    abserr = resasc* min (0.1D+01,(0.2D+03*abserr/resasc)**1.5D+00)
  if(resabs.gt.uflow/(0.5D+02*epmach)) abserr = max &
    ((epmach*0.5D+02)*resabs,abserr)

  return
end
subroutine dqk15i(f,boun,inf,a,b,result,abserr,resabs,resasc)

!*****************************************************************************80
!
!! DQK15I applies a 15 point Gauss-Kronrod quadrature on an infinite interval.
!
!
!     the abscissae and weights are supplied for the interval
!     (-1,1).  because of symmetry only the positive abscissae and
!     their corresponding weights are given.
!
!     xgk    - abscissae of the 15-point kronrod rule
!              xgk(2), xgk(4), ... abscissae of the 7-point
!              gauss rule
!              xgk(1), xgk(3), ...  abscissae which are optimally
!              added to the 7-point gauss rule
!
!     wgk    - weights of the 15-point kronrod rule
!
!     wg     - weights of the 7-point gauss rule, corresponding
!              to the abscissae xgk(2), xgk(4), ...
!              wg(1), wg(3), ... are set to zero.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  the original (infinite integration range is mapped
!      onto the interval (0,1) and (a,b) is a part of (0,1).
!      it is the purpose to compute
!      i = integral of transformed integrand over (a,b),
!      j = integral of abs(transformed integrand) over (a,b).
!
!  Parameters:
!
!      on entry
!        f      - real ( kind = 8 )
!                 fuction subprogram defining the integrand
!                 function f(x). the actual name for f needs to be
!                 declared e x t e r n a l in the calling program.
!
!        boun   - real ( kind = 8 )
!                 finite bound of original integration
!                 range (set to zero if inf = +2)
!
!        inf    - integer ( kind = 4 )
!                 if inf = -1, the original interval is
!                             (-infinity,bound),
!                 if inf = +1, the original interval is
!                             (bound,+infinity),
!                 if inf = +2, the original interval is
!                             (-infinity,+infinity) and
!                 the integral is computed as the sum of two
!                 integrals, one over (-infinity,0) and one over
!                 (0,+infinity).
!
!        a      - real ( kind = 8 )
!                 lower limit for integration over subrange
!                 of (0,1)
!
!        b      - real ( kind = 8 )
!                 upper limit for integration over subrange
!                 of (0,1)
!
!      on return
!        result - real ( kind = 8 )
!                 approximation to the integral i
!                 result is computed by applying the 15-point
!                 kronrod rule(resk) obtained by optimal addition
!                 of abscissae to the 7-point gauss rule(resg).
!
!        abserr - real ( kind = 8 )
!                 estimate of the modulus of the absolute error,
!                 which should equal or exceed abs(i-result)
!
!        resabs - real ( kind = 8 )
!                 approximation to the integral j
!
!        resasc - real ( kind = 8 )
!                 approximation to the integral of
!                 abs((transformed integrand)-i/(b-a)) over (a,b)
!
!  Local Parameters:
!
!     centr  - mid point of the interval
!     hlgth  - half-length of the interval
!     absc*  - abscissa
!     tabsc* - transformed abscissa
!     fval*  - function value
!     resg   - result of the 7-point gauss formula
!     resk   - result of the 15-point kronrod formula
!     reskh  - approximation to the mean value of the transformed
!              integrand over (a,b), i.e. to i/(b-a)
!
!     machine dependent constants
!
!     epmach is the largest relative spacing.
!     uflow is the smallest positive magnitude.
!
  implicit none

  real ( kind = 8 ) a,absc,absc1,absc2,abserr,b,boun,centr,dinf, &
    epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth, &
    resabs,resasc,resg,resk,reskh,result,tabsc1,tabsc2,uflow,wg,wgk, &
    xgk
  integer ( kind = 4 ) inf,j
  external f
  dimension fv1(7),fv2(7),xgk(8),wgk(8),wg(8)

  data wg(1) / 0.0d0 /
  data wg(2) / 0.129484966168869693270611432679082d0 /
  data wg(3) / 0.0d0 /
  data wg(4) / 0.279705391489276667901467771423780d0 /
  data wg(5) / 0.0d0 /
  data wg(6) / 0.381830050505118944950369775488975d0 /
  data wg(7) / 0.0d0 /
  data wg(8) / 0.417959183673469387755102040816327d0 /

  data xgk(1) / 0.991455371120812639206854697526329d0 /
  data xgk(2) / 0.949107912342758524526189684047851d0 /
  data xgk(3) / 0.864864423359769072789712788640926d0 /
  data xgk(4) / 0.741531185599394439863864773280788d0 /
  data xgk(5) / 0.586087235467691130294144838258730d0 /
  data xgk(6) / 0.405845151377397166906606412076961d0 /
  data xgk(7) / 0.207784955007898467600689403773245d0 /
  data xgk(8) / 0.000000000000000000000000000000000d0 /

  data wgk(1) / 0.022935322010529224963732008058970d0 /
  data wgk(2) / 0.063092092629978553290700663189204d0 /
  data wgk(3) / 0.104790010322250183839876322541518d0 /
  data wgk(4) / 0.140653259715525918745189590510238d0 /
  data wgk(5) / 0.169004726639267902826583426598550d0 /
  data wgk(6) / 0.190350578064785409913256402421014d0 /
  data wgk(7) / 0.204432940075298892414161999234649d0 /
  data wgk(8) / 0.209482141084727828012999174891714d0 /

  epmach = epsilon ( epmach )
  uflow = tiny ( uflow )
  dinf = min ( 1, inf )
  centr = 0.5D+00*(a+b)
  hlgth = 0.5D+00*(b-a)
  tabsc1 = boun+dinf*(0.1D+01-centr)/centr
  fval1 = f(tabsc1)
  if(inf.eq.2) fval1 = fval1+f(-tabsc1)
  fc = (fval1/centr)/centr
!
!  compute the 15-point kronrod approximation to
!  the integral, and estimate the error.
!
  resg = wg(8)*fc
  resk = wgk(8)*fc
  resabs =  abs ( resk)

  do j=1,7
    absc = hlgth*xgk(j)
    absc1 = centr-absc
    absc2 = centr+absc
    tabsc1 = boun+dinf*(0.1D+01-absc1)/absc1
    tabsc2 = boun+dinf*(0.1D+01-absc2)/absc2
    fval1 = f(tabsc1)
    fval2 = f(tabsc2)
    if(inf.eq.2) fval1 = fval1+f(-tabsc1)
    if(inf.eq.2) fval2 = fval2+f(-tabsc2)
    fval1 = (fval1/absc1)/absc1
    fval2 = (fval2/absc2)/absc2
    fv1(j) = fval1
    fv2(j) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(j)*fsum
    resabs = resabs+wgk(j)*( abs ( fval1)+ abs ( fval2))
  end do

  reskh = resk*0.5D+00
  resasc = wgk(8)* abs ( fc-reskh)

  do j=1,7
    resasc = resasc+wgk(j)*( abs ( fv1(j)-reskh)+ abs ( fv2(j)-reskh))
  end do

  result = resk*hlgth
  resasc = resasc*hlgth
  resabs = resabs*hlgth
  abserr =  abs ( (resk-resg)*hlgth)
  if(resasc.ne.0.0D+00.and.abserr.ne.0.d0) abserr = resasc* &
   min (0.1D+01,(0.2D+03*abserr/resasc)**1.5D+00)
  if(resabs.gt.uflow/(0.5D+02*epmach)) abserr = max &
   ((epmach*0.5D+02)*resabs,abserr)

  return
end
subroutine dqk15w(f,w,p1,p2,p3,p4,kp,a,b,result,abserr, resabs,resasc)

!*****************************************************************************80
!
!! DQK15W applies a 15 point Gauss-Kronrod rule for a weighted integrand.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  to compute i = integral of f*w over (a,b), with error
!                     estimate
!                 j = integral of abs(f*w) over (a,b)
!
!  Parameters:
!
!       on entry
!        f      - real ( kind = 8 )
!                 function subprogram defining the integrand
!                 function f(x). the actual name for f needs to be
!                 declared e x t e r n a l in the driver program.
!
!        w      - real ( kind = 8 )
!                 function subprogram defining the integrand
!                 weight function w(x). the actual name for w
!                 needs to be declared e x t e r n a l in the
!                 calling program.
!
!        p1, p2, p3, p4 - real ( kind = 8 )
!                 parameters in the weight function
!
!        kp     - integer ( kind = 4 )
!                 key for indicating the type of weight function
!
!        a      - real ( kind = 8 )
!                 lower limit of integration
!
!        b      - real ( kind = 8 )
!                 upper limit of integration
!
!      on return
!        result - real ( kind = 8 )
!                 approximation to the integral i
!                 result is computed by applying the 15-point
!                 kronrod rule (resk) obtained by optimal addition
!                 of abscissae to the 7-point gauss rule (resg).
!
!        abserr - real ( kind = 8 )
!                 estimate of the modulus of the absolute error,
!                 which should equal or exceed abs(i-result)
!
!        resabs - real ( kind = 8 )
!                 approximation to the integral of abs(f)
!
!        resasc - real ( kind = 8 )
!                 approximation to the integral of abs(f-i/(b-a))
!
!  Local Parameters:
!
!     the abscissae and weights are given for the interval (-1,1).
!     because of symmetry only the positive abscissae and their
!     corresponding weights are given.
!
!     xgk    - abscissae of the 15-point gauss-kronrod rule
!              xgk(2), xgk(4), ... abscissae of the 7-point
!              gauss rule
!              xgk(1), xgk(3), ... abscissae which are optimally
!              added to the 7-point gauss rule
!
!     wgk    - weights of the 15-point gauss-kronrod rule
!
!     wg     - weights of the 7-point gauss rule
!
!     centr  - mid point of the interval
!     hlgth  - half-length of the interval
!     absc*  - abscissa
!     fval*  - function value
!     resg   - result of the 7-point gauss formula
!     resk   - result of the 15-point kronrod formula
!     reskh  - approximation to the mean value of f*w over (a,b),
!              i.e. to i/(b-a)
!
!     machine dependent constants
!
!     epmach is the largest relative spacing.
!     uflow is the smallest positive magnitude.
!
  implicit none

  real ( kind = 8 ) a,absc,absc1,absc2,abserr,b,centr,dhlgth, &
    epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth, &
    p1,p2,p3,p4,resabs,resasc,resg,resk,reskh,result,uflow,w,wg,wgk, &
    xgk
  integer ( kind = 4 ) j,jtw,jtwm1,kp
  external f,w

  dimension fv1(7),fv2(7),xgk(8),wgk(8),wg(4)

  data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8)/ &
       0.9914553711208126D+00,     0.9491079123427585D+00, &
       0.8648644233597691D+00,     0.7415311855993944D+00, &
       0.5860872354676911D+00,     0.4058451513773972D+00, &
       0.2077849550078985D+00,     0.0000000000000000D+00/

  data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8)/ &
       0.2293532201052922D-01,     0.6309209262997855D-01, &
       0.1047900103222502D+00,     0.1406532597155259D+00, &
       0.1690047266392679D+00,     0.1903505780647854D+00, &
       0.2044329400752989D+00,     0.2094821410847278D+00/

  data wg(1),wg(2),wg(3),wg(4)/ &
       0.1294849661688697D+00,    0.2797053914892767D+00, &
       0.3818300505051889D+00,    0.4179591836734694D+00/

  epmach = epsilon ( epmach )
  uflow = tiny ( uflow )
  centr = 0.5D+00*(a+b)
  hlgth = 0.5D+00*(b-a)
  dhlgth =  abs ( hlgth)
!
!  compute the 15-point kronrod approximation to the
!  integral, and estimate the error.
!
  fc = f(centr)*w(centr,p1,p2,p3,p4,kp)
  resg = wg(4)*fc
  resk = wgk(8)*fc
  resabs =  abs ( resk)

  do j=1,3
    jtw = j*2
    absc = hlgth*xgk(jtw)
    absc1 = centr-absc
    absc2 = centr+absc
    fval1 = f(absc1)*w(absc1,p1,p2,p3,p4,kp)
    fval2 = f(absc2)*w(absc2,p1,p2,p3,p4,kp)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*( abs ( fval1)+ abs ( fval2))
  end do

  do j=1,4
    jtwm1 = j*2-1
    absc = hlgth*xgk(jtwm1)
    absc1 = centr-absc
    absc2 = centr+absc
    fval1 = f(absc1)*w(absc1,p1,p2,p3,p4,kp)
    fval2 = f(absc2)*w(absc2,p1,p2,p3,p4,kp)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*( abs ( fval1)+ abs ( fval2))
  end do

  reskh = resk*0.5D+00
  resasc = wgk(8)* abs ( fc-reskh)

  do j=1,7
    resasc = resasc+wgk(j)*( abs ( fv1(j)-reskh)+ abs ( fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr =  abs ( (resk-resg)*hlgth)
  if(resasc.ne.0.0D+00.and.abserr.ne.0.0D+00) &
    abserr = resasc* min (0.1D+01,(0.2D+03*abserr/resasc)**1.5D+00)
  if(resabs.gt.uflow/(0.5D+02*epmach)) abserr =  max ( (epmach* &
    0.5D+02)*resabs,abserr)

  return
end
subroutine dqk21(f,a,b,result,abserr,resabs,resasc)

!*****************************************************************************80
!
!! DQK21 carries out a 21 point Gauss-Kronrod quadrature rule.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  to compute i = integral of f over (a,b), with error
!                     estimate
!                 j = integral of abs(f) over (a,b)
!
!  Parameters:
!
!      on entry
!        f      - real ( kind = 8 )
!                 function subprogram defining the integrand
!                 function f(x). the actual name for f needs to be
!                 declared e x t e r n a l in the driver program.
!
!        a      - real ( kind = 8 )
!                 lower limit of integration
!
!        b      - real ( kind = 8 )
!                 upper limit of integration
!
!      on return
!        result - real ( kind = 8 )
!                 approximation to the integral i
!                 result is computed by applying the 21-point
!                 kronrod rule (resk) obtained by optimal addition
!                 of abscissae to the 10-point gauss rule (resg).
!
!        abserr - real ( kind = 8 )
!                 estimate of the modulus of the absolute error,
!                 which should not exceed abs(i-result)
!
!        resabs - real ( kind = 8 )
!                 approximation to the integral j
!
!        resasc - real ( kind = 8 )
!                 approximation to the integral of abs(f-i/(b-a))
!                 over (a,b)
!
!  Local Parameters:
!
!
!     the abscissae and weights are given for the interval (-1,1).
!     because of symmetry only the positive abscissae and their
!     corresponding weights are given.
!
!     xgk    - abscissae of the 21-point kronrod rule
!              xgk(2), xgk(4), ...  abscissae of the 10-point
!              gauss rule
!              xgk(1), xgk(3), ...  abscissae which are optimally
!              added to the 10-point gauss rule
!
!     wgk    - weights of the 21-point kronrod rule
!
!     wg     - weights of the 10-point gauss rule
!
!
! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.
!
!     centr  - mid point of the interval
!     hlgth  - half-length of the interval
!     absc   - abscissa
!     fval*  - function value
!     resg   - result of the 10-point gauss formula
!     resk   - result of the 21-point kronrod formula
!     reskh  - approximation to the mean value of f over (a,b),
!              i.e. to i/(b-a)
!
!
!     machine dependent constants
!
!     epmach is the largest relative spacing.
!     uflow is the smallest positive magnitude.
!
  implicit none

  real ( kind = 8 ) a,absc,abserr,b,centr,dhlgth, &
    epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc, &
    resg,resk,reskh,result,uflow,wg,wgk,xgk
  integer ( kind = 4 ) j,jtw,jtwm1
  external f
  dimension fv1(10),fv2(10),wg(5),wgk(11),xgk(11)

  data wg  (  1) / 0.066671344308688137593568809893332d0 /
  data wg  (  2) / 0.149451349150580593145776339657697d0 /
  data wg  (  3) / 0.219086362515982043995534934228163d0 /
  data wg  (  4) / 0.269266719309996355091226921569469d0 /
  data wg  (  5) / 0.295524224714752870173892994651338d0 /

  data xgk (  1) / 0.995657163025808080735527280689003d0 /
  data xgk (  2) / 0.973906528517171720077964012084452d0 /
  data xgk (  3) / 0.930157491355708226001207180059508d0 /
  data xgk (  4) / 0.865063366688984510732096688423493d0 /
  data xgk (  5) / 0.780817726586416897063717578345042d0 /
  data xgk (  6) / 0.679409568299024406234327365114874d0 /
  data xgk (  7) / 0.562757134668604683339000099272694d0 /
  data xgk (  8) / 0.433395394129247190799265943165784d0 /
  data xgk (  9) / 0.294392862701460198131126603103866d0 /
  data xgk ( 10) / 0.148874338981631210884826001129720d0 /
  data xgk ( 11) / 0.000000000000000000000000000000000d0 /

  data wgk (  1) / 0.011694638867371874278064396062192d0 /
  data wgk (  2) / 0.032558162307964727478818972459390d0 /
  data wgk (  3) / 0.054755896574351996031381300244580d0 /
  data wgk (  4) / 0.075039674810919952767043140916190d0 /
  data wgk (  5) / 0.093125454583697605535065465083366d0 /
  data wgk (  6) / 0.109387158802297641899210590325805d0 /
  data wgk (  7) / 0.123491976262065851077958109831074d0 /
  data wgk (  8) / 0.134709217311473325928054001771707d0 /
  data wgk (  9) / 0.142775938577060080797094273138717d0 /
  data wgk ( 10) / 0.147739104901338491374841515972068d0 /
  data wgk ( 11) / 0.149445554002916905664936468389821d0 /

  epmach = epsilon ( epmach )
  uflow = tiny ( uflow )
  centr = 0.5D+00*(a+b)
  hlgth = 0.5D+00*(b-a)
  dhlgth =  abs ( hlgth)
!
!  compute the 21-point kronrod approximation to
!  the integral, and estimate the absolute error.
!
  resg = 0.0D+00
  fc = f(centr)
  resk = wgk(11)*fc
  resabs =  abs ( resk)
  do j=1,5
    jtw = 2*j
    absc = hlgth*xgk(jtw)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*( abs ( fval1)+ abs ( fval2))
  end do

  do j = 1,5
    jtwm1 = 2*j-1
    absc = hlgth*xgk(jtwm1)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*( abs ( fval1)+ abs ( fval2))
  end do

  reskh = resk*0.5D+00
  resasc = wgk(11)* abs ( fc-reskh)

  do j=1,10
    resasc = resasc+wgk(j)*( abs ( fv1(j)-reskh)+ abs ( fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr =  abs ( (resk-resg)*hlgth)
  if(resasc.ne.0.0D+00.and.abserr.ne.0.0D+00) &
    abserr = resasc*min(0.1D+01,(0.2D+03*abserr/resasc)**1.5D+00)
  if(resabs.gt.uflow/(0.5D+02*epmach)) abserr = max &
    ((epmach*0.5D+02)*resabs,abserr)

  return
end
subroutine dqk31(f,a,b,result,abserr,resabs,resasc)

!*****************************************************************************80
!
!! DQK31 carries out a 31 point Gauss-Kronrod quadrature rule.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  to compute i = integral of f over (a,b) with error
!                     estimate
!                 j = integral of abs(f) over (a,b)
!
!  Parameters:
!
!      on entry
!        f      - real ( kind = 8 )
!                 function subprogram defining the integrand
!                 function f(x). the actual name for f needs to be
!                 declared e x t e r n a l in the calling program.
!
!        a      - real ( kind = 8 )
!                 lower limit of integration
!
!        b      - real ( kind = 8 )
!                 upper limit of integration
!
!      on return
!        result - real ( kind = 8 )
!                 approximation to the integral i
!                 result is computed by applying the 31-point
!                 gauss-kronrod rule (resk), obtained by optimal
!                 addition of abscissae to the 15-point gauss
!                 rule (resg).
!
!        abserr - double precison
!                 estimate of the modulus of the modulus,
!                 which should not exceed abs(i-result)
!
!        resabs - real ( kind = 8 )
!                 approximation to the integral j
!
!        resasc - real ( kind = 8 )
!                 approximation to the integral of abs(f-i/(b-a))
!                 over (a,b)
!
!  Local Parameters:
!
!
!     the abscissae and weights are given for the interval (-1,1).
!     because of symmetry only the positive abscissae and their
!     corresponding weights are given.
!
!     xgk    - abscissae of the 31-point kronrod rule
!              xgk(2), xgk(4), ...  abscissae of the 15-point
!              gauss rule
!              xgk(1), xgk(3), ...  abscissae which are optimally
!              added to the 15-point gauss rule
!
!     wgk    - weights of the 31-point kronrod rule
!
!     wg     - weights of the 15-point gauss rule
!
!
! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.
!
!     centr  - mid point of the interval
!     hlgth  - half-length of the interval
!     absc   - abscissa
!     fval*  - function value
!     resg   - result of the 15-point gauss formula
!     resk   - result of the 31-point kronrod formula
!     reskh  - approximation to the mean value of f over (a,b),
!              i.e. to i/(b-a)
!
!     machine dependent constants
!
!     epmach is the largest relative spacing.
!     uflow is the smallest positive magnitude.
!
  implicit none

  real ( kind = 8 ) a,absc,abserr,b,centr,dhlgth, &
    epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc, &
    resg,resk,reskh,result,uflow,wg,wgk,xgk
  integer ( kind = 4 ) j,jtw,jtwm1
  external f

  dimension fv1(15),fv2(15),xgk(16),wgk(16),wg(8)

  data wg  (  1) / 0.030753241996117268354628393577204d0 /
  data wg  (  2) / 0.070366047488108124709267416450667d0 /
  data wg  (  3) / 0.107159220467171935011869546685869d0 /
  data wg  (  4) / 0.139570677926154314447804794511028d0 /
  data wg  (  5) / 0.166269205816993933553200860481209d0 /
  data wg  (  6) / 0.186161000015562211026800561866423d0 /
  data wg  (  7) / 0.198431485327111576456118326443839d0 /
  data wg  (  8) / 0.202578241925561272880620199967519d0 /

  data xgk (  1) / 0.998002298693397060285172840152271d0 /
  data xgk (  2) / 0.987992518020485428489565718586613d0 /
  data xgk (  3) / 0.967739075679139134257347978784337d0 /
  data xgk (  4) / 0.937273392400705904307758947710209d0 /
  data xgk (  5) / 0.897264532344081900882509656454496d0 /
  data xgk (  6) / 0.848206583410427216200648320774217d0 /
  data xgk (  7) / 0.790418501442465932967649294817947d0 /
  data xgk (  8) / 0.724417731360170047416186054613938d0 /
  data xgk (  9) / 0.650996741297416970533735895313275d0 /
  data xgk ( 10) / 0.570972172608538847537226737253911d0 /
  data xgk ( 11) / 0.485081863640239680693655740232351d0 /
  data xgk ( 12) / 0.394151347077563369897207370981045d0 /
  data xgk ( 13) / 0.299180007153168812166780024266389d0 /
  data xgk ( 14) / 0.201194093997434522300628303394596d0 /
  data xgk ( 15) / 0.101142066918717499027074231447392d0 /
  data xgk ( 16) / 0.000000000000000000000000000000000d0 /

  data wgk (  1) / 0.005377479872923348987792051430128d0 /
  data wgk (  2) / 0.015007947329316122538374763075807d0 /
  data wgk (  3) / 0.025460847326715320186874001019653d0 /
  data wgk (  4) / 0.035346360791375846222037948478360d0 /
  data wgk (  5) / 0.044589751324764876608227299373280d0 /
  data wgk (  6) / 0.053481524690928087265343147239430d0 /
  data wgk (  7) / 0.062009567800670640285139230960803d0 /
  data wgk (  8) / 0.069854121318728258709520077099147d0 /
  data wgk (  9) / 0.076849680757720378894432777482659d0 /
  data wgk ( 10) / 0.083080502823133021038289247286104d0 /
  data wgk ( 11) / 0.088564443056211770647275443693774d0 /
  data wgk ( 12) / 0.093126598170825321225486872747346d0 /
  data wgk ( 13) / 0.096642726983623678505179907627589d0 /
  data wgk ( 14) / 0.099173598721791959332393173484603d0 /
  data wgk ( 15) / 0.100769845523875595044946662617570d0 /
  data wgk ( 16) / 0.101330007014791549017374792767493d0 /

  epmach = epsilon ( epmach )
  uflow = tiny ( uflow )
  centr = 0.5D+00*(a+b)
  hlgth = 0.5D+00*(b-a)
  dhlgth =  abs ( hlgth)
!
!  compute the 31-point kronrod approximation to
!  the integral, and estimate the absolute error.
!
  fc = f(centr)
  resg = wg(8)*fc
  resk = wgk(16)*fc
  resabs =  abs ( resk)

  do j=1,7
    jtw = j*2
    absc = hlgth*xgk(jtw)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*( abs ( fval1)+ abs ( fval2))
  end do

  do j = 1,8
    jtwm1 = j*2-1
    absc = hlgth*xgk(jtwm1)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*( abs ( fval1)+ abs ( fval2))
  end do

  reskh = resk*0.5D+00
  resasc = wgk(16)* abs ( fc-reskh)

  do j=1,15
    resasc = resasc+wgk(j)*( abs ( fv1(j)-reskh)+ abs ( fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr =  abs ( (resk-resg)*hlgth)
  if(resasc.ne.0.0D+00.and.abserr.ne.0.0D+00) &
    abserr = resasc* min (0.1D+01,(0.2D+03*abserr/resasc)**1.5D+00)
  if(resabs.gt.uflow/(0.5D+02*epmach)) abserr = max &
    ((epmach*0.5D+02)*resabs,abserr)

  return
end
subroutine dqk41 ( f, a, b, result, abserr, resabs, resasc )

!*****************************************************************************80
!
!! DQK41 carries out a 41 point Gauss-Kronrod quadrature rule.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  to compute i = integral of f over (a,b), with error
!                     estimate
!                 j = integral of abs(f) over (a,b)
!
!  Parameters:
!
!      on entry
!        f      - real ( kind = 8 )
!                 function subprogram defining the integrand
!                 function f(x). the actual name for f needs to be
!                 declared e x t e r n a l in the calling program.
!
!        a      - real ( kind = 8 )
!                 lower limit of integration
!
!        b      - real ( kind = 8 )
!                 upper limit of integration
!
!      on return
!        result - real ( kind = 8 )
!                 approximation to the integral i
!                 result is computed by applying the 41-point
!                 gauss-kronrod rule (resk) obtained by optimal
!                 addition of abscissae to the 20-point gauss
!                 rule (resg).
!
!        abserr - real ( kind = 8 )
!                 estimate of the modulus of the absolute error,
!                 which should not exceed abs(i-result)
!
!        resabs - real ( kind = 8 )
!                 approximation to the integral j
!
!        resasc - real ( kind = 8 )
!                 approximation to the integal of abs(f-i/(b-a))
!                 over (a,b)
!
!  Local Parameters:
!
!
!     the abscissae and weights are given for the interval (-1,1).
!     because of symmetry only the positive abscissae and their
!     corresponding weights are given.
!
!     xgk    - abscissae of the 41-point gauss-kronrod rule
!              xgk(2), xgk(4), ...  abscissae of the 20-point
!              gauss rule
!              xgk(1), xgk(3), ...  abscissae which are optimally
!              added to the 20-point gauss rule
!
!     wgk    - weights of the 41-point gauss-kronrod rule
!
!     wg     - weights of the 20-point gauss rule
!
!
! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.
!
!     centr  - mid point of the interval
!     hlgth  - half-length of the interval
!     absc   - abscissa
!     fval*  - function value
!     resg   - result of the 20-point gauss formula
!     resk   - result of the 41-point kronrod formula
!     reskh  - approximation to mean value of f over (a,b), i.e.
!              to i/(b-a)
!
!     machine dependent constants
!
!     epmach is the largest relative spacing.
!     uflow is the smallest positive magnitude.
!
  implicit none

  real ( kind = 8 ) a,absc,abserr,b,centr,dhlgth, &
    epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc, &
    resg,resk,reskh,result,uflow,wg,wgk,xgk
  integer ( kind = 4 ) j,jtw,jtwm1
  external f

  dimension fv1(20),fv2(20),xgk(21),wgk(21),wg(10)

  data wg  (  1) / 0.017614007139152118311861962351853d0 /
  data wg  (  2) / 0.040601429800386941331039952274932d0 /
  data wg  (  3) / 0.062672048334109063569506535187042d0 /
  data wg  (  4) / 0.083276741576704748724758143222046d0 /
  data wg  (  5) / 0.101930119817240435036750135480350d0 /
  data wg  (  6) / 0.118194531961518417312377377711382d0 /
  data wg  (  7) / 0.131688638449176626898494499748163d0 /
  data wg  (  8) / 0.142096109318382051329298325067165d0 /
  data wg  (  9) / 0.149172986472603746787828737001969d0 /
  data wg  ( 10) / 0.152753387130725850698084331955098d0 /

  data xgk (  1) / 0.998859031588277663838315576545863d0 /
  data xgk (  2) / 0.993128599185094924786122388471320d0 /
  data xgk (  3) / 0.981507877450250259193342994720217d0 /
  data xgk (  4) / 0.963971927277913791267666131197277d0 /
  data xgk (  5) / 0.940822633831754753519982722212443d0 /
  data xgk (  6) / 0.912234428251325905867752441203298d0 /
  data xgk (  7) / 0.878276811252281976077442995113078d0 /
  data xgk (  8) / 0.839116971822218823394529061701521d0 /
  data xgk (  9) / 0.795041428837551198350638833272788d0 /
  data xgk ( 10) / 0.746331906460150792614305070355642d0 /
  data xgk ( 11) / 0.693237656334751384805490711845932d0 /
  data xgk ( 12) / 0.636053680726515025452836696226286d0 /
  data xgk ( 13) / 0.575140446819710315342946036586425d0 /
  data xgk ( 14) / 0.510867001950827098004364050955251d0 /
  data xgk ( 15) / 0.443593175238725103199992213492640d0 /
  data xgk ( 16) / 0.373706088715419560672548177024927d0 /
  data xgk ( 17) / 0.301627868114913004320555356858592d0 /
  data xgk ( 18) / 0.227785851141645078080496195368575d0 /
  data xgk ( 19) / 0.152605465240922675505220241022678d0 /
  data xgk ( 20) / 0.076526521133497333754640409398838d0 /
  data xgk ( 21) / 0.000000000000000000000000000000000d0 /

  data wgk (  1) / 0.003073583718520531501218293246031d0 /
  data wgk (  2) / 0.008600269855642942198661787950102d0 /
  data wgk (  3) / 0.014626169256971252983787960308868d0 /
  data wgk (  4) / 0.020388373461266523598010231432755d0 /
  data wgk (  5) / 0.025882133604951158834505067096153d0 /
  data wgk (  6) / 0.031287306777032798958543119323801d0 /
  data wgk (  7) / 0.036600169758200798030557240707211d0 /
  data wgk (  8) / 0.041668873327973686263788305936895d0 /
  data wgk (  9) / 0.046434821867497674720231880926108d0 /
  data wgk ( 10) / 0.050944573923728691932707670050345d0 /
  data wgk ( 11) / 0.055195105348285994744832372419777d0 /
  data wgk ( 12) / 0.059111400880639572374967220648594d0 /
  data wgk ( 13) / 0.062653237554781168025870122174255d0 /
  data wgk ( 14) / 0.065834597133618422111563556969398d0 /
  data wgk ( 15) / 0.068648672928521619345623411885368d0 /
  data wgk ( 16) / 0.071054423553444068305790361723210d0 /
  data wgk ( 17) / 0.073030690332786667495189417658913d0 /
  data wgk ( 18) / 0.074582875400499188986581418362488d0 /
  data wgk ( 19) / 0.075704497684556674659542775376617d0 /
  data wgk ( 20) / 0.076377867672080736705502835038061d0 /
  data wgk ( 21) / 0.076600711917999656445049901530102d0 /

  epmach = epsilon ( epmach )
  uflow = tiny ( uflow )
  centr = 0.5D+00*(a+b)
  hlgth = 0.5D+00*(b-a)
  dhlgth =  abs ( hlgth)
!
!  compute the 41-point gauss-kronrod approximation to
!  the integral, and estimate the absolute error.
!
  resg = 0.0D+00
  fc = f(centr)
  resk = wgk(21)*fc
  resabs =  abs ( resk)

  do j=1,10
    jtw = j*2
    absc = hlgth*xgk(jtw)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*( abs ( fval1)+ abs ( fval2))
  end do

  do j = 1,10
    jtwm1 = j*2-1
    absc = hlgth*xgk(jtwm1)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*( abs ( fval1)+ abs ( fval2))
  end do

  reskh = resk*0.5D+00
  resasc = wgk(21)* abs ( fc-reskh)

  do j=1,20
    resasc = resasc+wgk(j)*( abs ( fv1(j)-reskh)+ abs ( fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr =  abs ( (resk-resg)*hlgth)
  if(resasc.ne.0.0D+00.and.abserr.ne.0.D+00) &
    abserr = resasc* min (0.1D+01,(0.2D+03*abserr/resasc)**1.5D+00)
  if(resabs.gt.uflow/(0.5D+02*epmach)) abserr = max &
    ((epmach*0.5D+02)*resabs,abserr)

  return
end
subroutine dqk51(f,a,b,result,abserr,resabs,resasc)

!*****************************************************************************80
!
!! DQK51 carries out a 51 point Gauss-Kronrod quadrature rule.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  to compute i = integral of f over (a,b) with error
!                     estimate
!                 j = integral of abs(f) over (a,b)
!
!  Parameters:
!
!      on entry
!        f      - real ( kind = 8 )
!                 function defining the integrand
!                 function f(x). the actual name for f needs to be
!                 declared e x t e r n a l in the calling program.
!
!        a      - real ( kind = 8 )
!                 lower limit of integration
!
!        b      - real ( kind = 8 )
!                 upper limit of integration
!
!      on return
!        result - real ( kind = 8 )
!                 approximation to the integral i
!                 result is computed by applying the 51-point
!                 kronrod rule (resk) obtained by optimal addition
!                 of abscissae to the 25-point gauss rule (resg).
!
!        abserr - real ( kind = 8 )
!                 estimate of the modulus of the absolute error,
!                 which should not exceed abs(i-result)
!
!        resabs - real ( kind = 8 )
!                 approximation to the integral j
!
!        resasc - real ( kind = 8 )
!                 approximation to the integral of abs(f-i/(b-a))
!                 over (a,b)
!
!  Local Parameters:
!
!     the abscissae and weights are given for the interval (-1,1).
!     because of symmetry only the positive abscissae and their
!     corresponding weights are given.
!
!     xgk    - abscissae of the 51-point kronrod rule
!              xgk(2), xgk(4), ...  abscissae of the 25-point
!              gauss rule
!              xgk(1), xgk(3), ...  abscissae which are optimally
!              added to the 25-point gauss rule
!
!     wgk    - weights of the 51-point kronrod rule
!
!     wg     - weights of the 25-point gauss rule
!
! gauss quadrature weights and kronron quadrature abscissae and weights
! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
! bell labs, nov. 1981.
!
!     centr  - mid point of the interval
!     hlgth  - half-length of the interval
!     absc   - abscissa
!     fval*  - function value
!     resg   - result of the 25-point gauss formula
!     resk   - result of the 51-point kronrod formula
!     reskh  - approximation to the mean value of f over (a,b),
!              i.e. to i/(b-a)
!
!     machine dependent constants
!
!     epmach is the largest relative spacing.
!     uflow is the smallest positive magnitude.
!
  implicit none

  real ( kind = 8 ) a,absc,abserr,b,centr,dhlgth, &
    epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc, &
    resg,resk,reskh,result,uflow,wg,wgk,xgk
  integer ( kind = 4 ) j,jtw,jtwm1
  external f

  dimension fv1(25),fv2(25),xgk(26),wgk(26),wg(13)

  data wg  (  1) / 0.011393798501026287947902964113235d0 /
  data wg  (  2) / 0.026354986615032137261901815295299d0 /
  data wg  (  3) / 0.040939156701306312655623487711646d0 /
  data wg  (  4) / 0.054904695975835191925936891540473d0 /
  data wg  (  5) / 0.068038333812356917207187185656708d0 /
  data wg  (  6) / 0.080140700335001018013234959669111d0 /
  data wg  (  7) / 0.091028261982963649811497220702892d0 /
  data wg  (  8) / 0.100535949067050644202206890392686d0 /
  data wg  (  9) / 0.108519624474263653116093957050117d0 /
  data wg  ( 10) / 0.114858259145711648339325545869556d0 /
  data wg  ( 11) / 0.119455763535784772228178126512901d0 /
  data wg  ( 12) / 0.122242442990310041688959518945852d0 /
  data wg  ( 13) / 0.123176053726715451203902873079050d0 /

  data xgk (  1) / 0.999262104992609834193457486540341d0 /
  data xgk (  2) / 0.995556969790498097908784946893902d0 /
  data xgk (  3) / 0.988035794534077247637331014577406d0 /
  data xgk (  4) / 0.976663921459517511498315386479594d0 /
  data xgk (  5) / 0.961614986425842512418130033660167d0 /
  data xgk (  6) / 0.942974571228974339414011169658471d0 /
  data xgk (  7) / 0.920747115281701561746346084546331d0 /
  data xgk (  8) / 0.894991997878275368851042006782805d0 /
  data xgk (  9) / 0.865847065293275595448996969588340d0 /
  data xgk ( 10) / 0.833442628760834001421021108693570d0 /
  data xgk ( 11) / 0.797873797998500059410410904994307d0 /
  data xgk ( 12) / 0.759259263037357630577282865204361d0 /
  data xgk ( 13) / 0.717766406813084388186654079773298d0 /
  data xgk ( 14) / 0.673566368473468364485120633247622d0 /
  data xgk ( 15) / 0.626810099010317412788122681624518d0 /
  data xgk ( 16) / 0.577662930241222967723689841612654d0 /
  data xgk ( 17) / 0.526325284334719182599623778158010d0 /
  data xgk ( 18) / 0.473002731445714960522182115009192d0 /
  data xgk ( 19) / 0.417885382193037748851814394594572d0 /
  data xgk ( 20) / 0.361172305809387837735821730127641d0 /
  data xgk ( 21) / 0.303089538931107830167478909980339d0 /
  data xgk ( 22) / 0.243866883720988432045190362797452d0 /
  data xgk ( 23) / 0.183718939421048892015969888759528d0 /
  data xgk ( 24) / 0.122864692610710396387359818808037d0 /
  data xgk ( 25) / 0.061544483005685078886546392366797d0 /
  data xgk ( 26) / 0.000000000000000000000000000000000d0 /

  data wgk (  1) / 0.001987383892330315926507851882843d0 /
  data wgk (  2) / 0.005561932135356713758040236901066d0 /
  data wgk (  3) / 0.009473973386174151607207710523655d0 /
  data wgk (  4) / 0.013236229195571674813656405846976d0 /
  data wgk (  5) / 0.016847817709128298231516667536336d0 /
  data wgk (  6) / 0.020435371145882835456568292235939d0 /
  data wgk (  7) / 0.024009945606953216220092489164881d0 /
  data wgk (  8) / 0.027475317587851737802948455517811d0 /
  data wgk (  9) / 0.030792300167387488891109020215229d0 /
  data wgk ( 10) / 0.034002130274329337836748795229551d0 /
  data wgk ( 11) / 0.037116271483415543560330625367620d0 /
  data wgk ( 12) / 0.040083825504032382074839284467076d0 /
  data wgk ( 13) / 0.042872845020170049476895792439495d0 /
  data wgk ( 14) / 0.045502913049921788909870584752660d0 /
  data wgk ( 15) / 0.047982537138836713906392255756915d0 /
  data wgk ( 16) / 0.050277679080715671963325259433440d0 /
  data wgk ( 17) / 0.052362885806407475864366712137873d0 /
  data wgk ( 18) / 0.054251129888545490144543370459876d0 /
  data wgk ( 19) / 0.055950811220412317308240686382747d0 /
  data wgk ( 20) / 0.057437116361567832853582693939506d0 /
  data wgk ( 21) / 0.058689680022394207961974175856788d0 /
  data wgk ( 22) / 0.059720340324174059979099291932562d0 /
  data wgk ( 23) / 0.060539455376045862945360267517565d0 /
  data wgk ( 24) / 0.061128509717053048305859030416293d0 /
  data wgk ( 25) / 0.061471189871425316661544131965264d0 /
  data wgk ( 26) / 0.061580818067832935078759824240066d0 /

  epmach = epsilon ( epmach )
  uflow = tiny ( uflow )
  centr = 0.5D+00*(a+b)
  hlgth = 0.5D+00*(b-a)
  dhlgth =  abs ( hlgth)
!
!  compute the 51-point kronrod approximation to
!  the integral, and estimate the absolute error.
!
  fc = f(centr)
  resg = wg(13)*fc
  resk = wgk(26)*fc
  resabs =  abs ( resk)

  do j=1,12
    jtw = j*2
    absc = hlgth*xgk(jtw)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*( abs ( fval1)+ abs ( fval2))
  end do

  do j = 1,13
    jtwm1 = j*2-1
    absc = hlgth*xgk(jtwm1)
    fval1 = f(centr-absc)
    fval2 = f(centr+absc)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*( abs ( fval1)+ abs ( fval2))
  end do

  reskh = resk*0.5D+00
  resasc = wgk(26)* abs ( fc-reskh)

  do j=1,25
    resasc = resasc+wgk(j)*( abs ( fv1(j)-reskh)+ abs ( fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr =  abs ( (resk-resg)*hlgth)
  if(resasc.ne.0.0D+00.and.abserr.ne.0.0D+00) &
    abserr = resasc* min (0.1D+01,(0.2D+03*abserr/resasc)**1.5D+00)
  if(resabs.gt.uflow/(0.5D+02*epmach)) abserr = max &
    ((epmach*0.5D+02)*resabs,abserr)

  return
end
subroutine dqk61(f,a,b,result,abserr,resabs,resasc)

!*****************************************************************************80
!
!! DQK61 carries out a 61 point Gauss-Kronrod quadrature rule.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  to compute i = integral of f over (a,b) with error
!                     estimate
!                 j = integral of  abs ( f) over (a,b)
!
!  Parameters:
!
!   on entry
!     f      - real ( kind = 8 )
!              function subprogram defining the integrand
!              function f(x). the actual name for f needs to be
!              declared e x t e r n a l in the calling program.
!
!     a      - real ( kind = 8 )
!              lower limit of integration
!
!     b      - real ( kind = 8 )
!              upper limit of integration
!
!   on return
!     result - real ( kind = 8 )
!              approximation to the integral i
!              result is computed by applying the 61-point
!              kronrod rule (resk) obtained by optimal addition of
!              abscissae to the 30-point gauss rule (resg).
!
!     abserr - real ( kind = 8 )
!              estimate of the modulus of the absolute error,
!              which should equal or exceed  abs ( i-result)
!
!     resabs - real ( kind = 8 )
!              approximation to the integral j
!
!     resasc - real ( kind = 8 )
!              approximation to the integral of  abs ( f-i/(b-a))
!
!  Local Parameters:
!
!     the abscissae and weights are given for the
!     interval (-1,1). because of symmetry only the positive
!     abscissae and their corresponding weights are given.
!
!     xgk   - abscissae of the 61-point kronrod rule
!             xgk(2), xgk(4)  ... abscissae of the 30-point
!             gauss rule
!             xgk(1), xgk(3)  ... optimally added abscissae
!             to the 30-point gauss rule
!
!     wgk   - weights of the 61-point kronrod rule
!
!     wg    - weigths of the 30-point gauss rule
!
!
!   gauss quadrature weights and kronron quadrature abscissae and weights
!   as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
!   bell labs, nov. 1981.
!
!     centr  - mid point of the interval
!     hlgth  - half-length of the interval
!     dabsc  - abscissa
!     fval*  - function value
!     resg   - result of the 30-point gauss rule
!     resk   - result of the 61-point kronrod rule
!     reskh  - approximation to the mean value of f
!              over (a,b), i.e. to i/(b-a)
!
!     machine dependent constants
!
!     epmach is the largest relative spacing.
!     uflow is the smallest positive magnitude.
!
  implicit none

  real ( kind = 8 ) a,dabsc,abserr,b,centr,dhlgth, &
    epmach,f,fc,fsum,fval1,fval2,fv1,fv2,hlgth,resabs,resasc, &
    resg,resk,reskh,result,uflow,wg,wgk,xgk
  integer ( kind = 4 ) j,jtw,jtwm1
  external f

  dimension fv1(30),fv2(30),xgk(31),wgk(31),wg(15)

  data wg  (  1) / 0.007968192496166605615465883474674d0 /
  data wg  (  2) / 0.018466468311090959142302131912047d0 /
  data wg  (  3) / 0.028784707883323369349719179611292d0 /
  data wg  (  4) / 0.038799192569627049596801936446348d0 /
  data wg  (  5) / 0.048402672830594052902938140422808d0 /
  data wg  (  6) / 0.057493156217619066481721689402056d0 /
  data wg  (  7) / 0.065974229882180495128128515115962d0 /
  data wg  (  8) / 0.073755974737705206268243850022191d0 /
  data wg  (  9) / 0.080755895229420215354694938460530d0 /
  data wg  ( 10) / 0.086899787201082979802387530715126d0 /
  data wg  ( 11) / 0.092122522237786128717632707087619d0 /
  data wg  ( 12) / 0.096368737174644259639468626351810d0 /
  data wg  ( 13) / 0.099593420586795267062780282103569d0 /
  data wg  ( 14) / 0.101762389748405504596428952168554d0 /
  data wg  ( 15) / 0.102852652893558840341285636705415d0 /

  data xgk (  1) / 0.999484410050490637571325895705811d0 /
  data xgk (  2) / 0.996893484074649540271630050918695d0 /
  data xgk (  3) / 0.991630996870404594858628366109486d0 /
  data xgk (  4) / 0.983668123279747209970032581605663d0 /
  data xgk (  5) / 0.973116322501126268374693868423707d0 /
  data xgk (  6) / 0.960021864968307512216871025581798d0 /
  data xgk (  7) / 0.944374444748559979415831324037439d0 /
  data xgk (  8) / 0.926200047429274325879324277080474d0 /
  data xgk (  9) / 0.905573307699907798546522558925958d0 /
  data xgk ( 10) / 0.882560535792052681543116462530226d0 /
  data xgk ( 11) / 0.857205233546061098958658510658944d0 /
  data xgk ( 12) / 0.829565762382768397442898119732502d0 /
  data xgk ( 13) / 0.799727835821839083013668942322683d0 /
  data xgk ( 14) / 0.767777432104826194917977340974503d0 /
  data xgk ( 15) / 0.733790062453226804726171131369528d0 /
  data xgk ( 16) / 0.697850494793315796932292388026640d0 /
  data xgk ( 17) / 0.660061064126626961370053668149271d0 /
  data xgk ( 18) / 0.620526182989242861140477556431189d0 /
  data xgk ( 19) / 0.579345235826361691756024932172540d0 /
  data xgk ( 20) / 0.536624148142019899264169793311073d0 /
  data xgk ( 21) / 0.492480467861778574993693061207709d0 /
  data xgk ( 22) / 0.447033769538089176780609900322854d0 /
  data xgk ( 23) / 0.400401254830394392535476211542661d0 /
  data xgk ( 24) / 0.352704725530878113471037207089374d0 /
  data xgk ( 25) / 0.304073202273625077372677107199257d0 /
  data xgk ( 26) / 0.254636926167889846439805129817805d0 /
  data xgk ( 27) / 0.204525116682309891438957671002025d0 /
  data xgk ( 28) / 0.153869913608583546963794672743256d0 /
  data xgk ( 29) / 0.102806937966737030147096751318001d0 /
  data xgk ( 30) / 0.051471842555317695833025213166723d0 /
  data xgk ( 31) / 0.000000000000000000000000000000000d0 /

  data wgk (  1) / 0.001389013698677007624551591226760d0 /
  data wgk (  2) / 0.003890461127099884051267201844516d0 /
  data wgk (  3) / 0.006630703915931292173319826369750d0 /
  data wgk (  4) / 0.009273279659517763428441146892024d0 /
  data wgk (  5) / 0.011823015253496341742232898853251d0 /
  data wgk (  6) / 0.014369729507045804812451432443580d0 /
  data wgk (  7) / 0.016920889189053272627572289420322d0 /
  data wgk (  8) / 0.019414141193942381173408951050128d0 /
  data wgk (  9) / 0.021828035821609192297167485738339d0 /
  data wgk ( 10) / 0.024191162078080601365686370725232d0 /
  data wgk ( 11) / 0.026509954882333101610601709335075d0 /
  data wgk ( 12) / 0.028754048765041292843978785354334d0 /
  data wgk ( 13) / 0.030907257562387762472884252943092d0 /
  data wgk ( 14) / 0.032981447057483726031814191016854d0 /
  data wgk ( 15) / 0.034979338028060024137499670731468d0 /
  data wgk ( 16) / 0.036882364651821229223911065617136d0 /
  data wgk ( 17) / 0.038678945624727592950348651532281d0 /
  data wgk ( 18) / 0.040374538951535959111995279752468d0 /
  data wgk ( 19) / 0.041969810215164246147147541285970d0 /
  data wgk ( 20) / 0.043452539701356069316831728117073d0 /
  data wgk ( 21) / 0.044814800133162663192355551616723d0 /
  data wgk ( 22) / 0.046059238271006988116271735559374d0 /
  data wgk ( 23) / 0.047185546569299153945261478181099d0 /
  data wgk ( 24) / 0.048185861757087129140779492298305d0 /
  data wgk ( 25) / 0.049055434555029778887528165367238d0 /
  data wgk ( 26) / 0.049795683427074206357811569379942d0 /
  data wgk ( 27) / 0.050405921402782346840893085653585d0 /
  data wgk ( 28) / 0.050881795898749606492297473049805d0 /
  data wgk ( 29) / 0.051221547849258772170656282604944d0 /
  data wgk ( 30) / 0.051426128537459025933862879215781d0 /
  data wgk ( 31) / 0.051494729429451567558340433647099d0 /

  epmach = epsilon ( epmach )
  uflow = tiny ( uflow )
  centr = 0.5D+00*(b+a)
  hlgth = 0.5D+00*(b-a)
  dhlgth =  abs ( hlgth)
!
!  compute the 61-point kronrod approximation to the
!  integral, and estimate the absolute error.
!
  resg = 0.0D+00
  fc = f(centr)
  resk = wgk(31)*fc
  resabs =  abs ( resk)

  do j=1,15
    jtw = j*2
    dabsc = hlgth*xgk(jtw)
    fval1 = f(centr-dabsc)
    fval2 = f(centr+dabsc)
    fv1(jtw) = fval1
    fv2(jtw) = fval2
    fsum = fval1+fval2
    resg = resg+wg(j)*fsum
    resk = resk+wgk(jtw)*fsum
    resabs = resabs+wgk(jtw)*( abs ( fval1)+ abs ( fval2))
  end do

  do j=1,15
    jtwm1 = j*2-1
    dabsc = hlgth*xgk(jtwm1)
    fval1 = f(centr-dabsc)
    fval2 = f(centr+dabsc)
    fv1(jtwm1) = fval1
    fv2(jtwm1) = fval2
    fsum = fval1+fval2
    resk = resk+wgk(jtwm1)*fsum
    resabs = resabs+wgk(jtwm1)*( abs ( fval1)+ abs ( fval2))
  end do

  reskh = resk*0.5D+00
  resasc = wgk(31)* abs ( fc-reskh)

  do j=1,30
    resasc = resasc+wgk(j)*( abs ( fv1(j)-reskh)+ abs ( fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr =  abs ( (resk-resg)*hlgth)
  if(resasc.ne.0.0D+00.and.abserr.ne.0.0D+00) &
    abserr = resasc* min (0.1D+01,(0.2D+03*abserr/resasc)**1.5D+00)
  if(resabs.gt.uflow/(0.5D+02*epmach)) abserr = max &
    ((epmach*0.5D+02)*resabs,abserr)

  return
end
subroutine dqmomo(alfa,beta,ri,rj,rg,rh,integr)

!*****************************************************************************80
!
!! DQMOMO computes modified Chebyshev moments.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  this routine computes modified chebsyshev moments. the k-th
!      modified chebyshev moment is defined as the integral over
!      (-1,1) of w(x)*t(k,x), where t(k,x) is the chebyshev
!      polynomial of degree k.
!
!  Parameters:
!
!     alfa   - real ( kind = 8 )
!              parameter in the weight function w(x), alfa.gt.(-1)
!
!     beta   - real ( kind = 8 )
!              parameter in the weight function w(x), beta.gt.(-1)
!
!     ri     - real ( kind = 8 )
!              vector of dimension 25
!              ri(k) is the integral over (-1,1) of
!              (1+x)**alfa*t(k-1,x), k = 1, ..., 25.
!
!     rj     - real ( kind = 8 )
!              vector of dimension 25
!              rj(k) is the integral over (-1,1) of
!              (1-x)**beta*t(k-1,x), k = 1, ..., 25.
!
!     rg     - real ( kind = 8 )
!              vector of dimension 25
!              rg(k) is the integral over (-1,1) of
!              (1+x)**alfa*log((1+x)/2)*t(k-1,x), k = 1, ..., 25.
!
!     rh     - real ( kind = 8 )
!              vector of dimension 25
!              rh(k) is the integral over (-1,1) of
!              (1-x)**beta*log((1-x)/2)*t(k-1,x), k = 1, ..., 25.
!
!     integr - integer ( kind = 4 )
!              input parameter indicating the modified
!              moments to be computed
!              integr = 1 compute ri, rj
!                     = 2 compute ri, rj, rg
!                     = 3 compute ri, rj, rh
!                     = 4 compute ri, rj, rg, rh
!
  implicit none

  real ( kind = 8 ) alfa,alfp1,alfp2,an,anm1,beta,betp1,betp2,ralf, &
    rbet,rg,rh,ri,rj
  integer ( kind = 4 ) i,im1,integr
  dimension rg(25),rh(25),ri(25),rj(25)

  alfp1 = alfa+0.1D+01
  betp1 = beta+0.1D+01
  alfp2 = alfa+0.2D+01
  betp2 = beta+0.2D+01
  ralf = 0.2D+01**alfp1
  rbet = 0.2D+01**betp1
!
!  compute ri, rj using a forward recurrence relation.
!
  ri(1) = ralf/alfp1
  rj(1) = rbet/betp1
  ri(2) = ri(1)*alfa/alfp2
  rj(2) = rj(1)*beta/betp2
  an = 0.2D+01
  anm1 = 0.1D+01

  do i=3,25
    ri(i) = -(ralf+an*(an-alfp2)*ri(i-1))/(anm1*(an+alfp1))
    rj(i) = -(rbet+an*(an-betp2)*rj(i-1))/(anm1*(an+betp1))
    anm1 = an
    an = an+0.1D+01
  end do

  if(integr.eq.1) go to 70
  if(integr.eq.3) go to 40
!
!  compute rg using a forward recurrence relation.
!
  rg(1) = -ri(1)/alfp1
  rg(2) = -(ralf+ralf)/(alfp2*alfp2)-rg(1)
  an = 0.2D+01
  anm1 = 0.1D+01
  im1 = 2

  do i=3,25
    rg(i) = -(an*(an-alfp2)*rg(im1)-an*ri(im1)+anm1*ri(i))/ &
    (anm1*(an+alfp1))
    anm1 = an
    an = an+0.1D+01
    im1 = i
  end do

  if(integr.eq.2) go to 70
!
!  compute rh using a forward recurrence relation.
!
   40 rh(1) = -rj(1)/betp1
  rh(2) = -(rbet+rbet)/(betp2*betp2)-rh(1)
  an = 0.2D+01
  anm1 = 0.1D+01
  im1 = 2

  do i=3,25
    rh(i) = -(an*(an-betp2)*rh(im1)-an*rj(im1)+ &
    anm1*rj(i))/(anm1*(an+betp1))
    anm1 = an
    an = an+0.1D+01
    im1 = i
  end do

  do i=2,25,2
    rh(i) = -rh(i)
  end do

   70 continue

  do i=2,25,2
    rj(i) = -rj(i)
  end do

   90 continue

  return
end
subroutine dqng ( f, a, b, epsabs, epsrel, result, abserr, neval, ier )

!*****************************************************************************80
!
!! DQNG estimates an integral, using non-adaptive integration.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  the routine calculates an approximation result to a
!      given definite integral i = integral of f over (a,b),
!      hopefully satisfying following claim for accuracy
!      abs(i-result).le.max(epsabs,epsrel*abs(i)).
!
!  Parameters:
!
!     f      - real ( kind = 8 )
!              function subprogram defining the integrand function
!              f(x). the actual name for f needs to be declared
!              e x t e r n a l in the driver program.
!
!     a      - real ( kind = 8 )
!              lower limit of integration
!
!     b      - real ( kind = 8 )
!              upper limit of integration
!
!     epsabs - real ( kind = 8 )
!              absolute accuracy requested
!     epsrel - real ( kind = 8 )
!              relative accuracy requested
!              if  epsabs.le.0
!              and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
!              the routine will end with ier = 6.
!
!   on return
!     result - real ( kind = 8 )
!              approximation to the integral i
!              result is obtained by applying the 21-point
!              gauss-kronrod rule (res21) obtained by optimal
!              addition of abscissae to the 10-point gauss rule
!              (res10), or by applying the 43-point rule (res43)
!              obtained by optimal addition of abscissae to the
!              21-point gauss-kronrod rule, or by applying the
!              87-point rule (res87) obtained by optimal addition
!              of abscissae to the 43-point rule.
!
!     abserr - real ( kind = 8 )
!              estimate of the modulus of the absolute error,
!              which should equal or exceed abs(i-result)
!
!     neval  - integer ( kind = 4 )
!              number of integrand evaluations
!
!     ier    - ier = 0 normal and reliable termination of the
!                      routine. it is assumed that the requested
!                      accuracy has been achieved.
!              ier.gt.0 abnormal termination of the routine. it is
!                      assumed that the requested accuracy has
!                      not been achieved.
!     error messages
!              ier = 1 the maximum number of steps has been
!                      executed. the integral is probably too
!                      difficult to be calculated by dqng.
!                  = 6 the input is invalid, because
!                      epsabs.le.0 and
!                      epsrel.lt.max(50*rel.mach.acc.,0.5d-28).
!                      result, abserr and neval are set to zero.
!
!  Local Parameters:
!
!     the data statements contain the
!     abscissae and weights of the integration rules used.
!
!     x1      abscissae common to the 10-, 21-, 43- and 87-
!             point rule
!     x2      abscissae common to the 21-, 43- and 87-point rule
!     x3      abscissae common to the 43- and 87-point rule
!     x4      abscissae of the 87-point rule
!     w10     weights of the 10-point formula
!     w21a    weights of the 21-point formula for abscissae x1
!     w21b    weights of the 21-point formula for abscissae x2
!     w43a    weights of the 43-point formula for abscissae x1, x3
!     w43b    weights of the 43-point formula for abscissae x3
!     w87a    weights of the 87-point formula for abscissae x1,
!             x2, x3
!     w87b    weights of the 87-point formula for abscissae x4
!
!
! gauss-kronrod-patterson quadrature coefficients for use in
! quadpack routine qng.  these coefficients were calculated with
! 101 decimal digit arithmetic by l. w. fullerton, bell labs, nov 1981.
!
!     centr  - mid point of the integration interval
!     hlgth  - half-length of the integration interval
!     fcentr - function value at mid point
!     absc   - abscissa
!     fval   - function value
!     savfun - array of function values which have already been
!              computed
!     res10  - 10-point gauss result
!     res21  - 21-point kronrod result
!     res43  - 43-point result
!     res87  - 87-point result
!     resabs - approximation to the integral of abs(f)
!     resasc - approximation to the integral of abs(f-i/(b-a))
!
!     machine dependent constants
!
!     epmach is the largest relative spacing.
!     uflow is the smallest positive magnitude.
!
  implicit none

  real ( kind = 8 ) a,absc,abserr,b,centr,dhlgth, &
    epmach,epsabs,epsrel,f,fcentr,fval,fval1,fval2,fv1,fv2, &
    fv3,fv4,hlgth,result,res10,res21,res43,res87,resabs,resasc, &
    reskh,savfun,uflow,w10,w21a,w21b,w43a,w43b,w87a,w87b,x1,x2,x3,x4
  integer ( kind = 4 ) ier,ipx,k,l,neval
  external f
  dimension fv1(5),fv2(5),fv3(5),fv4(5),x1(5),x2(5),x3(11),x4(22), &
    w10(5),w21a(5),w21b(6),w43a(10),w43b(12),w87a(21),w87b(23), &
    savfun(21)

  data x1    (  1) / 0.973906528517171720077964012084452d0 /
  data x1    (  2) / 0.865063366688984510732096688423493d0 /
  data x1    (  3) / 0.679409568299024406234327365114874d0 /
  data x1    (  4) / 0.433395394129247190799265943165784d0 /
  data x1    (  5) / 0.148874338981631210884826001129720d0 /
  data w10   (  1) / 0.066671344308688137593568809893332d0 /
  data w10   (  2) / 0.149451349150580593145776339657697d0 /
  data w10   (  3) / 0.219086362515982043995534934228163d0 /
  data w10   (  4) / 0.269266719309996355091226921569469d0 /
  data w10   (  5) / 0.295524224714752870173892994651338d0 /

  data x2    (  1) / 0.995657163025808080735527280689003d0 /
  data x2    (  2) / 0.930157491355708226001207180059508d0 /
  data x2    (  3) / 0.780817726586416897063717578345042d0 /
  data x2    (  4) / 0.562757134668604683339000099272694d0 /
  data x2    (  5) / 0.294392862701460198131126603103866d0 /
  data w21a  (  1) / 0.032558162307964727478818972459390d0 /
  data w21a  (  2) / 0.075039674810919952767043140916190d0 /
  data w21a  (  3) / 0.109387158802297641899210590325805d0 /
  data w21a  (  4) / 0.134709217311473325928054001771707d0 /
  data w21a  (  5) / 0.147739104901338491374841515972068d0 /
  data w21b  (  1) / 0.011694638867371874278064396062192d0 /
  data w21b  (  2) / 0.054755896574351996031381300244580d0 /
  data w21b  (  3) / 0.093125454583697605535065465083366d0 /
  data w21b  (  4) / 0.123491976262065851077958109831074d0 /
  data w21b  (  5) / 0.142775938577060080797094273138717d0 /
  data w21b  (  6) / 0.149445554002916905664936468389821d0 /
!
  data x3    (  1) / 0.999333360901932081394099323919911d0 /
  data x3    (  2) / 0.987433402908088869795961478381209d0 /
  data x3    (  3) / 0.954807934814266299257919200290473d0 /
  data x3    (  4) / 0.900148695748328293625099494069092d0 /
  data x3    (  5) / 0.825198314983114150847066732588520d0 /
  data x3    (  6) / 0.732148388989304982612354848755461d0 /
  data x3    (  7) / 0.622847970537725238641159120344323d0 /
  data x3    (  8) / 0.499479574071056499952214885499755d0 /
  data x3    (  9) / 0.364901661346580768043989548502644d0 /
  data x3    ( 10) / 0.222254919776601296498260928066212d0 /
  data x3    ( 11) / 0.074650617461383322043914435796506d0 /
  data w43a  (  1) / 0.016296734289666564924281974617663d0 /
  data w43a  (  2) / 0.037522876120869501461613795898115d0 /
  data w43a  (  3) / 0.054694902058255442147212685465005d0 /
  data w43a  (  4) / 0.067355414609478086075553166302174d0 /
  data w43a  (  5) / 0.073870199632393953432140695251367d0 /
  data w43a  (  6) / 0.005768556059769796184184327908655d0 /
  data w43a  (  7) / 0.027371890593248842081276069289151d0 /
  data w43a  (  8) / 0.046560826910428830743339154433824d0 /
  data w43a  (  9) / 0.061744995201442564496240336030883d0 /
  data w43a  ( 10) / 0.071387267268693397768559114425516d0 /
  data w43b  (  1) / 0.001844477640212414100389106552965d0 /
  data w43b  (  2) / 0.010798689585891651740465406741293d0 /
  data w43b  (  3) / 0.021895363867795428102523123075149d0 /
  data w43b  (  4) / 0.032597463975345689443882222526137d0 /
  data w43b  (  5) / 0.042163137935191811847627924327955d0 /
  data w43b  (  6) / 0.050741939600184577780189020092084d0 /
  data w43b  (  7) / 0.058379395542619248375475369330206d0 /
  data w43b  (  8) / 0.064746404951445885544689259517511d0 /
  data w43b  (  9) / 0.069566197912356484528633315038405d0 /
  data w43b  ( 10) / 0.072824441471833208150939535192842d0 /
  data w43b  ( 11) / 0.074507751014175118273571813842889d0 /
  data w43b  ( 12) / 0.074722147517403005594425168280423d0 /

  data x4    (  1) / 0.999902977262729234490529830591582d0 /
  data x4    (  2) / 0.997989895986678745427496322365960d0 /
  data x4    (  3) / 0.992175497860687222808523352251425d0 /
  data x4    (  4) / 0.981358163572712773571916941623894d0 /
  data x4    (  5) / 0.965057623858384619128284110607926d0 /
  data x4    (  6) / 0.943167613133670596816416634507426d0 /
  data x4    (  7) / 0.915806414685507209591826430720050d0 /
  data x4    (  8) / 0.883221657771316501372117548744163d0 /
  data x4    (  9) / 0.845710748462415666605902011504855d0 /
  data x4    ( 10) / 0.803557658035230982788739474980964d0 /
  data x4    ( 11) / 0.757005730685495558328942793432020d0 /
  data x4    ( 12) / 0.706273209787321819824094274740840d0 /
  data x4    ( 13) / 0.651589466501177922534422205016736d0 /
  data x4    ( 14) / 0.593223374057961088875273770349144d0 /
  data x4    ( 15) / 0.531493605970831932285268948562671d0 /
  data x4    ( 16) / 0.466763623042022844871966781659270d0 /
  data x4    ( 17) / 0.399424847859218804732101665817923d0 /
  data x4    ( 18) / 0.329874877106188288265053371824597d0 /
  data x4    ( 19) / 0.258503559202161551802280975429025d0 /
  data x4    ( 20) / 0.185695396568346652015917141167606d0 /
  data x4    ( 21) / 0.111842213179907468172398359241362d0 /
  data x4    ( 22) / 0.037352123394619870814998165437704d0 /
  data w87a  (  1) / 0.008148377384149172900002878448190d0 /
  data w87a  (  2) / 0.018761438201562822243935059003794d0 /
  data w87a  (  3) / 0.027347451050052286161582829741283d0 /
  data w87a  (  4) / 0.033677707311637930046581056957588d0 /
  data w87a  (  5) / 0.036935099820427907614589586742499d0 /
  data w87a  (  6) / 0.002884872430211530501334156248695d0 /
  data w87a  (  7) / 0.013685946022712701888950035273128d0 /
  data w87a  (  8) / 0.023280413502888311123409291030404d0 /
  data w87a  (  9) / 0.030872497611713358675466394126442d0 /
  data w87a  ( 10) / 0.035693633639418770719351355457044d0 /
  data w87a  ( 11) / 0.000915283345202241360843392549948d0 /
  data w87a  ( 12) / 0.005399280219300471367738743391053d0 /
  data w87a  ( 13) / 0.010947679601118931134327826856808d0 /
  data w87a  ( 14) / 0.016298731696787335262665703223280d0 /
  data w87a  ( 15) / 0.021081568889203835112433060188190d0 /
  data w87a  ( 16) / 0.025370969769253827243467999831710d0 /
  data w87a  ( 17) / 0.029189697756475752501446154084920d0 /
  data w87a  ( 18) / 0.032373202467202789685788194889595d0 /
  data w87a  ( 19) / 0.034783098950365142750781997949596d0 /
  data w87a  ( 20) / 0.036412220731351787562801163687577d0 /
  data w87a  ( 21) / 0.037253875503047708539592001191226d0 /
  data w87b  (  1) / 0.000274145563762072350016527092881d0 /
  data w87b  (  2) / 0.001807124155057942948341311753254d0 /
  data w87b  (  3) / 0.004096869282759164864458070683480d0 /
  data w87b  (  4) / 0.006758290051847378699816577897424d0 /
  data w87b  (  5) / 0.009549957672201646536053581325377d0 /
  data w87b  (  6) / 0.012329447652244853694626639963780d0 /
  data w87b  (  7) / 0.015010447346388952376697286041943d0 /
  data w87b  (  8) / 0.017548967986243191099665352925900d0 /
  data w87b  (  9) / 0.019938037786440888202278192730714d0 /
  data w87b  ( 10) / 0.022194935961012286796332102959499d0 /
  data w87b  ( 11) / 0.024339147126000805470360647041454d0 /
  data w87b  ( 12) / 0.026374505414839207241503786552615d0 /
  data w87b  ( 13) / 0.028286910788771200659968002987960d0 /
  data w87b  ( 14) / 0.030052581128092695322521110347341d0 /
  data w87b  ( 15) / 0.031646751371439929404586051078883d0 /
  data w87b  ( 16) / 0.033050413419978503290785944862689d0 /
  data w87b  ( 17) / 0.034255099704226061787082821046821d0 /
  data w87b  ( 18) / 0.035262412660156681033782717998428d0 /
  data w87b  ( 19) / 0.036076989622888701185500318003895d0 /
  data w87b  ( 20) / 0.036698604498456094498018047441094d0 /
  data w87b  ( 21) / 0.037120549269832576114119958413599d0 /
  data w87b  ( 22) / 0.037334228751935040321235449094698d0 /
  data w87b  ( 23) / 0.037361073762679023410321241766599d0 /

  epmach = epsilon ( epmach )
  uflow = tiny ( uflow )
!
!  test on validity of parameters
!
  result = 0.0D+00
  abserr = 0.0D+00
  neval = 0
  ier = 6
  if(epsabs.le.0.0D+00.and.epsrel.lt. max ( 0.5D+02*epmach,0.5d-28)) &
    go to 80
  hlgth = 0.5D+00*(b-a)
  dhlgth =  abs ( hlgth)
  centr = 0.5D+00*(b+a)
  fcentr = f(centr)
  neval = 21
  ier = 1
!
!  compute the integral using the 10- and 21-point formula.
!
  do 70 l = 1,3

    go to (5,25,45),l

    5 res10 = 0.0D+00
    res21 = w21b(6)*fcentr
    resabs = w21b(6)* abs ( fcentr)

    do k=1,5
      absc = hlgth*x1(k)
      fval1 = f(centr+absc)
      fval2 = f(centr-absc)
      fval = fval1+fval2
      res10 = res10+w10(k)*fval
      res21 = res21+w21a(k)*fval
      resabs = resabs+w21a(k)*( abs ( fval1)+ abs ( fval2))
      savfun(k) = fval
      fv1(k) = fval1
      fv2(k) = fval2
    end do

    ipx = 5

    do k=1,5
      ipx = ipx+1
      absc = hlgth*x2(k)
      fval1 = f(centr+absc)
      fval2 = f(centr-absc)
      fval = fval1+fval2
      res21 = res21+w21b(k)*fval
      resabs = resabs+w21b(k)*( abs ( fval1)+ abs ( fval2))
      savfun(ipx) = fval
      fv3(k) = fval1
      fv4(k) = fval2
    end do
!
!  test for convergence.
!
    result = res21*hlgth
    resabs = resabs*dhlgth
    reskh = 0.5D+00*res21
    resasc = w21b(6)* abs ( fcentr-reskh)

    do k = 1,5
      resasc = resasc+w21a(k)*( abs ( fv1(k)-reskh)+ abs ( fv2(k)-reskh)) &
                      +w21b(k)*( abs ( fv3(k)-reskh)+ abs ( fv4(k)-reskh))
    end do

    abserr =  abs ( (res21-res10)*hlgth)
    resasc = resasc*dhlgth
    go to 65
!
!  compute the integral using the 43-point formula.
!
25  res43 = w43b(12)*fcentr
    neval = 43

    do k=1,10
      res43 = res43+savfun(k)*w43a(k)
    end do

    do k=1,11
      ipx = ipx+1
      absc = hlgth*x3(k)
      fval = f(absc+centr)+f(centr-absc)
      res43 = res43+fval*w43b(k)
      savfun(ipx) = fval
    end do
!
!  test for convergence.
!
    result = res43*hlgth
    abserr =  abs ( (res43-res21)*hlgth)
    go to 65
!
!  compute the integral using the 87-point formula.
!
45  res87 = w87b(23)*fcentr
    neval = 87

    do k=1,21
      res87 = res87+savfun(k)*w87a(k)
    end do

    do k=1,22
      absc = hlgth*x4(k)
      res87 = res87+w87b(k)*(f(absc+centr)+f(centr-absc))
    end do

    result = res87*hlgth
    abserr =  abs ( (res87-res43)*hlgth)

65  continue

    if(resasc.ne.0.0D+00.and.abserr.ne.0.0D+00) then
      abserr = resasc* min (0.1D+01,(0.2D+03*abserr/resasc)**1.5D+00)
    end if

    if (resabs.gt.uflow/(0.5D+02*epmach)) then
      abserr = max ((epmach*0.5D+02)*resabs,abserr)
    end if

    if (abserr.le. max ( epsabs,epsrel* abs ( result))) then
      ier = 0
      return
    end if

70  continue

   80 call xerror('abnormal return from dqng ',26,ier,0)
  999 continue

  return
end
subroutine dqpsrt ( limit, last, maxerr, ermax, elist, iord, nrmax )

!*****************************************************************************80
!
!! DQPSRT maintains the order of a list of local error estimates.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
!***purpose  this routine maintains the descending ordering in the
!      list of the local error estimated resulting from the
!      interval subdivision process. at each call two error
!      estimates are inserted using the sequential search
!      method, top-down for the largest error estimate and
!      bottom-up for the smallest error estimate.
!
!  Parameters:
!
!        limit  - integer ( kind = 4 )
!                 maximum number of error estimates the list
!                 can contain
!
!        last   - integer ( kind = 4 )
!                 number of error estimates currently in the list
!
!        maxerr - integer ( kind = 4 )
!                 maxerr points to the nrmax-th largest error
!                 estimate currently in the list
!
!        ermax  - real ( kind = 8 )
!                 nrmax-th largest error estimate
!                 ermax = elist(maxerr)
!
!        elist  - real ( kind = 8 )
!                 vector of dimension last containing
!                 the error estimates
!
!        iord   - integer ( kind = 4 )
!                 vector of dimension last, the first k elements
!                 of which contain pointers to the error
!                 estimates, such that
!                 elist(iord(1)),...,  elist(iord(k))
!                 form a decreasing sequence, with
!                 k = last if last.le.(limit/2+2), and
!                 k = limit+1-last otherwise
!
!        nrmax  - integer ( kind = 4 )
!                 maxerr = iord(nrmax)
!
  implicit none

  real ( kind = 8 ) elist,ermax,errmax,errmin
  integer ( kind = 4 ) i,ibeg,ido,iord,isucc,j,jbnd,jupbn,k,last, &
    lim
  integer ( kind = 4 ) limit
  integer ( kind = 4 ) maxerr
  integer ( kind = 4 ) nrmax
  dimension elist(last),iord(last)
!
!  check whether the list contains more than
!  two error estimates.
!
  if(last.gt.2) go to 10
  iord(1) = 1
  iord(2) = 2
  go to 90
!
!  this part of the routine is only executed if, due to a
!  difficult integrand, subdivision increased the error
!  estimate. in the normal case the insert procedure should
!  start after the nrmax-th largest error estimate.
!
   10 errmax = elist(maxerr)

  ido = nrmax-1
  do i = 1,ido
    isucc = iord(nrmax-1)
    if(errmax.le.elist(isucc)) go to 30
    iord(nrmax) = isucc
    nrmax = nrmax-1
  end do
!
!  compute the number of elements in the list to be maintained
!  in descending order. this number depends on the number of
!  subdivisions still allowed.
!
   30 jupbn = last
  if(last.gt.(limit/2+2)) jupbn = limit+3-last
  errmin = elist(last)
!
!  insert errmax by traversing the list top-down,
!  starting comparison from the element elist(iord(nrmax+1)).
!
  jbnd = jupbn-1
  ibeg = nrmax+1

  do i=ibeg,jbnd
    isucc = iord(i)
    if(errmax.ge.elist(isucc)) go to 60
    iord(i-1) = isucc
  end do

  iord(jbnd) = maxerr
  iord(jupbn) = last
  go to 90
!
!  insert errmin by traversing the list bottom-up.
!
   60 iord(i-1) = maxerr
  k = jbnd

  do j=i,jbnd
    isucc = iord(k)
    if(errmin.lt.elist(isucc)) go to 80
    iord(k+1) = isucc
    k = k-1
  end do

  iord(i) = last
  go to 90
   80 iord(k+1) = last
!
!     set maxerr and ermax.
!
   90 maxerr = iord(nrmax)
  ermax = elist(maxerr)

  return
end
function dqwgtc ( x, c, p2, p3, p4, kp )

!*****************************************************************************80
!
!! DQWGTC defines the weight function used by DQC25C.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
  implicit none

  real ( kind = 8 ) dqwgtc
  real ( kind = 8 ) c,p2,p3,p4,x
  integer ( kind = 4 ) kp

  dqwgtc = 0.1D+01 / ( x - c )

  return
end
function dqwgtf(x,omega,p2,p3,p4,integr)

!*****************************************************************************80
!
!! DQWGTF defines the weight functions used by DQC25F.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
  implicit none

  real ( kind = 8 ) dqwgtf
  real ( kind = 8 ) dcos,dsin,omega,omx,p2,p3,p4,x
  integer ( kind = 4 ) integr

  omx = omega * x

  if ( integr == 1 ) then
    dqwgtf = cos ( omx )
  else
    dqwgtf = sin ( omx )
  end if

  return
end
function dqwgts ( x, a, b, alfa, beta, integr )

!*****************************************************************************80
!
!! DQWGTS defines the weight functions used by DQC25S.
!
!  Modified:
!
!    11 September 2015
!
!  Author:
!
!    Robert Piessens, Elise de Doncker
!
  implicit none

  real dqwgts
  real ( kind = 8 ) a,alfa,b,beta,bmx,x,xma
  integer ( kind = 4 ) integr

  xma = x - a
  bmx = b - x
  dqwgts = xma ** alfa * bmx ** beta
  go to (40,10,20,30),integr
   10 dqwgts = dqwgts* log ( xma )
  go to 40
   20 dqwgts = dqwgts* log ( bmx )
  go to 40
   30 dqwgts = dqwgts* log ( xma ) * log ( bmx )
   40 continue

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
  ampm = 'Noon'
    else
  ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
  ampm = 'PM'
    else if ( h == 12 ) then
  if ( n == 0 .and. s == 0 ) then
    ampm = 'Midnight'
  else
    ampm = 'AM'
  end if
    end if
  end if

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine xerror ( xmess, nmess, nerr, level )

!*****************************************************************************80
!
!! XERROR replaces the SLATEC XERROR routine.
!
!  Modified:
!
!    12 September 2015
!
  implicit none

  integer ( kind = 4 ) level
  integer ( kind = 4 ) nerr
  integer ( kind = 4 ) nmess
  character ( len = * ) xmess

  if ( 1 <= LEVEL ) then
    WRITE ( *,'(1X,A)') XMESS(1:NMESS)
    WRITE ( *,'('' ERROR NUMBER = '',I5,'', MESSAGE LEVEL = '',I5)') &
        NERR,LEVEL
  end if

  return
end

