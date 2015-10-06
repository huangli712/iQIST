!!!-----------------------------------------------------------------------
!!! project : CSML (Common Service Modules Library)
!!! program : leja       module
!!! source  : m_leja.f90
!!! type    : module
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 11/20/2011 by li huang (created)
!!!           08/17/2015 by li huang (last modified)
!!! purpose : calculate the product of time evolution operator and propagated
!!!           state using Newton interpolation and real leja points method.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

  module leja
     implicit none

!-------------------------------------------------------------------------
!::: declare global parameters (integer type)                          :::
!-------------------------------------------------------------------------

! dp: number precision, double precision for reals
     integer, private, parameter :: dp = kind(1.0d0)

! mystd: device descriptor, console output
     integer, private, parameter :: mystd = 6

!-------------------------------------------------------------------------
!::: declare global parameters (real type)                             :::
!-------------------------------------------------------------------------

! 0.0 in double precision form
     real(dp), private, parameter :: zero = 0.0_dp

! 1.0 in double precision form
     real(dp), private, parameter ::  one = 1.0_dp

! 2.0 in double precision form
     real(dp), private, parameter ::  two = 2.0_dp

! 0.5 in double precision form
     real(dp), private, parameter :: half = 0.5_dp

! $\epsilon$ in double precision form
     real(dp), private, parameter :: eps6 = 1.0E-06
     real(dp), private, parameter :: eps8 = 1.0E-08
     real(dp), private, parameter :: epst = 1.0E-10
     real(dp), private, parameter :: epss = 1.0E-12

! golden section number
     real(dp), private, parameter :: gamma= 0.61803398874989490_dp

!-------------------------------------------------------------------------
!::: declare common variables for real leja points algorithm           :::
!-------------------------------------------------------------------------

! status flag
     integer, private :: istat

! dimension of h matrix
     integer, private, save  :: ncfgs = 0

! maximum allowable number of non-zero elements in h matrix
     integer, private, save  :: nhmat = 0

! number of real leja points
     integer, private, save  :: nleja = 64

! effective number of real leja points: lleja = nleja - 1
     integer, private, save  :: lleja = 63

! type of real leja points, normal (1) or symmetrized (2)
     integer, private, save  :: ltype = 1

! number of consecutive error estimates to be averaged in order to filter
! error oscillations
     integer, private, save  :: lerrs = 5

! counter, to record the invoked times of leja_dsymv()
     real(dp), private, save :: lsymv = zero

! counter, degree of newton interpolation
     real(dp), private, save :: ldeny = zero

! counter, matrix-vector multiplication
     real(dp), private, save :: lmult = zero

! extrema of the real points in the Gersghorin’s circles of h matrix
     real(dp), private, save :: lrad  = zero
     real(dp), private, save :: rrad  = zero

! approximation of the real focal interval of a "minimal" ellipse which
! contains the numerical range of the underlying h matrix, and ceta is
! the center and gamm the capacity (length/4) of the interval [lrad,rrad]
     real(dp), private, save :: ceta  = zero
     real(dp), private, save :: gamm  = zero

!-------------------------------------------------------------------------
!::: declare common arrays for real leja points algorithm              :::
!-------------------------------------------------------------------------

! real leja points \xi
     real(dp), private, allocatable, save :: lvec(:)

! divided differences for \phi function
     real(dp), private, allocatable, save :: dvec(:)

! auxiliary vector: gvec = ceta + gamm * \xi
     real(dp), private, allocatable, save :: gvec(:)

! error vector, used to estimate the erroneous boundary
     real(dp), private, allocatable, save :: verr(:)

! auxiliary matrix: lmat(i,j) = one / (\xi_i - \xi_j)
     real(dp), private, allocatable, save :: lmat(:,:)

! lop_h, lop_jh, lop_ih, the h matrix in compressed sparse row format
     integer, private, allocatable, save  :: lop_ih(:)
     integer, private, allocatable, save  :: lop_jh(:)
     real(dp), private, allocatable, save :: lop_h(:)

!-------------------------------------------------------------------------
!::: declare accessibility for module routines                         :::
!-------------------------------------------------------------------------

     public  :: leja_setup_param

     public  :: leja_setup_array
     public  :: leja_final_array

     public  :: leja_build_spmat

     public  :: leja_trace_count
     public  :: leja_reset_count

     public  :: leja_dsymv

     private :: leja_make_lejas
     private :: leja_make_gersh

     private :: getmax
     private :: getflj
     private :: csrelm
     private :: newton
     private :: divdif
     private :: dsspmv

  contains ! encapsulated functionality

!>>> setup key parameters for real leja points algorithm
  subroutine leja_setup_param(nc, nh, nl)
     implicit none

! external arguments
! dimension of h matrix
     integer, intent(in) :: nc

! maximum allowable number of non-zero elements in h matrix
     integer, intent(in) :: nh

! number of real leja points
     integer, intent(in) :: nl

! setup common variables
     ncfgs = nc
     nhmat = nh
     nleja = nl

! check nleja, when ltype == 2, nleja must be odd number
     if ( ltype == 2 .and. mod(nleja,2) == 0 ) then
         nleja = nleja - 1
     endif ! back if ( ltype == 2 .and. mod(nleja,2) == 0 ) block

! if nleja is too small, terminal code immediately
     if ( nleja < 10 ) then
         write(mystd,'(a)') 'leja: the number of leja points should be larger than 10'
         STOP
     endif ! back if ( nleja < 10 ) block

! evaluate lleja
     lleja = nleja - 1

     return
  end subroutine leja_setup_param

!>>> allocate memory for module variables
  subroutine leja_setup_array()
     implicit none

! allocate memory
     allocate(lvec(0:nleja-1), stat=istat)
     allocate(dvec(0:nleja-1), stat=istat)
     allocate(gvec(0:nleja-1), stat=istat)
     allocate(verr(0:nleja-1), stat=istat)

     allocate(lmat(0:nleja-1,0:nleja-1), stat=istat)

     allocate(lop_ih(ncfgs+1), stat=istat)
     allocate(lop_jh(nhmat),   stat=istat)
     allocate(lop_h(nhmat),    stat=istat)

! check status
     if ( istat /= 0 ) then
         write(mystd,'(a)') 'leja: can not allocate enough memory'
         STOP
     endif ! back if ( istat /= 0 ) block

! initialize them
     lvec   = zero
     dvec   = zero
     gvec   = zero
     verr   = zero

     lmat   = zero

     lop_ih = 0
     lop_jh = 0
     lop_h  = zero

     return
  end subroutine leja_setup_array

!>>> deallocate memory for module variables
  subroutine leja_final_array()
     implicit none

! deallocate memory
     if ( allocated(lvec)   ) deallocate(lvec  )
     if ( allocated(dvec)   ) deallocate(dvec  )
     if ( allocated(gvec)   ) deallocate(gvec  )
     if ( allocated(verr)   ) deallocate(verr  )

     if ( allocated(lmat)   ) deallocate(lmat  )

     if ( allocated(lop_ih) ) deallocate(lop_ih)
     if ( allocated(lop_jh) ) deallocate(lop_jh)
     if ( allocated(lop_h)  ) deallocate(lop_h )

     return
  end subroutine leja_final_array

!>>> setup h matrix, which is in sparse matrix style
  subroutine leja_build_spmat(nc, nh, Th, Tjh, Tih)
     implicit none

! external arguments
! dimension of h matrix
     integer, intent(in)  :: nc

! maximum allowable number of non-zero elements in h matrix
     integer, intent(in)  :: nh

! Th, Tjh, Tih, input matrix in compressed sparse row format
     integer, intent(in)  :: Tih(nc+1)
     integer, intent(in)  :: Tjh(nh)
     real(dp), intent(in) :: Th(nh)

! check the dimension for safety reason
     if ( nc /= ncfgs .or. nh /= nhmat ) then
         write(mystd,'(a)') 'leja: the dimension of input matrix is not correct'
         STOP
     endif ! back if ( nc /= ncfgs .or. nh /= nhmat ) block

! setup h matrix, its elements should not be changed during the newton
! interpolation process
     lop_h  = Th
     lop_jh = Tjh
     lop_ih = Tih

! determine the extrema of the real points in the Gersghorin's circles 
! of h matrix
     call leja_make_gersh()

! as a byproduct, we calulate the real leja points
     call leja_make_lejas()

     return
  end subroutine leja_build_spmat

!>>> to dump the runtime counter information for real leja points algorithm
  subroutine leja_trace_count()
     implicit none

! dump the statistics to screen
     write(mystd,'(4X,a)') 'real leja points algorithm statistics:'
     write(mystd,'(4X,a,3i10)')   'count:', int( lsymv ) , int( ldeny ) , int( lmult )
     write(mystd,'(4X,a,3f10.5)') 'ratio:', lsymv / lsymv, ldeny / lsymv, lmult / lsymv

! do not forget to reset the counter
     call leja_reset_count()

     return
  end subroutine leja_trace_count

!>>> to reset the runtime counter for real leja points algorithm
  subroutine leja_reset_count()
     implicit none

     lsymv = zero
     ldeny = zero
     lmult = zero

     return
  end subroutine leja_reset_count

!>>> to calculate w = exp( dt h ) v using newton interpolation and real
! leja points method, in which h is a symmetrix matrix
  subroutine leja_dsymv(dt, vvec, wvec)
     implicit none

! external arguments
! \delta\tau, interval for imaginary time
     real(dp), intent(in)  :: dt

! initial propagated state
     real(dp), intent(in)  :: vvec(ncfgs)

! final propagated state
     real(dp), intent(out) :: wvec(ncfgs)

! local variables
! degree for newton interpolation
     integer  :: deny

! iteration number
     integer  :: iter

! \sigma parameter used to detect fast convergence
     real(dp) :: sigm

! time-related variables for time marching scheme
! prev: previous time steps
! curr: current time steps
! rest: the rest time interval, when rest is zero, then simulation is over
     real(dp) :: prev
     real(dp) :: curr
     real(dp) :: rest

! auxiliary vectors
     real(dp) :: pvec(ncfgs)
     real(dp) :: qvec(ncfgs)

! initialize iteration number
     iter = 0

! initialize vectors
     pvec = zero
     wvec = vvec

! initialize time interval
     prev = zero
     curr = min( dt, real(lleja) / (3.0_dp * gamm) )
     rest = dt

     LEJA_TIME_LOOP: do while ( rest > zero .and. iter < 100 )

! increase iteration number
         iter = iter + 1

! when the current and the previous substeps (time interval) are different
! the divided differences (d0 , . . . , dM ) need to be (re)computed
         if ( curr /= prev ) then
             call divdif(curr); prev = curr
         endif ! back if ( curr /= prev ) block

! computes q = p_{m}(H)w \sim \phi(tH)w using newton interpolation and real
! leja points, where w = Hp + v
         call newton(deny, wvec, qvec)

! adjust new time step (time interval) and then update p and w vectors
         if ( deny >= lleja ) then
! real leja points is not enough to get convergence, time interval should
! be decreased
             curr = curr / two
         else
! update p vector: p = p + curr * q
             pvec = pvec + curr * qvec

! evaluate rest of time interval
             rest = rest - curr

! evaluate new time steps
             if ( rest > zero ) then
! update w vector: w = Hp + v
                 call dsspmv(pvec, vvec, wvec); lmult = lmult + one
                 sigm = curr * gamm / deny
                 if ( sigm > one ) then
                     curr = min(sigm * curr, rest, real(lleja) / gamm)
                 else
                     curr = min(curr, rest)
                 endif ! back if ( sigm > one ) block
             endif ! back if ( rest > zero ) block
         endif ! back if ( deny >= lleja ) block
     enddo LEJA_TIME_LOOP ! over do while loop

! final results: w = Hp + v \sim exp (dt H)v
     call dsspmv(pvec, vvec, wvec); lmult = lmult + one

! deal with the statistics counter
     lsymv = lsymv + one; ldeny = ldeny + deny

! check the iteration number
     if ( iter >= 100 ) then
         write(mystd,'(a)') 'leja: it is hard to achieve convergence'
         STOP
     endif ! back if ( iter >= 100 ) block

     return
  end subroutine leja_dsymv

!>>> evaluate real leja points
  subroutine leja_make_lejas()
     implicit none

! local variables
! loop variables
     integer :: i
     integer :: k

     real(dp) :: xmax
     real(dp) :: fmax

     real(dp) :: fmaxlc
     real(dp) :: xmaxlc

     real(dp) :: x0, x1

     real(dp) :: xsaver(0:nleja)
     real(dp) :: xorder(0:nleja)

! setup initial conditions
     if ( ltype == 1 ) then
         xsaver(0) = +two
         xsaver(1) = -two
         xsaver(2) = zero
         xorder(0) = -two
         xorder(1) = zero
         xorder(2) = +two
     else
         xsaver(0) = zero
         xsaver(1) = -two
         xsaver(2) = +two
         xorder(0) = -two
         xorder(1) = zero
         xorder(2) = +two
     endif ! back if ( ltype == 1 ) block

! generate real leja points
     do i=3,lleja,ltype
         fmax = zero
         xmax = zero
         call dlasrt('I', i, xorder, istat)
! now find the maximum over each subinterval
         do k=1,(i-1) / ltype
             x0 = xorder(k-1)
             x1 = xorder(k)
             call getmax(x0, x1, xmaxlc, fmaxlc, xorder, i)
             if ( fmaxlc >= fmax ) then
                 fmax = fmaxlc
                 xmax = xmaxlc
             endif ! back if ( fmax1c >= fmax ) block
         enddo ! over k={1,(i-1) / ltype} loop
! the next is overwritten if type = 2
         xsaver(i) = +xmax; xsaver(i+1) = -xmax
         xorder(i) = +xmax; xorder(i+1) = -xmax
     enddo ! over i={3,lleja} loop

! copy final leja points from xsaver to lvec
     do i=0,lleja
         lvec(i) = xsaver(i)
     enddo ! over i={0,lleja} loop

! build gvec vector
     do i=0,lleja
         gvec(i) = ceta + gamm * lvec(i)
     enddo ! over i={0,lleja} loop

! build lmat matrix
     do i=0,lleja
         do k=0,lleja
             if ( i /= k ) then
                 lmat(i,k) = one / ( lvec(i) - lvec(k) )
             endif ! back if ( i /= k ) block
         enddo ! over k={0,lleja} loop
     enddo ! over i={0,lleja} loop

     return
  end subroutine leja_make_lejas

!>>> estimate the numerical distribution range of eigenvalues for the 
! underlying matrix by the Gersghorin's theorem
  subroutine leja_make_gersh()
     implicit none

! local variables
! loop variables
     integer  :: i
     integer  :: j

! radius for Gersghorin's discs
     real(dp) :: disc(ncfgs)

! left and right boundaries for Gersghorin's discs
     real(dp) :: lv(ncfgs)
     real(dp) :: rv(ncfgs)

! initialize disc
     disc = zero

! evaluate disc by the Gersghorin's theorem
     do i=1,ncfgs
         do j=1,ncfgs
             if ( i /= j ) then
                 disc(i) = disc(i) + abs( csrelm(i, j)  )
             endif ! back if ( i /= j ) block
         enddo ! over j={1,ncfgs} loop
     enddo ! over i={1,ncfgs} loop

! evaluate lv and rv
! for the ith circle, csrelm(i,i) is the center, and disc(i) is the radius
     do i=1,ncfgs
         lv(i) = csrelm(i, i) - disc(i)
         rv(i) = csrelm(i, i) + disc(i)
     enddo ! over i={1,ncfgs} loop

! evaluate the extrema among lv and rv
     lrad = minval(lv); rrad = maxval(rv)

! special tricks to avoid numerical instability
     if ( abs(lrad) == abs(rrad) ) then
         lrad = lrad - one - 0.810310_dp
         rrad = rrad + one + 0.951049_dp
     endif ! back if ( abs(lrad) == abs(rrad) ) block

! evaluate ceta and gamm parameters
     ceta = ( lrad + rrad ) / two
     gamm = ( rrad - lrad ) / two / two

     return
  end subroutine leja_make_gersh

!>>> subroutine used by leja_make_lejas() subroutine
  subroutine getmax(x0, x1, xmax, fmax, x, n)
     implicit none

! external arguments
     integer, intent(in)     :: n
     real(dp), intent(inout) :: x0, x1
     real(dp), intent(inout) :: xmax
     real(dp), intent(inout) :: fmax
     real(dp), intent(in)    :: x(0:n-1)

! local variables
     real(dp) :: dist
     real(dp) :: a,  b
     real(dp) :: fa, fb

     dist = gamma * ( x1 - x0 )
     a = x1 - dist; b = x0 + dist
     fa = getflj(n, x, a); fb = getflj(n, x, b)

     do while ( x1 - x0 >= epss )
         dist = gamma * dist
         if ( fa <= fb ) then
             x0 = a
             a = b; b = x0 + dist
             fa = fb; fb = getflj(n, x, b)
         else
             x1 = b
             b = a; a = x1 - dist
             fb = fa; fa = getflj(n, x, a)
         endif ! back if ( fa <= fb ) block
     enddo ! over do while loop

     fmax = max(fa, fb)
     xmax = ( a + b + sign(fa - fb, one) * ( a - b ) ) / two

     return
  end subroutine getmax

!>>> function used by leja_make_lejas() subroutine
  real(dp) &
  function getflj(n, x, t) result(value)
     implicit none

! external arguments
     integer, intent(in)  :: n
     real(dp), intent(in) :: t
     real(dp), intent(in) :: x(0:n-1)

     value = abs( product( x - t ) )

     return
  end function getflj

!>>> this function returns the element H(i,j) of matrix H
  real(dp) &
  function csrelm(i, j) result(elm)
     implicit none

! external arguments
! the row index of the element sought
     integer, intent(in) :: i

! the column index of the element sought
     integer, intent(in) :: j

! local variables
! loop index
     integer :: k

! memory address of H(i,j)
     integer :: addr

! initialization
     addr = 0
     elm = zero

! scan the row - exit as soon as a(i,j) is found
     do k=lop_ih(i),lop_ih(i+1)-1
         if ( lop_jh(k) == j ) then
             addr = k; EXIT
         endif ! back if ( lop_jh(k) == j ) block
     enddo ! over k={lop_ih(i),lop_ih(i+1)-1} loop

! the required element is contained in sparse matrix
     if ( addr /= 0 ) then
         elm = lop_h(addr)
     endif ! back if ( addr /= 0 ) block

     return
  end function csrelm

!>>> computes q = p_{m}(H)w using newton interpolation and real leja points
  subroutine newton(deny, wvec, qvec)
     implicit none

! external arguments
! degree for newton interpolation
     integer, intent(out)  :: deny

! initial propagated state
     real(dp), intent(in)  :: wvec(ncfgs)

! final propagated state
     real(dp), intent(out) :: qvec(ncfgs)

! local variables
! used to estimate numerical error
     real(dp) :: beta_nrm
     real(dp) :: curr_err

! dummy vectors
     real(dp) :: uvec(ncfgs)
     real(dp) :: zvec(ncfgs)

     deny = 0

     uvec = wvec
     zvec = zero
     qvec = dvec(0) * wvec

     beta_nrm = sqrt( dot_product(wvec, wvec) )
     curr_err = beta_nrm
     verr (0) = beta_nrm * dvec(0)
     do while ( curr_err > beta_nrm * eps6 .and. deny < lleja )
         call dsspmv(one/gamm * uvec, -gvec(deny)/gamm * uvec, uvec); lmult = lmult + one
         deny = deny + 1
         verr(deny) = abs( dvec(deny) ) * sqrt( dot_product(uvec, uvec) )
         qvec = qvec + dvec(deny) * uvec
         if ( deny >= lerrs - 1 ) then
             curr_err = sum( verr(deny - lerrs + 1:deny) ) / real(lerrs)
         else
             curr_err = sum( verr(0:deny) ) / real(deny + 1)
         endif ! back if ( deny > = lerrs - 1 ) block
     enddo ! over do while loop

     return
  end subroutine newton

!>>> to return the divided differences for the newton interpolation of
! \phi(tH) up to degree M
  subroutine divdif(curr)
     implicit none

! external arguments
! current time steps
     real(dp), intent(in) :: curr

! local variables
! loop index
     integer :: i
     integer :: j

! dummy variables
     real(dp) :: x

     do i=0,lleja
         x = curr * gvec(i)
         dvec(i) = ( exp( x ) - one ) / x
         if ( abs(x) < epss ) dvec(i) = one
     enddo ! over i={0,lleja} loop

     do i=0,lleja
         do j=1,i
             dvec(i) = ( dvec(i) - dvec(j-1) ) * lmat(i,j-1)
         enddo ! over j={1,i} loop
     enddo ! over i={0,lleja} loop

     return
  end subroutine divdif

!>>> perform sparse matrix-vector operation: v = H x + y
  subroutine dsspmv(x, y, v)
     implicit none

! external arguments
! input vector
     real(dp), intent(in)  :: x(ncfgs)
     real(dp), intent(in)  :: y(ncfgs)

! output vector
     real(dp), intent(out) :: v(ncfgs)

! local variables
! loop index
     integer :: i
     integer :: j

     v = y
     do i=1,ncfgs
         do j=lop_ih(i),lop_ih(i+1)-1
             v(i) = v(i) + lop_h(j) * x( lop_jh(j) )
         enddo ! over j={lop_ih(i),lop_ih(i+1)-1} loop
     enddo ! over i={1,ncfgs} loop

     return
  end subroutine dsspmv

  end module leja
