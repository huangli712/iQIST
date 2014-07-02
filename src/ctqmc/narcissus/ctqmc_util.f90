!-------------------------------------------------------------------------
! project : narcissus
! program : ctqmc_dmat_inv
!           ctqmc_zmat_inv
!           ctqmc_dmat_det
!           ctqmc_zmat_det
!           ctqmc_make_uumat
!           ctqmc_make_state
!           ctqmc_make_sbess
!           ctqmc_time_sorter
!           ctqmc_time_qsorter
!           ctqmc_time_qscorer
!           ctqmc_time_builder
!           ctqmc_time_analyzer
! source  : ctqmc_util.f90
! type    : functions & subroutines
! author  : li huang (email:huangli712@yahoo.com.cn)
! history : 10/01/2008 by li huang
!           02/08/2009 by li huang
!           09/23/2009 by li huang
!           09/26/2009 by li huang
!           11/17/2009 by li huang
!           11/21/2009 by li huang
!           12/18/2009 by li huang
!           12/22/2009 by li huang
!           12/29/2009 by li huang
!           01/12/2010 by li huang
!           02/27/2010 by li huang
!           06/08/2010 by li huang
!           06/22/2010 by li huang
! purpose : to provide utility functions and subroutines for hybridization
!           expansion version continuous time quantum Monte Carlo (CTQMC)
!           quantum impurity solver
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

# define prefix '>>> used time: '

# define iolst1 prefix , mday, ' d ', mhou, ' h ', mmin, ' m in this iteration.'
# define iolst2 prefix , mhou, ' h ', mmin, ' m in this iteration.'
# define iolst3 prefix , mmin, ' m ', msec, ' s in this iteration.'
# define iolst4 prefix , msec, ' s in this iteration.'

# define iolst5 prefix , nday, ' d ', nhou, ' h ', nmin, ' m in total iteration.'
# define iolst6 prefix , nhou, ' h ', nmin, ' m in total iteration.'
# define iolst7 prefix , nmin, ' m ', nsec, ' s in total iteration.'
# define iolst8 prefix , nsec, ' s in total iteration.'

!>>> invert real(dp) matrix using lapack subroutines
  subroutine ctqmc_dmat_inv(ndim, dmat)
     use constants, only : dp

     implicit none

! external arguments
! dimension of dmat matrix
     integer, intent(in) :: ndim

! object matrix, on entry, it contains the original matrix, on exit,
! it is destroyed and replaced with the inversed matrix
     real(dp), intent(inout) :: dmat(ndim,ndim)

! local variables
! error flag
     integer  :: ierror

! working arrays for lapack subroutines
     integer  :: ipiv(ndim)
     real(dp) :: work(ndim)

! computes the LU factorization of a general m-by-n matrix, need lapack
! package, dgetrf subroutine
     call dgetrf(ndim, ndim, dmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         call ctqmc_print_error('ctqmc_dmat_inv','error in lapack subroutine dgetrf')
     endif

! computes the inverse of an LU-factored general matrix, need lapack
! package, dgetri subroutine
     call dgetri(ndim, dmat, ndim, ipiv, work, ndim, ierror)
     if ( ierror /= 0 ) then
         call ctqmc_print_error('ctqmc_dmat_inv','error in lapack subroutine dgetri')
     endif

     return
  end subroutine ctqmc_dmat_inv

!>>> invert complex(dp) matrix using lapack subroutines
  subroutine ctqmc_zmat_inv(ndim, zmat)
     use constants, only : dp

     implicit none

! external arguments
! dimension of zmat matrix
     integer, intent(in) :: ndim

! object matrix, on entry, it contains the original matrix, on exit,
! it is destroyed and replaced with the inversed matrix
     complex(dp), intent(inout) :: zmat(ndim,ndim)

! local variables
! error flag
     integer     :: ierror

! working arrays for lapack subroutines
     integer     :: ipiv(ndim)
     complex(dp) :: work(ndim)

! computes the LU factorization of a general m-by-n matrix, need lapack
! package, zgetrf subroutine
     call zgetrf(ndim, ndim, zmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         call ctqmc_print_error('ctqmc_zmat_inv','error in lapack subroutine zgetrf')
     endif

! computes the inverse of an LU-factored general matrix, need lapack
! package, zgetri subroutine
     call zgetri(ndim, zmat, ndim, ipiv, work, ndim, ierror)
     if ( ierror /= 0 ) then
         call ctqmc_print_error('ctqmc_zmat_inv','error in lapack subroutine zgetri')
     endif

     return
  end subroutine ctqmc_zmat_inv

!>>> calculate the determinant of a real(dp) matrix
  subroutine ctqmc_dmat_det(ndim, dmat, ddet)
     use constants, only : dp, one, cone

     implicit none

! external arguments
! dimension of dmat matrix
     integer, intent(in) :: ndim

! determinant of dmat matrix
     real(dp), intent(out) :: ddet

! object matrix, on entry, it contains the original matrix, on exit,
! it is destroyed and replaced with the L and U matrix
     real(dp), intent(inout) :: dmat(ndim,ndim)

! local variables
! loop index
     integer  :: i

! error flag
     integer  :: ierror

! size of working array work
     integer  :: lwork

! working arrays for lapack subroutines: dgetrf
     integer  :: ipiv(ndim)

! working arrays for lapack subroutines: dgeev
     real(dp) :: work(4*ndim)

! real and imaginary parts of the computed eigenvalues
     real(dp) :: wi(ndim)
     real(dp) :: wr(ndim)

! left and right eigenvectors
     real(dp) :: vl(ndim,ndim)
     real(dp) :: vr(ndim,ndim)

! dummy arrays, used to save dmat
     real(dp) :: amat(ndim,ndim)

! used to calculate determinant
     complex(dp) :: cres

! setup lwork
     lwork = 4 * ndim

! copy dmat to amat at first
     amat = dmat

!-------------------------------------------------------------------------
! method A: preferred method
!-------------------------------------------------------------------------
! computes the LU factorization of a general m-by-n matrix, need lapack
! package, dgetrf subroutine
     call dgetrf(ndim, ndim, dmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         call ctqmc_print_exception('ctqmc_dmat_det','error in lapack subroutine dgetrf')
     endif

! calculate determinant
     ddet = one
     do i=1,ndim
         if ( ipiv(i) == i ) then
             ddet = ddet * ( +dmat(i,i) )
         else
             ddet = ddet * ( -dmat(i,i) )
         endif
     enddo ! over i={1,ndim} loop

! everything is ok!
     if ( ierror == 0 ) RETURN

!-------------------------------------------------------------------------
! method B: as a backup
!-------------------------------------------------------------------------
! diagonalize amat to obtain its eigenvalues: wr and wi
     call dgeev('N', 'N', ndim, amat, ndim, wr, wi, vl, ndim, vr, ndim, work, lwork, ierror)
     if ( ierror /= 0 ) then
         call ctqmc_print_error('ctqmc_dmat_det','error in lapack subroutine dgeev')
     endif

! evaluate the final determinant
     cres = cone
     do i=1,ndim
         cres = cres * dcmplx( wr(i), wi(i) )
     enddo ! over i={1,ndim} loop
     ddet = cres

     return
  end subroutine ctqmc_dmat_det

!>>> calculate the determinant of a complex(dp) matrix
  subroutine ctqmc_zmat_det(ndim, zmat, zdet)
     use constants, only : dp, cone

     implicit none

! external arguments
! dimension of zmat matrix
     integer, intent(in) :: ndim

! determinant of zmat matrix
     real(dp), intent(out) :: zdet

! object matrix, on entry, it contains the original matrix, on exit,
! it is destroyed and replaced with the L and U matrix
     real(dp), intent(inout) :: zmat(ndim,ndim)

! local variables
! loop index
     integer :: i

! error flag
     integer :: ierror

! working arrays for lapack subroutines
     integer :: ipiv(ndim)

! computes the LU factorization of a general m-by-n matrix, need lapack
! package, zgetrf subroutine
     call zgetrf(ndim, ndim, zmat, ndim, ipiv, ierror)
     if ( ierror /= 0 ) then
         call ctqmc_print_error('ctqmc_zmat_det','error in lapack subroutine zgetrf')
     endif

! calculate determinant
     zdet = cone
     do i=1,ndim
         if ( ipiv(i) == i ) then
             zdet = zdet * ( +zmat(i,i) )
         else
             zdet = zdet * ( -zmat(i,i) )
         endif
     enddo ! over i={1,ndim} loop

     return
  end subroutine ctqmc_zmat_det

!>>> to build general U interaction matrix: uumat, using my own style
! note: do not support spin-flip and pair-hopping term so far
! note: only Uc and Jz are need, the other Coulomb interaction parameters
! are used as backup
  subroutine ctqmc_make_uumat(uumat)
     use constants
     use control

     implicit none

! external arguments
! Coulomb interaction matrix
     real(dp), intent(out) :: uumat(norbs, norbs)

! local variables
! loop index
     integer  :: i
     integer  :: j

     integer  :: k
     integer  :: m

! control flag, whether the shift for chemical potential is loaded already
     integer, save :: touch = 0

! Coulomb interaction shift introduced by dynamical screening effect
     real(dp) :: shift

! dummy u vector
     real(dp) :: ut(nband*(norbs-1))

! evaluate Coulomb interaction shift
     select case ( isscr )

         case (1)
             shift = zero               ! normal model, recover azalea

         case (2)
             shift = two * lc * lc / wc ! holstein-hubbard model

         case (3)
             shift = two * lc * lc / wc ! plasmon pole model

         case (4)
             shift = two * lc * wc      ! ohmic model

         case (99)
             shift = lc                 ! realistic materials

     end select

! initialize it
     uumat = zero

! calculate it
     k = 0
     do i=1,norbs-1
         do j=i+1,norbs
             k = k + 1
             if ( i <= nband .and. j > nband ) then
                 m = j - nband
                 if ( m == i ) then
                     ut(k) = Uc
                 else
                     ut(k) = Uc - 2.0_dp * Jz
                 endif
             else
                 ut(k) = Uc - 3.0_dp * Jz
             endif

             uumat(i,j) = ut(k) - shift
             uumat(j,i) = ut(k) - shift
         enddo ! over j={i+1,norbs} loop
     enddo ! over i={1,norbs-1} loop

! shift chemical potential as a byproduct
     if ( touch == 0 ) then
         mune = mune - shift / two
         touch = 1 ! shut down the chemical potential shift
     endif

     return
  end subroutine ctqmc_make_uumat

!>>> convert current atomic state array into a decimal number (state index)
  subroutine ctqmc_make_state(norbs, pstat, state)
     implicit none

! external arguments
! index of atomic state
     integer, intent(out) :: pstat

! number of orbitals
     integer, intent(in)  :: norbs

! atomic state array
     integer, intent(in)  :: state(norbs)

! local variables
! loop index
     integer :: i

! init pstat
     pstat = 1

! evaluate pstat, for example, 0101 = 0*2^0 + 1*2^1 + 0*2^2 + 1*2^3 = 10
     do i=1,norbs
         if ( state(i) > 0 ) pstat = pstat + ishft(1, i-1)
     enddo ! over i={1,norbs} loop

     return
  end subroutine ctqmc_make_state

!>>> computes the spherical Bessel functions of the first kind, j_l(x),
! for argument x and l=0,1,\ldots,l_{max}. The recursion relation:
!     j_{l+1}(x)=\frac{2l+1}{x}j_l(x)-j_{l-1}(x)
! is used either downwards for x < l or upwards for x >= l. for x << 1,
! the asymtotic form is used:
!     j_l(x) \approx \frac{x^l}{(2l+1)!!}
! this procedure is numerically stable and accurate to near this machine
! precision for l <= 50
  subroutine ctqmc_make_sbess(lmax, x, jl)
     use constants

     implicit none

! external arguments
! maximum order of spherical Bessel function
     integer, intent(in)   :: lmax

! real argument
     real(dp), intent(in)  :: x

! array of returned values
     real(dp), intent(out) :: jl(0:lmax)

! local parameters
! staring value for l above lmax (suitable for lmax < 50)
     integer, parameter  :: lst  = 25

! rescale limit
     real(dp), parameter :: rsc  = 1.0D100
     real(dp), parameter :: rsci = one / rsc

! local variables
! loop index
     integer  :: l

! real(dp) dummy variables
     real(dp) :: xi, jt
     real(dp) :: j0, j1
     real(dp) :: t1, t2

! check the range of input variables
     if ( lmax < 0 .or. lmax > 50 ) then
         call ctqmc_print_error('ctqmc_make_sbess','lmax is out of range')
     endif

     if ( x < zero .or. x > 1.0E5 ) then
         call ctqmc_print_error('ctqmc_make_sbess','x is out of range')
     endif

! treat x << 1
     if ( x < eps8 ) then
         jl(0) = one
         t1 = one; t2 = one
         do l=1,lmax
             t1 = t1 / (two * l + one)
             t2 = t2 * x
             jl(l) = t2 * t1
         enddo ! over l={1,lmax} loop
         RETURN
     endif ! back if ( x < eps8 ) block

     xi = one / x

! for x < lmax recurse down
     if ( x < lmax ) then
         if ( lmax == 0 ) then
             jl(0) = sin(x) / x
             RETURN
         endif ! back if ( lmax == 0 ) block

! start from truly random numbers
         j0 = 0.6370354636449841609d0 * rsci
         j1 = 0.3532702964695481204d0 * rsci
         do l=lmax+lst,lmax+1,-1
             jt = j0 * (two * l + one) * xi - j1
             j1 = j0
             j0 = jt
! check for overflow
             if ( abs(j0) > rsc ) then
! rescale
                 jt = jt * rsci
                 j1 = j1 * rsci
                 j0 = j0 * rsci
             endif ! back if ( abs(j0) > rsc ) block
         enddo ! over l={lmax+lst,lmax+1} loop

         do l=lmax,0,-1
             jt = j0 * (two * l + one) * xi - j1
             j1 = j0
             j0 = jt
! check for overflow
             if ( abs(j0) > rsc ) then
! rescale
                 jt = jt * rsci
                 j1 = j1 * rsci
                 j0 = j0 * rsci
                 jl(l+1:lmax) = jl(l+1:lmax) * rsci
             endif ! back if ( abs(j0) > rsc ) block
             jl(l) = j1
         enddo ! over l={lmax,0} loop
! rescaling constant
         t1 = one / ( ( jl(0) - x * jl(1) ) * cos(x) + x * jl(0) * sin(x) )
         jl = t1 * jl
     else
! for large x recurse up
         jl(0) = sin(x) * xi
         if ( lmax == 0 ) RETURN
         jl(1) = ( jl(0) - cos(x) ) * xi
         if ( lmax == 1 ) RETURN
         j0 = jl(0)
         j1 = jl(1)
         do l=2,lmax
             jt = (two * l - one ) * j1 * xi - j0
             j0 = j1
             j1 = jt
             jl(l) = j1
         enddo ! over l={2,lmax} loop
     endif ! back if ( x < lmax ) block

     return
  end subroutine ctqmc_make_sbess

!>>> using bubble sort algorithm to sort a real dataset, the slowest algorithm
  subroutine ctqmc_time_sorter(nsize, list)
     use constants, only : dp

     implicit none

! external arguments
! grab the number of values from the calling code
     integer, intent(in) :: nsize

! dataset to be sorted
     real(dp), intent(inout) :: list(nsize)

! local variables
! dataset index
     integer  :: i = 0
     integer  :: j = 0

! dummy variables
     real(dp) :: swap

! basically we just loop through every element to compare it against
! every other element
! this loop increments i which is our starting point for the comparison
     sort_loop1: do i=nsize,1,-1
! this loop increments j which is the ending point for the comparison
         sort_loop2: do j=1,i-1
! swap the two elements here
             exchange: if ( list(j) > list(j+1) ) then
                 swap = list(j)
                 list(j) = list(j+1)
                 list(j+1) = swap
             endif exchange
         enddo sort_loop2 ! over j={1,i-1} loop
     enddo sort_loop1 ! over i={nsize,1,-1} loop

     return
  end subroutine ctqmc_time_sorter

!>>> sets up for the quick sort recursive method
  subroutine ctqmc_time_qsorter(nsize, list)
     use constants, only : dp

     implicit none

! external arguments
! grab the number of values from the calling code
     integer, intent(in) :: nsize

! dataset to be sorted
     real(dp), intent(inout) :: list(nsize)

! kicks off the recursive process
     call ctqmc_time_qscorer(1, nsize, nsize, list)

     return
  end subroutine ctqmc_time_qsorter

!>>> this is the actually recursive portion of the quicksort algorithm
  recursive &
  subroutine ctqmc_time_qscorer(pstart, pend, nsize, list)
     use constants, only : dp

     implicit none

! external arguments
! start point
     integer, intent(in) :: pstart

! end point
     integer, intent(in) :: pend

! size of array
     integer, intent(in) :: nsize

! dataset to be sorted
     real(dp), intent(inout) :: list(nsize)

! local variables
! used to find out list(left) > kaux and list(right) < kaux
     integer  :: left, right

! used to record list(pstart)
     real(dp) :: kaux

! used to swap data
     real(dp) :: taux

! setup left and right
     left = pstart
     right = pend + 1

! only in right > left, the data is to be sorted
     if ( right > left ) then

! record list(pstart) at first
         kaux = list(pstart)

         do while ( .true. )

! find out where list(left) < kaux
             do while ( .true. )
                 left = left + 1
                 if ( list(left)  > kaux .or. left  >= pend   ) EXIT
             enddo ! over do while loop

! find out where list(right) > kaux
             do while ( .true. )
                 right = right - 1
                 if ( list(right) < kaux .or. right <= pstart ) EXIT
             enddo ! over do while loop

! we should ensure right is larger than left
             if ( right <= left ) EXIT

! exchange data between list(left) and list(right)
             taux = list(left)
             list(left) = list(right)
             list(right) = taux

         enddo ! over do while loop

! exchange data between list(pstart) and list(right)
        list(pstart) = list(right)
        list(right) = kaux

! sort data from pstart to right-1
        call ctqmc_time_qscorer(pstart, right-1, nsize, list)

! sort data from right+1 to pend
        call ctqmc_time_qscorer(right+1, pend, nsize, list)

     endif ! back if ( right > left ) block

     return
  end subroutine ctqmc_time_qscorer

!>>> returns a string containing date and time in human-readable format
  subroutine ctqmc_time_builder(date_time_string)
     implicit none

! external arguments
! output date and time
     character (len = 20), intent(out) :: date_time_string

! local variables
! used to extract data from a standard fortran call: date_and_time()
     integer :: date_time(8)

! string for current date
     character (len = 12) :: cdate

! string for current time
     character (len = 08) :: ctime

! month array
     character (len = 03) :: months(12)

! init the month array
     months( 1) = 'Jan'; months( 2) = 'Feb'; months( 3) = 'Mar'
     months( 4) = 'Apr'; months( 5) = 'May'; months( 6) = 'Jun'
     months( 7) = 'Jul'; months( 8) = 'Aug'; months( 9) = 'Sep'
     months(10) = 'Oct'; months(11) = 'Nov'; months(12) = 'Dec'

! obtain current date and time
     call date_and_time(values = date_time)

! convert date and time from integer to string
     write (cdate,'(1X,a3,1X,i2,1X,i4)') months(date_time(2)), date_time(3), date_time(1)
     write (ctime, '(i2,":",i2,":",i2)') date_time(5), date_time(6), date_time(7)

! build final output string by concating them
     date_time_string = ctime // cdate

     return
  end subroutine ctqmc_time_builder

!>>> used to print the iteration timing information about continuous time
! quantum Monte Carlo quantum impurity solver.
  subroutine ctqmc_time_analyzer(time_iter, time_niter)
     use constants

     implicit none

! external arguments
! time used in this iteration
     real(dp), intent(in) :: time_iter

! time used in total iteration
     real(dp), intent(in) :: time_niter

! local variables
! number of days
     integer  :: mday, nday

! number of hours
     integer  :: mhou, nhou

! number of minutes
     integer  :: mmin, nmin

! number of seconds
     real(dp) :: msec, nsec

     mday = time_iter / 86400
     msec = time_iter - 86400 * mday
     mhou = msec / 3600
     msec = msec - 3600 * mhou
     mmin = msec / 60
     msec = msec - 60 * mmin

     nday = time_niter / 86400
     nsec = time_niter - 86400 * nday
     nhou = nsec / 3600
     nsec = nsec - 3600 * nhou
     nmin = nsec / 60
     nsec = nsec - 60 * nmin

     if      ( mday > 0 ) then
         write(mystd,'(4X,3(a,i2),a)')     iolst1

     else if ( mhou > 0 ) then
         write(mystd,'(4X,2(a,i2),a)')     iolst2

     else if ( mmin > 0 ) then
         write(mystd,'(4X,a,i2,a,f5.2,a)') iolst3

     else
         write(mystd,'(4X,a,f5.2,a)')      iolst4

     endif ! back if ( mday > 0 ) block

     if      ( nday > 0 ) then
         write(mystd,'(4X,3(a,i2),a)')     iolst5

     else if ( nhou > 0 ) then
         write(mystd,'(4X,2(a,i2),a)')     iolst6

     else if ( nmin > 0 ) then
         write(mystd,'(4X,a,i2,a,f5.2,a)') iolst7

     else
         write(mystd,'(4X,a,f5.2,a)')      iolst8

     endif ! back if ( nday > 0 ) block

     return
  end subroutine ctqmc_time_analyzer
