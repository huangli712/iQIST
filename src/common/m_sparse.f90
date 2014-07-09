!!!-----------------------------------------------------------------------
!!! project : CSML (Common Service Modules Library)
!!! program : sparse
!!!           sparse@sparse_format_csrdns
!!!           sparse@sparse_format_csrdns_
!!!           sparse@sparse_format_dnscsr
!!!           sparse@sparse_format_dnscsr_
!!!           sparse@sparse_format_unicsr
!!!           sparse@sparse_format_unicsr_
!!!           sparse@sparse_matrix_copyer
!!!           sparse@sparse_matrix_copyer_
!!!           sparse@sparse_matrix_getter
!!!           sparse@sparse_matrix_getter_
!!!           sparse@sparse_matmul_amuvec
!!!           sparse@sparse_matmul_amuvec_
!!!           sparse@sparse_matmul_amumat
!!!           sparse@sparse_matmul_amumat_
!!!           sparse@sparse_matmul_amudia
!!!           sparse@sparse_matmul_amudia_
!!!           sparse@sparse_matmul_diamua
!!!           sparse@sparse_matmul_diamua_
!!! source  : m_sparse.f90
!!! type    : module
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 02/01/2010 by li huang
!!!           07/09/2014 by li huang
!!! purpose : the purpose of this module is to implement important sparse
!!!           matrix/vector operations, including matrix multiplication,
!!!           format convertion, etc. the internal format of sparse matrix
!!!           used in this module is CSR (compressed sparse row) format.
!!! status  : unstable
!!! comment : only support real(dp) and complex(dp) data types
!!!-----------------------------------------------------------------------

  module sparse
     implicit none

!-------------------------------------------------------------------------
!::: declare global parameters                                         :::
!-------------------------------------------------------------------------

! dp: number precision, double precision for real and complex number
     integer, private, parameter :: dp = kind(1.0d0)

! mystd: device descriptor, console output
     integer, private, parameter :: mystd = 6

!-------------------------------------------------------------------------
!::: declare accessibility for module routines                         :::
!-------------------------------------------------------------------------

! CSR -> DNS
     private :: sparse_format_csrdns  ! real(dp) version
     private :: sparse_format_csrdns_ ! complex(dp) version

! DNS -> CSR
     private :: sparse_format_dnscsr  ! real(dp) version
     private :: sparse_format_dnscsr_ ! complex(dp) version

! UNI -> CSR
     private :: sparse_format_unicsr  ! real(dp) version
     private :: sparse_format_unicsr_ ! complex(dp) version

! CSR -> CSR
     private :: sparse_matrix_copyer  ! real(dp) version
     private :: sparse_matrix_copyer_ ! complex(dp) version

! CSR -> extractor
     private :: sparse_matrix_getter  ! real(dp) version
     private :: sparse_matrix_getter_ ! complex(dp) version

! CSR X VEC
     private :: sparse_matmul_amuvec  ! real(dp) version
     private :: sparse_matmul_amuvec_ ! complex(dp) version

! CSR X CSR
     private :: sparse_matmul_amumat  ! real(dp) version
     private :: sparse_matmul_amumat_ ! complex(dp) version

! CSR X DIA
     private :: sparse_matmul_amudia  ! real(dp) version
     private :: sparse_matmul_amudia_ ! complex(dp) version

! DIA X CSR
     private :: sparse_matmul_diamua  ! real(dp) version
     private :: sparse_matmul_diamua_ ! complex(dp) version

!-------------------------------------------------------------------------
!::: declare interface and module procedure                            :::
!-------------------------------------------------------------------------

     public :: sparse_csr_to_dns
     interface sparse_csr_to_dns
         module procedure sparse_format_csrdns
         module procedure sparse_format_csrdns_
     end interface sparse_csr_to_dns

     public :: sparse_dns_to_csr
     interface sparse_dns_to_csr
         module procedure sparse_format_dnscsr
         module procedure sparse_format_dnscsr_
     end interface sparse_dns_to_csr

     public :: sparse_uni_to_csr
     interface sparse_uni_to_csr
         module procedure sparse_format_unicsr
         module procedure sparse_format_unicsr_
     end interface sparse_uni_to_csr

     public :: sparse_csr_cp_csr
     interface sparse_csr_cp_csr
         module procedure sparse_matrix_copyer
         module procedure sparse_matrix_copyer_
     end interface sparse_csr_cp_csr

     public :: sparse_csr_cp_elm
     interface sparse_csr_cp_elm
         module procedure sparse_matrix_getter
         module procedure sparse_matrix_getter_
     end interface sparse_csr_cp_elm

     public :: sparse_csr_mv_vec
     interface sparse_csr_mv_vec
         module procedure sparse_matmul_amuvec
         module procedure sparse_matmul_amuvec_
     end interface sparse_csr_mv_vec

     public :: sparse_csr_mm_csr
     interface sparse_csr_mm_csr
         module procedure sparse_matmul_amumat
         module procedure sparse_matmul_amumat_
     end interface sparse_csr_mm_csr

     public :: sparse_csr_mm_dia
     interface sparse_csr_mm_dia
         module procedure sparse_matmul_amudia
         module procedure sparse_matmul_amudia_
     end interface sparse_csr_mm_dia

     public :: sparse_dia_mm_csr
     interface sparse_dia_mm_csr
         module procedure sparse_matmul_diamua
         module procedure sparse_matmul_diamua_
     end interface sparse_dia_mm_csr

  contains ! encapsulated functionality

!>>> converts a row-stored sparse matrix into a densely stored one
  subroutine sparse_format_csrdns(nrow, ncol, nmax, a, ja, ia, dns)
     implicit none

! external arguments
! row dimension of dense matrix
     integer, intent(in)   :: nrow

! column dimension of dense matrix
     integer, intent(in)   :: ncol

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays a and ja
     integer, intent(in)   :: nmax

! a, ja, ia, input matrix in compressed sparse row format
     integer, intent(in)   :: ia(nrow+1)
     integer, intent(in)   :: ja(nmax)
     real(dp), intent(in)  :: a(nmax)

! array where to store dense matrix
     real(dp), intent(out) :: dns(nrow,ncol)

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k

! init dns matrix
     dns = 0.0_dp

! convert sparse matrix to dense matrix
     do i=1,nrow
         do k=ia(i),ia(i+1)-1
             j = ja(k)
             if ( j > ncol ) then
                 write(mystd,'(a)') 'sparse: error in sparse_format_csrdns'
                 STOP
             endif ! back if ( j > ncol ) block
             dns(i,j) = a(k)
         enddo ! over k={ia(i),ia(i+1)-1} loop
     enddo ! over i={1,nrow} loop

     return
  end subroutine sparse_format_csrdns

!>>> converts a row-stored sparse matrix into a densely stored one
  subroutine sparse_format_csrdns_(nrow, ncol, nmax, sa, ja, ia, dns)
     implicit none

! external arguments
! row dimension of dense matrix
     integer, intent(in)      :: nrow

! column dimension of dense matrix
     integer, intent(in)      :: ncol

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays sa and ja
     integer, intent(in)      :: nmax

! sa, ja, ia, input matrix in compressed sparse row format
     integer, intent(in)      :: ia(nrow+1)
     integer, intent(in)      :: ja(nmax)
     complex(dp), intent(in)  :: sa(nmax)

! array where to store dense matrix
     complex(dp), intent(out) :: dns(nrow,ncol)

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k

! init dns matrix
     dns = dcmplx(0.0_dp, 0.0_dp)

! convert sparse matrix to dense matrix
     do i=1,nrow
         do k=ia(i),ia(i+1)-1
             j = ja(k)
             if ( j > ncol ) then
                 write(mystd,'(a)') 'sparse: error in sparse_format_csrdns_'
                 STOP
             endif ! back if ( j > ncol ) block
             dns(i,j) = sa(k)
         enddo ! over k={ia(i),ia(i+1)-1} loop
     enddo ! over i={1,nrow} loop

     return
  end subroutine sparse_format_csrdns_

!>>> converts a densely stored matrix into a row orientied compactly sparse matrix
  subroutine sparse_format_dnscsr(nrow, ncol, nmax, dns, a, ja, ia)
     implicit none

! external arguments
! row dimension of dense matrix
     integer, intent(in)   :: nrow

! column dimension of dense matrix
     integer, intent(in)   :: ncol

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays a and ja
     integer, intent(in)   :: nmax

! input densely stored matrix
     real(dp), intent(in)  :: dns(nrow,ncol)

! a, ja, ia, output matrix in compressed sparse row format
     integer, intent(out)  :: ia(nrow+1)
     integer, intent(out)  :: ja(nmax)
     real(dp), intent(out) :: a(nmax)

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k

! init sparse matrix
     a = 0.0_dp
     ia = 0
     ja = 0

     k = 1
     ia(1) = 1
     do i=1,nrow
         do j=1,ncol
             if ( dns(i,j) == 0.0_dp ) CYCLE
             ja(k) = j
             a(k) = dns(i,j)
             k = k + 1
             if ( k > nmax ) then
                 write(mystd,'(a)') 'sparse: error in sparse_format_dnscsr'
                 STOP
             endif ! back if ( k > nmax ) block
         enddo ! over j={1,ncol} loop
         ia(i+1) = k
     enddo ! over i={1,nrow} loop

     return
  end subroutine sparse_format_dnscsr

!>>> converts a densely stored matrix into a row orientied compactly sparse matrix
  subroutine sparse_format_dnscsr_(nrow, ncol, nmax, dns, sa, ja, ia)
     implicit none

! external arguments
! row dimension of dense matrix
     integer, intent(in)      :: nrow

! column dimension of dense matrix
     integer, intent(in)      :: ncol

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays sa and ja
     integer, intent(in)      :: nmax

! input densely stored matrix
     complex(dp), intent(in)  :: dns(nrow,ncol)

! sa, ja, ia, output matrix in compressed sparse row format
     integer, intent(out)     :: ia(nrow+1)
     integer, intent(out)     :: ja(nmax)
     complex(dp), intent(out) :: sa(nmax)

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k

! init sparse matrix
     sa = dcmplx(0.0_dp, 0.0_dp)
     ia = 0
     ja = 0

     k = 1
     ia(1) = 1
     do i=1,nrow
         do j=1,ncol
             if ( real( dns(i,j) ) == 0.0_dp .and. aimag( dns(i,j) ) == 0.0_dp ) CYCLE
             ja(k) = j
             sa(k) = dns(i,j)
             k = k + 1
             if ( k > nmax ) then
                 write(mystd,'(a)') 'sparse: error in sparse_format_dnscsr_'
                 STOP
             endif ! back if ( k > nmax ) block
         enddo ! over j={1,ncol} loop
         ia(i+1) = k
     enddo ! over i={1,nrow} loop

     return
  end subroutine sparse_format_dnscsr_

!>>> converts a densely stored identity matrix into a row orientied compactly sparse matrix
  subroutine sparse_format_unicsr(nrow, nmax, a, ja, ia)
     implicit none

! external arguments
! row dimension of dense matrix
     integer, intent(in)   :: nrow

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays a and ja
     integer, intent(in)   :: nmax

! a, ja, ia, output matrix in compressed sparse row format
     integer, intent(out)  :: ia(nrow+1)
     integer, intent(out)  :: ja(nmax)
     real(dp), intent(out) :: a(nmax)

! local variables
! loop index
     integer :: i

     if ( nrow > nmax ) then
         write(mystd,'(a)') 'sparse: error in sparse_format_unicsr'
         STOP
     endif ! back if ( nrow > nmax ) block

     do i=1,nrow
         ia(i) = i
         ja(i) = i
         a(i) = 1.0_dp
     enddo ! over i={1,nrow} loop

     ia(nrow+1) = nrow + 1
     ja(nrow+1:nmax) = 0
     a(nrow+1:nmax) = 0.0_dp

     return
  end subroutine sparse_format_unicsr

!>>> converts a densely stored identity matrix into a row orientied compactly sparse matrix
  subroutine sparse_format_unicsr_(nrow, nmax, sa, ja, ia)
     implicit none

! external arguments
! row dimension of dense matrix
     integer, intent(in)      :: nrow

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays sa and ja
     integer, intent(in)      :: nmax

! sa, ja, ia, output matrix in compressed sparse row format
     integer, intent(out)     :: ia(nrow+1)
     integer, intent(out)     :: ja(nmax)
     complex(dp), intent(out) :: sa(nmax)

! local variables
! loop index
     integer :: i

     if ( nrow > nmax ) then
         write(mystd,'(a)') 'sparse: error in sparse_format_unicsr_'
         STOP
     endif ! back if ( nrow > nmax ) block

     do i=1,nrow
         ia(i) = i
         ja(i) = i
         sa(i) = dcmplx(1.0_dp, 0.0_dp)
     enddo ! over i={1,nrow} loop

     ia(nrow+1) = nrow + 1
     ja(nrow+1:nmax) = 0
     do i=nrow+1,nmax
         sa(i) = dcmplx(0.0_dp, 0.0_dp)
     enddo ! over i={nrow+1,nmax} loop

     return
  end subroutine sparse_format_unicsr_

!>>> copy data between two row orientied compactly sparse matrices
  subroutine sparse_matrix_copyer(nrow, nmax, a, ja, ia, b, jb, ib)
     implicit none

! external arguments
! row dimension of dense matrix
     integer, intent(in)   :: nrow

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays a and ja
     integer, intent(in)   :: nmax

! a, ja, ia, input matrix in compressed sparse row format
     integer, intent(in)   :: ia(nrow+1)
     integer, intent(in)   :: ja(nmax)
     real(dp), intent(in)  :: a(nmax)

! b, jb, ib, output matrix in compressed sparse row format
     integer, intent(out)  :: ib(nrow+1)
     integer, intent(out)  :: jb(nmax)
     real(dp), intent(out) :: b(nmax)

! local variables
! loop index
     integer :: i

     do i=1,nrow+1
         ib(i) = ia(i)
     enddo ! over i={1,nrow+1} loop

     do i=ia(1),ia(nrow+1)-1
         jb(i) = ja(i)
     enddo ! over i={ia(1),ia(nrow+1)-1} loop

     do i=ia(1),ia(nrow+1)-1
         b(i) = a(i)
     enddo ! over i={ia(1),ia(nrow+1)-1} loop

     return
  end subroutine sparse_matrix_copyer

!>>> copy data between two row orientied compactly sparse matrices
  subroutine sparse_matrix_copyer_(nrow, nmax, sa, ja, ia, sb, jb, ib)
     implicit none

! external arguments
! row dimension of dense matrix
     integer, intent(in)      :: nrow

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays sa and ja
     integer, intent(in)      :: nmax

! sa, ja, ia, input matrix in compressed sparse row format
     integer, intent(in)      :: ia(nrow+1)
     integer, intent(in)      :: ja(nmax)
     complex(dp), intent(in)  :: sa(nmax)

! sb, jb, ib, output matrix in compressed sparse row format
     integer, intent(out)     :: ib(nrow+1)
     integer, intent(out)     :: jb(nmax)
     complex(dp), intent(out) :: sb(nmax)

! local variables
! loop index
     integer :: i

     do i=1,nrow+1
         ib(i) = ia(i)
     enddo ! over i={1,nrow+1} loop

     do i=ia(1),ia(nrow+1)-1
         jb(i) = ja(i)
     enddo ! over i={ia(1),ia(nrow+1)-1} loop

     do i=ia(1),ia(nrow+1)-1
         sb(i) = sa(i)
     enddo ! over i={ia(1),ia(nrow+1)-1} loop

     return
  end subroutine sparse_matrix_copyer_

!>>> this function returns the element a(i,j) of matrix a
  real(dp) &
  function sparse_matrix_getter(i, j, nrow, nmax, a, ja, ia) result(elm)
     implicit none

! external arguments
! the row index of the element sought
     integer, intent(in)  :: i

! the column index of the element sought
     integer, intent(in)  :: j

! row dimension of dense matrix
     integer, intent(in)  :: nrow

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays a and ja
     integer, intent(in)  :: nmax

! a, ja, ia, input matrix in compressed sparse row format
     integer, intent(in)  :: ia(nrow+1)
     integer, intent(in)  :: ja(nmax)
     real(dp), intent(in) :: a(nmax)

! local variables
! loop index
     integer :: k

! memory address of a(i,j)
     integer :: addr

! initialization
     addr = 0
     elm = 0.0_dp

! scan the row - exit as soon as a(i,j) is found
     do k=ia(i),ia(i+1)-1
         if ( ja(k) == j ) then
             addr = k
             EXIT
         endif ! back if ( ja(k) == j ) block
     enddo ! over k={ia(i),ia(i+1)-1} loop

! the required element is contained in sparse matrix
     if ( addr /= 0 ) then
         elm = a(addr)
     endif ! back if ( addr /= 0 ) block

     return
  end function sparse_matrix_getter

!>>> this function returns the element sa(i,j) of matrix sa
  complex(dp) &
  function sparse_matrix_getter_(i, j, nrow, nmax, sa, ja, ia) result(elm)
     implicit none

! external arguments
! the row index of the element sought
     integer, intent(in)     :: i

! the column index of the element sought
     integer, intent(in)     :: j

! row dimension of dense matrix
     integer, intent(in)     :: nrow

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays sa and ja
     integer, intent(in)     :: nmax

! sa, ja, ia, input matrix in compressed sparse row format
     integer, intent(in)     :: ia(nrow+1)
     integer, intent(in)     :: ja(nmax)
     complex(dp), intent(in) :: sa(nmax)

! local variables
! loop index
     integer :: k

! memory address of sa(i,j)
     integer :: addr

! initialization
     addr = 0
     elm = dcmplx(0.0_dp, 0.0_dp)

! scan the row - exit as soon as sa(i,j) is found
     do k=ia(i),ia(i+1)-1
         if ( ja(k) == j ) then
             addr = k
             EXIT
         endif ! back if ( ja(k) == j ) block
     enddo ! over k={ia(i),ia(i+1)-1} loop

! the required element is contained in sparse matrix
     if ( addr /= 0 ) then
         elm = sa(addr)
     endif ! back if ( addr /= 0 ) block

     return
  end function sparse_matrix_getter_

!>>> multiplies a matrix by a vector using the dot product form
  subroutine sparse_matmul_amuvec(nrow, ncol, nmax, a, ja, ia, x, y)
     implicit none

! external arguments
! row dimension of dense matrix
     integer, intent(in)   :: nrow

! column dimension of dense matrix
     integer, intent(in)   :: ncol

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays a and ja
     integer, intent(in)   :: nmax

! a, ja, ia, input matrix in compressed sparse row format
     integer, intent(in)   :: ia(nrow+1)
     integer, intent(in)   :: ja(nmax)
     real(dp), intent(in)  :: a(nmax)

! vector, length equal to the column dimension of the dense matrix
     real(dp), intent(in)  :: x(ncol)

! vector, real array of length nrow, containing the product y = A . x
     real(dp), intent(out) :: y(nrow)

! local variables
! loop index
     integer :: i
     integer :: k

! zero out output vector
     y = 0.0_dp

! compute the inner product of row i with vector x
     do i=1,nrow
         do k=ia(i),ia(i+1)-1
             y(i) = y(i) + a(k) * x( ja(k) )
         enddo ! over k={ia(i),ia(i+1)-1} loop
     enddo ! over i={1,nrow} loop

     return
  end subroutine sparse_matmul_amuvec

!>>> multiplies a matrix by a vector using the dot product form
  subroutine sparse_matmul_amuvec_(nrow, ncol, nmax, sa, ja, ia, sx, sy)
     implicit none

! external arguments
! row dimension of dense matrix
     integer, intent(in)      :: nrow

! column dimension of dense matrix
     integer, intent(in)      :: ncol

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays sa and ja
     integer, intent(in)      :: nmax

! sa, ja, ia, input matrix in compressed sparse row format
     integer, intent(in)      :: ia(nrow+1)
     integer, intent(in)      :: ja(nmax)
     complex(dp), intent(in)  :: sa(nmax)

! vector, length equal to the column dimension of the dense matrix
     complex(dp), intent(in)  :: sx(ncol)

! vector, complex(dp) array of length nrow, containing the product y = A . x
     complex(dp), intent(out) :: sy(nrow)

! local variables
! loop index
     integer :: i
     integer :: k

! zero out output vector
     sy = dcmplx(0.0_dp, 0.0_dp)

! compute the inner product of row i with vector sx
     do i=1,nrow
         do k=ia(i),ia(i+1)-1
             sy(i) = sy(i) + sa(k) * sx( ja(k) )
         enddo ! over k={ia(i),ia(i+1)-1} loop
     enddo ! over i={1,nrow} loop

     return
  end subroutine sparse_matmul_amuvec_

!>>> performs the matrix by matrix product C = A * B
  subroutine sparse_matmul_amumat(nrow, ndim, ncol, nmax, a, ja, ia, b, jb, ib, c, jc, ic)
     implicit none

! external arguments
! the row dimension of matrix A = row dimension of matrix C
     integer, intent(in)   :: nrow

! the column dimension of matrix A = row dimension of matrix B
     integer, intent(in)   :: ndim

! the column dimension of matrix B = column dimension of matrix C
     integer, intent(in)   :: ncol

! the length of the arrays c and jc
! sparse_matmul_amumat() will stop if the result matrix C has a number of
! elements that exceeds nmax
     integer, intent(in)   :: nmax

! a, ja, ia, matrix A in compressed sparse row format
     integer, intent(in)   :: ia(nrow+1)
     integer, intent(in)   :: ja(nmax)
     real(dp), intent(in)  :: a(nmax)

! b, jb, ib, matrix B in compressed sparse row format
     integer, intent(in)   :: ib(ndim+1)
     integer, intent(in)   :: jb(nmax)
     real(dp), intent(in)  :: b(nmax)

! c, jc, ic, resulting matrix C in compressed sparse row format
     integer, intent(out)  :: ic(nrow+1)
     integer, intent(out)  :: jc(nmax)
     real(dp), intent(out) :: c(nmax)

! local variables
! loop index
     integer :: i, j, k

! loop index
     integer :: ka, kb

! dummy integer variables
     integer :: p, q

! integer work array of length equal to the number of columns in matrix B,
! which is an array that has nonzero value if the column index already
! exist, in which case the value is the index of that column
     integer :: iw(ncol)

! dummy real(dp) variables, used to improve the ratio of floating point
! operations to memory accesses
     real(dp) :: atmp, btmp

! init work array
     iw = 0

! init C sparse matrix
     ic(1) = 1

     q = 0
     do i=1,nrow
         do ka=ia(i),ia(i+1)-1
             j = ja(ka)
             atmp = a(ka)
             do kb=ib(j),ib(j+1)-1
                 k = jb(kb)
                 btmp = b(kb)

                 p = iw(k)
                 if ( p == 0 ) then
                     q = q + 1
                     iw(k) = q
                     jc(q) = k
                     c(q) = atmp * btmp
                 else
                     c(p) = c(p) + atmp * btmp
                 endif ! back if ( p == 0 ) block
             enddo ! over kb={ib(j),ib(j+1)-1} loop
         enddo ! over ka={ia(i),ia(i+1)-1} loop

! done this row i, so set work array to zero again
         do k=ic(i),q
             iw( jc( k ) ) = 0
         enddo ! over k={ic(i),q} loop
         ic(i+1) = q + 1
     enddo ! over i={1,nrow} loop

! check the number of nonzero elements
     if ( q > nmax ) then
         write(mystd,'(a)') 'sparse: error in sparse_format_amumat'
         STOP
     endif ! back if ( q > nmax ) block

     return
  end subroutine sparse_matmul_amumat

!>>> performs the matrix by matrix product C = A * B
  subroutine sparse_matmul_amumat_(nrow, ndim, ncol, nmax, sa, ja, ia, sb, jb, ib, sc, jc, ic)
     implicit none

! external arguments
! the row dimension of matrix A = row dimension of matrix C
     integer, intent(in)      :: nrow

! the column dimension of matrix A = row dimension of matrix B
     integer, intent(in)      :: ndim

! the column dimension of matrix B = column dimension of matrix C
     integer, intent(in)      :: ncol

! the length of the arrays sc and jc
! sparse_matmul_amumat() will stop if the result matrix C has a number of
! elements that exceeds nmax
     integer, intent(in)      :: nmax

! sa, ja, ia, matrix A in compressed sparse row format
     integer, intent(in)      :: ia(nrow+1)
     integer, intent(in)      :: ja(nmax)
     complex(dp), intent(in)  :: sa(nmax)

! sb, jb, ib, matrix B in compressed sparse row format
     integer, intent(in)      :: ib(ndim+1)
     integer, intent(in)      :: jb(nmax)
     complex(dp), intent(in)  :: sb(nmax)

! sc, jc, ic, resulting matrix C in compressed sparse row format
     integer, intent(out)     :: ic(nrow+1)
     integer, intent(out)     :: jc(nmax)
     complex(dp), intent(out) :: sc(nmax)

! local variables
! loop index
     integer :: i, j, k

! loop index
     integer :: ka, kb

! dummy integer variables
     integer :: p, q

! integer work array of length equal to the number of columns in matrix B,
! which is an array that has nonzero value if the column index already
! exist, in which case the value is the index of that column
     integer :: iw(ncol)

! dummy complex(dp) variables, used to improve the ratio of floating point
! operations to memory accesses
     complex(dp) :: atmp, btmp

! init work array
     iw = 0

! init C sparse matrix
     ic(1) = 1

     q = 0
     do i=1,nrow
         do ka=ia(i),ia(i+1)-1
             j = ja(ka)
             atmp = sa(ka)
             do kb=ib(j),ib(j+1)-1
                 k = jb(kb)
                 btmp = sb(kb)

                 p = iw(k)
                 if ( p == 0 ) then
                     q = q + 1
                     iw(k) = q
                     jc(q) = k
                     sc(q) = atmp * btmp
                 else
                     sc(p) = sc(p) + atmp * btmp
                 endif ! back if ( p == 0 ) block
             enddo ! over kb={ib(j),ib(j+1)-1} loop
         enddo ! over ka={ia(i),ia(i+1)-1} loop

! done this row i, so set work array to zero again
         do k=ic(i),q
             iw( jc( k ) ) = 0
         enddo ! over k={ic(i),q} loop
         ic(i+1) = q + 1
     enddo ! over i={1,nrow} loop

! check the number of nonzero elements
     if ( q > nmax ) then
         write(mystd,'(a)') 'sparse: error in sparse_format_amumat_'
         STOP
     endif ! back if ( q > nmax ) block

     return
  end subroutine sparse_matmul_amumat_

!>>> performs the matrix by matrix product B = A * Diag
  subroutine sparse_matmul_amudia(nrow, nmax, a, ja, ia, diag, b, jb, ib)
     implicit none

! external arguments
! the row dimension of dense matrix
     integer, intent(in)   :: nrow

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays b and jb
     integer, intent(in)   :: nmax

! a, ja, ia, matrix A in compressed sparse row format
     integer, intent(in)   :: ia(nrow+1)
     integer, intent(in)   :: ja(nmax)
     real(dp), intent(in)  :: a(nmax)

! diagonal matrix stored as a vector diag
     real(dp), intent(in)  :: diag(nrow)

! b, jb, ib, resulting matrix B in compressed sparse row format
     integer, intent(out)  :: ib(nrow+1)
     integer, intent(out)  :: jb(nmax)
     real(dp), intent(out) :: b(nmax)

! local variables
! loop index
     integer :: i
     integer :: k

! loop index
     integer :: k1
     integer :: k2

! init B sparse matrix
     b = 0.0_dp
     ib = 0
     jb = 0

! scale each element
     do i=1,nrow
         k1 = ia(i)
         k2 = ia(i+1) - 1
         do k=k1,k2
             b(k) = a(k) * diag( ja(k) )
         enddo ! over k={k1,k2} loop
     enddo ! over i={1,nrow} loop

     do i=1,nrow+1
         ib(i) = ia(i)
     enddo ! over i={1,nrow+1} loop

     do k=ia(1),ia(nrow+1)-1
         jb(k) = ja(k)
     enddo ! over k={ia(1),ia(nrow+1)-1} loop

     return
  end subroutine sparse_matmul_amudia

!>>> performs the matrix by matrix product B = A * Diag
  subroutine sparse_matmul_amudia_(nrow, nmax, sa, ja, ia, diag, sb, jb, ib)
     implicit none

! external arguments
! the row dimension of dense matrix
     integer, intent(in)      :: nrow

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays sb and jb
     integer, intent(in)      :: nmax

! sa, ja, ia, matrix A in compressed sparse row format
     integer, intent(in)      :: ia(nrow+1)
     integer, intent(in)      :: ja(nmax)
     complex(dp), intent(in)  :: sa(nmax)

! diagonal matrix stored as a vector diag
     complex(dp), intent(in)  :: diag(nrow)

! sb, jb, ib, resulting matrix B in compressed sparse row format
     integer, intent(out)     :: ib(nrow+1)
     integer, intent(out)     :: jb(nmax)
     complex(dp), intent(out) :: sb(nmax)

! local variables
! loop index
     integer :: i
     integer :: k

! loop index
     integer :: k1
     integer :: k2

! init B sparse matrix
     sb = dcmplx(0.0_dp, 0.0_dp)
     ib = 0
     jb = 0

! scale each element
     do i=1,nrow
         k1 = ia(i)
         k2 = ia(i+1) - 1
         do k=k1,k2
             sb(k) = sa(k) * diag( ja(k) )
         enddo ! over k={k1,k2} loop
     enddo ! over i={1,nrow} loop

     do i=1,nrow+1
         ib(i) = ia(i)
     enddo ! over i={1,nrow+1} loop

     do k=ia(1),ia(nrow+1)-1
         jb(k) = ja(k)
     enddo ! over k={ia(1),ia(nrow+1)-1} loop

     return
  end subroutine sparse_matmul_amudia_

!>>> performs the matrix by matrix product B = Diag * A
  subroutine sparse_matmul_diamua(nrow, nmax, diag, a, ja, ia, b, jb, ib)
     implicit none

! external arguments
! the row dimension of dense matrix
     integer, intent(in)   :: nrow

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays b and jb
     integer, intent(in)   :: nmax

! diagonal matrix stored as a vector diag
     real(dp), intent(in)  :: diag(nrow)

! a, ja, ia, matrix A in compressed sparse row format
     integer, intent(in)   :: ia(nrow+1)
     integer, intent(in)   :: ja(nmax)
     real(dp), intent(in)  :: a(nmax)

! b, jb, ib, resulting matrix B in compressed sparse row format
     integer, intent(out)  :: ib(nrow+1)
     integer, intent(out)  :: jb(nmax)
     real(dp), intent(out) :: b(nmax)

! local variables
! loop index
     integer :: i
     integer :: k

! loop index
     integer :: k1
     integer :: k2

! init B sparse matrix
     b = 0.0_dp
     ib = 0
     jb = 0

! normalize each row
     do i=1,nrow
         k1 = ia(i)
         k2 = ia(i+1) - 1
         do k=k1,k2
             b(k) = a(k) * diag(i)
         enddo ! over k={k1,k2} loop
     enddo ! over i={1,nrow} loop

     do i=1,nrow+1
         ib(i) = ia(i)
     enddo ! over i={1,nrow+1} loop

     do k=ia(1),ia(nrow+1)-1
         jb(k) = ja(k)
     enddo ! over k={ia(1),ia(nrow+1)-1} loop

     return
  end subroutine sparse_matmul_diamua

!>>> performs the matrix by matrix product B = Diag * A
  subroutine sparse_matmul_diamua_(nrow, nmax, diag, sa, ja, ia, sb, jb, ib)
     implicit none

! external arguments
! the row dimension of dense matrix
     integer, intent(in)      :: nrow

! maximum number of nonzero elements allowed
! this should be set to be the lengths of the arrays sb and jb
     integer, intent(in)      :: nmax

! diagonal matrix stored as a vector diag
     complex(dp), intent(in)  :: diag(nrow)

! sa, ja, ia, matrix A in compressed sparse row format
     integer, intent(in)      :: ia(nrow+1)
     integer, intent(in)      :: ja(nmax)
     complex(dp), intent(in)  :: sa(nmax)

! sb, jb, ib, resulting matrix B in compressed sparse row format
     integer, intent(out)     :: ib(nrow+1)
     integer, intent(out)     :: jb(nmax)
     complex(dp), intent(out) :: sb(nmax)

! local variables
! loop index
     integer :: i
     integer :: k

! loop index
     integer :: k1
     integer :: k2

! init B sparse matrix
     sb = dcmplx(0.0_dp, 0.0_dp)
     ib = 0
     jb = 0

! normalize each row
     do i=1,nrow
         k1 = ia(i)
         k2 = ia(i+1) - 1
         do k=k1,k2
             sb(k) = sa(k) * diag(i)
         enddo ! over k={k1,k2} loop
     enddo ! over i={1,nrow} loop

     do i=1,nrow+1
         ib(i) = ia(i)
     enddo ! over i={1,nrow+1} loop

     do k=ia(1),ia(nrow+1)-1
         jb(k) = ja(k)
     enddo ! over k={ia(1),ia(nrow+1)-1} loop

     return
  end subroutine sparse_matmul_diamua_

  end module sparse
