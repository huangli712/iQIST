!-------------------------------------------------------------------------
! project : jasmine
! program : atomic_make_basis_fullspace
!           state_pick
! source  : atomic_basis.f90
! type    : subroutines
! author  : yilin wang (email: qhwyl2006@126.com)
! history : 07/09/2014 by yilin wang
! purpose : make Fock basis 
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> make Fock basis for full Hilbert space, this subroutine is taken from
! Dr. LiangDu's (duleung@gmail.com) atomic program
subroutine atomic_make_basis_fullspace()
    use control, only: norbs, ncfgs
    use m_basis_fullspace  

    implicit none

    ! external variables
    integer, external :: state_pick

    ! local variables
    ! loop index
    integer :: i, j, k
    ! basis counter
    integer :: basis_count
    ! number of electrons for Fock state
    integer :: nelec

    ! first, allocate memory for basis related matrices
    call alloc_m_basis_fullspace()

    do i=0,norbs
        dim_sub_n(i) = state_pick(i, norbs)
    enddo 

    do i=0, norbs
        do j=0, 2**norbs-1
            nelec = 0
            do k=1,norbs
                if( btest(j, k-1) ) nelec = nelec + 1
            enddo 
            if ( nelec .eq. i ) then
                basis_count = basis_count + 1
                dec_basis(basis_count) = j
                index_basis(j) = basis_count
            endif 
        enddo 
    enddo 

    ! construct binary from of Fock basis
    do i=1,ncfgs
        do j=1,norbs
            if( btest(dec_basis(i), j-1) ) bin_basis(j, i) = 1
        enddo 
    enddo 

    ! dump atomic configurations to file "atom.basis.dat"
    call atomic_write_basis()

    return
end subroutine atomic_make_basis_fullspace

!>>> calculate combination algebra 
function state_pick(ntiny, nlarg) result(value)
    implicit none

    ! external variables
    integer, intent(in) :: ntiny
    integer, intent(in) :: nlarg

    ! local variables
    integer :: i
    ! auxiliary integer variable
    integer :: nlow
    ! numberator of the combination algebra
    real(8) :: numer
    ! denominator of the combination algebra
    real(8) :: denom
    ! result value of the combination algebra
    integer :: value

    ! transform the combination algebra
    nlow = min(ntiny, nlarg-ntiny)

    ! numerator in combination algebra
    numer = 1.0D0
    do i=nlarg-nlow+1,nlarg
       numer = numer * dble(i)
    enddo ! over i={nlarg-nlow+1,nlarg} loop

    ! denominator in combination algebra
    denom = 1.0D0
    do i=1,nlow
       denom = denom * dble(i)
    enddo ! over i={1,nlow} loop

    ! result value
    value = nint(numer / denom)

    return
end function state_pick


