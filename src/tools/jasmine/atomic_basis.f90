  subroutine atomic_make_basis_full()
     use constants
     use control
     use m_basis_fullspace

     implicit none

! local variables
! loop index
     integer :: i, j

! basis counter
     integer :: basis_count

! number of electrons for Fock state
     integer :: nelec

! first, allocate memory for basis related matrices
     call alloc_m_basis_fullspace()

     do i=0,norbs
         dim_sub_N(i) = state_pick(i, norbs)
     enddo 

     do i=0, norbs
         do j=0, 2**norbs-1
             call verif( j, norbs, nelec )
             if ( nelec .eq. i ) then
                 basis_count = basis_count + 1
                 dec_basis(basis_count) = j
                 index_basis(j) = basis_count
             endif 
         enddo 
     enddo 

! construct inverse binary code from a decimal number 
     do i=1,ncfgs
         do j=1,norbs
             if( btest(dec_basis(i), j-1) ) bin_basis(j, i) = 1
         enddo 
     enddo 

! dump atomic configurations to file "atom.basis.dat"
     call atomic_write_basis()

     return
  end subroutine atomic_make_basis_fullspace

!>>> electron number denoted by a decimal number
  subroutine verif( inum, norbs, nbits )
     implicit none

! external arguments
! decimal number
     integer, intent(in)  :: inum

! number of orbits
     integer, intent(in)  :: norbs

! number of total electrons
     integer, intent(out) :: nbits

! local variables
! local index over orbits
     integer :: ipos

! number of total electrons
     nbits = 0
     do ipos=1,norbs
         if( btest(inum, ipos-1) ) nbits = nbits + 1
     enddo ! over ipos={1,norbs} loop

     return
  end subroutine verif

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


