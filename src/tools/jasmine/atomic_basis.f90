!!!------------------------------------------------------------------------
!!! project : jasmine
!!! program : atomic_make_basis_fullspace
!!!           atomic_combination
!!! source  : atomic_basis.f90
!!! type    : subroutines
!!! author  : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!! purpose : make Fock basis 
!!! input   :
!!! output  :
!!! status  : unstable
!!! comment : 
!!!-------------------------------------------------------------------------

!!>>> make Fock basis for full Hilbert space, this subroutine is modified based on
! Dr. LiangDu's (duleung@gmail.com) atomic program
  subroutine atomic_make_basis_fullspace()
     use control,            only: norbs, ncfgs
     use m_basis_fullspace,  only: dim_sub_n, dec_basis, bin_basis, index_basis  
  
     implicit none
  
! external variables
     integer, external :: atomic_combination
  
! local variables
! loop index
     integer :: i, j, k

! basis counter
     integer :: basis_count

! number of electrons for Fock state
     integer :: nelec
  
! initialize them
     dim_sub_n = 0
     dec_basis = 0
     bin_basis = 0
     index_basis = 0
  
! it is a number of combination C_{norbs}^{i}
     do i=0,norbs
         dim_sub_n(i) = atomic_combination(i, norbs)
     enddo 
  
! construct decimal form and index of Fock basis
     basis_count = 0
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
  
! construct binary form of Fock basis
     do i=1,ncfgs
         do j=1,norbs
             if( btest(dec_basis(i), j-1) ) bin_basis(j, i) = 1
         enddo 
     enddo 
  
! dump Fock basis to file "atom.basis.dat" for reference
     call atomic_write_basis()
  
     return
  end subroutine atomic_make_basis_fullspace

!!>>> calculate combination algebra 
  function atomic_combination(ntiny, nlarg) result(value)
     use constants, only: dp

     implicit none
  
! external variables
! the small number
     integer, intent(in) :: ntiny

! the large number
     integer, intent(in) :: nlarg
 
! result value of the combination algebra
     integer :: value
  
! local variables
! loop index
     integer :: i

! auxiliary integer variable
     integer :: nlow

! numberator of the combination algebra
     real(dp) :: numer

! denominator of the combination algebra
     real(dp) :: denom

 
! find the minimum number 
     nlow = min(ntiny, nlarg-ntiny)
  
! numerator in combination algebra
     numer = 1.0_dp
     do i=nlarg-nlow+1,nlarg
        numer = numer * dble(i)
     enddo 
  
! denominator in combination algebra
     denom = 1.0_dp
     do i=1,nlow
        denom = denom * dble(i)
     enddo 
  
! result value
     value = nint(numer / denom)
  
     return
  end function atomic_combination
