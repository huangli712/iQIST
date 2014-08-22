!!!-------------------------------------------------------------------------
!!! project : jasmine
!!! program : atomic_mkoccu_fullspace
!!! source  : atomic_occu.f90
!!! type    : subroutines
!!! author  : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!!           08/22/2014 by yilin wang
!!! purpose : make occupancy of eigensates of atomic Hamiltonian 
!!! status  : unstable
!!! comment :
!!!-------------------------------------------------------------------------

!!>>> atomic_mkoccu_fullspace: make occupancy for atomic eigenstates, fullspace case
  subroutine atomic_mkoccu_fullspace()
     use constants, only: zero, one
     use control, only: ncfgs, norbs

     use m_basis_fullspace, only: bin_basis
     use m_glob_fullspace, only: occu_mat, hmat_eigvec
  
     implicit none
  
! local variables
! loop index over orbits
     integer :: iorb
  
! loop index over configurations
     integer :: ibas
  
     occu_mat = zero
     do ibas=1,ncfgs
         do iorb=1,norbs
             if (bin_basis(iorb, ibas) .eq. 1) then
                 occu_mat(ibas, ibas) = occu_mat(ibas, ibas) + one
             endif
         enddo 
     enddo 
  
     call atomic_tran_repr_real(ncfgs, occu_mat, hmat_eigvec)
  
     return
  end subroutine atomic_mkoccu_fullspace
