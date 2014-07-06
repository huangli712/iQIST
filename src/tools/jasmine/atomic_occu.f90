!>>> make occupancy for atomic eigenstates, fullspace case
  subroutine atomic_make_occumat_fullspace()
     use constants
     use control
     use m_basis_fullspace
     use m_mpmat_fullspace

     implicit none

! local variables
! loop index over orbits
     integer :: iorb

! loop index over configurations
     integer :: ibas

! initialize nnmat to be zero
     do ibas=1,ncfgs
         do iorb=1,norbs
             if (bin_basis(iorb, ibas) .eq. 1) then
                 mp_occu_mat(ibas, ibas) = mp_occu_nmat(ibas, ibas) + one
             endif
         enddo ! over iorb={1,norbs} loop
     enddo ! over ibas={1,ncfgs} loop

     return
  end subroutine atomic_make_occumat_fullspace

  subroutine atomic_build_nnmat(norbs, ncfgs, invcd, nnmat)
     use constants

     implicit none

! external arguments
! number of orbits
     integer, intent(in) :: norbs

! number of configurations
     integer, intent(in) :: ncfgs

! binary representation of Fock basis
     integer, intent(in) :: invcd(norbs, ncfgs)

! nn operator matrix in many body basis
     complex(dp), intent(out) :: nnmat(ncfgs, ncfgs)

! local variables
! loop index over orbits
     integer :: iorb

! loop index over configurations
     integer :: ibas

! initialize nnmat to be zero
     nnmat = czero

     do ibas=1,ncfgs
         do iorb=1,norbs
             if (invcd(iorb, ibas) .eq. 1) then
                 nnmat(ibas, ibas) = nnmat(ibas, ibas) + done
             endif
         enddo ! over iorb={1,norbs} loop
     enddo ! over ibas={1,ncfgs} loop

     return
  end subroutine atomic_build_nnmat

