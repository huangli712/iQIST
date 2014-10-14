!!!-------------------------------------------------------------------------
!!! project : jasmine
!!! program : atomic_driver_fullspace
!!!           atomic_driver_sectors
!!! source  : atomic_driver.f90
!!! type    : subroutines
!!! author  : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!!           08/13/2014 by yilin wang
!!!           08/22/2014 by yilin wang
!!! purpose : solve atomic problem for different CTQMC trace algorithms
!!! status  : unstable
!!! comment :
!!!-------------------------------------------------------------------------

!!>>> atomic_driver_fullspace: CTQMC direct matrices multiplications 
!!>>> trace algorithm, use full Hilbert space 
  subroutine atomic_driver_fullspace()
     use constants, only : dp, mystd
     use control, only : ncfgs

     use m_glob_fullspace, only : hmat, hmat_eigval, hmat_eigvec
     use m_glob_fullspace, only : alloc_m_glob_fullspace, dealloc_m_glob_fullspace
  
     implicit none
  
! local variables
! a temp matrix
     real(dp) :: tmp_mat(ncfgs, ncfgs)

! whether the Hamiltonian is real ? 
     logical :: lreal 
  
! allocate memory 
     write(mystd, "(2X,a)") "jasmine >>> allocate memory of global variables for fullspace case ..."
     write(mystd,*)
     call alloc_m_glob_fullspace()
  
! build atomic many particle Hamiltonian matrix
     write(mystd, "(2X,a)") "jasmine >>> make atomic many particle Hamiltonian ..."
     write(mystd,*)
     call atomic_mkhmat_fullspace()
  
! check whether the many particle Hamiltonian is real 
     write(mystd, "(2X,a)") "jasmine >>> check whether Hamiltonian is real or not ..."
     write(mystd,*)
     call atomic_check_realmat(ncfgs, hmat, lreal)
     if (lreal .eqv. .false.) then
         call s_print_error('atomic_driver_fullspace', 'hmat is not real !')
     else
         write(mystd, "(2X,a)") "jasmine >>> the atomic Hamiltonian is real"
         write(mystd,*)
     endif
  
! diagonalize hmat
     write(mystd, "(2X,a)") "jasmine >>> diagonalize the atomic Hamiltonian ..."
     write(mystd,*)
     tmp_mat = real(hmat)
     call s_eig_sy(ncfgs, ncfgs, tmp_mat, hmat_eigval, hmat_eigvec)
  
! build fmat
! first, build fmat of annihilation operators in Fock basis
! then, transform them to the eigen basis
     write(mystd, "(2X,a)") "jasmine >>> make fmat for annihilation fermion operators ... "
     write(mystd,*)
     call atomic_mkfmat_fullspace()

! build occupancy number
     write(mystd, "(2X,a)") "jasmine >>> make occupancy number of atomic eigenstates ... "
     write(mystd,*)
     call atomic_mkoccu_fullspace()    

! write eigenvalues of hmat to file 'atom.eigval.dat'
     write(mystd, "(2X,a)") "jasmine >>> write eigenvalue, eigenvector, and atom.cix to files ..."
     write(mystd,*)
     call atomic_write_eigval_fullspace()
  
! write eigenvectors of hmat to file 'atom.eigvec.dat'
     call atomic_write_eigvec_fullspace()
   
! write eigenvalue of hmat, occupany number of eigenstates and 
! fmat of annihilation fermion operators to file "atom.cix"
! this is for begonia, lavender codes of iQIST package
     call atomic_write_atomcix_fullspace()
  
! deallocate memory
     write(mystd, "(2X,a)") "jasmine >>> free memory of global variables for fullspace case ... "
     write(mystd,*)
     call dealloc_m_glob_fullspace()
  
     return
  end subroutine atomic_driver_fullspace
  
!!>>> atomic_driver_sectors: CTQMC trace algorithm: use good quantum numbers (GQNs)
  subroutine atomic_driver_sectors()
     use constants, only : mystd
     use m_glob_sectors, only : nsectors, sectors, dealloc_m_glob_sectors
  
     implicit none

! local variables
! loop index
     integer :: i

! whether the Hamiltonian is real ?
     logical :: lreal

! make all the sectors, allocate m_glob_sectors memory inside
     write(mystd, "(2X,a)") "jasmine >>> determine sectors by good quantum numbers (GQNs)... "
     write(mystd,*)
     call atomic_mksectors()
 
! make atomic Hamiltonian
     write(mystd, "(2X,a)") "jasmine >>> make atomic Hamiltonian for each sector ... "
     write(mystd,*)
     call atomic_mkhmat_sectors()
  
! check whether the many particle Hamiltonian is real 
     write(mystd, "(2X,a)") "jasmine >>> check whether Hamiltonian is real or not ..."
     write(mystd,*)
     do i=1, nsectors 
         call atomic_check_realmat(sectors(i)%ndim, sectors(i)%myham, lreal)
         if (lreal .eqv. .false.) then
             call s_print_error('atomic_solve_sectors', 'hmat is not real !')
         endif
     enddo
     write(mystd, "(2X,a)") "jasmine >>> the Hamiltonian is real"
     write(mystd,*)
  
! diagonalize Hamiltonian of each sector one by one
     write(mystd, "(2X,a)") "jasmine >>> diagonalize atomic Hamiltonian for each sector ... "
     write(mystd,*)
     call atomic_diaghmat_sectors()
  
! make fmat of both creation and annihilation operators for each sector
     write(mystd, "(2X,a)") "jasmine >>> make fmat for each sector ..."
     write(mystd,*)
     call atomic_mkfmat_sectors()
  
     write(mystd, "(2X,a)") "jasmine >>> write eigenvalue, eigenvector, and atom.cix to files ... "
     write(mystd,*)
! write eigenvalues to file 'atom.eigval.dat'
     call atomic_write_eigval_sectors()
  
! write eigenvectors to file 'atom.eigvec.dat'
     call atomic_write_eigvec_sectors()
  
! write information of sectors to file 'atom.cix'
     call atomic_write_atomcix_sectors()

! free memory
     write(mystd, "(2X,a)") "jasmine >>> free memory for sectors case ..."
     write(mystd,*)
     call dealloc_m_glob_sectors()
 
     return
  end subroutine atomic_driver_sectors
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
!!!----------------------------------------------------------------------------
!!! project : jasmine
!!! program : atomic_mkfmat_fullspace
!!!           atomic_mkfmat_sectors
!!!           atomic_rotate_fmat
!!!           atomic_make_construct
!!!           atomic_make_eliminate 
!!! source  : atomic_fmat.f90
!!! type    : subroutines
!!! author  : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!! purpose : make fmat
!!! status  : unstable
!!! comment : these subroutines are based on Dr. LiangDu's (duleung@gmail.com) 
!!!           atomic program
!!!----------------------------------------------------------------------------

!!>>> atomic_mkfmat_fullspace: make fmat for annihilation operators 
!!>>> for full space case
  subroutine atomic_mkfmat_fullspace()
     use control, only : norbs, ncfgs
     use m_glob_fullspace, only : anni_fmat, hmat_eigvec
     use m_basis_fullspace, only : dec_basis, index_basis
  
     implicit none
  
! local variables
! loop index
     integer :: i,j,k

! left Fock state
     integer :: left

! right Fock state
     integer :: right

! the sign change 
     integer :: isgn
  
     do i=1,norbs
         do j=1,ncfgs
             right = dec_basis(j)
             if (btest(right, i-1) .eqv. .true.) then
                call atomic_make_eliminate(i, right, left, isgn)
                k = index_basis(left)
                anni_fmat(k, j, i) = dble(isgn)
             endif
         enddo 
     enddo 
  
! rotate it to the atomic eigenvector basis
     do i=1, norbs
         call atomic_tran_repr_real(ncfgs, anni_fmat(:,:,i), hmat_eigvec)
     enddo
  
     return
  end subroutine atomic_mkfmat_fullspace

!!>>> atomic_mkfmat_sectors: build fmat for good quantum numbers (GQNs) algorithm
  subroutine atomic_mkfmat_sectors()
     use constants, only : zero
     use control, only : norbs

     use m_basis_fullspace, only : dec_basis, index_basis
     use m_glob_sectors, only : nsectors, sectors
     use m_sector, only : alloc_one_fmat
  
     implicit none
  
! local variables
! loop index 
     integer :: iorb
     integer :: ifermi
     integer :: isect, jsect
     integer :: ibas, jbas
     integer :: i

! sign change due to commute relation
     integer :: isgn

! auxiliary integer variables
     integer :: jold, jnew
  
! loop over all the sectors
     do isect=1, nsectors
! loop over all the orbitals
         do iorb=1,norbs
! loop over the creation and annihilation fermion operators
             do ifermi=0, 1 
                 jsect = sectors(isect)%next_sector(iorb, ifermi) 
                 if (jsect == -1) cycle
! allocate memory for fmat
                 sectors(isect)%myfmat(iorb, ifermi)%n = sectors(jsect)%ndim
                 sectors(isect)%myfmat(iorb, ifermi)%m = sectors(isect)%ndim
                 call alloc_one_fmat(sectors(isect)%myfmat(iorb, ifermi))
                 sectors(isect)%myfmat(iorb,ifermi)%item = zero
! build fmat
                 do jbas=1, sectors(isect)%ndim
                     jold = dec_basis(sectors(isect)%mybasis(jbas))
! for creation fermion operator
                     if (ifermi == 1 .and. ( btest(jold, iorb-1) .eqv. .false. )) then
                         call atomic_make_construct(iorb, jold, jnew, isgn)
! for annihilation fermion operator
                     elseif (ifermi == 0 .and. ( btest(jold, iorb-1) .eqv. .true. )) then
                         call atomic_make_eliminate(iorb, jold, jnew, isgn)
                     else
                         cycle
                     endif
                     ibas = index_basis(jnew)
                     do i=1, sectors(jsect)%ndim 
                         if (ibas == sectors(jsect)%mybasis(i)) then
                             ibas = i
                             sectors(isect)%myfmat(iorb, ifermi)%item(ibas, jbas) = dble(isgn)
                             exit
                         endif
                     enddo
                 enddo  ! over jbas={1, sectors(isect)%ndim} loop
! roate fmat to atomic eigenstates basis
                 call atomic_rotate_fmat(sectors(jsect)%ndim, sectors(isect)%ndim, sectors(jsect)%myeigvec, &
                     sectors(isect)%myfmat(iorb, ifermi)%item, sectors(isect)%myeigvec)
             enddo ! over ifermi={0,1} loop
         enddo ! over iorb={1, norbs} loop
     enddo ! over isect={1,nsectors} loop
  
     return
  end subroutine atomic_mkfmat_sectors

!!>>> atomic_rotate_fmat: rotate fmat from Fock basis to eigenstates basis
  subroutine atomic_rotate_fmat(ndimx, ndimy, amat, bmat, cmat)
     use constants, only: dp, zero, one
     implicit none
     
! external variables
     integer, intent(in) :: ndimx
     integer, intent(in) :: ndimy
     real(dp), intent(in) :: amat(ndimx, ndimx)
     real(dp), intent(inout) :: bmat(ndimx, ndimy)
     real(dp), intent(in) :: cmat(ndimy, ndimy)
  
! local variables
     real(dp) :: tmp_mat(ndimx, ndimy)
     real(dp) :: amat_t(ndimx, ndimx)
     real(dp) :: alpha
     real(dp) :: betta
  
     amat_t = transpose(amat)
     tmp_mat = zero
  
     alpha = one; betta = zero
     call dgemm('N', 'N', ndimx, ndimy, ndimy, &
                           alpha, bmat, ndimx, &
                                  cmat, ndimy, &
                        betta, tmp_mat, ndimx  )
  
     alpha = one; betta = zero
     call dgemm('N', 'N', ndimx, ndimy, ndimx, &
                         alpha, amat_t, ndimx, &
                               tmp_mat, ndimx, &
                           betta, bmat, ndimx  )
  
  
     return
  end subroutine atomic_rotate_fmat

!!>>> atomic_make_construct: create one electron on ipos 
!!>>> of |jold> to deduce |jnew>
  subroutine atomic_make_construct(ipos, jold, jnew, isgn)
     implicit none
  
! external argument
! position number (serial number of orbit)
     integer, intent(in) :: ipos
  
! old Fock state and new Fock state
     integer, intent(in ):: jold
     integer, intent(out):: jnew
  
! sgn due to anti-commute relation between fernions
     integer, intent(out):: isgn
  
! local variables
! loop index over orbit
     integer :: iorb
  
     if (btest(jold, ipos-1) .eqv. .true.) then
         call s_print_error("atomic_construct", "severe error happened")
     endif
  
     isgn = 0
     do iorb=1,ipos-1
        if (btest(jold, iorb-1)) isgn = isgn + 1
     enddo
     isgn = mod(isgn, 2)
  
     isgn = (-1)**isgn
     jnew = jold + 2**(ipos-1)
  
     return
  end subroutine atomic_make_construct

!!>>> atomic_make_eliminate: destroy one electron on ipos 
!!>>> of |jold> to deduce |jnew>
  subroutine atomic_make_eliminate(ipos, jold, jnew, isgn)
      implicit none
  
! external argument
! position number (serial number of orbit)
      integer, intent(in)  :: ipos
  
! old Fock state and new Fock state
      integer, intent(in ) :: jold
      integer, intent(out) :: jnew
  
! sgn due to anti-commute relation between fernions
      integer, intent(out) :: isgn
  
! local variables
! loop index
      integer :: iorb
  
      if (btest(jold, ipos-1) .eqv. .false.) then
          call s_print_error("atomic_eliminate", "severe error happened")
      endif 
  
      isgn = 0
      do iorb=1,ipos-1
          if (btest(jold, iorb-1)) isgn = isgn + 1
      enddo
      isgn = mod(isgn, 2)
  
      isgn = (-1)**isgn
      jnew = jold - 2**(ipos-1)
  
      return
  end subroutine atomic_make_eliminate

