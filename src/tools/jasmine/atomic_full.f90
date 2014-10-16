!!!-----------------------------------------------------------------------
!!! project : jasmine
!!! program : atomic_make_foccu
!!!           atomic_make_ffmat
!!!           atomic_make_fhmat 
!!! source  : atomic_full.f90
!!! type    : subroutines
!!! author  : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!!           08/22/2014 by yilin wang
!!! purpose : make occupancy of eigensates of atomic Hamiltonian 
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> atomic_make_foccu: make occupancy for atomic eigenstates, fullspace case
  subroutine atomic_make_foccu()
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
  end subroutine atomic_make_foccu

!!>>> atomic_make_ffmat: make fmat for annihilation operators 
!!>>> for full space case
  subroutine atomic_make_ffmat()
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
                call atomic_make_c(i, right, left, isgn)
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
  end subroutine atomic_make_ffmat

!!>>> atomic_make_fhmat: make atomic Hamiltonian for the full space
  subroutine atomic_make_fhmat()
     use constants, only : czero, epst
     use control, only : norbs, ncfgs

     use m_basis_fullspace, only : dec_basis, index_basis, bin_basis
     use m_spmat, only : eimpmat, cumat
     use m_glob_fullspace, only : hmat
  
     implicit none
  
! local variables
! loop index
     integer :: i
     integer :: ibas, jbas
     integer :: alpha, betta
     integer :: delta, gamma

! sign change due to fermion anti-commute relation
     integer :: sgn

! new basis state after fermion operators act
     integer :: knew

! binary code form of a state
     integer :: code(norbs)
  
! start to make Hamiltonian
! initialize hmat
     hmat = czero
  
! first, two fermion operator terms 
     do jbas=1, ncfgs
         alploop: do alpha=1,norbs
         betloop: do betta=1,norbs
  
             sgn = 0
             knew = dec_basis(jbas)
             code(1:norbs) = bin_basis(1:norbs, jbas)
           
             if ( abs(eimpmat(alpha, betta)) .lt. epst ) cycle
  
! simulate one annihilation operator
             if (code(betta) == 1) then
                 do i=1,betta-1
                     if (code(i) == 1) sgn = sgn + 1
                 enddo 
                 code(betta) = 0
  
! simulate one creation operator
                 if (code(alpha) == 0) then
                     do i=1,alpha-1
                         if (code(i) == 1) sgn = sgn + 1
                     enddo
                     code(alpha) = 1
  
! determine the row number and hamiltonian matrix elememt
                     knew = knew - 2**(betta-1)
                     knew = knew + 2**(alpha-1)
                     sgn  = mod(sgn, 2)
                     ibas = index_basis(knew)
                     if (ibas == 0) then
                         call s_print_error('atomic_mkhmat_fullspace', &
                                            'error while determining row1')
                     endif
  
                     hmat(ibas,jbas) = hmat(ibas,jbas) + eimpmat(alpha,betta) * (-1.0d0)**sgn 
  
                 endif ! back if (code(alpha) == 0) block
             endif ! back if (code(betta) == 1) block
  
         enddo betloop ! over betta={1,norbs} loop
         enddo alploop ! over alpha={1,norbs} loop
     enddo ! over jbas={1,ncfgs} loop
  
! four fermion operator terms (coulomb interaction)
     do jbas=1, ncfgs
         alphaloop : do alpha=1,norbs
         bettaloop : do betta=1,norbs
         deltaloop : do delta=1,norbs
         gammaloop : do gamma=1,norbs
  
             sgn  = 0
             knew = dec_basis(jbas)
             code(1:norbs) = bin_basis(1:norbs, jbas)
  
             if ((alpha .eq. betta) .or. (delta .eq. gamma)) cycle
             if ( abs(cumat(alpha,betta,delta,gamma)) .lt. epst ) cycle
  
! simulate two annihilation operators
             if ((code(delta) == 1) .and. (code(gamma) == 1)) then
                 do i=1,gamma-1
                     if(code(i) == 1) sgn = sgn + 1
                 enddo 
                 code(gamma) = 0
  
                 do i=1,delta-1
                     if(code(i) == 1) sgn = sgn + 1
                 enddo 
                 code(delta) = 0
  
! simulate two creation operator
                 if ((code(alpha) == 0) .and. (code(betta) == 0)) then
                     do i=1,betta-1
                         if(code(i) == 1) sgn = sgn + 1
                     enddo 
                     code(betta) = 1
  
                     do i=1,alpha-1
                         if(code(i) == 1) sgn = sgn + 1
                     enddo 
                     code(alpha) = 1
  
! determine the row number and hamiltonian matrix elememt
                     knew = knew - 2**(gamma-1) - 2**(delta-1)
                     knew = knew + 2**(betta-1) + 2**(alpha-1)
                     sgn  = mod(sgn, 2)
                     ibas = index_basis(knew)
                     if (ibas == 0) then
                         call s_print_error('atomic_mkhmat_fullspace', 'error while determining row3')
                     endif
  
                     hmat(ibas,jbas) = hmat(ibas,jbas) + cumat(alpha,betta,delta,gamma) * (-1.0d0)**sgn
  
                 endif ! back if ((code(delta) == 1) .and. (code(gamma) == 1)) block
             endif ! back if ((code(alpha) == 0) .and. (code(betta) == 0)) block
  
         enddo gammaloop ! over gamma={1,norbs-1} loop
         enddo deltaloop ! over delta={gamma+1,norbs} loop
         enddo bettaloop ! over betta={alpha+1,norbs} loop
         enddo alphaloop ! over alpha={1,norbs-1} loop
     enddo  
  
     return
  end subroutine atomic_make_fhmat
