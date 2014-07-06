!>>> subroutines used to make many particle matrices
!>>> make two fermion operator matrices
  subroutine atomic_make_two_fermion(norbs, ncfgs, eimp, dec_basis, bin_basis, index_basis, mp_mat)
     use constants
     
     implicit none

! external variables
! number of orbitals
     integer, intent(in) :: norbs

! number of configurations
     integer, intent(in) :: ncfgs

! the single-particle matrices
     complex(dp), intent(in) :: eimp(norbs, norbs)

! the decimal form of basis
     integer, intent(in) :: dec_basis(ncfgs)

! the binary form of basis
     integer, intent(in) :: bin_basis(norbs, ncfgs)

! the index of basis
     integer, intent(in) :: index_basis(0: ncfgs-1)

! the calculated many particle matrix
     complex(dp), intent(out) :: mp_mat(ncfgs, ncfgs)

! local variables
! loop index
     integer :: i,j

! loop index over orbits
     integer :: iorb

! sign change due to fermion anti-commute relation
     integer :: sgn

! new basis state after four fermion operation
     integer :: knew

! loop index over Fock basis
     integer :: ibas
     integer :: jbas

! index for general interaction matrix
     integer :: alpha
     integer :: betta
     integer :: delta
     integer :: gamma

! binary code representation of a state
     integer :: code(norbs)

! initialize mp_mat
     mp_mat = czero

! two fermion operator terms 
     do jbas=1, ncfgs
         alploop: do alpha=1,norbs
         betloop: do betta=1,norbs

             sgn = 0
             knew = dec_basis(jbas)
             code(1:norbs) = bin_basis(1:norbs, jbas)
             if ( abs(eimp(alpha, betta)) .lt. eps6 ) cycle

! simulate one eliminate operator
             if (code(betta) == 1) then
                 do i=1,betta-1
                     if (code(i) == 1) sgn = sgn + 1
                 enddo 
                 code(betta) = 0

! simulate one construct operator
                 if (code(alpha) == 0) then
                     do i=1,alpha-1
                         if (code(i) == 1) sgn = sgn + 1
                     enddo
                     code(alpha) = 1

! determine the column number and hamiltonian matrix elememt
                     knew = knew - 2**(betta-1)
                     knew = knew + 2**(alpha-1)
 
                     sgn  = mod(sgn, 2)
                     ibas = index_basis(knew)
                     if (ibas == 0) stop "error while determining row1"
                     mp_mat(ibas,jbas) = mp_mat(ibas,jbas) + eimp(alpha,betta) * (-1.0d0)**sgn 

                 endif ! back if (code(alpha) == 0) block
             endif ! back if (betta == 1) block

         enddo betloop ! over betta={1,norbs} loop
         enddo alploop ! over alpha={1,norbs} loop
     enddo ! over jbas={1,nbas} loop

     return
  end subroutine atomic_make_two_fermion

!>>> make four fermion operators
  subroutine atomic_make_four_fermion(norbs, ncfgs, cumat, dec_basis, bin_basis, index_basis, mp_mat)
     use constants
     
     implicit none

! external variables
! number of orbitals
     integer, intent(in) :: norbs

! number of configurations
     integer, intent(in) :: ncfgs

! the single-particle tensor
     complex(dp), intent(in) :: cumat(norbs, norbs, norbs, norbs)

! the decimal form of basis
     integer, intent(in) :: dec_basis(ncfgs)

! the binary form of basis
     integer, intent(in) :: bin_basis(norbs, ncfgs)

! the index of basis
     integer, intent(in) :: index_basis(0: ncfgs-1)

! the calculated many particle matrix
     complex(dp), intent(out) :: mp_mat(ncfgs, ncfgs)

! local variables
! loop index
     integer :: i,j

! loop index over orbits
     integer :: iorb

! sign change due to fermion anti-commute relation
     integer :: sgn

! new basis state after four fermion operation
     integer :: knew

! loop index over Fock basis
     integer :: ibas
     integer :: jbas

! index for general interaction matrix
     integer :: alpha
     integer :: betta
     integer :: delta
     integer :: gamma

! binary code representation of a state
     integer :: code(norbs)

! initialize mp_mat
     mp_mat = czero

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
            if ( abs(cumat(alpha,betta,delta,gamma)) .lt. eps6 ) cycle

! simulate two eliminate operator
            if ((code(delta) == 1) .and. (code(gamma) == 1)) then
                do i=1,gamma-1
                    if(code(i) == 1) sgn = sgn + 1
                enddo ! over i={1,gamma-1} loop
                code(gamma) = 0

                do i=1,delta-1
                    if(code(i) == 1) sgn = sgn + 1
                enddo ! over i={1,delta-1} loop
                code(delta) = 0

! simulate two construct operator
                if ((code(alpha) == 0) .and. (code(betta) == 0)) then
                    do i=1,betta-1
                        if(code(i) == 1) sgn = sgn + 1
                    enddo ! over i={1,betta-1} loop
                    code(betta) = 1

                    do i=1,alpha-1
                        if(code(i) == 1) sgn = sgn + 1
                    enddo ! over i={1,alpha-1} loop
                    code(alpha) = 1

! determine the column number and hamiltonian matrix elememt
                    knew = knew - 2**(gamma-1) - 2**(delta-1)
                    knew = knew + 2**(betta-1) + 2**(alpha-1)

                    sgn  = mod(sgn, 2)
                    ibas = index_basis(knew)
                    if (ibas == 0) stop "error while determining row3"
                    mp_mat(ibas,jbas) = mp_mat(ibas,jbas) + cumat(alpha,betta,delta,gamma) * (-1.0d0)**sgn

                endif ! back if ((code(delta) == 1) .and. (code(gamma) == 1)) block
            endif ! back if ((code(alpha) == 0) .and. (code(betta) == 0)) block

        enddo gammaloop ! over gamma={1,norbs-1} loop
        enddo deltaloop ! over delta={gamma+1,norbs} loop
        enddo bettaloop ! over betta={alpha+1,norbs} loop
        enddo alphaloop ! over alpha={1,norbs-1} loop

     enddo ! over jbas={1,nbas} loop

     return
  end subroutine atomic_make_four_fermion
