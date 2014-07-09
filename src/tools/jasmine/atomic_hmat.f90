!-------------------------------------------------------------------------
! project : jasmine
! program : atomic_mkhmat_fullspace
!         : atomic_mkhmat_sectors
! source  : atomic_hmat.f90
! type    : subroutines
! author  : yilin wang (email: qhwyl2006@126.com)
! history : 07/09/2014 by yilin wang
! purpose : make Hamltonian matrices
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> make atomic Hamiltonian for the full space
subroutine atomic_mkhmat_fullspace()
    use constants,         only: czero, epst
    use control,           only: norbs, ncfgs
    use m_basis_fullspace, only: dec_basis, index_basis, bin_basis
    use m_spmat,           only: eimpmat, cumat
    use m_glob_fullspace,  only: hmat

    implicit none

    ! local variables
    ! loop index
    integer :: i,j
    integer :: iorb
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
                    if (ibas == 0) stop "error while determining row1"

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
                    if (ibas == 0) stop "error while determining row3"

                    hmat(ibas,jbas) = hmat(ibas,jbas) + cumat(alpha,betta,delta,gamma) * (-1.0d0)**sgn

                endif ! back if ((code(delta) == 1) .and. (code(gamma) == 1)) block
            endif ! back if ((code(alpha) == 0) .and. (code(betta) == 0)) block

        enddo gammaloop ! over gamma={1,norbs-1} loop
        enddo deltaloop ! over delta={gamma+1,norbs} loop
        enddo bettaloop ! over betta={alpha+1,norbs} loop
        enddo alphaloop ! over alpha={1,norbs-1} loop
    enddo  

    return
end subroutine atomic_mkhmat_fullspace

subroutine atomic_mkhmat_sectors()
    use constants,         only: dp, czero, epst
    use control,           only: norbs, ncfgs
    use m_basis_fullspace, only: dec_basis, index_basis, bin_basis
    use m_spmat,           only: eimpmat, cumat
    use m_glob_sectors,    only: nsectors, sectors

    implicit none

    ! local variables
    ! loop index
    integer :: i,j
    integer :: isect
    integer :: iorb, jorb
    integer :: ibas, jbas
    integer :: alpha, betta
    integer :: delta, gamma
    ! sign change due to fermion anti-commute relation
    integer :: isgn
    ! new basis state after four fermion operation
    integer :: knew
    ! binary form of a Fock state
    integer :: code(norbs)
    ! whether in some sector
    logical :: insect
     
    do isect=1, nsectors
        sectors(isect)%myham = czero

        !---------------------------------------------------------------------------------------!
        ! two fermion operators
        do jbas=1,sectors(isect)%ndim

            alploop: do alpha=1,norbs
            betloop: do betta=1,norbs

                isgn = 0
                knew = dec_basis(sectors(isect)%mybasis(jbas))
                code(1:norbs) = bin_basis(1:norbs, sectors(isect)%mybasis(jbas))

                if ( abs(eimpmat(alpha, betta)) .lt. epst ) cycle

                ! simulate one annihilation operator
                if (code(betta) == 1) then
                    do i=1,betta-1
                        if (code(i) == 1) isgn = isgn + 1
                    enddo 
                    code(betta) = 0

                    ! simulate one creation operator
                    if (code(alpha) == 0) then
                        do i=1,alpha-1
                            if (code(i) == 1) isgn = isgn + 1
                        enddo
                        code(alpha) = 1

                        ! determine the row number and hamiltonian matrix elememt
                        knew = knew - 2**(betta-1)
                        knew = knew + 2**(alpha-1)
                        isgn  = mod(isgn, 2)
                        ibas = index_basis(knew)
                        if (ibas == 0) stop "error while determining row1"

                        insect = .false.
                        do i=1, sectors(isect)%ndim 
                            if (sectors(isect)%mybasis(i) == ibas) then
                                ibas = i
                                insect = .true.
                            endif
                        enddo

                        if (insect) then
                            sectors(isect)%myham(ibas,jbas) = sectors(isect)%myham(ibas,jbas) + &
                                                           eimpmat(alpha, betta) * (-1.0d0)**isgn 
                        endif

                    endif ! back if (code(alpha) == 0) block
                endif ! back if (betta == 1) block

            enddo betloop ! over betta={1,norbs} loop
            enddo alploop ! over alpha={1,norbs} loop
        enddo ! over jbas={1,sectors(isect)%ndim} loop
        !---------------------------------------------------------------------------------------!

        !---------------------------------------------------------------------------------------!
        ! four fermion operators
        do jbas=1,sectors(isect)%ndim
            alphaloop : do alpha=1,norbs
            bettaloop : do betta=1,norbs
            gammaloop : do gamma=1,norbs
            deltaloop : do delta=1,norbs

                isgn = 0
                knew = dec_basis(sectors(isect)%mybasis(jbas))
                code(1:norbs) = bin_basis(1:norbs, sectors(isect)%mybasis(jbas))

                ! very important if single particle basis has been rotated
                if ((alpha .eq. betta) .or. (delta .eq. gamma)) cycle
                if ( abs(cumat(alpha,betta,delta,gamma)) .lt. epst ) cycle

                ! simulate two annihilation operators
                if ((code(delta) == 1) .and. (code(gamma) == 1)) then
                    do i=1,gamma-1
                        if(code(i) == 1) isgn = isgn + 1
                    enddo 
                    code(gamma) = 0

                    do i=1,delta-1
                        if(code(i) == 1) isgn = isgn + 1
                    enddo 
                    code(delta) = 0

                    ! simulate two creation operators
                    if ((code(alpha) == 0) .and. (code(betta) == 0)) then
                        do i=1,betta-1
                            if(code(i) == 1) isgn = isgn + 1
                        enddo 
                        code(betta) = 1

                        do i=1,alpha-1
                            if(code(i) == 1) isgn = isgn + 1
                        enddo
                        code(alpha) = 1

                        ! determine the row number and hamiltonian matrix elememt
                        knew = knew - 2**(gamma-1) - 2**(delta-1)
                        knew = knew + 2**(betta-1) + 2**(alpha-1)
                        ibas = index_basis(knew)
                        isgn = mod(isgn, 2)

                        insect = .false.
                        do i=1, sectors(isect)%ndim 
                            if (sectors(isect)%mybasis(i) == ibas) then
                                ibas = i
                                insect = .true.
                            endif
                        enddo

                        if (insect) then
                            sectors(isect)%myham(ibas,jbas) = sectors(isect)%myham(ibas,jbas) + &
                                                 cumat(alpha,betta,delta,gamma) * (-1.0d0)**isgn
                        endif

                    endif ! back if ((code(delta) == 1) .and. (code(gamma) == 1)) block
                endif ! back if ((code(alpha) == 0) .and. (code(betta) == 0)) block

            enddo deltaloop ! over delta={gamma+1,norbs} loop
            enddo gammaloop ! over gamma={1,norbs-1} loop
            enddo bettaloop ! over betta={alpha+1,norbs} loop
            enddo alphaloop ! over alpha={1,norbs-1} loop
        enddo ! over jbas={1,sectors(isect)%ndim} loop
        !---------------------------------------------------------------------------------------!

    enddo ! over i={1, nsectors}

    return
end subroutine atomic_mkhmat_sectors

subroutine atomic_diag_hmat_sectors()
    use m_glob_sectors

    implicit none

    ! local variables
    integer :: i
   
    do i=1, nsectors
        call diag_one_sector(sectors(i)%ndim, sectors(i)%myham, &
                               sectors(i)%myeigval, sectors(i)%myeigvec)
    enddo

    return
end subroutine atomic_diag_hmat_sectors

subroutine diag_one_sector(ndim, amat, eval, evec)
     use constants, only: dp
     implicit none

     ! external variables
     ! the order of the matrix amat
     integer, intent(in) :: ndim

     ! original real symmetric matrix to compute eigenval and eigenvector
     complex(dp), intent(in) :: amat(ndim, ndim)

     ! if info = 0, the eigenvalues in ascending order.
     real(dp), intent(out) :: eval(ndim)

     ! if info = 0, orthonormal eigenvectors of the matrix A
     real(dp), intent(out) :: evec(ndim, ndim)

     ! local variables
     integer :: i, j
     real(dp) :: hmat(ndim, ndim)

     do i=1, ndim
         do j=1, ndim
             hmat(i,j) = real(amat(i, j))
         enddo
     enddo

     call dmat_dsyev(ndim, ndim, hmat, eval, evec)

     return
end subroutine diag_one_sector 
