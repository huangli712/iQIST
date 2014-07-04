!=========================================================================!
! project : strawberry
! program : atomic_build_Hmtrx
! history : 09/29/2011
! author  : xidai and duliang (email:duleung@gmail.com)
! purpose : construct atomic hamiltonian matrix
! comment : 
!=========================================================================!
subroutine atomic_make_hmtrx()
    use constants
    use control
    use mod_global 
    implicit none

    ! local variables
    ! loop index
    integer :: isub
    integer :: iorb
    integer :: jorb
    integer :: i, j

    ! sign change due to fermion anti-commute relation
    integer :: isgn

    ! loop index over Fock basis
    integer :: ibas, jbas

    ! new basis state after four fermion operation
    integer :: knew

    ! index for general interaction matrix
    integer :: alpha, betta
    integer :: delta, gamma

    ! binary code representation of a state
    integer :: code(norbs)

    ! auxiliary complex(dp) variables
    complex(dp) :: ztmpa

    ! whether in this subspace
    logical :: insub
     
    do isub=1, nsubs
        !=========================================================================!
        ! two fermion operator terms (crystalline electric field)
        !=========================================================================!
        subspaces(isub)%myham = czero

        do jbas=1,subspaces(isub)%ndim
            alploop: do alpha=1,norbs
            betloop: do betta=1,norbs

                isgn = 0
                knew = basis(subspaces(isub)%mybasis(jbas))
                code(1:norbs) = invcd(1:norbs, subspaces(isub)%mybasis(jbas))
                if ( abs(cemat(alpha, betta)) .lt. eps6 ) cycle

                ! simulate one eliminate operator
                if (code(betta) == 1) then
                    do i=1,betta-1
                        if (code(i) == 1) isgn = isgn + 1
                    enddo 
                    code(betta) = 0

                    ! simulate one construct operator
                    if (code(alpha) == 0) then
                        do i=1,alpha-1
                            if (code(i) == 1) isgn = isgn + 1
                        enddo
                        code(alpha) = 1

                        ! determine the column number and hamiltonian matrix elememt
                        knew = knew - 2**(betta-1)
                        knew = knew + 2**(alpha-1)
 
                        isgn  = mod(isgn, 2)
                        ibas = invsn(knew)
                        insub = .false.
                        do i=1, subspaces(isub)%ndim 
                            if (subspaces(isub)%mybasis(i) == ibas) then
                                ibas = i
                                insub = .true.
                            endif
                        enddo
                        if (ibas == 0) stop "error while determining row1"
                        if (insub) then
                            subspaces(isub)%myham(ibas,jbas) = subspaces(isub)%myham(ibas,jbas) + &
                                                           cemat(alpha,betta) * (-1.0d0)**isgn 
                        endif

                    endif ! back if (code(alpha) == 0) block
                endif ! back if (betta == 1) block

            enddo betloop ! over betta={1,norbs} loop
            enddo alploop ! over alpha={1,norbs} loop
        enddo ! over jbas={1,nbas} loop

        ! four fermion terms in local Hamiltonian
        do jbas=1,subspaces(isub)%ndim
           alphaloop : do alpha=1,norbs
           bettaloop : do betta=1,norbs
           gammaloop : do gamma=1,norbs
           deltaloop : do delta=1,norbs

               isgn = 0
               knew = basis(subspaces(isub)%mybasis(jbas))
               code(1:norbs) = invcd(1:norbs, subspaces(isub)%mybasis(jbas))
               !# very important if single particle basis rotated
               if ((alpha .eq. betta) .or. (delta .eq. gamma)) cycle
               if ( abs(cumat(alpha,betta,delta,gamma)) .lt. eps6 ) cycle

               ! simulate two eliminate operator
               if ((code(delta) == 1) .and. (code(gamma) == 1)) then
                   do i=1,gamma-1
                       if(code(i) == 1) isgn = isgn + 1
                   enddo ! over i={1,gamma-1} loop
                   code(gamma) = 0

                   do i=1,delta-1
                       if(code(i) == 1) isgn = isgn + 1
                   enddo ! over i={1,delta-1} loop
                   code(delta) = 0

                   ! simulate two construct operator
                   if ((code(alpha) == 0) .and. (code(betta) == 0)) then
                       do i=1,betta-1
                           if(code(i) == 1) isgn = isgn + 1
                       enddo ! over i={1,betta-1} loop
                       code(betta) = 1

                       do i=1,alpha-1
                           if(code(i) == 1) isgn = isgn + 1
                       enddo ! over i={1,alpha-1} loop
                       code(alpha) = 1

                       ! determine the column number and hamiltonian matrix elememt
                       knew = knew - 2**(gamma-1) - 2**(delta-1)
                       knew = knew + 2**(betta-1) + 2**(alpha-1)

                       ibas = invsn(knew)
                       isgn = mod(isgn, 2)
                       insub = .false.
                       do i=1, subspaces(isub)%ndim 
                           if (subspaces(isub)%mybasis(i) == ibas) then
                               ibas = i
                               insub = .true.
                           endif
                       enddo
                       if (insub) then
                           subspaces(isub)%myham(ibas,jbas) = subspaces(isub)%myham(ibas,jbas) + &
                                                cumat(alpha,betta,delta,gamma) * (-1.0d0)**isgn
                       endif

                       ! simplly check the fermion anti-commute relation
                       if ( isgn /= 0 ) then
                       !$  stop "something wrong in atomic_make_hmat, pls check carefull"
                       endif ! back if ( sgn /= 0 ) block

                   endif ! back if ((code(delta) == 1) .and. (code(gamma) == 1)) block
               endif ! back if ((code(alpha) == 0) .and. (code(betta) == 0)) block

           enddo deltaloop ! over delta={gamma+1,norbs} loop
           enddo gammaloop ! over gamma={1,norbs-1} loop
           enddo bettaloop ! over betta={alpha+1,norbs} loop
           enddo alphaloop ! over alpha={1,norbs-1} loop
        enddo ! over jbas={1,ncfgs} loop

    enddo ! over i={1, nusbs}

    open(mytst, file='test-ham.dat', status='unknown')
    do isub=1, nsubs
        do jbas=1,subspaces(isub)%ndim
            do ibas=1,subspaces(isub)%ndim
                if (abs(subspaces(isub)%myham(ibas, jbas)) .lt. eps6) cycle
                    write(mytst, '(3i8, 2f17.10)') isub, subspaces(isub)%mybasis(ibas), &
                         subspaces(isub)%mybasis(jbas), subspaces(isub)%myham(ibas, jbas)
            enddo ! over ibas={1,ncfgs} loop
        enddo ! over jbas={1,ncfgs} loop
    enddo
    close(mytst)
    return
end subroutine atomic_make_hmtrx

subroutine atomic_diag_hmtrx()
    use constants
    use control
    use mod_global 
    implicit none

    ! local variables
    integer :: isub
   
    do isub=1, nsubs
        call diag_one_subspace(subspaces(isub)%ndim, subspaces(isub)%myham, &
                               subspaces(isub)%myeigval, subspaces(isub)%myeigvec)
    enddo

    return
end subroutine atomic_diag_hmtrx

subroutine diag_one_subspace(ndim, amat, eval, evec)
     use constants
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
end subroutine diag_one_subspace 


