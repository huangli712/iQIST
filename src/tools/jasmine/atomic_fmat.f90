!-------------------------------------------------------------------------
! project : jasmine
! program : atomic_make_annifmat_fullspace
!         : atomic_make_fmat_sectors
!         : rotate_fmat
!         : atomic_construct
!         : atomic_eliminate 
! source  : atomic_fmat.f90
! type    : subroutines
! author  : yilin wang (email: qhwyl2006@126.com)
! history : 07/09/2014 by yilin wang
! purpose : make fmat
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> make fmat for annihilation operators for full space
subroutine atomic_make_annifmat_fullspace()
    use control, only: norbs, ncfgs
    use m_glob_fullspace, only: anni_fmat, hmat_eigvec
    use m_basis_fullspace, only: dec_basis, index_basis

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
               call atomic_eliminate(i, right, left, isgn)
               k = index_basis(left)
               anni_fmat(k, j, i) = dble(isgn)
            endif
        enddo 
    enddo 

    ! rotate it to the atomic eigenvector basis
    do i=1, norbs
        call atomic_tran_represent_real(ncfgs, anni_fmat(:,:,i), hmat_eigvec)
    enddo

    return
end subroutine atomic_make_annifmat_fullspace

!>>> build fmat for good quantum algorithm
subroutine atomic_make_fmat_sectors()
     use control, only: norbs
     use m_basis_fullspace, only: dec_basis, index_basis
     use m_glob_sectors 
     use m_sector
     implicit none

     ! local variables
     ! loop index ober orbits
     integer :: iorb
     integer :: ifermi
     integer :: isect, jsect
     integer :: ibas, jbas
     integer :: i
     ! sign change due to commute relation
     integer :: isgn
     ! auxiliary integer variables
     integer :: jold, jnew

     do isect=1, nsectors
         do iorb=1,norbs
             do ifermi=0, 1 
                 jsect = sectors(isect)%next_sector(iorb, ifermi) 
                 if (jsect == -1) cycle
                 ! allocate memory for fmat
                 sectors(isect)%myfmat(iorb, ifermi)%n = sectors(jsect)%ndim
                 sectors(isect)%myfmat(iorb, ifermi)%m = sectors(isect)%ndim
                 call alloc_one_fmat(sectors(isect)%myfmat(iorb, ifermi))

                 ! build fmat
                 do jbas=1, sectors(isect)%ndim
                     jold = dec_basis(sectors(isect)%mybasis(jbas))
                     if (ifermi == 1 .and. ( btest(jold, iorb-1) .eqv. .false. )) then
                         call atomic_construct(iorb, jold, jnew, isgn)
                     elseif (ifermi == 0 .and. ( btest(jold, iorb-1) .eqv. .true. )) then
                         call atomic_eliminate(iorb, jold, jnew, isgn)
                     endif

                     ibas = index_basis(jnew)
                     do i=1, sectors(jsect)%ndim 
                         if (ibas == sectors(jsect)%mybasis(i)) then
                             ibas = i
                             exit
                         endif
                     enddo
                     sectors(isect)%myfmat(iorb, ifermi)%item(ibas, jbas) = dble(isgn)
                 enddo 
                 call rotate_fmat(sectors(jsect)%ndim, sectors(isect)%ndim, sectors(jsect)%myeigvec, &
                     sectors(isect)%myfmat(iorb, ifermi)%item, sectors(isect)%myeigvec)

             enddo ! over ifermi={0,1} loop
         enddo ! over iorb={1, norbs} loop
     enddo ! over isect={1,nsectors} loop

     return
end subroutine atomic_make_fmat_sectors

! calculate A^T * B * C
subroutine rotate_fmat(ndimx, ndimy, amat, bmat, cmat)
    use constants, only: dp, zero
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

    amat_t = transpose(amat)
    tmp_mat = zero
    call dmat_dgemm(ndimx, ndimy, ndimy, bmat, cmat,  tmp_mat)
    call dmat_dgemm(ndimx, ndimx, ndimy, amat_t,  tmp_mat, bmat)

    return
end subroutine rotate_fmat

!>>> create one electron on ipos of |jold) to deduce |jnew)
subroutine atomic_construct(ipos, jold, jnew, isgn)
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
        stop "severe error happened in atomic_construct"
    endif ! back if (btest(jold, ipos-1) .eqv. .true.) block

    isgn = 0
    do iorb=1,ipos-1
       if (btest(jold, iorb-1)) isgn = isgn + 1
    enddo ! over i={1,ipos-1} loop
    isgn = mod(isgn, 2)

    isgn = (-1)**isgn
    jnew = jold + 2**(ipos-1)

    return
end subroutine atomic_construct

!>>> destroy one electron on ipos of |jold) to deduce |jnew)
subroutine atomic_eliminate(ipos, jold, jnew, isgn)
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
    ! loop index over orbit
    integer :: iorb

    if (btest(jold, ipos-1) .eqv. .false.) then
        stop "severe error happened in atomic_eliminate"
    endif ! back if (btest(jold, ipos-1) .eqv. .false.) block

    isgn = 0
    do iorb=1,ipos-1
        if (btest(jold, iorb-1)) isgn = isgn + 1
    enddo ! over i={1,ipos-1} loop
    isgn = mod(isgn, 2)

    isgn = (-1)**isgn
    jnew = jold - 2**(ipos-1)

    return
end subroutine atomic_eliminate

