!>>> build fmat for full space
  subroutine atomic_make_annifmat_fullspace()
     use constants
     use control
     use m_mpmat_fullspace
     use m_basis_fullspace

     implicit none

! local variables
! loop index
     integer :: i,j,k

! left Fock state
     integer :: left

! right Fock state
     integer :: right

     do i=1,norbs
         do j=1,ncfgs
             right = dec_basis(j)
             if (btest(right, i-1) .eqv. .true.) then
                call atomic_eliminate(i, right, left, isgn)
                k = index_basis(left)
                mp_anni_mat(k, j, i) = dble(isgn)
             endif
         enddo 
     enddo 

! rotate it to the atomic eigenvector basis
     do i=1, norbs
         call atomic_tran_represent_real(ncfgs, mp_anni_mat(:,:,i), mp_hmat_eigvec)
     enddo

     return
  end subroutine atomic_make_annifmat_fullspace

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


!>>> build matrix for construct operator on Fock basis
subroutine build_cfmat_fock_njz()
     use constants
     use control
     use mod_fmat
     use mod_global
     implicit none

     ! local variables
     ! loop index ober orbits
     integer :: iorb
     ! sign change due to commute relation
     integer :: isgn
     ! auxiliary integer variables
     integer :: jold, jnew
     ! loop index over configurations
     integer :: isub, jsub
     integer :: ibas, jbas
     integer :: i

     ! setup matrix for construct operator
     do iorb=1,norbs
         do isub=1, nsubs
             if (c_towhich(isub, iorb) == -1) cycle
             c_fmat(isub, iorb)%iorb = iorb
             c_fmat(isub, iorb)%x = isub
             c_fmat(isub, iorb)%y = c_towhich(isub,iorb)
             c_fmat(isub, iorb)%ndimx = subspaces(isub)%ndim
             c_fmat(isub, iorb)%ndimy = subspaces(c_towhich(isub,iorb))%ndim

             call alloc_one_fmat(c_fmat(isub, iorb))
             c_fmat(isub, iorb)%elem = zero

             ! build elements of c_fmat(isub, iorb)
             do jbas=1, subspaces(isub)%ndim
                 jold = basis(subspaces(isub)%mybasis(jbas))
                 if (btest(jold, iorb-1) .eqv. .false.) then
                     call atomic_construct(iorb, jold, jnew, isgn)
                     ibas = invsn(jnew)
                     ! check which basis in c_towhich(isub,iorb) subspace
                     do i=1, subspaces(c_towhich(isub,iorb))%ndim 
                         if (ibas == subspaces(c_towhich(isub,iorb))%mybasis(i)) then
                             ibas = i
                         endif
                     enddo
                     c_fmat(isub, iorb)%elem(ibas, jbas) = dble(isgn)
                 endif
             enddo
             ! end building elements of c_fmat(isub, iorb)
         enddo ! over isub={1, nusbs} loop
     enddo ! over iorb={1,norbs} loop

     return
end subroutine build_cfmat_fock_njz

subroutine cfmat_fock2eigen_basis_njz()
     use constants
     use control
     use mod_fmat
     use mod_global
     implicit none

     ! local variables
     ! loop index ober orbits
     integer :: iorb
     ! loop index over configurations
     integer :: isub, jsub
     ! dimension of isub and jsub
     integer :: ndimx, ndimy
     integer :: i

     do iorb=1, norbs
         do isub=1, nsubs
             jsub = c_towhich(isub, iorb)
             ! if x --> y subspaces exists
             if (jsub /= -1) then
                 ndimx = c_fmat(isub, iorb)%ndimx 
                 ndimy = c_fmat(isub, iorb)%ndimy 
                 call rotate_fmat(ndimy, ndimx, subspaces(jsub)%myeigvec, &
                     c_fmat(isub,iorb)%elem, subspaces(isub)%myeigvec)
             endif
         enddo
     enddo

     return
end subroutine cfmat_fock2eigen_basis_njz

! calculate A^T * B * C
subroutine rotate_fmat(ndimx, ndimy, amat, bmat, cmat)
    use constants
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

subroutine atomic_build_cfmat_njz()
    use mod_global 
    implicit none

    ! first, allocate memory of cfmat
    call alloc_glob_fmat()

    ! then, build fmat on fock basis
    call build_cfmat_fock_njz()

    ! finally, rotate it to eigenvalue basis
    call cfmat_fock2eigen_basis_njz()

    return
end subroutine atomic_build_cfmat_njz

subroutine atomic_dump_fmtrx(norbs, ncfgs, eval, nval, fmat)
     use constants

     implicit none

! number of orbits
     integer, intent(in) :: norbs

! number of configurations
     integer, intent(in) :: ncfgs

     real(dp), intent(in) :: eval(ncfgs)
     real(dp), intent(in) :: nval(ncfgs)

! matrix form of construct operator
     complex(dp), intent(in) :: fmat(ncfgs, ncfgs, norbs)

! local variables
     integer :: i, j, k
     integer :: icnt(5)

     real(dp) :: sxyz 
     sxyz = 0.0d0

! open data file: atom.cix
     open(mytmp, file='atom.cix', form='formatted', status='unknown')

! write eigenvalues
     write(mytmp,'(a)') '# eigenvalues: index | energy | occupy | spin'
     do i=1,ncfgs
         write(mytmp,'(i5,3f22.15)') i, eval(i), nval(i), sxyz
     enddo ! over i={1,ncfgs} loop

! write f matrix element
     write(mytmp,'(a)') '# f matrix element: alpha | beta | orbital | fmat'
     do i=1,norbs
         do j=1,ncfgs
             do k=1,ncfgs
                 write(mytmp,'(3i5,2f22.15)') k, j, i, (fmat(k,j,i))
             enddo ! over k={1,ncfgs} loop
         enddo ! over j={1,ncfgs} loop
     enddo ! over i={1,norbs} loop

! close input file
     close(mytmp)

     icnt = 0
     do j=1,ncfgs
         do k=1,ncfgs
             if (abs(fmat(j,k,1)) .lt. 1.0d-9) then
                 icnt(1) = icnt(1) + 1
             else if (abs(fmat(j,k,1)) .lt. 1.0d-3) then
                 icnt(2) = icnt(2) + 1
             else if (abs(fmat(j,k,1)) .lt. 1.0d-2) then
                 icnt(3) = icnt(3) + 1
             else if (abs(fmat(j,k,1)) .lt. 1.0d-1) then
                 icnt(4) = icnt(4) + 1
             else if (abs(fmat(j,k,1)) .lt. 1.0d0) then
                 icnt(5) = icnt(5) + 1
             endif
        enddo
   enddo
   print*, "< 1.0d-9: ", icnt(1)
   print*, "< 1.0d-3: ", icnt(2)
   print*, "< 1.0d-2: ", icnt(3)
   print*, "< 1.0d-1: ", icnt(4)
   print*, "< 1.0d+0: ", icnt(5)
   print*, sum(icnt)

     return
end subroutine atomic_dump_fmtrx
