
  subroutine atomic_find_subspace(sind, find, sector_size, sector_basis)
     use control, only : ncfgs

     implicit none

     integer, intent(out) :: sind
     integer, intent(in) :: find
     integer, intent(in) :: sector_size(ncfgs)
     integer, intent(in) :: sector_basis(ncfgs,ncfgs)

     integer :: m
     integer :: n

     sind = 0

     m1loop: do m=1,ncfgs
         do n=1,sector_size(m)
             if ( sector_basis(m,n) == find ) then
                 sind = m
                 EXIT m1loop
             endif
         enddo
     enddo m1loop

     call s_assert(sind /= 0)

     return
  end subroutine atomic_find_subspace

  subroutine atomic_test_ad()
     use constants, only : dp
     use constants, only : mystd

     use control, only : nband, norbs, ncfgs

     use m_fock, only : hmat
     use m_fock, only : bin_basis, dec_basis, ind_basis

     implicit none

     integer, parameter :: max_mapping = 200000

     integer :: i
     integer :: j
     integer :: m
     integer :: jold, jnew, iorb
     integer :: ia, ib, isgn
     integer :: iup, idn
     integer :: HA, HL, HU
     integer :: N, Sz, Jz, Ap

     integer, allocatable :: sector_size(:), sector_size_(:)
     integer, allocatable :: sector_basis(:,:), sector_basis_(:,:)
     integer, allocatable :: sector_N(:), sector_Sz(:), sector_Jz(:), sector_Ap(:)
     integer, allocatable :: Mup(:,:)
     integer, allocatable :: Mdn(:,:)

     allocate(sector_size(ncfgs))
     allocate(sector_size_(ncfgs))
     allocate(sector_basis(ncfgs,ncfgs))
     allocate(sector_basis_(ncfgs,ncfgs))
     allocate(Mup(max_mapping,2))
     allocate(Mdn(max_mapping,2))
     allocate(sector_N(ncfgs))
     allocate(sector_Sz(ncfgs))
     allocate(sector_Jz(ncfgs))
     allocate(sector_Ap(ncfgs))

     ! initialization
     sector_size = 1
     sector_basis = 0
     do i=1,ncfgs
         sector_basis(i,1) = i
     enddo
     sector_N = 0
     sector_Sz = 0
     sector_Jz = 0
     sector_Ap = 0

     ! phase 1
     do i=1,ncfgs
         do j=1,ncfgs
             if ( abs(hmat(i,j)) > 0.0_dp ) then

                 call atomic_find_subspace(ia, i, sector_size, sector_basis)
                 call atomic_find_subspace(ib, j, sector_size, sector_basis)

                 if ( ia /= ib ) then
                     call s_assert(sector_size(ia) >= 1)
                     call s_assert(sector_size(ib) >= 1)
                     !
                     do m=1,sector_size(ib)
                         sector_basis(ia, sector_size(ia) + m) = sector_basis(ib,m)
                     enddo
                     !
                     sector_size(ia) = sector_size(ia) + sector_size(ib)
                     sector_size(ib) = 0
                     sector_basis(ib,:) = 0
                 endif

             endif
         enddo
     enddo

     sector_size_ = 0
     sector_basis_ = 0
     m = 0
     do i=1,ncfgs
         if ( sector_size(i) > 0 ) then
             m = m + 1
             sector_size_(m) = sector_size(i)
             sector_basis_(m,:) = sector_basis(i,:)
         endif
     enddo
     sector_size = sector_size_
     sector_basis = sector_basis_

     ! print subspaces
     m = 0
     do i=1,ncfgs
         if ( sector_size(i) > 0 ) then
             m = m + 1
             write(mystd,'(a,i6)') 'subspace -> ', m
             write(mystd,'(a,i6)') 'size :', sector_size(i)
             write(mystd,'(a)') 'basis :'
             do j=1,sector_size(i)
                 write(mystd,'(14i1)') bin_basis(:,sector_basis(i,j))
             enddo
             print *
         endif
     enddo


     DO iorb = 1,norbs

     ! create mapping
     iup = 0
     idn = 0
     Mup = 0
     Mdn = 0

     do i=1,ncfgs
         call atomic_find_subspace(ia, i, sector_size, sector_basis)

         ! c^+
         if ( bin_basis(iorb,i) == 0 ) then
             jold = dec_basis(i)
             call atomic_make_cdagger(iorb, jold, jnew, isgn)
             j = ind_basis(jnew)

             call atomic_find_subspace(ib, j, sector_size, sector_basis)

             iup = iup + 1
             call s_assert(iup <= max_mapping)
             Mup(iup,1) = ia
             Mup(iup,2) = ib
         endif

         ! c
         if ( bin_basis(iorb,i) == 1 ) then
             jold = dec_basis(i)
             call atomic_make_c(iorb, jold, jnew, isgn)
             j = ind_basis(jnew)

             call atomic_find_subspace(ib, j, sector_size, sector_basis)

             idn = idn + 1
             call s_assert(idn <= max_mapping)
             Mdn(idn,1) = ia
             Mdn(idn,2) = ib
         endif
     enddo

     print *, '# orb: ', iorb
     print *, 'number of Mup:', iup
     print *, 'number of Mdn:', idn

     !
     !do i=1,iup
     !    print *, i, Mup(i,:), Mdn(i,:)
     !enddo

     do i=1,max_mapping
         HA = Mup(i,1)
         HL = Mup(i,1)
         HU = Mup(i,2)

         if ( HL /= 0 .and. HU /= 0 ) then
         call zigzag(1, HA, HL, HU, sector_size, sector_basis, Mup, Mdn, max_mapping)
         endif
     enddo

     enddo

     ! print subspaces
     m = 0
     do i=1,ncfgs
         if ( sector_size(i) > 0 ) then
             m = m + 1
             write(mystd,'(a,i6)') 'subspace -> ', m
             write(mystd,'(a,i6)') 'size :', sector_size(i)
             write(mystd,'(a)') 'basis :'
             do j=1,sector_size(i)
                 write(mystd,'(i,2X,14i1)') j, bin_basis(:,sector_basis(i,j))
             enddo
             call atomic_sector_N(N, sector_size(i), sector_basis(i,:))
             call atomic_sector_Sz(Sz, sector_size(i), sector_basis(i,:))
             !call atomic_sector_Jz(Jz, sector_size(i), sector_basis(i,:))
             Ap = 1
             do j=1,m-1
                 if ( ( sector_N(j) == N ) .and. ( sector_Sz(j) == Sz ) ) then
                     Ap = Ap + 1
                 endif
             enddo
             write(mystd, '(a, i3)') 'N :', N
             write(mystd, '(a, i3)') 'Sz:', Sz
             !write(mystd, '(a, i3)') 'Jz:', Jz
             write(mystd, '(a, i3)') 'AP:', Ap
             sector_N(m) = N
             sector_Sz(m) = Sz
             !sector_Jz(m) = Jz
             sector_Ap(m) = Ap
             print *
         endif
     enddo

     STOP

     return
  end subroutine atomic_test_ad

  subroutine atomic_sector_N(GQN_N, sector_size, sector_basis)
     use control, only : norbs, ncfgs
     use m_fock, only : bin_basis

     implicit none

     integer, intent(in) :: sector_size
     integer, intent(in) :: sector_basis(ncfgs)
     integer, intent(out) :: GQN_N

     integer :: i
     integer :: basis(norbs)
     integer :: N

     GQN_N = 999
     do i=1,sector_size
         basis = bin_basis(:,sector_basis(i))
         N = sum(basis)
         if ( i == 1 ) then
             GQN_N = N
         else
             if ( N /= GQN_N ) then
                 STOP "wrong in GQN(N)"
             endif
         endif
     enddo

     return
  end subroutine atomic_sector_N

  subroutine atomic_sector_Sz(GQN_Sz, sector_size, sector_basis)
     use control, only : nband, norbs, ncfgs
     use m_fock, only : bin_basis

     implicit none

     integer, intent(in) :: sector_size
     integer, intent(in) :: sector_basis(ncfgs)
     integer, intent(out) :: GQN_Sz

     integer :: i
     integer :: basis(norbs)
     integer :: Sz

     GQN_Sz = 999
     do i=1,sector_size
         basis = bin_basis(:,sector_basis(i))
         Sz = sum(basis(1:nband)) - sum(basis(nband+1:norbs))
         if ( i == 1 ) then
             GQN_Sz = Sz
         else
             if ( Sz /= GQN_Sz ) then
                 STOP "wrong in GQN(Sz)"
             endif
         endif
     enddo

     return
  end subroutine atomic_sector_Sz

  subroutine atomic_sector_Jz(GQN_Jz, sector_size, sector_basis)
     use control, only : nband, norbs, ncfgs
     use m_fock, only : bin_basis

     implicit none

     integer, intent(in) :: sector_size
     integer, intent(in) :: sector_basis(ncfgs)
     integer, intent(out) :: GQN_Jz

     integer :: i, k
     integer :: basis(norbs)
     integer :: Jz
     integer :: good_jz(norbs)

     call atomic_make_gjz(good_jz)

     GQN_Jz = 999
     do i=1,sector_size
         basis = bin_basis(:,sector_basis(i))

         Jz = 0
         do k=1,norbs
             Jz = Jz + good_jz(k) * basis(k)
         enddo ! over k={1,norbs} loop

         if ( i == 1 ) then
             GQN_Jz = Jz
         else
             if ( Jz /= GQN_Jz ) then
                 STOP "wrong in GQN(Jz)"
             endif
         endif
     enddo

     return
  end subroutine atomic_sector_Jz

recursive &
  subroutine zigzag(up_or_down, HA, HL, HU, sector_size, sector_basis, Mup, Mdn, max_mapping)
     use control, only : ncfgs

     implicit none

     integer, intent(in) :: up_or_down
     integer, intent(in) :: HA
     integer, intent(in) :: HL
     integer, intent(in) :: HU
     integer, intent(in) :: max_mapping
     integer, intent(inout) :: sector_size(ncfgs)
     integer, intent(inout) :: sector_basis(ncfgs,ncfgs)
     integer, intent(inout) :: Mup(max_mapping,2)
     integer, intent(inout) :: Mdn(max_mapping,2)

     integer :: i, j
     integer :: HB

     if ( up_or_down == 1 ) then
         do i=1,max_mapping
             if ( Mup(i,1) == HA ) then
                 HB = Mup(i,2)
                 Mup(i,:) = 0

                 if ( HB /= HU ) then
                     print *, 'merge up:', HA, HB, HU
                     do j=1,sector_size(HB)
                         sector_basis(HU, sector_size(HU) + j) = sector_basis(HB,j)
                     enddo
                     sector_size(HU) = sector_size(HU) + sector_size(HB)
                     sector_size(HB) = 0
                 endif

                 call zigzag(2, HB, HL, HU, sector_size, sector_basis, Mup, Mdn, max_mapping)
             endif
         enddo 
     else
         do i=1,max_mapping
             if ( Mdn(i,1) == HA ) then
                 HB = Mdn(i,2)
                 Mdn(i,:) = 0

                 if ( HB /= HL ) then
                     print *, 'merge dn:', HA, HB, HL
                     do j=1,sector_size(HB)
                         sector_basis(HL, sector_size(HL) + j) = sector_basis(HB,j)
                     enddo
                     sector_size(HL) = sector_size(HL) + sector_size(HB)
                     sector_size(HB) = 0
                 endif

                 call zigzag(1, HB, HL, HU, sector_size, sector_basis, Mup, Mdn, max_mapping)
             endif
         enddo 
     endif

     return
  end subroutine zigzag