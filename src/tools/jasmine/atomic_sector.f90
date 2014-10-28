!!!-----------------------------------------------------------------------
!!! project : jasmine
!!! program : atomic_make_sfmat
!!!           atomic_make_shmat
!!!           atomic_diag_shmat
!!!           atomic_make_sectors
!!! source  : atomic_sector.f90
!!! type    : subroutines
!!! author  : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!!           10/27/2014 by li huang
!!! purpose :
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> atomic_make_sfmat: build F-matrix (fmat) for good quantum numbers
!!>>> (GQNs) algorithm
  subroutine atomic_make_sfmat()
     use constants, only : zero

     use control, only : norbs
     use m_full, only : dec_basis, index_basis
     use m_sector, only : nsectors, sectors
     use m_sector, only : alloc_one_fmat

     implicit none

! local variables
! loop index
     integer :: i
     integer :: iorb
     integer :: ifermi
     integer :: isect
     integer :: jsect
     integer :: ibas
     integer :: jbas

! sign change due to commute relation
     integer :: isgn

! auxiliary integer variables
     integer :: jold
     integer :: jnew

! loop over all the sectors
     do isect=1,nsectors
! loop over all the orbitals
         do iorb=1,norbs
! loop over the creation and annihilation fermion operators
             do ifermi=0,1
                 jsect = sectors(isect)%next(iorb,ifermi)
                 if ( jsect == -1 ) CYCLE
! allocate memory for fmat
                 sectors(isect)%fmat(iorb,ifermi)%n = sectors(jsect)%ndim
                 sectors(isect)%fmat(iorb,ifermi)%m = sectors(isect)%ndim
                 call alloc_one_fmat(sectors(isect)%fmat(iorb,ifermi))
                 sectors(isect)%fmat(iorb,ifermi)%val = zero
! build fmat
                 do jbas=1,sectors(isect)%ndim
                     jold = dec_basis(sectors(isect)%basis(jbas))
! for creation fermion operator
                     if ( ifermi == 1 .and. ( btest(jold, iorb-1) .eqv. .false. ) ) then
                         call atomic_make_cdagger(iorb, jold, jnew, isgn)
! for annihilation fermion operator
                     elseif ( ifermi == 0 .and. ( btest(jold, iorb-1) .eqv. .true. ) ) then
                         call atomic_make_c(iorb, jold, jnew, isgn)
                     else
                         CYCLE
                     endif ! back if ( ifermi == 1 .and. ( btest(jold, iorb-1) .eqv. .false. ) ) block
                     ibas = index_basis(jnew)
                     do i=1,sectors(jsect)%ndim
                         if ( ibas == sectors(jsect)%basis(i) ) then
                             ibas = i
                             sectors(isect)%fmat(iorb, ifermi)%val(ibas,jbas) = dble(isgn)
                             EXIT
                         endif ! back if ( ibas == sectors(jsect)%basis(i) ) block
                     enddo ! over i={1,sectors(jsect)%ndim} loop
                 enddo  ! over jbas={1, sectors(isect)%ndim} loop
! roate fmat to atomic eigenstates basis
                 call atomic_tran_fmat(sectors(jsect)%ndim, &
                                       sectors(isect)%ndim, &
                                       sectors(jsect)%evec, &
                                       sectors(isect)%fmat(iorb,ifermi)%val, &
                                       sectors(isect)%evec)
             enddo ! over ifermi={0,1} loop
         enddo ! over iorb={1,norbs} loop
     enddo ! over isect={1,nsectors} loop

     return
  end subroutine atomic_make_sfmat

!!>>> atomic_make_shmat: make Hamiltonian for each sector one by one
  subroutine atomic_make_shmat()
     use constants, only : dp, one, epst, czero

     use control, only : norbs, ncfgs
     use m_full, only : dec_basis, index_basis, bin_basis
     use m_spmat, only : emat, umat
     use m_sector, only : nsectors, sectors

     implicit none

! local variables
! loop index
     integer :: i
     integer :: isect
     integer :: ibas
     integer :: jbas
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

     do isect=1,nsectors
         sectors(isect)%hmat = czero

! two fermions term
!-------------------------------------------------------------------------
         do jbas=1,sectors(isect)%ndim
             alploop: do alpha=1,norbs
                 betloop: do betta=1,norbs

                     isgn = 0
                     knew = dec_basis(sectors(isect)%basis(jbas))
                     code(1:norbs) = bin_basis(1:norbs, sectors(isect)%basis(jbas))

                     if ( abs(emat(alpha, betta)) < epst ) CYCLE

! simulate one annihilation operator
                     if ( code(betta) == 1 ) then
                         do i=1,betta-1
                             if ( code(i) == 1 ) isgn = isgn + 1
                         enddo ! over i={1,betta-1} loop
                         code(betta) = 0

! simulate one creation operator
                         if ( code(alpha) == 0 ) then
                             do i=1,alpha-1
                                 if ( code(i) == 1 ) isgn = isgn + 1
                             enddo ! over i={1,alpha-1} loop
                             code(alpha) = 1

! determine the row number and hamiltonian matrix elememt
                             knew = knew - 2**(betta-1)
                             knew = knew + 2**(alpha-1)
                             isgn = mod(isgn, 2)
                             ibas = index_basis(knew)
                             if ( ibas == 0 ) then
                                 call s_print_error('atomic_make_shmat','error while determining row')
                             endif ! back if ( ibas == 0 ) block

                             insect = .false.
                             do i=1,sectors(isect)%ndim
                                 if ( sectors(isect)%basis(i) == ibas ) then
                                     ibas = i; insect = .true.
                                 endif ! back if ( sectors(isect)%basis(i) == ibas ) block
                             enddo ! over i={1,sectors(isect)%ndim} loop

                             if ( insect .eqv. .true. ) then
                                 sectors(isect)%hmat(ibas,jbas) = sectors(isect)%hmat(ibas,jbas) + emat(alpha,betta) * (-one)**isgn
                             endif ! back if ( insect .eqv. .true. ) block
                         endif ! back if ( code(alpha) == 0 ) block
                     endif ! back if ( code(betta_ == 1 ) block

                 enddo betloop ! over betta={1,norbs} loop
             enddo alploop ! over alpha={1,norbs} loop
         enddo ! over jbas={1,sectors(isect)%ndim} loop

! four fermions term
!-------------------------------------------------------------------------
         do jbas=1,sectors(isect)%ndim
             alphaloop : do alpha=1,norbs
                 bettaloop : do betta=1,norbs
                     gammaloop : do gamma=1,norbs
                         deltaloop : do delta=1,norbs

                             isgn = 0
                             knew = dec_basis(sectors(isect)%basis(jbas))
                             code(1:norbs) = bin_basis(1:norbs, sectors(isect)%basis(jbas))

! very important if single particle basis has been rotated
                             if ( ( alpha == betta ) .or. ( delta == gamma ) ) CYCLE
                             if ( abs(umat(alpha,betta,delta,gamma)) < epst ) CYCLE

! simulate two annihilation operators
                             if ( ( code(delta) == 1 ) .and. ( code(gamma) == 1 ) ) then
                                 do i=1,gamma-1
                                     if ( code(i) == 1 ) isgn = isgn + 1
                                 enddo ! over i={1,gamma-1} loop
                                 code(gamma) = 0
                                 do i=1,delta-1
                                     if(code(i) == 1) isgn = isgn + 1
                                 enddo ! over i={1,delta-1} loop
                                 code(delta) = 0

! simulate two creation operators
                                 if ( ( code(alpha) == 0 ) .and. ( code(betta) == 0 ) ) then
                                     do i=1,betta-1
                                         if (code(i) == 1) isgn = isgn + 1
                                     enddo ! over i={1,betta-1} loop
                                     code(betta) = 1
                                     do i=1,alpha-1
                                         if (code(i) == 1) isgn = isgn + 1
                                     enddo ! over i={1,alpha-1} loop
                                     code(alpha) = 1

! determine the row number and hamiltonian matrix elememt
                                     knew = knew - 2**(gamma-1) - 2**(delta-1)
                                     knew = knew + 2**(betta-1) + 2**(alpha-1)
                                     ibas = index_basis(knew)
                                     isgn = mod(isgn, 2)
                                     if ( ibas == 0 ) then
                                         call s_print_error('atomic_make_shmat','error while determining row')
                                     endif ! back if ( ibas == 0 ) block

                                     insect = .false.
                                     do i=1,sectors(isect)%ndim
                                         if ( sectors(isect)%basis(i) == ibas ) then
                                             ibas = i; insect = .true.
                                         endif ! back if ( sectors(isect)%basis(i) == ibas ) block
                                     enddo ! over i={1,sectors(isect)%ndim} loop

                                     if ( insect .eqv. .true. ) then
                                         sectors(isect)%hmat(ibas,jbas) = sectors(isect)%hmat(ibas,jbas) + umat(alpha,betta,delta,gamma) * (-one)**isgn
                                     endif ! back if ( insect .eqv. .true. ) block
                                 endif ! back if ( ( code(alpha) == 0 ) .and. ( code(betta) == 0 ) ) block
                             endif ! back if ( ( code(delta) == 1 ) .and. ( code(gamma) == 1 ) ) block

                         enddo deltaloop ! over delta={gamma+1,norbs} loop
                     enddo gammaloop ! over gamma={1,norbs-1} loop
                 enddo bettaloop ! over betta={alpha+1,norbs} loop
             enddo alphaloop ! over alpha={1,norbs-1} loop
         enddo ! over jbas={1,sectors(isect)%ndim} loop

     enddo ! over i={1,nsectors} loop

     return
  end subroutine atomic_make_shmat

!!>>> atomic_diag_shmat: diagonalize the Hamiltonian for each sector
  subroutine atomic_diag_shmat()
     use constants, only : dp

     use m_sector, only : nsectors, sectors

     implicit none

! local variables
! loop index
     integer :: i

! dummy hamiltonian matrix
     real(dp), allocatable :: hmat(:,:)

     do i=1,nsectors
         allocate( hmat(sectors(i)%ndim, sectors(i)%ndim) )
         hmat = real( sectors(i)%hmat )
         call s_eig_sy(sectors(i)%ndim, sectors(i)%ndim, hmat, sectors(i)%eval, sectors(i)%evec)
         deallocate( hmat )
     enddo ! over i={1,nsectors} loop

     return
  end subroutine atomic_diag_shmat

!!>>> atomic_make_sectors: determine all the sectors for good quantum numbers
!!>>> a sector consists of some many particle Fock states labeled by
!!>>> good quantum numbers
  subroutine atomic_make_sectors()
     use constants, only : zero

     use control, only : ictqmc
     use control, only : nband, norbs, ncfgs
     use m_full, only : dim_sub_n, bin_basis
     use m_sector, only : max_dim_sect, ave_dim_sect
     use m_sector, only : nsectors, sectors
     use m_sector, only : alloc_one_sector
     use m_sector, only : alloc_m_sector

     implicit none

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k
     integer :: l

! total electrons
     integer :: my_ntot

! Sz value
     integer :: my_sz

! Jz value
     integer :: my_jz

! PS value
     integer :: my_ps

! a counter
     integer :: counter

! index of Fock basis
     integer :: ibasis

! number of sectors
     integer :: nsect

! which sector point to
     integer :: which_sect

! can point to next sector
     logical :: can

! the sz and jz values for each orbital
     integer :: orb_good_sz(norbs)
     integer :: orb_good_jz(norbs)

! a temp binary form of Fock basis
     integer :: tmp_basis(norbs)

! good quantum number N, Sz, PS for each Fock state
     integer :: fock_good_ntot(ncfgs)
     integer :: fock_good_sz(ncfgs)
     integer :: fock_good_ps(ncfgs)
     integer :: fock_good_jz(ncfgs)

! good quantum number N, Sz, PS for each sector
     integer :: sect_good_ntot(ncfgs)
     integer :: sect_good_sz(ncfgs)
     integer :: sect_good_ps(ncfgs)
     integer :: sect_good_jz(ncfgs)

! dimension of each sector
     integer :: ndims(ncfgs)

! sector basis index
     integer :: sector_basis(ncfgs,ncfgs)

! make good_sz and good_jz
!-------------------------------------------------------------------------
     orb_good_sz = 0
     orb_good_jz = 0
     call atomic_make_gsz(orb_good_sz)

! jz only valid for nband==3, 5, 7
     if ( nband == 3 .or. nband == 5 .or. nband == 7 ) then
         call atomic_make_gjz(orb_good_jz)
     endif ! back if ( nband == 3 .or. nband == 5 .or. nband == 7 ) block

! build good quantum numbers for each Fock state
!-------------------------------------------------------------------------
     counter = 0
     fock_good_ntot = 0
     fock_good_sz = 0
     fock_good_ps = 0
     fock_good_jz = 0
! loop over all number of total electrons
     do i=0,norbs
! loop over each state
         do j=1,dim_sub_n(i)
             counter = counter + 1
! build N
             fock_good_ntot(counter) = i
! build Sz
             my_sz = 0
             do k=1,norbs
                 my_sz = my_sz + orb_good_sz(k) * bin_basis(k, counter)
             enddo ! over k={1,norbs} loop
             fock_good_sz(counter) = my_sz
! build Jz
              my_jz = 0
              do k=1,norbs
                  my_jz = my_jz + orb_good_jz(k) * bin_basis(k, counter)
              enddo ! over k={1,norbs} loop
              fock_good_jz(counter) = my_jz
! build PS number
             do k=1,nband
                 fock_good_ps(counter) = fock_good_ps(counter) + (2**k) * (bin_basis(2*k-1,counter) - bin_basis(2*k,counter))**2
             enddo ! over k={1,nband} loop
         enddo ! over j={1,dim_sub_n(i)} loop
     enddo ! over i={0,norbs} loop

! loop over all the Fock states to determine sectors
!-------------------------------------------------------------------------
     nsect = 0
     ndims = 0
     sector_basis = 0
     do i=1,ncfgs
         my_ntot = fock_good_ntot(i)
         if ( ictqmc == 3 .or. ictqmc == 4 ) then
             my_sz = fock_good_sz(i)
         endif ! back if ( ictqmc == 3 .or. ictqmc == 4 ) block
         if ( ictqmc == 4 ) then
             my_ps = fock_good_ps(i)
         endif ! back if ( ictqmc == 4 ) block
         if ( ictqmc == 5 ) then
             my_jz = fock_good_jz(i)
         endif ! back if ( ictqmc == 5 ) block

         if ( nsect == 0 ) then
             sect_good_ntot(1) = my_ntot
             if ( ictqmc == 3 .or. ictqmc == 4 ) then
                 sect_good_sz(1) = my_sz
             endif ! back if ( ictqmc == 3 .or. ictqmc == 4 ) block
             if (ictqmc == 4) then
                 sect_good_ps(1) = my_ps
             endif ! back if ( ictqmc == 4 ) block
             if ( ictqmc == 5 ) then
                 sect_good_jz(1) = my_jz
             endif ! back if ( ictqmc == 5 ) block

             nsect = nsect + 1
             ndims(1) = ndims(1) + 1
             sector_basis(ndims(1),1) = i
         else
! loop over the exists sectors
             which_sect = -1
             do j=1,nsect
! compare two sectors
                 select case (ictqmc)
                     case (2)
                         if ( sect_good_ntot(j) == my_ntot ) then
                             which_sect = j; EXIT
                         endif ! back if ( sect_good_ntot(j) == my_ntot ) block
                     case (3)
                         if ( sect_good_ntot(j) == my_ntot ) then
                             if ( sect_good_sz(j) == my_sz ) then
                                 which_sect = j; EXIT
                             endif ! back if ( sect_good_sz(j) == my_sz ) block
                         endif ! back if ( sect_good_ntot(j) == my_ntot ) block
                     case (4)
                         if ( sect_good_ntot(j) == my_ntot ) then
                             if ( sect_good_sz(j) == my_sz ) then
                                 if ( sect_good_ps(j) == my_ps) then
                                     which_sect = j; EXIT
                                 endif ! back if ( sect_good_ps(j) == my_ps) block
                             endif ! back if ( sect_good_sz(j) == my_sz ) block
                         endif ! back if ( sect_good_ntot(j) == my_ntot ) block
                     case (5)
                         if ( sect_good_ntot(j) == my_ntot ) then
                             if ( sect_good_jz(j) == my_jz ) then
                                 which_sect = j; EXIT
                             endif ! back if ( sect_good_jz(j) == my_jz ) block
                         endif ! back if ( sect_good_ntot(j) == my_ntot ) block
                 end select
             enddo ! over j={1,nsect} loop
! new sector
             if ( which_sect == -1 ) then
                 nsect = nsect + 1
                 sect_good_ntot(nsect) = my_ntot
                 if ( ictqmc == 3 .or. ictqmc == 4 ) then
                     sect_good_sz(nsect) = my_sz
                 endif ! back if ( ictqmc == 3 .or. ictqmc == 4 ) block
                 if ( ictqmc == 4 ) then
                     sect_good_ps(nsect) = my_ps
                 endif ! back if ( ictqmc == 4 ) block
                 if ( ictqmc == 5 ) then
                     sect_good_jz(nsect) = my_jz
                 endif ! back if ( ictqmc == 5 ) block
                 ndims(nsect) = ndims(nsect) + 1
                 sector_basis(ndims(nsect),nsect) = i
! old sector
             else
                 ndims(which_sect) = ndims(which_sect) + 1
                 sector_basis(ndims(which_sect), which_sect) = i
             endif ! back if ( which_sect == -1 ) block
         endif ! back if ( nsect == 0 ) block
     enddo ! over i={1,ncfgs} loop

! after we know how many sectors and the dimension of each sector are there,
! we can allocate memory for global variables for sectors
!-------------------------------------------------------------------------
     max_dim_sect = 0
     ave_dim_sect = zero
     nsectors = nsect
     call alloc_m_sector()
! now we will build each sector
     counter = 1
     do i=1,nsect
         sectors(i)%ndim = ndims(i)
         sectors(i)%nele = sect_good_ntot(i)
         sectors(i)%nops = norbs
         sectors(i)%istart = counter
         counter = counter + ndims(i)
! allocate memory for each sector
         call alloc_one_sector( sectors(i) )
! set basis for each sector
         do j=1,ndims(i)
             sectors(i)%basis(j) = sector_basis(j,i)
         enddo ! over j={1,ndims(i)} loop
     enddo ! over i={1,nsect} loop

! make next sector index
!-------------------------------------------------------------------------
! loop over all the sectors
     do i=1,nsectors
! loop over all the orbtials
         do j=1,norbs
! loop over creation and annihilation fermion operators
             do k=0,1
                 which_sect = -1
! we should check each state in this sector
                 can = .false.
                 do l=1,sectors(i)%ndim
                     ibasis = sectors(i)%basis(l)
! for creation fermion operator
                     if ( k == 1 .and. bin_basis(j,ibasis) == 0 ) then
                         tmp_basis = bin_basis(:, ibasis)
                         can = .true.
                         EXIT
! for annihilation fermion operator
                     elseif ( k==0 .and. bin_basis(j, ibasis) == 1 ) then
                         tmp_basis = bin_basis(:, ibasis)
                         can = .true.
                         EXIT
                     endif ! back if ( k == 1 .and. bin_basis(j,ibasis) == 0 ) block
                 enddo ! over l={1,sectors(i)%ndim} loop

                 if ( can == .true. ) then
                     select case (ictqmc)
                         case (2)
                             if ( k == 1 ) then
                                 my_ntot = sect_good_ntot(i) + 1
                             else
                                 my_ntot = sect_good_ntot(i) - 1
                             endif ! back if ( k == 1 ) block
! loop over all sectors to see which sector it will point to
                             do l=1,nsectors
                                 if ( sect_good_ntot(l) == my_ntot ) then
                                     which_sect = l; EXIT
                                 endif ! back if ( sect_good_ntot(l) == my_ntot ) block
                             enddo ! over l={1,nsectors} loop

                         case (3)
                             if ( k == 1 ) then
                                 my_ntot = sect_good_ntot(i) + 1
                                 my_sz = sect_good_sz(i) + orb_good_sz(j)
                             else
                                 my_ntot = sect_good_ntot(i) - 1
                                 my_sz = sect_good_sz(i) - orb_good_sz(j)
                             endif ! back if ( k == 1 ) block
! loop over all sectors to see which sector it will point to
                             do l=1,nsectors
                                 if ( sect_good_ntot(l) == my_ntot ) then
                                     if ( sect_good_sz(l) == my_sz ) then
                                         which_sect = l; EXIT
                                     endif ! back if ( sect_good_sz(l) == my_sz ) block
                                 endif ! back if ( sect_good_ntot(l) == my_ntot ) block
                             enddo ! over l={1,nsectors} loop

                         case (4)
                             if ( k == 1 ) then
                                 my_ntot = sect_good_ntot(i) + 1
                                 my_sz = sect_good_sz(i) + orb_good_sz(j)
                                 tmp_basis(j) = 1
                             else
                                 my_ntot = sect_good_ntot(i) - 1
                                 my_sz   = sect_good_sz(i) - orb_good_sz(j)
                                 tmp_basis(j) = 0
                             endif ! back if ( k == 1 ) block
! calculate new PS number
                             my_ps = 0
                             do l=1,nband
                                 my_ps = my_ps + (2**l) * ( tmp_basis(2*l-1) - tmp_basis(2*l) )**2
                             enddo ! over l={1,nband} loop
! loop over all sectors to see which sector it will point to
                             do l=1,nsectors
                                 if ( sect_good_ntot(l) == my_ntot ) then
                                     if ( sect_good_sz(l) == my_sz ) then
                                         if ( sect_good_ps(l) == my_ps) then
                                             which_sect = l; EXIT
                                         endif ! back if ( sect_good_ps(l) == my_ps) block
                                     endif ! back if ( sect_good_sz(l) == my_sz ) block
                                 endif ! back if ( sect_good_ntot(l) == my_ntot ) block
                             enddo ! over l={1,nsectors} loop

                         case (5)
                             if ( k == 1 ) then
                                 my_ntot = sect_good_ntot(i) + 1
                                 my_jz = sect_good_jz(i) + orb_good_jz(j)
                             else
                                 my_ntot = sect_good_ntot(i) - 1
                                 my_jz = sect_good_jz(i) - orb_good_jz(j)
                             endif ! back if ( k == 1 ) block
! loop over all sectors to see which sector it will point to
                             do l=1,nsectors
                                 if ( sect_good_ntot(l) == my_ntot ) then
                                     if ( sect_good_jz(l) == my_jz ) then
                                         which_sect = l; EXIT
                                     endif ! back if ( sect_good_jz(l) == my_jz ) block
                                 endif ! back if ( sect_good_ntot(l) == my_ntot ) block
                             enddo ! over l={1,nsectors} loop
                     end select ! back select case (ictqmc) block
                 endif  ! back if ( can == .true. ) block
                 sectors(i)%next(j,k) = which_sect
             enddo ! over k={0,1} loop
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,nsectors} loop

! calculate the maximum and average dimensions of sectors
!-------------------------------------------------------------------------
     max_dim_sect = 0
     counter = 0
     do i=1,nsectors
         if (sectors(i)%ndim > max_dim_sect) max_dim_sect = sectors(i)%ndim
         counter = counter + sectors(i)%ndim
     enddo ! over i={1,nsectors} loop
     ave_dim_sect = real(counter) / real(nsectors)

! dump sector information for reference
!-------------------------------------------------------------------------
     call atomic_dump_sector(sect_good_ntot, sect_good_sz, sect_good_ps, sect_good_jz)

     return
  end subroutine atomic_make_sectors
