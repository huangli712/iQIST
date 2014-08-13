!!!-------------------------------------------------------------------------
!!! project : jasmine
!!! program : atomic_mksectors_n
!!!         : atomic_mksectors_nsz
!!!         : atomic_mksectors_nszps
!!!         : atomic_mksectors_njz
!!!         : make_good_sz
!!!         : make_good_jz
!!! source  : atomic_mksector.f90
!!! type    : subroutines
!!! author  : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!! purpose : make sectors by using good quantum numbers
!!! input   :
!!! output  :
!!! status  : unstable
!!! comment :
!!!-------------------------------------------------------------------------

!!>>> determine all the sectors for good quantum N case
!! a sector consists of some many particle Fock states labeled by 
!! good quantum number N 
  subroutine atomic_mksectors_n()
     use constants
     use control
     use m_basis_fullspace
     use m_sector
     use m_glob_sectors
  
     implicit none
  
! which sector points to
     integer :: which_sect

! fock basis index
     integer :: ibasis

! a counter 
     integer :: counter

! total electrons
     integer :: myntot

! can point to next sector 
     logical :: can  

! loop index
     integer :: i,j,k,l
  
!----------------------------------------------------------------
! allocate memory for global variables of sectors
     nsectors = norbs + 1
     max_dim_sect = 0
     ave_dim_sect = zero
     call alloc_m_glob_sectors()

! now, build each sector
     counter = 1
     do i=1, nsectors
         sectors(i)%ndim = dim_sub_n(i-1)
         sectors(i)%nelectron = i-1
         sectors(i)%nops = norbs
         sectors(i)%istart = counter 
! allocate memory for each sector
         call alloc_one_sector( sectors(i) )  
! set basis for each sector
         do j=1, sectors(i)%ndim
             sectors(i)%mybasis(j) = counter + j - 1 
         enddo
         counter = counter + dim_sub_n(i-1)
     enddo
!----------------------------------------------------------------
  
!----------------------------------------------------------------
! make next_sector index
! loop for all the sectors
     do i=1, nsectors   
! loop for all the orbitals
         do j=1, norbs              
! loop for creation and annihilation fermion operators
             do k=0,1                   
                 which_sect = -1
! we should lookup each basis in this sector
                 can = .false.
                 do l=1, sectors(i)%ndim
                     ibasis = sectors(i)%mybasis(l)
! for creation fermion operator
                     if (k==1 .and. bin_basis(j,ibasis) == 0) then
                         can = .true.
                         exit
! for annihilation fermion operator
                     elseif (k==0 .and. bin_basis(j, ibasis) == 1) then
                         can = .true. 
                         exit
                     endif 
                 enddo 
  
                 if (can == .true.) then
                     if (k==1) then
                         myntot = sectors(i)%nelectron + 1
                     else
                         myntot = sectors(i)%nelectron - 1
                     endif
! loop over all sectors to see which sector it will point to 
                     if (myntot >=0 .and. myntot <= norbs) then
                         which_sect = myntot + 1
                     endif 
                 endif  ! back to if (can == .true.) block
                 sectors(i)%next_sector(j,k) = which_sect 
             enddo ! over k={0,1} loop
         enddo ! over j={1,norbs} loop
     enddo ! over i={1, nsectors} loop
!----------------------------------------------------------------
  
!----------------------------------------------------------------
! dump sector information for reference
! calculate the maximum and average dimensions of sectors
     max_dim_sect = 0
     counter = 0 
     do i=1, nsectors
         if (sectors(i)%ndim > max_dim_sect) max_dim_sect = sectors(i)%ndim
         counter = counter + sectors(i)%ndim
     enddo
     ave_dim_sect = real(counter) / real(nsectors)
  
     open(mytmp, file='atom.sector.dat')
     write(mytmp, '(a,I10)')    '#number_sectors : ', nsectors
     write(mytmp, '(a,I10)')    '#max_dim_sectors: ', max_dim_sect
     write(mytmp, '(a,F16.8)')  '#ave_dim_sectors: ', ave_dim_sect
     write(mytmp, '(a)') '#      i | electron(i) |     ndim(i) |           j |   fock_basis(j,i) |  '
     do i=1, nsectors
         do j=1, sectors(i)%ndim
             write(mytmp,'(I10,4X,I10,4X,I10,4X,I10,8X, 14I1)') i, sectors(i)%nelectron, &
                                   sectors(i)%ndim, j, bin_basis(:, sectors(i)%mybasis(j)) 
         enddo
     enddo
     close(mytmp)
!----------------------------------------------------------------
  
     return
  end subroutine atomic_mksectors_n

!!>>> determine all the sectors for good quantum N, Sz case
!! a sector consists of some many particle Fock states labeled by 
!! good quantum number N, Sz 
  subroutine atomic_mksectors_nsz()
     use constants
     use control
     use m_basis_fullspace
     use m_sector
     use m_glob_sectors
  
     implicit none
  
! local variables
! the maximum number of sectors
     integer :: max_nsect

! the maximum dimension of each sector
     integer :: max_ndim

! the sz value for each single particle orbital
     integer :: orb_good_sz(norbs)

! good quantum number N, Sz for each Fock state
     integer :: fock_good_ntot(ncfgs)
     integer :: fock_good_sz(ncfgs)

! good quantum number N, Sz for each sector
     integer, allocatable :: sect_good_ntot(:)
     integer, allocatable :: sect_good_sz(:)

! dimension of each sector
     integer, allocatable :: ndims(:)

! sector basis index
     integer, allocatable :: sector_basis(:,:)

! number of sectors
     integer :: nsect

! which sector points to
     integer :: which_sect

! total electrons
     integer :: myntot

! Sz value
     integer :: mysz

! a counter
     integer :: counter

! fock basis
     integer :: ibasis

! loop index
     integer :: i,j,k,l

! can point to next sector 
     logical :: can  
  
     max_nsect = ncfgs
     max_ndim  = ncfgs
! allocate memory
     allocate(sect_good_ntot(max_nsect))
     allocate(sect_good_sz(max_nsect))
     allocate(ndims(max_nsect))
     allocate(sector_basis(max_ndim, max_nsect))
  
!----------------------------------------------------------------
! make good_sz
     call make_good_sz(orb_good_sz)

! build good quantum numbers for each Fock state
     counter = 0
! loop over all number of total electrons
     do i=0, norbs
         do j=1, dim_sub_n(i)
             counter = counter + 1
             mysz = 0
             do k=1, norbs
                 mysz = mysz + orb_good_sz(k) * bin_basis(k, counter) 
             enddo
             fock_good_ntot(counter) = i
             fock_good_sz(counter)   = mysz
         enddo  
     enddo
!----------------------------------------------------------------
  
!----------------------------------------------------------------
! loop over all the Fock states to determine sectors
     nsect = 0
     ndims = 0
     sector_basis = 0
     do i=1, ncfgs    
         myntot = fock_good_ntot(i)
         mysz   = fock_good_sz(i)
         if (nsect==0 ) then
             sect_good_ntot(1) = myntot
             sect_good_sz(1)   = mysz
             nsect = nsect + 1
             ndims(1) = ndims(1) + 1 
             sector_basis(ndims(1),1) = i
         else
! loop over the exists sectors
             which_sect = -1
             do j=1, nsect
! compare two sectors
                 if (sect_good_ntot(j) == myntot .and. sect_good_sz(j) == mysz) then
                     which_sect = j
                     EXIT
                 endif
             enddo 
! new sector
             if( which_sect == -1 ) then
                 nsect = nsect + 1
                 sect_good_ntot(nsect) = myntot
                 sect_good_sz(nsect)   = mysz
                 ndims(nsect) = ndims(nsect) + 1
                 sector_basis(ndims(nsect), nsect) = i
! oldsector
             else
                 ndims(which_sect) = ndims(which_sect) + 1 
                 sector_basis(ndims(which_sect), which_sect) = i
             endif
         endif ! back to if (nsect == 0) then block 
     enddo  ! over i={1,ncfgs} loop 
!----------------------------------------------------------------
  
!----------------------------------------------------------------
! after we know how many sectors and the dimension of each sector,
! we can allocate memory for global variables for sectors
     max_dim_sect = 0
     ave_dim_sect = zero
     nsectors = nsect
     call alloc_m_glob_sectors()

! now we will build each sector
     counter = 1
     do i=1, nsect
         sectors(i)%ndim = ndims(i)
         sectors(i)%nelectron = sect_good_ntot(i)
         sectors(i)%nops = norbs
         sectors(i)%istart = counter 
         counter = counter + ndims(i)
! allocate memory for each sector
         call alloc_one_sector( sectors(i) )  
! set Fock basis for each sector
         do j=1, ndims(i)
             sectors(i)%mybasis(j) = sector_basis(j,i) 
         enddo
     enddo
!----------------------------------------------------------------
  
!----------------------------------------------------------------
! make next_sector index
! loop over all sectors
     do i=1, nsectors
! loop over all the orbitals
         do j=1, norbs 
! loop over creation and annihilation fermion operators
             do k=0,1 
                 which_sect = -1
! we should check each state in this sector
                 can = .false.
                 do l=1, sectors(i)%ndim
                     ibasis = sectors(i)%mybasis(l)
! for creation fermion operator
                     if (k==1 .and. bin_basis(j,ibasis) == 0) then
                         can = .true.
                         exit
! for annihilation fermion operator
                     elseif (k==0 .and. bin_basis(j, ibasis) == 1) then
                         can = .true. 
                         exit
                     endif 
                 enddo 
  
                 if (can == .true.) then
                     if (k==1) then
                         myntot = sect_good_ntot(i) + 1
                         mysz   = sect_good_sz(i) + orb_good_sz(j)
                     else
                         myntot = sect_good_ntot(i) - 1
                         mysz   = sect_good_sz(i) - orb_good_sz(j)
                     endif
! loop over all sectors to see which sector it will points to 
                     do l=1, nsectors
                         if (sect_good_ntot(l) == myntot .and. sect_good_sz(l) == mysz) then
                             which_sect = l
                             exit 
                         endif 
                     enddo 
                 endif  ! back to if (can == .true.) block
                 sectors(i)%next_sector(j,k) = which_sect 
             enddo ! over k={0,1} loop
         enddo ! over j={1,norbs} loop
     enddo ! over i={1, nsectors} loop
!----------------------------------------------------------------
  
!----------------------------------------------------------------
! dump sector information for reference
! calculate the maximum and average dimensions of sectors
     max_dim_sect = 0
     counter = 0 
     do i=1, nsectors
         if (sectors(i)%ndim > max_dim_sect) max_dim_sect = sectors(i)%ndim
         counter = counter + sectors(i)%ndim
     enddo
     ave_dim_sect = real(counter) / real(nsectors)
  
     open(mytmp, file='atom.sector.dat')
     write(mytmp, '(a,I10)')    '#number_sectors : ', nsectors
     write(mytmp, '(a,I10)')    '#max_dim_sectors: ', max_dim_sect
     write(mytmp, '(a,F16.8)')  '#ave_dim_sectors: ', ave_dim_sect
     write(mytmp, '(a)') '#      i | electron(i) |       Sz(i) |     ndim(i) |           j |   fock_basis(j,i) |  '
     do i=1, nsectors
         do j=1, sectors(i)%ndim
             write(mytmp,'(I10,4X,I10,4X,I10,4X,I10,4X,I10,8X,14I1)') i, sect_good_ntot(i),&
                  sect_good_sz(i), sectors(i)%ndim, j, bin_basis(:, sectors(i)%mybasis(j)) 
         enddo
     enddo
     close(mytmp)
!----------------------------------------------------------------
  
! free memeory
     if (allocated(sect_good_ntot)) deallocate(sect_good_ntot)
     if (allocated(sect_good_sz))   deallocate(sect_good_sz)
     if (allocated(ndims))          deallocate(ndims)
     if (allocated(sector_basis))   deallocate(sector_basis)
  
     return
  end subroutine atomic_mksectors_nsz

!!>>> determine all the sectors for good quantum N, Sz, PS case
!! a sector consists of some many particle Fock states labeled by 
!! good quantum number N, Sz, PS
  subroutine atomic_mksectors_nszps()
     use constants
     use control
     use m_basis_fullspace
     use m_sector
     use m_glob_sectors
  
     implicit none
  
! local variables
! the maximum number of sectors
     integer :: max_nsect

! the maximum dimension of each sector
     integer :: max_ndim

! the sz value for each orbital
     integer :: orb_good_sz(norbs)

! good quantum number N, Sz, PS for each Fock state
     integer :: fock_good_ntot(ncfgs)
     integer :: fock_good_sz(ncfgs)
     integer :: fock_good_ps(ncfgs)

! good quantum number N, Sz, PS for each sector
     integer, allocatable :: sect_good_ntot(:)
     integer, allocatable :: sect_good_sz(:)
     integer, allocatable :: sect_good_ps(:)

! dimension of each sector
     integer, allocatable :: ndims(:)

! sector basis index
     integer, allocatable :: sector_basis(:,:)

! number of sectors
     integer :: nsect

! which sector point to
     integer :: which_sect

! a temp binary form of Fock basis
     integer :: tmp_basis(norbs)

! total electrons
     integer :: myntot

! Sz value
     integer :: mysz

! PS value
     integer :: myps

! a counter
     integer :: counter

! index of Fock basis
     integer :: ibasis

! loop index
     integer :: i,j,k,l

! can point to next sector
     logical :: can  
  
     max_nsect = ncfgs
     max_ndim = ncfgs
! allocate memory
     allocate(sect_good_ntot(max_nsect))
     allocate(sect_good_sz(max_nsect))
     allocate(sect_good_ps(max_nsect))
     allocate(ndims(max_nsect))
     allocate(sector_basis(max_ndim, max_nsect))
  
!----------------------------------------------------------------
! make good_sz
     call make_good_sz(orb_good_sz)
! build good quantum numbers for each Fock state
     counter = 0
     fock_good_ps = 0
! loop over all number of total electrons
     do i=0, norbs
! loop over each state 
         do j=1, dim_sub_n(i)
             counter = counter + 1
             mysz = 0
             do k=1, norbs
                 mysz = mysz + orb_good_sz(k) * bin_basis(k, counter) 
             enddo
             fock_good_ntot(counter) = i
             fock_good_sz(counter) = mysz
! build PS number
             do k=1, nband
                 fock_good_ps(counter) = fock_good_ps(counter) + &
                 2**k * (bin_basis(2*k-1,counter) - bin_basis(2*k,counter))**2
             enddo
         enddo  
     enddo
!----------------------------------------------------------------
  
!----------------------------------------------------------------
! loop over all the Fock states to determine sectors
     nsect = 0
     ndims = 0
     sector_basis = 0
     do i=1, ncfgs    
         myntot = fock_good_ntot(i)
         mysz   = fock_good_sz(i)
         myps   = fock_good_ps(i)
         if (nsect==0) then
             sect_good_ntot(1) = myntot
             sect_good_sz(1)   = mysz
             sect_good_ps(1)   = myps
             nsect = nsect + 1
             ndims(1) = ndims(1) + 1 
             sector_basis(ndims(1),1) = i
         else
! loop over the exists sectors
             which_sect = -1
             do j=1, nsect
! compare two sectors
                 if ( sect_good_ntot(j) == myntot .and. sect_good_sz(j) == mysz &
                     .and. sect_good_ps(j) == myps) then
                     which_sect = j
                     EXIT
                 endif
             enddo 
! new sector
             if( which_sect == -1 ) then
                 nsect = nsect + 1
                 sect_good_ntot(nsect) = myntot
                 sect_good_sz(nsect)   = mysz
                 sect_good_ps(nsect)   = myps
                 ndims(nsect) = ndims(nsect) + 1
                 sector_basis(ndims(nsect), nsect) = i
! old sector
             else
                 ndims(which_sect) = ndims(which_sect) + 1 
                 sector_basis(ndims(which_sect), which_sect) = i
             endif
         endif ! back to if (nsect == 0) then block 
     enddo ! over i={1,ncfgs} loop
!----------------------------------------------------------------
  
!----------------------------------------------------------------
! after we know how many sectors and the dimension of each sector,
! we can allocate memory for global variables for sectors
     max_dim_sect = 0
     ave_dim_sect = zero
     nsectors = nsect
     call alloc_m_glob_sectors()
! now we will build each sector
     counter = 1
     do i=1, nsect
         sectors(i)%ndim = ndims(i)
         sectors(i)%nelectron = sect_good_ntot(i)
         sectors(i)%nops = norbs
         sectors(i)%istart = counter 
         counter = counter + ndims(i)
! allocate memory for each sector 
         call alloc_one_sector( sectors(i) )  
! set basis for each sector
         do j=1, ndims(i)
             sectors(i)%mybasis(j) = sector_basis(j,i) 
         enddo
     enddo
!----------------------------------------------------------------
  
!----------------------------------------------------------------
! make next_sector index
! loop over all the sectors
     do i=1, nsectors
! loop over all the orbtials
         do j=1, norbs 
! loop over creation and annihilation fermion operators
             do k=0,1 
                 which_sect = -1
! we should check each state in this sector
                 can = .false.
                 do l=1, sectors(i)%ndim
                     ibasis = sectors(i)%mybasis(l)
! for creation fermion operator
                     if (k==1 .and. bin_basis(j,ibasis) == 0) then
                         tmp_basis = bin_basis(:, ibasis)
                         can = .true.
                         exit
! for annihilation fermion operator
                     elseif (k==0 .and. bin_basis(j, ibasis) == 1) then
                         tmp_basis = bin_basis(:, ibasis)
                         can = .true. 
                         exit
                     endif 
                 enddo 
  
                 if (can == .true.) then
                     if (k==1) then
                         myntot = sect_good_ntot(i) + 1
                         mysz   = sect_good_sz(i) + orb_good_sz(j)
                         tmp_basis(j) = 1
                     else
                         myntot = sect_good_ntot(i) - 1
                         mysz   = sect_good_sz(i) - orb_good_sz(j)
                         tmp_basis(j) = 0
                     endif
! calculate new PS number
                     myps = 0
                     do l=1, nband
                         myps = myps + 2**l * ( tmp_basis(2*l-1) - tmp_basis(2*l) )**2
                     enddo
! loop over all sectors to see which sector it will point to 
                     do l=1, nsectors
                         if (sect_good_ntot(l) == myntot .and. sect_good_sz(l) == mysz &
                             .and. sect_good_ps(l) == myps) then
                             which_sect = l
                             exit 
                         endif 
                     enddo 
                 endif  ! back to if (can == .true.) block
                 sectors(i)%next_sector(j,k) = which_sect 
             enddo ! over k={0,1} loop
         enddo ! over j={1,norbs} loop
     enddo ! over i={1, nsectors} loop
!----------------------------------------------------------------
  
!----------------------------------------------------------------
! dump sector information for reference
! calculate the maximum and average dimensions of sectors
     max_dim_sect = 0
     counter = 0 
     do i=1, nsectors
         if (sectors(i)%ndim > max_dim_sect) max_dim_sect = sectors(i)%ndim
         counter = counter + sectors(i)%ndim
     enddo
     ave_dim_sect = real(counter) / real(nsectors)
  
     open(mytmp, file='atom.sector.dat')
     write(mytmp, '(a,I10)')    '#number_sectors : ', nsectors
     write(mytmp, '(a,I10)')    '#max_dim_sectors: ', max_dim_sect
     write(mytmp, '(a,F16.8)')  '#ave_dim_sectors: ', ave_dim_sect
     write(mytmp, '(a)') '#      i | electron(i) |       Sz(i) |       PS(i) |     nd&
                         im(i) |           j |    fock_basis(j,i) |  '
     do i=1, nsectors
         do j=1, sectors(i)%ndim
             write(mytmp,'(I10,4X,I10,4X,I10,4X,I10,4X,I10,4X,I10,8X,14I1)') i, sect_good_ntot(i), &
           sect_good_sz(i), sect_good_ps(i), sectors(i)%ndim, j, bin_basis(:, sectors(i)%mybasis(j)) 
         enddo
     enddo
     close(mytmp)
!----------------------------------------------------------------
  
! free memeory
     if (allocated(sect_good_ntot)) deallocate(sect_good_ntot)
     if (allocated(sect_good_sz))   deallocate(sect_good_sz)
     if (allocated(ndims))          deallocate(ndims)
     if (allocated(sector_basis))   deallocate(sector_basis)
  
     return
  end subroutine atomic_mksectors_nszps
  
!!>>> determine all the sectors for good quantum N, Jz case
!! a sector consists of some many particle Fock states labeled by 
!! good quantum number N, Jz 
  subroutine atomic_mksectors_njz()
      use constants
      use control
      use m_basis_fullspace
      use m_sector
      use m_glob_sectors
  
      implicit none
  
! local variables
! the maximum number of sectors
      integer :: max_nsect

! the maximum dimension of each sector
      integer :: max_ndim

! the jz value for each |j2,jz> single particle orbital
      integer :: orb_good_jz(norbs)

! good quantum number N, Jz for each Fock state
      integer :: fock_good_ntot(ncfgs)
      integer :: fock_good_jz(ncfgs)

! good quantum number N, Jz for each sector
      integer, allocatable :: sect_good_ntot(:)
      integer, allocatable :: sect_good_jz(:)

! dimension of each sector
      integer, allocatable :: ndims(:)

! sector basis index
      integer, allocatable :: sector_basis(:,:)

! number of sectors
      integer :: nsect

! which sector point to
      integer :: which_sect

! total electrons
      integer :: myntot

! Jz value
      integer :: myjz

! a counter
      integer :: counter

! index of Fock basis
      integer :: ibasis

! loop index
      integer :: i,j,k,l

! can point to next sector
      logical :: can  
  
  
      max_nsect = ncfgs
      max_ndim = ncfgs
! allocate memory
      allocate(sect_good_ntot(max_nsect))
      allocate(sect_good_jz(max_nsect))
      allocate(ndims(max_nsect))
      allocate(sector_basis(max_ndim, max_nsect))
  
!----------------------------------------------------------------
! make good_jz
      call make_good_jz(orb_good_jz)
! build good quantum numbers for each Fock state
      counter = 0
! loop over all the number of total electrons
      do i=0, norbs
          do j=1, dim_sub_n(i)
              counter = counter + 1
              myjz = 0
              do k=1, norbs
                  myjz = myjz + orb_good_jz(k) * bin_basis(k, counter) 
              enddo
              fock_good_ntot(counter) = i
              fock_good_jz(counter) = myjz
          enddo  
      enddo
!----------------------------------------------------------------

!----------------------------------------------------------------
! loop over all the Fock states to determine sectors
      nsect = 0
      ndims = 0
      sector_basis = 0
      do i=1, ncfgs    
          myntot = fock_good_ntot(i)
          myjz   = fock_good_jz(i)
          if (nsect==0) then
              sect_good_ntot(1) = myntot
              sect_good_jz(1)   = myjz
              nsect = nsect + 1
              ndims(1) = ndims(1) + 1 
              sector_basis(ndims(1),1) = i
          else
! loop over the exists sectors
              which_sect = -1
              do j=1, nsect
! compare two sectors
                  if ( sect_good_ntot(j) == myntot .and. sect_good_jz(j) == myjz) then
                      which_sect = j
                      EXIT
                  endif
              enddo 
! new sector
              if( which_sect == -1 ) then
                  nsect = nsect + 1
                  sect_good_ntot(nsect) = myntot
                  sect_good_jz(nsect)   = myjz
                  ndims(nsect) = ndims(nsect) + 1
                  sector_basis(ndims(nsect), nsect) = i
! old sector
              else
                  ndims(which_sect) = ndims(which_sect) + 1 
                  sector_basis(ndims(which_sect), which_sect) = i
              endif
          endif ! back to if (nsect == 0) then block 
      enddo ! over i={1,ncfgs} loop
!----------------------------------------------------------------

!----------------------------------------------------------------
! after we know how many sectors and the dimension of each sector,
! we can allocate memory for global variables for sectors
      max_dim_sect = 0
      ave_dim_sect = zero
      nsectors = nsect
      call alloc_m_glob_sectors()
! now we will build each sector
      counter = 1
      do i=1, nsect
          sectors(i)%ndim = ndims(i)
          sectors(i)%nelectron = sect_good_ntot(i)
          sectors(i)%nops = norbs
          sectors(i)%istart = counter 
          counter = counter + ndims(i)
! allocate memory for each sector
          call alloc_one_sector( sectors(i) )  
! set basis for each sector
          do j=1, ndims(i)
              sectors(i)%mybasis(j) = sector_basis(j,i) 
          enddo
      enddo
!----------------------------------------------------------------

!----------------------------------------------------------------
! make next_sector index
! loop over all the sectors
      do i=1, nsectors
! loop over all the orbitals
          do j=1, norbs 
! loop over creation and annihilation fermion operators
              do k=0,1 
                  which_sect = -1
! we should check each state in this sector
                  can = .false.
                  do l=1, sectors(i)%ndim
                      ibasis = sectors(i)%mybasis(l)
! for creation fermion operator
                      if (k==1 .and. bin_basis(j,ibasis) == 0 ) then
                          can = .true.
                          exit
! for annihilation fermion operator
                      elseif (k==0 .and. bin_basis(j, ibasis) == 1) then
                          can = .true. 
                          exit
                      endif 
                  enddo 
  
                  if (can == .true.) then
                      if (k==1) then
                          myntot = sect_good_ntot(i) + 1
                          myjz   = sect_good_jz(i) + orb_good_jz(j)
                      else
                          myntot = sect_good_ntot(i) - 1
                          myjz   = sect_good_jz(i) - orb_good_jz(j)
                      endif
! loop over all sectors to see which sector it will point to 
                      do l=1, nsectors
                          if (sect_good_ntot(l) == myntot .and. sect_good_jz(l) == myjz) then
                              which_sect = l
                              exit 
                          endif 
                      enddo 
                  endif  ! back to if (can == .true.) block
                  sectors(i)%next_sector(j,k) = which_sect 
              enddo ! over k={0,1} loop
          enddo ! over j={1,norbs} loop
      enddo ! over i={1, nsectors} loop
!----------------------------------------------------------------

!----------------------------------------------------------------
! dump sector information for reference
! calculate the maximum and average dimensions of sectors
      max_dim_sect = 0
      counter = 0 
      do i=1, nsectors
          if (sectors(i)%ndim > max_dim_sect) max_dim_sect = sectors(i)%ndim
          counter = counter + sectors(i)%ndim
      enddo
      ave_dim_sect = counter / real(nsectors)
  
      open(mytmp, file='atom.sector.dat')
      write(mytmp, '(a,I10)')    '#number_sectors : ', nsectors
      write(mytmp, '(a,I10)')    '#max_dim_sectors: ', max_dim_sect
      write(mytmp, '(a,F16.8)')  '#ave_dim_sectors: ', ave_dim_sect
      write(mytmp, '(a)') '#      i | electron(i) |       Jz(i) |     ndim(i) |           j |   fock_basis(j,i) |  '
      do i=1, nsectors
          do j=1, sectors(i)%ndim
              write(mytmp,'(I10,4X,I10,4X,I10,4X,I10,4X,I10,8X,14I1)') i, sect_good_ntot(i),&
                   sect_good_jz(i), sectors(i)%ndim, j, bin_basis(:, sectors(i)%mybasis(j)) 
          enddo
      enddo
      close(mytmp)
!----------------------------------------------------------------

! free memeory
      if (allocated(sect_good_ntot)) deallocate(sect_good_ntot)
      if (allocated(sect_good_jz))   deallocate(sect_good_jz)
      if (allocated(ndims))          deallocate(ndims)
      if (allocated(sector_basis))   deallocate(sector_basis)
  
      return
  end subroutine atomic_mksectors_njz
  
!!>>> make sz for each orbital
  subroutine make_good_sz(good_sz)
     use control, only: norbs
  
     implicit none
  
! external variables
     integer, intent(out) :: good_sz(norbs)

! local variables
     integer :: i
  
     do i=1, norbs
         if (mod(i,2) /= 0 ) then
             good_sz(i) = 1
         else
             good_sz(i) = -1
         endif
     enddo
  
     return
  end subroutine make_good_sz
  
!>>> make jz for each orbital
  subroutine make_good_jz(good_jz)
     use control, only: nband, norbs
  
     implicit none
  
! external variables
     integer, intent(out) :: good_jz(norbs)
  
     if (nband == 3) then
! j=1/2
         good_jz(1) = -1
         good_jz(2) =  1
! j=3/2
         good_jz(3) = -3
         good_jz(4) = -1
         good_jz(5) =  1
         good_jz(6) =  3
     elseif (nband == 5) then
! j=3/2
         good_jz(1) = -3
         good_jz(2) = -1
         good_jz(3) =  1
         good_jz(4) =  3
! j=5/2
         good_jz(5) = -5
         good_jz(6) = -3
         good_jz(7) = -1
         good_jz(8) =  1
         good_jz(9) =  3
         good_jz(10)=  5
     elseif (nband == 7) then
! j=5/2
         good_jz(1) = -5
         good_jz(2) = -3
         good_jz(3) = -1
         good_jz(4) =  1
         good_jz(5) =  3
         good_jz(6) =  5
! j=7/2
         good_jz(7) = -7
         good_jz(8) = -5
         good_jz(9) = -3
         good_jz(10)= -1
         good_jz(11)=  1
         good_jz(12)=  3
         good_jz(13)=  5
         good_jz(14)=  7
     else
         call s_print_error('make_good_jz', 'not implemented for this norbs value !')
     endif
  
     return
  end subroutine make_good_jz
