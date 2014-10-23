!!!-----------------------------------------------------------------------
!!! project : jasmine
!!! program : atomic_dump_fock
!!!           atomic_dump_natural
!!!           atomic_dump_emat
!!!           atomic_dump_feigval
!!!           atomic_dump_feigvec
!!!           atomic_dump_fcix
!!!           atomic_dump_seigval
!!!           atomic_dump_seigvec
!!!           atomic_dump_scix
!!!           atomic_dump_sector
!!! source  : atomic_dump.f90
!!! type    : subroutines
!!! author  : yilin wang (email: qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang
!!!           08/22/2014 by yilin wang
!!!           10/20/2014 by li huang
!!! purpose : write output files
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> atomic_dump_fock: write Fock basis to file 'atom.basis.dat'
  subroutine atomic_dump_fock()
     use constants, only : mytmp

     use control, only : ncfgs
     use m_full, only : dec_basis, index_basis, bin_basis

     implicit none

! local variables
     integer :: i

! open file 'atom.basis.dat' to write
     open(mytmp, file='atom.basis.dat', form='formatted', status='unknown')

! write the header
     write(mytmp,'(a)') '#      i |  decimal(i) |    index(i) |      binary(i) |'

! write the data
     do i=1,ncfgs
         write(mytmp,"(I10,4X)",advance="no") i
         write(mytmp,"(I10,4X)",advance="no") dec_basis(i)
         write(mytmp,"(I10,8X)",advance="no") index_basis(dec_basis(i))
         write(mytmp,"(14I1)") bin_basis(:,i)
     enddo ! over i={1,ncfgs} loop

! close data file
     close(mytmp)

     return
  end subroutine atomic_dump_fock

!!>>> atomic_dump_natural: write the transformation matrix from the
!!>>> original basis to natural basis
  subroutine atomic_dump_natural(info)
     use constants, only : dp, mytmp

     use control, only : norbs
     use m_spmat, only : tmat

! external arguments
! message to be displayed in the file header
     character(len=*), intent(in) :: info

! local variables
! loop index
     integer :: i
     integer :: j

! open file 'atom.natural.dat' to write
     open(mytmp, file='atom.natural.dat', form='formatted', status='unknown')

! write the header
     write(mytmp,'(a)') info
     write(mytmp,'(a)') '#      i |       j |    umat_real(i,j) |    umat_imag(i,j) |'

! write the data
     do i=1,norbs
         do j=1,norbs
             write(mytmp,'(2I10,2F20.10)') j, i, tmat(j,i)
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine atomic_dump_natural

!!>>> atomic_dump_emat: write on-site impurity energy on natural basis
  subroutine atomic_dump_emat()
     use constants, only : mytmp

     use control, only : isoc
     use control, only : nband, norbs
     use m_spmat, only : emat

     implicit none

! local variables
! loop index
     integer :: i

! auxiliary integer variable used to convert the spin sequence
     integer :: s_order

! open file 'atom.eimp.dat' to write
     open(mytmp, file='atom.eimp.dat', form='formatted', status='unknown')

! write the data
     do i=1,norbs
         if ( isoc == 0 ) then
             if ( i <= nband ) then
                 s_order = 2*i-1
             else
                 s_order = 2*(i-nband)
             endif ! back if ( i <= nband ) block
         else
             s_order = i
         endif ! back if ( isoc == 0 ) block
         write(mytmp,'(I10,F20.10)') i, real(emat(s_order,s_order))
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine atomic_dump_emat

!!>>> atomic_dump_feigval: write eigenvalue for full Hilbert space to
!!>>> file 'atom.eigval.dat'
  subroutine atomic_dump_feigval()
     use constants, only : mytmp

     use control, only : ncfgs
     use m_full, only : eval, occu

     implicit none

! local variables
! loop index
     integer :: i

! open file 'atom.eigval.dat' to write
     open(mytmp, file='atom.eigval.dat', form='formatted', status='unknown')

! write the header
     write(mytmp,'(a)') '#       i |     eigenvalue(i) |      occupancy(i) |'

! write the data
     do i=1,ncfgs
         write(mytmp,"(I10,2F20.10)") i, eval(i), occu(i,i)
     enddo ! over i={1,ncfgs} loop

! close data file
     close(mytmp)

     return
  end subroutine atomic_dump_feigval

!!>>> atomic_dump_feigvec: write eigenvector for full Hilbert space to
!!>>> file 'atom.eigvec.dat'
  subroutine atomic_dump_feigvec()
     use constants, only : eps6, mytmp

     use control, only : ncfgs
     use m_full, only : bin_basis
     use m_full, only : evec

     implicit none

! local variables
! loop index
     integer :: i
     integer :: j

! open file 'atom.eigvec.dat' to write
     open(mytmp, file='atom.eigvec.dat', form='formatted', status='unknown')

! write the header
     write(mytmp,'(a)') '#      i |       j |  eigenvector(i,j) |       fock_basis(i) |'

! write the data
     do i=1,ncfgs
         do j=1,ncfgs
             if ( abs(evec(j,i)) < eps6 ) CYCLE
             write(mytmp,"(2I10,F20.10,8X,14I1)") j, i, evec(j,i), bin_basis(:,j)
         enddo ! over j={1,ncfgs} loop
     enddo ! over i={1,ncfgs} loop

! close data file
     close(mytmp)

     return
  end subroutine atomic_dump_feigvec

!!>>> atomic_dump_fcix: write atom.cix file which are only compatible with
!!>>> the BEGONIA and LAVENDER components
  subroutine atomic_dump_fcix()
     use constants, only : zero, mytmp

     use control, only : isoc
     use control, only : nband, norbs, ncfgs
     use m_full, only : eval, occu, fmat

     implicit none

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k

! auxiliary integer variable used to convert the spin sequence
     integer :: s_order

! open file 'atom.cix' to write
     open(mytmp, file='atom.cix', form='formatted', status='unknown')

! write eigenvalues
     write(mytmp,'(a)') '# eigenvalues: index | energy | occupy | spin'
     do i=1,ncfgs
         write(mytmp,'(I10,3F20.10)') i, eval(i), occu(i,i), zero
     enddo ! over i={1,ncfgs} loop

! write F-matrix
! for non-soc case, the spin order of CTQMC is like up, up, up, dn, dn, dn
! but the spin order of this program is up, dn, up, dn, up, dn. So we have
! to adjust it here. However for soc case, it doesn't matter
     write(mytmp,'(a)') '# f matrix element: alpha | beta | orbital | fmat'
     do i=1,norbs
         if ( isoc == 0 ) then
             if ( i <= nband ) then
                 s_order = 2*i-1
             else
                 s_order = 2*(i-nband)
             endif ! back if ( i <= nband ) block
         else
             s_order = i
         endif ! back if ( isoc == 0 ) block

         do j=1,ncfgs
             do k=1,ncfgs
                 write(mytmp,'(3I10,F20.10)') k, j, i, fmat(k,j,s_order)
             enddo ! back k={1,ncfgs} loop
         enddo ! over j={1,ncfgs} loop
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine atomic_dump_fcix

!!>>> atomic_dump_seigval: write eigenvalue of sectors
!!>>> to file 'atom.eigval.dat'
  subroutine atomic_dump_seigval()
     use constants, only : mytmp

     use m_sector, only : nsectors
     use m_sector, only : sectors

     implicit none

! local variables
! loop index
     integer :: i
     integer :: j

! the counter
     integer :: counter

! open file 'atom.eigval.dat' to write
     open(mytmp, file='atom.eigval.dat', form='formatted', status='unknown')

! write the header
     write(mytmp, '(a)') "#      i |     sect(i) | electron(i) |           j |   eigenvalue(j,i) |"

! write the data
     counter = 0
     do i=1,nsectors
         do j=1,sectors(i)%ndim
             counter = counter + 1
             write(mytmp,"(I10,4X)",advance="no") counter
             write(mytmp,"(I10,4X)",advance="no") i
             write(mytmp,"(I10,4X)",advance="no") sectors(i)%nele
             write(mytmp,"(I10)",advance="no") j
             write(mytmp,"(F20.10)") sectors(i)%eigval(j)
         enddo ! over i={1,nsectors} loop
     enddo ! over j={1,sectors(i)%ndim} loop

! close data file
     close(mytmp)

     return
  end subroutine atomic_dump_seigval

!!>>> atomic_dump_seigvec: write eigenvector of sectors
!!>>> to file 'atom.eigval.dat'
  subroutine atomic_dump_seigvec()
     use constants, only : eps6, mytmp

     use m_full, only : bin_basis
     use m_sector, only : nsectors
     use m_sector, only : sectors

     implicit none

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k

! the counter
     integer :: counter

! open file 'atom.eigvec.dat' to write
     open(mytmp, file="atom.eigvec.dat", form='formatted', status='unknown')

! write the header
     write(mytmp, '(a)') '#      i |       j |       k |       eigvec(j,k) |       fock_basis(j) |'

! write the data
     counter = 0
     do i=1,nsectors
         do j=1,sectors(i)%ndim
             do k=1,sectors(i)%ndim
                 if ( abs( sectors(i)%eigvec(k,j) ) < eps6 ) CYCLE
                 write(mytmp,"(I10)",advance="no") i
                 write(mytmp,"(I10)",advance="no") k + counter
                 write(mytmp,"(I10)",advance="no") j + counter
                 write(mytmp,"(F20.10)",advance="no") sectors(i)%eigvec(k,j)
                 write(mytmp,"(8X,14I1)") bin_basis(:,sectors(i)%basis(k))
             enddo ! over k={1,sectors(i)%ndim} loop
         enddo ! over j={1,sectors(i)%ndim} loop
         counter = counter + sectors(i)%ndim
     enddo ! over i={1,nsectors} loop

! close data file
     close(mytmp)

     return
  end subroutine atomic_dump_seigvec

!!>>> atomic_dump_scix: write atom.cix file which are only compatible with
!!>>> the PANSY and MANJUSHAKA components
  subroutine atomic_dump_scix()
     use constants, only : mytmp

     use control, only : icu, isoc
     use control, only : Uc, Uv, Jz, Js, Jp
     use control, only : Ud, Jh
     use m_sector, only : nsectors, max_dim_sect, ave_dim_sect
     use m_sector, only : sectors

     implicit none

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k

! auxiliary integer variable used to convert the spin sequence
     integer :: s_order

! open 'atom.cix' to write
     open(mytmp, file='atom.cix', form='formatted', status='unknown')

! write header
     write(mytmp,"(a)") "# This file is generated by JASMINE, can be ONLY used by PANSY and MANJUSHAKA CTQMC codes."
     write(mytmp,"(a)") "# WARNING: please DON'T change this file manually !"
     if ( icu == 1 ) then
         write(mytmp,"(a)") "# Kanamori-parameterized Coulomb interaction: "
         write(mytmp,"(5(a,F9.5,2X))") "# Uc:", Uc, "Uv:", Uv, "Jz:", Jz, "Js:", Js, "Jp:", Jp
     else
         write(mytmp,"(a)") "# Slater-parameterized Coulomb interaction: "
         write(mytmp,"(2(a,F9.5,2X))") "# Ud:", Ud, "Jh:", Jh
     endif ! back if ( icu == 1 ) block
     write(mytmp,"(a)") "#================================================================================"

! write number of sectors
     write(mytmp,"(a)") "# NUMBER OF SECTORS | MAXIMUM DIMENSION OF SECTORS | AVERAGE DIMENSION OF SECTORS"
     write(mytmp,"(I10,15X,I10,15X,F20.10)") nsectors, max_dim_sect, ave_dim_sect

! write dimension, total electrons, next_sector, eigenvalue of each sector
     do i=1,nsectors
         write(mytmp,"(a)") "# SECT_INFO: INDEX  |  NDIM  |  NELEC  |   NOPS  |  ISTART"
         write(mytmp,"(6X,5I10)") i, sectors(i)%ndim, sectors(i)%nele, sectors(i)%nops, sectors(i)%istart

! write next sector
         write(mytmp, "(4X,a)") "# NEXT_SECTOR    F     F^{\dagger}"
         do j=1,sectors(i)%nops
! adjust the orbital order for CTQMC, up, up, up, dn, dn, dn
             if ( isoc == 0 ) then
                 if (j <= sectors(i)%nops / 2) then
                     s_order = 2*j-1
                 else
                     s_order = 2*(j - sectors(i)%nops / 2)
                 endif ! back if (j <= sectors(i)%nops / 2) block
             else
                 s_order = j
             endif ! back if ( isoc == 0 ) block
             write(mytmp, "(2X, 3I10)") j, sectors(i)%next(s_order,0), sectors(i)%next(s_order,1)
         enddo ! over j={1,sectors(i)%nops} loop

! write eigeanvalue
         write(mytmp, "(4X,a)") "# EIGENVALUES"
         do j=1,sectors(i)%ndim
             write(mytmp, "(2X,I10, F20.10)") j, sectors(i)%eigval(j)
         enddo ! over j={1,sectors(i)%ndim} loop
     enddo ! over i={1,nsectors} loop

! close data file
     close(mytmp)

! open 'atom.fmat' to write
     open(mytmp, file='atom.fmat', form='unformatted', status='unknown')

! write the data
     do i=1,nsectors
         do j=1,sectors(i)%nops
             if ( isoc == 0 ) then
                 if (j <= sectors(i)%nops / 2) then
                     s_order = 2*j-1
                 else
                     s_order = 2*(j-sectors(i)%nops / 2)
                 endif ! back if (j <= sectors(i)%nops / 2) block
             else
                 s_order = j
             endif ! back if ( isoc == 0 ) block
             do k=0,1
                 if ( sectors(i)%next(s_order,k) == -1 ) CYCLE
                 write(mytmp) sectors(i)%fmat(s_order,k)%item(:,:)
             enddo  ! over k={0,1} loop
         enddo ! over j={1,sectors(i)%nops} loop
     enddo  ! over i={1,nsect} loop

! close data file
     close(mytmp)

     return
  end subroutine atomic_dump_scix

!!>>> atomic_dump_sector: write out the configuration for each sector
!!>>> to file atom.sector.dat
  subroutine atomic_dump_sector(sect_good_ntot, sect_good_sz, sect_good_ps, sect_good_jz)
     use constants, only : mytmp

     use control, only : ictqmc
     use control, only : ncfgs
     use m_full, only : bin_basis
     use m_sector, only : nsectors, max_dim_sect, ave_dim_sect
     use m_sector, only : sectors

     implicit none

! external arguments
! good quantum number: N
     integer, intent(in) :: sect_good_ntot(ncfgs)

! good quantum number: Sz
     integer, intent(in) :: sect_good_sz(ncfgs)

! good quantum number: PS
     integer, intent(in) :: sect_good_ps(ncfgs)

! good quantum number: Jz
     integer, intent(in) :: sect_good_jz(ncfgs)

! local variables
! loop index
     integer :: i
     integer :: j

! open 'atom.sector.dat' to write
     open(mytmp, file='atom.sector.dat', form='formatted', status='unknown')

! write header
     write(mytmp, '(a,I10)')    '#number_sectors : ', nsectors
     write(mytmp, '(a,I10)')    '#max_dim_sectors: ', max_dim_sect
     write(mytmp, '(a,F16.8)')  '#ave_dim_sectors: ', ave_dim_sect

! write the data
     select case (ictqmc)
         case (1)
             call s_print_error('atomic_dump_sector','this function is not implemented')

         case (2)
             write(mytmp,'(a)',advance='no') '#      i |'
             write(mytmp,'(a)',advance='no') ' electron(i) '
             write(mytmp,'(a)') '|     ndim(i) |           j |   fock_basis(j,i) |'
             do i=1,nsectors
                 do j=1,sectors(i)%ndim
                     write(mytmp,'(I10,4X)',advance='no') i
                     write(mytmp,'(I10,4X)',advance='no') sectors(i)%nele
                     write(mytmp,'(I10,4X)',advance='no') sectors(i)%ndim
                     write(mytmp,'(I10,8X)',advance='no') j
                     write(mytmp,'(14I1)') bin_basis(:, sectors(i)%basis(j))
                 enddo ! over j={1,sectors(i)%ndim} loop
             enddo ! over i={1,nsectors} loop

         case(3)
             write(mytmp,'(a)',advance='no') '#      i |'
             write(mytmp,'(a)',advance='no') ' electron(i) |       Sz(i)' 
             write(mytmp,'(a)') '|     ndim(i) |           j |   fock_basis(j,i) |'
             do i=1,nsectors
                 do j=1,sectors(i)%ndim
                     write(mytmp,'(I10,4X)',advance='no') i
                     write(mytmp,'(I10,4X)',advance='no') sect_good_ntot(i)
                     write(mytmp,'(I10,4X)',advance='no') sect_good_sz(i)
                     write(mytmp,'(I10,4X)',advance='no') sectors(i)%ndim
                     write(mytmp,'(I10,8X)',advance='no') j
                     write(mytmp,'(14I1)') bin_basis(:, sectors(i)%basis(j))
                 enddo ! over j={1,sectors(i)%ndim} loop
             enddo ! over i={1,nsectors} loop

         case(4)
             write(mytmp,'(a)',advance='no') '#      i |'
             write(mytmp,'(a)',advance='no') ' electron(i) |       Sz(i) |       PS(i) '
             write(mytmp,'(a)') '|     ndim(i) |           j |    fock_basis(j,i) |'
             do i=1,nsectors
                 do j=1,sectors(i)%ndim
                     write(mytmp,'(I10,4X)',advance='no') i
                     write(mytmp,'(I10,4X)',advance='no') sect_good_ntot(i)
                     write(mytmp,'(I10,4X)',advance='no') sect_good_sz(i)
                     write(mytmp,'(I10,4X)',advance='no') sect_good_ps(i)
                     write(mytmp,'(I10,4X)',advance='no') sectors(i)%ndim
                     write(mytmp,'(I10,8X)',advance='no') j
                     write(mytmp,'(14I1)') bin_basis(:, sectors(i)%basis(j))
                 enddo ! over j={1,sectors(i)%ndim} loop
             enddo ! over i={1,nsectors} loop

          case(5)
              write(mytmp,'(a)',advance='no') '#      i |'
              write(mytmp,'(a)',advance='no') ' electron(i) |       Jz(i) '
              write(mytmp,'(a)') '|     ndim(i) |           j |   fock_basis(j,i) |'
              do i=1,nsectors
                  do j=1,sectors(i)%ndim
                     write(mytmp,'(I10,4X)',advance='no') i
                     write(mytmp,'(I10,4X)',advance='no') sect_good_ntot(i)
                     write(mytmp,'(I10,4X)',advance='no') sect_good_jz(i)
                     write(mytmp,'(I10,4X)',advance='no') sectors(i)%ndim
                     write(mytmp,'(I10,8X)',advance='no') j
                     write(mytmp,'(14I1)') bin_basis(:, sectors(i)%basis(j))
                  enddo ! over j={1,sectors(i)%ndim} loop
              enddo ! over i={1,nsectors} loop
     end select ! back select case (ictqmc) block

! close data file
     close(mytmp)

     return
   end subroutine atomic_dump_sector
