!!!-----------------------------------------------------------------------
!!! project : iqist @ jasmine
!!! program : atomic_dump_fock
!!!           atomic_dump_tmat
!!!           atomic_dump_emat
!!!           atomic_dump_umat
!!!           atomic_dump_feigval
!!!           atomic_dump_feigvec
!!!           atomic_dump_fcix
!!!           atomic_dump_seigval
!!!           atomic_dump_seigvec
!!!           atomic_dump_scix
!!!           atomic_dump_sector
!!! source  : atomic_dump.f90
!!! type    : subroutines
!!! author  : yilin wang (email:qhwyl2006@126.com)
!!! history : 07/09/2014 by yilin wang (created)
!!!           01/23/2024 by li huang (last modified)
!!! purpose : write some essential arrays and data structures to files
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub atomic_dump_fock
!!
!! write Fock basis to file atom.fock.dat
!!
  subroutine atomic_dump_fock()
     use constants, only : mytmp

     use control, only : ncfgs

     use m_fock, only : bin_basis
     use m_fock, only : dec_basis
     use m_fock, only : ind_basis

     implicit none

!! local variables
     ! loop index over Fock states
     integer :: i

     ! used to draw a dashed line
     character (len=1) :: dash(75)

!! [body

     ! setup dash
     dash = '-'

     ! open file atom.fock.dat to write
     open(mytmp, file='atom.fock.dat', form='formatted', status='unknown')

     ! write the header
     write(mytmp,'(75a1)') dash ! dashed line
     write(mytmp,'(a)') '# i | decimal | index | binary'
     write(mytmp,'(75a1)') dash ! dashed line

     ! write the data
     do i=1,ncfgs
         write(mytmp,'(i6)', advance='no') i
         write(mytmp,'(i6)', advance='no') dec_basis(i)
         write(mytmp,'(i6)', advance='no') ind_basis(dec_basis(i))
         write(mytmp,'(4X,*(i1))') bin_basis(:,i)
     enddo ! over i={1,ncfgs} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine atomic_dump_fock

!!
!! @sub atomic_dump_tmat
!!
!! write the transformation matrix from the original basis to the
!! natural basis (eigenstates of H_{CFS} + H_{SOC})
!!
  subroutine atomic_dump_tmat()
     use constants, only : mytmp

     use control, only : norbs

     use m_spmat, only : tmat

!! local variables
     ! loop index
     integer :: i
     integer :: j

     ! used to draw a dashed line
     character (len=1) :: dash(75)

!! [body

     ! setup dash
     dash = '-'

     ! open file atom.tmat.dat to write
     open(mytmp, file='atom.tmat.dat', form='formatted', status='unknown')

     ! write the header
     write(mytmp,'(75a1)') dash ! dashed line
     write(mytmp,'(a)') '# i | j | tmat_real | tmat_imag'
     write(mytmp,'(75a1)') dash ! dashed line

     ! write the data
     do i=1,norbs
         do j=1,norbs
             write(mytmp,'(2i6,2f16.8)') i, j, tmat(i,j)
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine atomic_dump_tmat

!!
!! @sub atomic_dump_emat
!!
!! write onsite impurity energy on natural basis
!! it should be a diagonal matrix
!!
  subroutine atomic_dump_emat()
     use constants, only : mytmp

     use control, only : isoc
     use control, only : nband, norbs

     use m_spmat, only : emat

     implicit none

!! local variables
     ! loop index
     integer :: i

     ! auxiliary integer variable used to convert the spin sequence
     integer :: s_order

     ! used to draw a dashed line
     character (len=1) :: dash(75)

!! [body

     ! setup dash
     dash = '-'

     ! open file atom.emat.dat to write
     open(mytmp, file='atom.emat.dat', form='formatted', status='unknown')

     ! write the header
     write(mytmp,'(75a1)') dash ! dashed line
     write(mytmp,'(a)') '# i | emat_real | emat_imag'
     write(mytmp,'(75a1)') dash ! dashed line

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
         write(mytmp,'(i6,2f16.8)') i, emat(s_order,s_order)
     enddo ! over i={1,norbs} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine atomic_dump_emat

!!
!! @sub atomic_dump_umat
!!
!! write onsite Coulomb interaction tensor
!!
  subroutine atomic_dump_umat()
     use constants, only : dp
     use constants, only : zero, two
     use constants, only : epst
     use constants, only : mytmp

     use control, only : icu
     use control, only : nband, norbs

     use m_spmat, only : umat

     implicit none

!! local variables
     ! loop index
     integer  :: i
     integer  :: j
     integer  :: k
     integer  :: l

     ! rank-2 Coulomb interaction tensor
     real(dp) :: umat_t(norbs,norbs)

     ! used to draw a dashed line
     character (len=1) :: dash(75)

!! [body

     ! setup dash
     dash = '-'

     ! open file atom.umat.dat to write
     open(mytmp, file='atom.umat.dat', form='formatted', status='unknown')

     ! write the header
     write(mytmp,'(75a1)') dash ! dashed line
     write(mytmp,'(a)') '# i | j | k | l | umat_real | umat_imag'
     write(mytmp,'(75a1)') dash ! dashed line

     ! write the data, only the non-zero elements are outputed
     ! note: we do not change the spin sequence here
     do i=1,norbs
         do j=1,norbs
             do k=1,norbs
                 do l=1,norbs
                     if ( abs( umat(i,j,k,l) ) > epst ) then
                         write(mytmp,'(4i6,2f16.8)') i, j, k, l, umat(i,j,k,l)
                     endif ! back if ( abs( umat(i,j,k,l) ) > epst ) block
                 enddo ! over l={1,norbs} loop
             enddo ! over k={1,norbs} loop
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

     ! close data file
     close(mytmp)

     ! initialize rank-2 Coulomb interaction tensor
     umat_t = zero

     ! Kanamori type
     if ( icu == 1 .or. icu == 3 ) then
         do i=1,norbs
             do j=i+1,norbs
                 umat_t(i,j) = real( umat(i,j,j,i) )
                 umat_t(j,i) = umat_t(i,j)
             enddo ! over j={i+1,norbs} loop
         enddo ! over i={1,norbs} loop
     ! Slater type
     elseif ( icu == 2 ) then
         do i=1,norbs
             do j=i+1,norbs
                 if ( mod(i,2) == mod(j,2) ) then
                     umat_t(i,j) = two * real(umat(i,j,j,i) - umat(i,j,i,j))
                 else
                     umat_t(i,j) = two * real(umat(i,j,j,i))
                 endif  ! back if ( mod(i,2) == mod(j,2) ) block
                 umat_t(j,i) = umat_t(i,j)
             enddo ! over j={i+1,norbs} loop
         enddo ! over i={1,norbs} loop
     endif ! back if ( icu == 1 .or. icu == 3 ) block

     ! open file solver.umat.in to write
     ! this file is used as input for the the other ctqmc code
     open(mytmp, file='solver.umat.in', form='formatted', status='unknown')

     ! write the data, all of the elements are written
     ! note: we have to change the spin sequence here
     do i=1,norbs
         if ( i <= nband ) then
             k = 2*i-1
         else
             k = 2*(i-nband)
         endif ! back if ( i <= nband ) block

         do j=1,norbs
             if ( j <= nband ) then
                 l = 2*j-1
             else
                 l = 2*(j-nband)
             endif ! back if ( j <= nband ) block

             write(mytmp,'(2i6,f16.8)') i, j, umat_t(k,l)
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine atomic_dump_umat

!!
!! @sub atomic_dump_feigval
!!
!! write eigenvalues in full Hilbert space to the file atom.eigval.dat
!!
  subroutine atomic_dump_feigval()
     use constants, only : mytmp

     use control, only : ncfgs

     use m_fock, only : eval

     implicit none

!! local variables
     ! loop index
     integer :: i

     ! used to draw a dashed line
     character (len=1) :: dash(75)

!! [body

     ! setup dash
     dash = '-'

     ! open file atom.eigval.dat to write
     open(mytmp, file='atom.eigval.dat', form='formatted', status='unknown')

     ! write the header
     write(mytmp,'(75a1)') dash ! dashed line
     write(mytmp,'(a)') '# i | eigenvalues'
     write(mytmp,'(75a1)') dash ! dashed line

     ! write the data
     do i=1,ncfgs
         write(mytmp,'(i6,f16.8)') i, eval(i)
     enddo ! over i={1,ncfgs} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine atomic_dump_feigval

!!
!! @sub atomic_dump_feigvec
!!
!! write eigenvectors in full Hilbert space to the file atom.eigvec.dat
!!
  subroutine atomic_dump_feigvec()
     use constants, only : eps6
     use constants, only : mytmp

     use control, only : ncfgs

     use m_fock, only : evec

     implicit none

!! local variables
     ! loop index
     integer :: i
     integer :: j

     ! used to draw a dashed line
     character (len=1) :: dash(75)

!! [body

     ! setup dash
     dash = '-'

     ! open file atom.eigvec.dat to write
     open(mytmp, file='atom.eigvec.dat', form='formatted', status='unknown')

     ! write the header
     write(mytmp,'(75a1)') dash ! dashed line
     write(mytmp,'(a)') '# i | j | eigenvectors'
     write(mytmp,'(75a1)') dash ! dashed line

     ! write the data
     do i=1,ncfgs
         do j=1,ncfgs
             if ( abs( evec(j,i) ) > eps6 ) then
                 write(mytmp,'(2i6,f16.8)') j, i, evec(j,i)
             endif ! back if ( abs( evec(j,i) ) > eps6 ) block
         enddo ! over j={1,ncfgs} loop
     enddo ! over i={1,ncfgs} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine atomic_dump_feigvec

!!
!! @sub atomic_dump_fcix
!!
!! write atom.cix file for the ctqmc solver. the file format is designed
!! for the begonia and lavender codes. however, the two codes are already
!! deprecated. we retain this subroutine only for reference. 
!!
  subroutine atomic_dump_fcix()
     use constants, only : epst
     use constants, only : mytmp

     use control, only : cname
     use control, only : ictqmc, icu, icf, isoc
     use control, only : nband, nspin, norbs, ncfgs
     use control, only : nmini, nmaxi
     use control, only : Uc, Uv, Js, Jp, Jz
     use control, only : Ud, Jh
     use control, only : mune, lambda

     use version, only : V_MAIL 

     use m_fock, only : eval, evec
     use m_fock, only : occu, spin
     use m_fock, only : fmat

     implicit none

!! local variables
     ! loop index
     integer :: i
     integer :: j
     integer :: k

     ! auxiliary integer variable used to convert the spin sequence
     integer :: s_order

     ! used to draw a dashed line
     character (len=1) :: dash(75)

     ! string for current date and time
     character (len = 20) :: date_time_string

!! [body

     ! setup dash
     dash = '-'

     ! obtain current date and time
     call s_time_builder(date_time_string)

     ! open file atom.cix to write
     open(mytmp, file='atom.cix', form='formatted', status='unknown')

     ! write the header
     write(mytmp,'(a)') '# WARNING : DO NOT MODIFY THIS FILE MANUALLY!'
     write(mytmp,'(a)') '# File    : atom.cix'
     write(mytmp,'(a)') '# Format  : v1.3, designed for BEGONIA and LAVENDER'
     write(mytmp,'(a)') '# Built   : by '//cname//' code at '//date_time_string
     write(mytmp,'(a)') '# Support : any problem, please contact me: '//V_MAIL
     write(mytmp,*)
     write(mytmp,*)

     ! write configurations
     write(mytmp,'(75a1)') dash ! dashed line
     write(mytmp,'(a)') '# PARAMETERS:'
     write(mytmp,'(75a1)') dash ! dashed line
     write(mytmp,'(4i8,12X,a)') 1, icu, icf, isoc, 'VER ICU ICF ISOC'
     write(mytmp,'(4i8,12X,a)') nband, nspin, norbs, ncfgs, 'NBAND NSPIN NORBS NCFGS'
     write(mytmp,'(2i8,28X,a)') nmini, nmaxi, 'NMINI NMAXI'
     write(mytmp,'(5f8.4,04X,a)') Uc, Uv, Jz, Js, Jp, 'Uc Uv Jz Js Jp'
     write(mytmp,'(2f8.4,28X,a)') Ud, Jh, 'Ud Jh'
     write(mytmp,'(2f8.4,28X,a)') mune, lambda, 'mune lambda'

     ! write eigenvalues
     write(mytmp,'(75a1)') dash ! dashed line
     write(mytmp,'(a)') '# EIGENVALUES: INDEX | ENERGY | OCCUPY | SPIN'
     write(mytmp,'(75a1)') dash ! dashed line
     do i=1,ncfgs
         write(mytmp,'(i10,3f20.10)') i, eval(i), occu(i,i), spin(i,i)
     enddo ! over i={1,ncfgs} loop

     ! write eigenvectors
     ! only for the camellia code
     if ( ictqmc == 0 ) then
         write(mytmp,'(75a1)') dash ! dashed line
         write(mytmp,'(a)') '# EIGENVECTORS: ALPHA | BETA | EVEC'
         write(mytmp,'(75a1)') dash ! dashed line
         do i=1,ncfgs
             do j=1,ncfgs
                 if ( abs( evec(i,j) ) > epst ) then
                     write(mytmp,'(2i10,f20.10)') i, j, evec(i,j)
                 endif ! back if ( abs( evec(i,j) ) > epst ) block
             enddo ! over j={1,ncfgs} loop
         enddo ! over i={1,ncfgs} loop
     endif ! back if ( ictqmc == 0 ) block

     ! write annihilation operator matrix in atomic eigenbasis
     !
     ! for non-soc case, the spin order of ctqmc is like
     !     up, up, up, dn, dn, dn
     !
     ! but the spin order of this program is
     !     up, dn, up, dn, up, dn
     !
     ! so we have to adjust it here. however if the spin-orbit coupling
     ! is enabled, the spin order doesn't matter.
     write(mytmp,'(75a1)') dash ! dashed line
     write(mytmp,'(a)') '# F MATRIX ELEMENT: ALPHA | BETA | FLAVOR | FMAT'
     write(mytmp,'(75a1)') dash ! dashed line
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
         !
         do j=1,ncfgs
             do k=1,ncfgs
                 if ( abs( fmat(k,j,s_order) ) > epst ) then
                     write(mytmp,'(3i10,f20.10)') k, j, i, fmat(k,j,s_order)
                 endif ! back if ( abs( fmat(k,j,s_order) ) > epst ) block
             enddo ! back k={1,ncfgs} loop
         enddo ! over j={1,ncfgs} loop
     enddo ! over i={1,norbs} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine atomic_dump_fcix

!!
!! @sub atomic_dump_seigval
!!
!! write eigenvalues in all sectors to the file atom.eigval.dat
!!
  subroutine atomic_dump_seigval()
     use constants, only : mytmp

     use m_sector, only : nsectors
     use m_sector, only : sectors

     implicit none

!! local variables
     ! loop index
     integer :: i
     integer :: j

     ! the counter
     integer :: counter

     ! used to draw a dashed line
     character (len=1) :: dash(75)

!! [body

     ! setup dash
     dash = '-'

     ! open file atom.eigval.dat to write
     open(mytmp, file='atom.eigval.dat', form='formatted', status='unknown')

     ! write the header
     write(mytmp,'(75a1)') dash ! dashed line
     write(mytmp,'(a)') '# i | sector | j | eigenvalues'
     write(mytmp,'(75a1)') dash ! dashed line

     ! write the data
     counter = 0
     do i=1,nsectors
         do j=1,sectors(i)%ndim
             counter = counter + 1
             write(mytmp,'(3i6,f16.8)') counter, i, j, sectors(i)%eval(j)
         enddo ! over i={1,nsectors} loop
     enddo ! over j={1,sectors(i)%ndim} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine atomic_dump_seigval

!!
!! @sub atomic_dump_seigvec
!!
!! write eigenvectors in all sectors to the file atom.eigvec.dat
!!
  subroutine atomic_dump_seigvec()
     use constants, only : eps6
     use constants, only : mytmp

     use m_sector, only : nsectors
     use m_sector, only : sectors

     implicit none

!! local variables
     ! loop index
     integer :: i
     integer :: j
     integer :: k

     ! used to draw a dashed line
     character (len=1) :: dash(75)

!! [body

     ! setup dash
     dash = '-'

     ! open file atom.eigvec.dat to write
     open(mytmp, file='atom.eigvec.dat', form='formatted', status='unknown')

     ! write the header
     write(mytmp,'(75a1)') dash ! dashed line
     write(mytmp,'(a)') '# sector | i | j | eigenvectors'
     write(mytmp,'(75a1)') dash ! dashed line

     ! write the data
     do i=1,nsectors
         do j=1,sectors(i)%ndim
             do k=1,sectors(i)%ndim
                 if ( abs( sectors(i)%evec(k,j) ) > eps6 ) then
                     write(mytmp,'(3i6,f16.8)') i, k, j, sectors(i)%evec(k,j)
                 endif ! back if ( abs( sectors(i)%evec(j,k) ) > eps6 ) block
             enddo ! over k={1,sectors(i)%ndim} loop
         enddo ! over j={1,sectors(i)%ndim} loop
     enddo ! over i={1,nsectors} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine atomic_dump_seigvec

!!
!! @sub atomic_dump_scix
!!
!! write atom.cix file for the ctqmc solver. the file format is designed
!! for the manjushaka code.
!!
  subroutine atomic_dump_scix()
     use constants, only : epst
     use constants, only : mytmp

     use control, only : cname
     use control, only : icu, icf, isoc
     use control, only : nband, nspin, norbs, ncfgs
     use control, only : nmini, nmaxi
     use control, only : Uc, Uv, Jz, Js, Jp
     use control, only : Ud, Jh
     use control, only : mune, lambda

     use version, only : V_MAIL 

     use m_sector, only : nsectors
     use m_sector, only : sectors
     use m_sector, only : max_dim_sect
     use m_sector, only : ave_dim_sect

     implicit none

!! local variables
     ! loop index
     integer :: i
     integer :: j
     integer :: k
     integer :: m
     integer :: n

     ! counter for the non-zero matrix elements
     integer :: counter

     ! auxiliary integer variable used to convert the spin sequence
     integer :: s_order

     ! used to draw a dashed line
     character (len=1) :: dash(75)

     ! string for current date and time
     character (len = 20) :: date_time_string

!! [body

     ! setup dash
     dash = '-'

     ! obtain current date and time
     call s_time_builder(date_time_string)

     ! open atom.cix to write
     open(mytmp, file='atom.cix', form='formatted', status='unknown')

     ! write header
     write(mytmp,'(a)') '# WARNING : DO NOT MODIFY THIS FILE MANUALLY!'
     write(mytmp,'(a)') '# File    : atom.cix'
     write(mytmp,'(a)') '# Format  : v2.3, designed for MANJUSHAKA'
     write(mytmp,'(a)') '# Built   : by '//cname//' code at '//date_time_string
     write(mytmp,'(a)') '# Support : any problem, please contact me: '//V_MAIL
     write(mytmp,*)
     write(mytmp,*)

     ! write configurations
     write(mytmp,'(75a1)') dash ! dashed line
     write(mytmp,'(a)') '# PARAMETERS:'
     write(mytmp,'(75a1)') dash ! dashed line
     write(mytmp,'(4i8,12X,a)') 2, icu, icf, isoc, 'VER ICU ICF ISOC'
     write(mytmp,'(4i8,12X,a)') nband, nspin, norbs, ncfgs, 'NBAND NSPIN NORBS NCFGS'
     write(mytmp,'(2i8,28X,a)') nmini, nmaxi, 'NMINI NMAXI'
     write(mytmp,'(5f8.4,04X,a)') Uc, Uv, Jz, Js, Jp, 'Uc Uv Jz Js Jp'
     write(mytmp,'(2f8.4,28X,a)') Ud, Jh, 'Ud Jh'
     write(mytmp,'(2f8.4,28X,a)') mune, lambda, 'mune lambda'

     ! write summary of sectors
     write(mytmp,'(75a1)') dash ! dashed line
     write(mytmp,'(a)') '# SECTORS:'
     write(mytmp,'(75a1)') dash ! dashed line
     write(mytmp,'(a)') '# SUMMARY: NSECTORS | MAX_DIM_SECT | AVE_DIM_SECT'
     write(mytmp,'(5X,2(i10,2X),f20.10)') nsectors, max_dim_sect, ave_dim_sect

     ! write dimension, total electrons, next sector, eigenvalue of each sector
     do i=1,nsectors
         write(mytmp,'(a)') '# SECT_INFO: INDEX | NDIM | NOPS | ISTART | NE | SZ | JZ | PS'
         write(mytmp,'(12X,8i6)') i, sectors(i)%ndim,   &
                                     sectors(i)%nops,   &
                                     sectors(i)%istart, &
                                     sectors(i)%nele,   &
                                     sectors(i)%sz,     &
                                     sectors(i)%jz,     &
                                     sectors(i)%ps

         ! write next sector
         write(mytmp,'(4X,a)') '# NEXT SECTOR    F     F^{\DAGGER}'
         do j=1,sectors(i)%nops
             ! adjust the orbital order for ctqmc, up, up, up, dn, dn, dn
             if ( isoc == 0 ) then
                 if (j <= sectors(i)%nops / 2) then
                     s_order = 2*j-1
                 else
                     s_order = 2*(j - sectors(i)%nops / 2)
                 endif ! back if (j <= sectors(i)%nops / 2) block
             else
                 s_order = j
             endif ! back if ( isoc == 0 ) block
             write(mytmp,'(2X,3i10)') j, sectors(i)%next(s_order,0), sectors(i)%next(s_order,1)
         enddo ! over j={1,sectors(i)%nops} loop

         ! write eigeanvalue
         write(mytmp,'(4X,a)') '# EIGENVALUES'
         do j=1,sectors(i)%ndim
             write(mytmp,'(2X,i10,f20.10)') j, sectors(i)%eval(j)
         enddo ! over j={1,sectors(i)%ndim} loop
     enddo ! over i={1,nsectors} loop

     ! write F-matrix of each sector
     write(mytmp,'(75a1)') dash ! dashed line
     write(mytmp,'(a)') '# F-MATRIX:'
     write(mytmp,'(75a1)') dash ! dashed line

     ! write the data
     do i=1,nsectors
         do j=1,sectors(i)%nops
             ! adjust the orbital order for ctqmc, up, up, up, dn, dn, dn
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
                 write(mytmp,'(a)') '# SECTOR | FLAVOR | DAGGER |      N |      M | SPARSE'
                 write(mytmp,'(2X,6(i6,3X))') i, j, k, &
                                              sectors(i)%fmat(s_order,k)%n, &
                                              sectors(i)%fmat(s_order,k)%m, &
                                              count( abs(sectors(i)%fmat(s_order,k)%val) > epst )
                 counter = 0
                 do n=1,sectors(i)%fmat(s_order,k)%n
                     do m=1,sectors(i)%fmat(s_order,k)%m
                         if ( abs( sectors(i)%fmat(s_order,k)%val(n,m) ) > epst ) then
                             counter = counter + 1
                             write(mytmp,'(2i6,f20.10)') n, m, sectors(i)%fmat(s_order,k)%val(n,m)
                         endif ! back if ( abs( sectors(i)%fmat(s_order,k)%val(n,m) ) > epst ) block
                     enddo ! over m={1,sectors(i)%fmat(s_order,k)%m} loop
                 enddo ! over n={1,sectors(i)%fmat(s_order,k)%n} loop
                 call s_assert( counter == count( abs(sectors(i)%fmat(s_order,k)%val) > epst ) )
             enddo  ! over k={0,1} loop
         enddo ! over j={1,sectors(i)%nops} loop
     enddo  ! over i={1,nsect} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine atomic_dump_scix

!!
!! @sub atomic_dump_sector
!!
!! write out the configuration for each sector to file atom.sector.dat
  subroutine atomic_dump_sector(sect_good_ntot, sect_good_sz, sect_good_ps, sect_good_jz)
     use constants, only : mytmp

     use control, only : ictqmc
     use control, only : ncfgs

     use m_fock, only : bin_basis

     use m_sector, only : nsectors
     use m_sector, only : sectors

     use m_sector, only : max_dim_sect
     use m_sector, only : ave_dim_sect

     implicit none

!! external arguments
     ! good quantum number: N
     integer, intent(in) :: sect_good_ntot(ncfgs)

     ! good quantum number: Sz
     integer, intent(in) :: sect_good_sz(ncfgs)

     ! good quantum number: PS
     integer, intent(in) :: sect_good_ps(ncfgs)

     ! good quantum number: Jz
     integer, intent(in) :: sect_good_jz(ncfgs)

!! local variables
     ! loop index
     integer :: i
     integer :: j

     ! used to draw a dashed line
     character (len=1) :: dash(75)

!! [body

     ! setup dash
     dash = '-'

     ! open 'atom.sector.dat' to write
     open(mytmp, file='atom.sector.dat', form='formatted', status='unknown')

     ! write header
     write(mytmp,'(a,i10)')   '# number_sectors : ', nsectors
     write(mytmp,'(a,i10)')   '# max_dim_sectors: ', max_dim_sect
     write(mytmp,'(a,f10.5)') '# ave_dim_sectors: ', ave_dim_sect
     write(mytmp,'(75a1)') dash ! dashed line

     ! write the data
     select case (ictqmc)
         case (1)
             call s_print_error('atomic_dump_sector','this case is not implemented')

         case (2)
             write(mytmp,'(a)') '# i | Ntot | Ndim | j | Fock'
             write(mytmp,'(75a1)') dash ! dashed line
             do i=1,nsectors
                 do j=1,sectors(i)%ndim
                     write(mytmp,'(i4,2X)',advance='no') i
                     write(mytmp,'(i4,2X)',advance='no') sectors(i)%nele
                     write(mytmp,'(i4,2X)',advance='no') sectors(i)%ndim
                     write(mytmp,'(i4,2X)',advance='no') j
                     write(mytmp,'(14i1)') bin_basis(:,sectors(i)%basis(j))
                 enddo ! over j={1,sectors(i)%ndim} loop
             enddo ! over i={1,nsectors} loop

         case(3)
             write(mytmp,'(a)') '# i | Ntot | Sz | Ndim | j | Fock'
             write(mytmp,'(75a1)') dash ! dashed line
             do i=1,nsectors
                 do j=1,sectors(i)%ndim
                     write(mytmp,'(i4,2X)',advance='no') i
                     write(mytmp,'(i4,2X)',advance='no') sect_good_ntot(i)
                     write(mytmp,'(i4,2X)',advance='no') sect_good_sz(i)
                     write(mytmp,'(i4,2X)',advance='no') sectors(i)%ndim
                     write(mytmp,'(i4,2X)',advance='no') j
                     write(mytmp,'(14i1)') bin_basis(:,sectors(i)%basis(j))
                 enddo ! over j={1,sectors(i)%ndim} loop
             enddo ! over i={1,nsectors} loop

         case(4)
             write(mytmp,'(a)') '# i | Ntot | Sz | PS | Ndim | j | Fock'
             write(mytmp,'(75a1)') dash ! dashed line
             do i=1,nsectors
                 do j=1,sectors(i)%ndim
                     write(mytmp,'(i4,2X)',advance='no') i
                     write(mytmp,'(i4,2X)',advance='no') sect_good_ntot(i)
                     write(mytmp,'(i4,2X)',advance='no') sect_good_sz(i)
                     write(mytmp,'(i4,2X)',advance='no') sect_good_ps(i)
                     write(mytmp,'(i4,2X)',advance='no') sectors(i)%ndim
                     write(mytmp,'(i4,2X)',advance='no') j
                     write(mytmp,'(14i1)') bin_basis(:,sectors(i)%basis(j))
                 enddo ! over j={1,sectors(i)%ndim} loop
             enddo ! over i={1,nsectors} loop

          case(5)
              write(mytmp,'(a)') '# i | Ntot | Jz | Ndim | j | Fock'
              write(mytmp,'(75a1)') dash ! dashed line
              do i=1,nsectors
                  do j=1,sectors(i)%ndim
                     write(mytmp,'(i4,2X)',advance='no') i
                     write(mytmp,'(i4,2X)',advance='no') sect_good_ntot(i)
                     write(mytmp,'(i4,2X)',advance='no') sect_good_jz(i)
                     write(mytmp,'(i4,2X)',advance='no') sectors(i)%ndim
                     write(mytmp,'(i4,2X)',advance='no') j
                     write(mytmp,'(14i1)') bin_basis(:,sectors(i)%basis(j))
                  enddo ! over j={1,sectors(i)%ndim} loop
              enddo ! over i={1,nsectors} loop
     end select ! back select case (ictqmc) block

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine atomic_dump_sector
