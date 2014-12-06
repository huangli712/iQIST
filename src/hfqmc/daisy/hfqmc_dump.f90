!!!-----------------------------------------------------------------------
!!! project : daisy
!!! program : hfqmc_dump_gtau
!!!           hfqmc_dump_wtau
!!!           hfqmc_dump_gbin
!!!           hfqmc_dump_grnf
!!!           hfqmc_dump_wssf
!!!           hfqmc_dump_sigf
!!!           hfqmc_dump_nmat
!!!           hfqmc_dump_quas
!!! source  : hfqmc_dump.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 12/23/2009 by li huang
!!!           03/08/2010 by li huang
!!!           12/06/2014 by li huang
!!! purpose : To dump key observables produced by the Hirsch-Fye quantum
!!!           Monte Carlo (HFQMC) quantum impurity solver and dynamical
!!!           mean field theory (DMFT) self-consistent engine to files
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> hfqmc_dump_gtau: write out impurity green's function in imaginary
!!>>> time space
  subroutine hfqmc_dump_gtau(tmesh, gtau)
     use constants, only : dp, one, mytmp

     use control, only : nband, norbs
     use control, only : ntime
     use control, only : beta

     implicit none

! external arguments
! imaginary time mesh
     real(dp), intent(in) :: tmesh(ntime)

! impurity green's function
     real(dp), intent(in) :: gtau(ntime,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: solver.green.dat
     open(mytmp, file='solver.green.dat', form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,ntime
             write(mytmp,'(2i6,3f12.6)') i, j, tmesh(j), -gtau(j,i), -gtau(j,i+nband)
         enddo ! over j={1,ntime} loop
         write(mytmp,'(2i6,3f12.6)') i, ntime+1, beta, gtau(1,i)-one, gtau(1,i+nband)-one
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine hfqmc_dump_gtau

!>>> write out bath weiss's function in imaginary time space
  subroutine hfqmc_dump_wtau(tmesh, wtau)
     use constants
     use control

     implicit none

! external arguments
! imaginary time mesh
     real(dp), intent(in) :: tmesh(ntime)

! bath weiss's function
     real(dp), intent(in) :: wtau(ntime,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: solver.weiss.dat
     open(mytmp, file='solver.weiss.dat', form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,ntime
             write(mytmp,'(2i5,3f12.6)') i, j, tmesh(j), -wtau(j,i), -wtau(j,i+nband)
         enddo ! over j={1,ntime} loop
         write(mytmp,'(2i5,3f12.6)') i, ntime+1, beta, wtau(1,i)-one, wtau(1,i+nband)-one
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine hfqmc_dump_wtau

!>>> write out impurity green's function in imaginary time space (binning mode)
  subroutine hfqmc_dump_gbin(ibin, tmesh, gtau)
     use constants
     use control

     implicit none

! external arguments
! current bin index, integer representation
     integer, intent(in)  :: ibin

! imaginary time mesh
     real(dp), intent(in) :: tmesh(ntime)

! impurity green's function
     real(dp), intent(in) :: gtau(ntime,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! current bin index, string representation
     character(len=10) :: sbin

! open data file: solver.green.bin.x
     write(sbin,'(i10)') ibin ! convert ibin to sbin
     open(mytmp, file='solver.green.bin.'//trim(adjustl(sbin)), form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,ntime
             write(mytmp,'(2i5,3f12.6)') i, j, tmesh(j), -gtau(j,i), -gtau(j,i+nband)
         enddo ! over j={1,ntime} loop
         write(mytmp,'(2i5,3f12.6)') i, ntime+1, beta, gtau(1,i)-one, gtau(1,i+nband)-one
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine hfqmc_dump_gbin

!>>> write out impurity green's function in matsubara frequency space
  subroutine hfqmc_dump_grnf(rmesh, grnf)
     use constants
     use control

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in) :: rmesh(mfreq)

! impurity green's function
     complex(dp), intent(in) :: grnf(mfreq,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: solver.grn.dat
     open(mytmp, file='solver.grn.dat', form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,mfreq
             write(mytmp,'(i5,5f16.8)') i, rmesh(j), &
                                    real(grnf(j,i)), &
                                   aimag(grnf(j,i)), &
                              real(grnf(j,i+nband)), &
                             aimag(grnf(j,i+nband))
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine hfqmc_dump_grnf

!>>> write out bath weiss's function in matsubara frequency space
  subroutine hfqmc_dump_wssf(rmesh, wssf)
     use constants
     use control

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in) :: rmesh(mfreq)

! bath weiss's function
     complex(dp), intent(in) :: wssf(mfreq,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: solver.wss.dat
     open(mytmp, file='solver.wss.dat', form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,mfreq
             write(mytmp,'(i5,5f16.8)') i, rmesh(j), &
                                    real(wssf(j,i)), &
                                   aimag(wssf(j,i)), &
                              real(wssf(j,i+nband)), &
                             aimag(wssf(j,i+nband))
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine hfqmc_dump_wssf

!>>> write out self-energy function in matsubara frequency space
  subroutine hfqmc_dump_sigf(rmesh, sigf)
     use constants
     use control

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in) :: rmesh(mfreq)

! self-energy function
     complex(dp), intent(in) :: sigf(mfreq,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: solver.sgm.dat
     open(mytmp, file='solver.sgm.dat', form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,mfreq
             write(mytmp,'(i5,5f16.8)') i, rmesh(j), &
                                    real(sigf(j,i)), &
                                   aimag(sigf(j,i)), &
                              real(sigf(j,i+nband)), &
                             aimag(sigf(j,i+nband))
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine hfqmc_dump_sigf

!>>> write out the occupation matrix and double occupation matrix
  subroutine hfqmc_dump_nmat(nmat, nnmat)
     use constants
     use control

     implicit none

! external arguments
! occupation matrix data
     real(dp), intent(in) :: nmat(norbs)

! double occupation matrix data
     real(dp), intent(in) :: nnmat(norbs,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: solver.nmat.dat
     open(mytmp, file='solver.nmat.dat', form='formatted', status='unknown')

! write it
     write(mytmp,'(a)') '  < n_i >   data:'
     do i=1,norbs
         write(mytmp,'(i5,f12.6)') i, nmat(i)
     enddo ! over i={1,norbs} loop
     write(mytmp,'(a5,f12.6)') 'sup', sum( nmat(1:nband) )
     write(mytmp,'(a5,f12.6)') 'sdn', sum( nmat(nband+1:norbs) )
     write(mytmp,'(a5,f12.6)') 'sum', sum( nmat(1:norbs) )

     write(mytmp,'(a)') '< n_i n_j > data:'
     do i=1,norbs
         do j=1,norbs
             write(mytmp,'(2i5,f12.6)') i, j, nnmat(i,j)
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine hfqmc_dump_nmat

!>>> write out the quasiparticle weight
  subroutine hfqmc_dump_quas(quas)
     use constants
     use control

     implicit none

! external arguments
! quasiparticle weight
     real(dp), intent(in) :: quas(norbs)

! local variables
! loop index
     integer :: i

! open data file: solver.quas.dat
     open(mytmp, file='solver.quas.dat', form='formatted', status='unknown')

! write it
     do i=1,norbs
         write(mytmp,'(i5,f12.6)') i, quas(i)
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine hfqmc_dump_quas
