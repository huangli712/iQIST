!!!-----------------------------------------------------------------------
!!! project : camellia
!!! program : ctqmc_dump_gtau
!!!           ctqmc_dump_wtau
!!!           ctqmc_dump_htau
!!!           ctqmc_dump_grnf
!!!           ctqmc_dump_wssf
!!!           ctqmc_dump_hybf
!!!           ctqmc_dump_sigf
!!!           ctqmc_dump_hub1
!!!           ctqmc_dump_hist
!!!           ctqmc_dump_prob
!!!           ctqmc_dump_nmat
!!!           ctqmc_dump_kmat
!!!           ctqmc_dump_lmat
!!!           ctqmc_dump_twop
!!!           ctqmc_dump_pair
!!! source  : ctqmc_dump.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           08/17/2015 by li huang (last modified)
!!! purpose : dump key observables produced by the hybridization expansion
!!!           version continuous time quantum Monte Carlo (CTQMC) quantum
!!!           impurity solver and dynamical mean field theory (DMFT) self
!!!           -consistent engine to disk files
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> dump data on imaginary time axis                                 <<<
!!========================================================================

!!>>> ctqmc_dump_gtau: write out impurity green's function in imaginary
!!>>> time space
  subroutine ctqmc_dump_gtau(tmesh, gtau, gerr)
     use constants, only : dp, mytmp

     use control, only : norbs
     use control, only : ntime

     implicit none

! external arguments
! imaginary time mesh
     real(dp), intent(in) :: tmesh(ntime)

! impurity green's function
     real(dp), intent(in) :: gtau(ntime,norbs,norbs)
     real(dp), intent(in) :: gerr(ntime,norbs,norbs)

! local variables
! loop index
     integer  :: i
     integer  :: j

! scaled impurity green's function
     real(dp) :: gaux(ntime,norbs,norbs)
     real(dp) :: gtmp(ntime,norbs,norbs)

! evaluate gaux and gtmp at first
     call ctqmc_make_gtau(tmesh, gtau, gaux)
     call ctqmc_make_gtau(tmesh, gerr, gtmp)

! open data file: solver.green.dat
     open(mytmp, file='solver.green.dat', form='formatted', status='unknown')

! write it
     do i=1,norbs
         do j=1,ntime
             write(mytmp,'(2i6,3f12.6)') i, j, tmesh(j), gaux(j,i,i), gtmp(j,i,i)
         enddo ! over j={1,ntime} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_gtau

!!>>> ctqmc_dump_wtau: write out bath weiss's function in imaginary
!!>>> time space
  subroutine ctqmc_dump_wtau(tmesh, wtau)
     use constants, only : dp, zero, mytmp

     use control, only : norbs
     use control, only : ntime

     implicit none

! external arguments
! imaginary time mesh
     real(dp), intent(in) :: tmesh(ntime)

! bath weiss's function
     real(dp), intent(in) :: wtau(ntime,norbs,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: solver.weiss.dat
     open(mytmp, file='solver.weiss.dat', form='formatted', status='unknown')

! write it
     do i=1,norbs
         do j=1,ntime
             write(mytmp,'(2i6,3f12.6)') i, j, tmesh(j), wtau(j,i,i), zero
         enddo ! over j={1,ntime} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_wtau

!!>>> ctqmc_dump_htau: write out hybridization function in imaginary
!!>>> time space
  subroutine ctqmc_dump_htau(tmesh, htau)
     use constants, only : dp, zero, mytmp

     use control, only : norbs
     use control, only : ntime

     implicit none

! external arguments
! imaginary time mesh
     real(dp), intent(in) :: tmesh(ntime)

! hybridization function
     real(dp), intent(in) :: htau(ntime,norbs,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: solver.hybri.dat
     open(mytmp, file='solver.hybri.dat', form='formatted', status='unknown')

! write it
     do i=1,norbs
         do j=1,ntime
             write(mytmp,'(2i6,3f12.6)') i, j, tmesh(j), htau(j,i,i), zero
         enddo ! over j={1,ntime} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_htau

!!========================================================================
!!>>> dump data on matsubara frequency axis                            <<<
!!========================================================================

!!>>> ctqmc_dump_grnf: write out impurity green's function in matsubara
!!>>> frequency space
  subroutine ctqmc_dump_grnf(rmesh, grnf, gerr)
     use constants, only : dp, mytmp

     use control, only : norbs
     use control, only : mfreq

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(mfreq)

! impurity green's function
     complex(dp), intent(in) :: grnf(mfreq,norbs,norbs)
     complex(dp), intent(in) :: gerr(mfreq,norbs,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: solver.grn.dat
     open(mytmp, file='solver.grn.dat', form='formatted', status='unknown')

! write it
     do i=1,norbs
         do j=1,mfreq
             write(mytmp,'(i6,5f16.8)') i, rmesh(j), &
              real(grnf(j,i,i)), aimag(grnf(j,i,i)), &
              real(gerr(j,i,i)), aimag(gerr(j,i,i))
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_grnf

!!>>> ctqmc_dump_wssf: write out bath weiss's function in matsubara
!!>>> frequency space
  subroutine ctqmc_dump_wssf(rmesh, wssf)
     use constants, only : dp, zero, mytmp

     use control, only : norbs
     use control, only : mfreq

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(mfreq)

! bath weiss's function
     complex(dp), intent(in) :: wssf(mfreq,norbs,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: solver.wss.dat
     open(mytmp, file='solver.wss.dat', form='formatted', status='unknown')

! write it
     do i=1,norbs
         do j=1,mfreq
             write(mytmp,'(i6,5f16.8)') i, rmesh(j), &
              real(wssf(j,i,i)), aimag(wssf(j,i,i)), &
                                         zero, zero
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_wssf

!!>>> ctqmc_dump_hybf: write out hybridization function in matsubara
!!>>> frequency space
  subroutine ctqmc_dump_hybf(rmesh, hybf)
     use constants, only : dp, zero, mytmp

     use control, only : norbs
     use control, only : mfreq

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(mfreq)

! hybridization function
     complex(dp), intent(in) :: hybf(mfreq,norbs,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: solver.hyb.dat
     open(mytmp, file='solver.hyb.dat', form='formatted', status='unknown')

! write it
     do i=1,norbs
         do j=1,mfreq
             write(mytmp,'(i6,5f16.8)') i, rmesh(j), &
              real(hybf(j,i,i)), aimag(hybf(j,i,i)), &
                                         zero, zero
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_hybf

!!>>> ctqmc_dump_sigf: write out self-energy function in matsubara
!!>>> frequency space
  subroutine ctqmc_dump_sigf(rmesh, sigf)
     use constants, only : dp, zero, mytmp

     use control, only : norbs
     use control, only : mfreq

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(mfreq)

! self-energy function
     complex(dp), intent(in) :: sigf(mfreq,norbs,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: solver.sgm.dat
     open(mytmp, file='solver.sgm.dat', form='formatted', status='unknown')

! write it
     do i=1,norbs
         do j=1,mfreq
             write(mytmp,'(i6,5f16.8)') i, rmesh(j), &
              real(sigf(j,i,i)), aimag(sigf(j,i,i)), &
                                         zero, zero
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_sigf

!!>>> ctqmc_dump_hub1: write out impurity green's function and self-energy
!!>>> function obtained by hubbard-I approximation in matsubara frequency
!!>>> space
  subroutine ctqmc_dump_hub1(rmesh, ghub, shub)
     use constants, only : dp, mytmp

     use control, only : norbs
     use control, only : mfreq

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(mfreq)

! impurity green's function by hubbard-I approximation
     complex(dp), intent(in) :: ghub(mfreq,norbs)

! self-energy function by hubbard-I approximation
     complex(dp), intent(in) :: shub(mfreq,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: solver.hub.dat
     open(mytmp, file='solver.hub.dat', form='formatted', status='unknown')

! write it
     do i=1,norbs
         do j=1,mfreq
             write(mytmp,'(i6,5f16.8)') i, rmesh(j), &
                  real(ghub(j,i)), aimag(ghub(j,i)), &
                  real(shub(j,i)), aimag(shub(j,i))
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_hub1

!!========================================================================
!!>>> dump data of physical observables                                <<<
!!========================================================================

!!>>> ctqmc_dump_hist: write out the Monte Carlo sampling histogram for
!!>>> perturbation expansion series
  subroutine ctqmc_dump_hist(hist, herr)
     use constants, only : dp, mytmp

     use control, only : mkink

     implicit none

! external arguments
! histogram data
     real(dp), intent(in) :: hist(mkink)
     real(dp), intent(in) :: herr(mkink)

! local variables
! loop index
     integer  :: i

! scaled histogram data
     real(dp) :: haux(mkink)
     real(dp) :: htmp(mkink)

! evaluate haux and htmp at first
     haux = hist / sum(hist)
     htmp = herr / sum(hist)

! open data file: solver.hist.dat
     open(mytmp, file='solver.hist.dat', form='formatted', status='unknown')

! write it
     write(mytmp,'(a)') '# histogram: order | count | percent'
     do i=1,mkink
         write(mytmp,'(i6,i12,2f12.6)') i, int( hist(i) ), haux(i), htmp(i)
     enddo ! over i={1,mkink} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_hist
