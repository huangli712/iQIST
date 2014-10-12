!!!-----------------------------------------------------------------------
!!! project : narcissus
!!! program : ctqmc_dump_gtau
!!!           ctqmc_dump_wtau
!!!           ctqmc_dump_htau
!!!           ctqmc_dump_ktau
!!!           ctqmc_dump_gbin
!!!           ctqmc_dump_grnf
!!!           ctqmc_dump_wssf
!!!           ctqmc_dump_hybf
!!!           ctqmc_dump_sigf
!!!           ctqmc_dump_hub1
!!!           ctqmc_dump_hist
!!!           ctqmc_dump_prob
!!!           ctqmc_dump_nmat
!!!           ctqmc_dump_schi
!!!           ctqmc_dump_ochi
!!!           ctqmc_dump_twop
!!!           ctqmc_dump_vrtx
!!!           ctqmc_dump_pair
!!! source  : ctqmc_dump.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 09/16/2009 by li huang
!!!           08/23/2010 by li huang
!!!           09/12/2014 by li huang
!!!           10/13/2014 by li huang
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
  subroutine ctqmc_dump_gtau(tmesh, gtau)
     use constants, only : dp, mytmp

     use control, only : nband, norbs
     use control, only : ntime

     implicit none

! external arguments
! imaginary time mesh
     real(dp), intent(in) :: tmesh(ntime)

! impurity green's function
     real(dp), intent(in) :: gtau(ntime,norbs,norbs)

! local variables
! loop index
     integer  :: i
     integer  :: j

! scaled impurity green's function
     real(dp) :: gaux(ntime,norbs,norbs)

! evaluate gaux first
     call ctqmc_make_gtau(tmesh, gtau, gaux)

! open data file: solver.green.dat
     open(mytmp, file='solver.green.dat', form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,ntime
             write(mytmp,'(2i6,3f12.6)') i, j, tmesh(j), gaux(j,i,i), gaux(j,i+nband,i+nband)
         enddo ! over j={1,ntime} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_gtau

!!>>> ctqmc_dump_wtau: write out bath weiss's function in imaginary
!!>>> time space
  subroutine ctqmc_dump_wtau(tmesh, wtau)
     use constants, only : dp, mytmp

     use control, only : nband, norbs
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
     do i=1,nband
         do j=1,ntime
             write(mytmp,'(2i6,3f12.6)') i, j, tmesh(j), wtau(j,i,i), wtau(j,i+nband,i+nband)
         enddo ! over j={1,ntime} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_wtau

!!>>> ctqmc_dump_htau: write out hybridization function in imaginary
!!>>> time space
  subroutine ctqmc_dump_htau(tmesh, htau)
     use constants, only : dp, mytmp

     use control, only : nband, norbs
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
     do i=1,nband
         do j=1,ntime
             write(mytmp,'(2i6,3f12.6)') i, j, tmesh(j), htau(j,i,i), htau(j,i+nband,i+nband)
         enddo ! over j={1,ntime} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_htau

!!>>> ctqmc_dump_ktau: write out kernel function in imaginary time space
  subroutine ctqmc_dump_ktau(tmesh, ktau)
     use constants, only : dp, mytmp

     use control, only : ntime

     implicit none

! external arguments
! imaginary time mesh
     real(dp), intent(in) :: tmesh(ntime)

! kernel function
     real(dp), intent(in) :: ktau(ntime)

! local variables
! loop index
     integer :: i

! open data file: solver.kernel.dat
     open(mytmp, file='solver.kernel.dat', form='formatted', status='unknown')

! write it
     do i=1,ntime
         write(mytmp,'(i6,2f12.6)') i, tmesh(i), ktau(i)
     enddo ! over i={1,ntime} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_ktau

!!>>> ctqmc_dump_gbin: write out impurity green's function in imaginary
!!>>> time space (generated in binning mode)
  subroutine ctqmc_dump_gbin(ibin, tmesh, gtau)
     use constants, only : dp, mytmp

     use control, only : nband, norbs
     use control, only : ntime

     implicit none

! external arguments
! current bin index, integer representation
     integer, intent(in)  :: ibin

! imaginary time mesh
     real(dp), intent(in) :: tmesh(ntime)

! impurity green's function
     real(dp), intent(in) :: gtau(ntime,norbs,norbs)

! local variables
! loop index
     integer  :: i
     integer  :: j

! scaled impurity green's function
     real(dp) :: gaux(ntime,norbs,norbs)

! current bin index, string representation
     character(len=10) :: sbin

! evaluate gaux first
     call ctqmc_make_gtau(tmesh, gtau, gaux)

! open data file: solver.green.bin.x
     write(sbin,'(i10)') ibin ! convert ibin to sbin
     open(mytmp, file='solver.green.bin.'//trim(adjustl(sbin)), form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,ntime
             write(mytmp,'(2i6,3f12.6)') i, j, tmesh(j), gaux(j,i,i), gaux(j,i+nband,i+nband)
         enddo ! over j={1,ntime} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_gbin

!!========================================================================
!!>>> dump data on matsubara frequency axis                            <<<
!!========================================================================

!!>>> ctqmc_dump_grnf: write out impurity green's function in matsubara
!!>>> frequency space
  subroutine ctqmc_dump_grnf(rmesh, grnf)
     use constants, only : dp, mytmp

     use control, only : nband, norbs
     use control, only : mfreq

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(mfreq)

! impurity green's function
     complex(dp), intent(in) :: grnf(mfreq,norbs,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: solver.grn.dat
     open(mytmp, file='solver.grn.dat', form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,mfreq
             write(mytmp,'(i6,5f16.8)') i, rmesh(j), &
                                  real(grnf(j,i,i)), &
                                 aimag(grnf(j,i,i)), &
                      real(grnf(j,i+nband,i+nband)), &
                     aimag(grnf(j,i+nband,i+nband))
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_grnf

!!>>> ctqmc_dump_wssf: write out bath weiss's function in matsubara
!!>>> frequency space
  subroutine ctqmc_dump_wssf(rmesh, wssf)
     use constants, only : dp, mytmp

     use control, only : nband, norbs
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
     do i=1,nband
         do j=1,mfreq
             write(mytmp,'(i6,5f16.8)') i, rmesh(j), &
                                  real(wssf(j,i,i)), &
                                 aimag(wssf(j,i,i)), &
                      real(wssf(j,i+nband,i+nband)), &
                     aimag(wssf(j,i+nband,i+nband))
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_wssf

!!>>> ctqmc_dump_hybf: write out hybridization function in matsubara
!!>>> frequency space
  subroutine ctqmc_dump_hybf(rmesh, hybf)
     use constants, only : dp, mytmp

     use control, only : nband, norbs
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
     do i=1,nband
         do j=1,mfreq
             write(mytmp,'(i6,5f16.8)') i, rmesh(j), &
                                  real(hybf(j,i,i)), &
                                 aimag(hybf(j,i,i)), &
                      real(hybf(j,i+nband,i+nband)), &
                     aimag(hybf(j,i+nband,i+nband))
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_hybf

!!>>> ctqmc_dump_sigf: write out self-energy function in matsubara
!!>>> frequency space
  subroutine ctqmc_dump_sigf(rmesh, sigf)
     use constants, only : dp, mytmp

     use control, only : nband, norbs
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
     do i=1,nband
         do j=1,mfreq
             write(mytmp,'(i6,5f16.8)') i, rmesh(j), &
                                  real(sigf(j,i,i)), &
                                 aimag(sigf(j,i,i)), &
                      real(sigf(j,i+nband,i+nband)), &
                     aimag(sigf(j,i+nband,i+nband))
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

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
                                    real(ghub(j,i)), &
                                   aimag(ghub(j,i)), &
                                    real(shub(j,i)), &
                                   aimag(shub(j,i))
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
  subroutine ctqmc_dump_hist(hist)
     use constants, only : dp, mytmp

     use control, only : mkink

     implicit none

! external arguments
! histogram data
     real(dp), intent(in) :: hist(mkink)

! local variables
! loop index
     integer  :: i

! scaled histogram data
     real(dp) :: haux(mkink)

! evaluate haux at first
     haux = hist / sum(hist)

! open data file: solver.hist.dat
     open(mytmp, file='solver.hist.dat', form='formatted', status='unknown')

! write it
     write(mytmp,'(a)') '# histogram: order | count | percent'
     do i=1,mkink
         write(mytmp,'(i6,i12,f12.6)') i, int( hist(i) ), haux(i)
     enddo ! over i={1,mkink} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_hist

!!>>> ctqmc_dump_prob: write out the probability of eigenstates of local
!!>>> hamiltonian matrix
  subroutine ctqmc_dump_prob(prob)
     use constants, only : dp, zero, half, mytmp

     use control, only : nband, norbs, ncfgs

     implicit none

! external arguments
! probability data of eigenstates
     real(dp), intent(in) :: prob(ncfgs)

! local variables
! loop index
     integer  :: i
     integer  :: j

! occupation number of eigenstates
     integer  :: noccs(ncfgs)

! net spin of eigenstates
     integer  :: soccs(ncfgs)

! atomic basis sets
     integer  :: basis(ncfgs,norbs)

! probability of occupation number distribution
     real(dp) :: oprob(0:norbs)

! probability of net spin distribution
     real(dp) :: sprob(-nband:nband)

! build atomic basis set, we do not order them according to their
! occupation numbers
     do i=1,ncfgs
         do j=1,norbs
             if ( btest(i-1,j-1) .eqv. .true. ) then
                 basis(i,j) = 1
             else
                 basis(i,j) = 0
             endif ! back if ( btest(i-1,j-1) .eqv. .true. ) block
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,ncfgs} loop

! build occupation numbers for atomic basis set
     do i=1,ncfgs
         noccs(i) = sum( basis(i,:) )
     enddo ! over i={1,ncfgs} loop

! build net spin for eigenstates
     do i=1,ncfgs
         soccs(i) = ( sum( basis(i,1:nband) ) - sum( basis(i,nband+1:norbs) ) )
     enddo ! over i={1,ncfgs} loop

! evaluate oprob
     oprob = zero
     do i=1,ncfgs
         j = noccs(i)
         oprob(j) = oprob(j) + prob(i)
     enddo ! over i={1,ncfgs} loop

! evaluate sprob
     sprob = zero
     do i=1,ncfgs
         j = soccs(i)
         sprob(j) = sprob(j) + prob(i)
     enddo ! over i={1,ncfgs} loop

! open data file: solver.prob.dat
     open(mytmp, file='solver.prob.dat', form='formatted', status='unknown')

! write it
     write(mytmp,'(a)') '# state probability: index | prob | occupy | spin'
     do i=1,ncfgs
         write(mytmp,'(i6,3f12.6)') i, prob(i), real(noccs(i)), real(soccs(i)) * half
     enddo ! over i={1,ncfgs} loop

     write(mytmp,'(a)') '# orbital probability: index | occupy | prob'
     do i=0,norbs
         write(mytmp,'(i6,2f12.6)') i+1, real(i), oprob(i)
     enddo ! over i={0,norbs} loop
     write(mytmp,'(a6,12X,f12.6)') 'sum', sum(oprob)

     write(mytmp,'(a)') '# spin probability: index | spin | prob'
     do i=-nband,nband
         write(mytmp,'(i6,2f12.6)') i+nband+1, i*half, sprob(i)
     enddo ! over i={-nband,nband} loop
     write(mytmp,'(a6,12X,f12.6)') 'sum', sum(sprob)

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_prob

!!>>> ctqmc_dump_nmat: write out the occupation matrix and double
!!>>> occupation matrix
  subroutine ctqmc_dump_nmat(nmat, nnmat)
     use constants, only : dp, mytmp

     use control, only : nband, norbs

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
     write(mytmp,'(a)') '#   < n_i >   data:'
     do i=1,norbs
         write(mytmp,'(i6,f12.6)') i, nmat(i)
     enddo ! over i={1,norbs} loop
     write(mytmp,'(a6,f12.6)') 'sup', sum( nmat(1:nband) )
     write(mytmp,'(a6,f12.6)') 'sdn', sum( nmat(nband+1:norbs) )
     write(mytmp,'(a6,f12.6)') 'sum', sum( nmat(1:norbs) )

     write(mytmp,'(a)') '# < n_i n_j > data:'
     do i=1,norbs
         do j=1,norbs
             write(mytmp,'(2i6,f12.6)') i, j, nnmat(i,j)
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_nmat

!!>>> ctqmc_dump_schi: write out the spin-spin correlation function
  subroutine ctqmc_dump_schi(schi, sschi)
     use constants, only : dp, mytmp

     use control, only : isvrt
     use control, only : nband
     use control, only : ntime
     use context, only : tmesh

     implicit none

! external arguments
! spin-spin correlation function data, < Sz(0) Sz(\tau) >, totally-averaged
     real(dp), intent(in) :: schi(ntime)

! spin-spin correlation function data, < Sz(0) Sz(\tau) >, orbital-resolved
     real(dp), intent(in) :: sschi(ntime,nband)

! local variables
! loop index
     integer :: i
     integer :: j

! check if we need to dump spin-spin correlation function data
     if ( .not. btest(isvrt, 1) ) RETURN

! open data file: solver.schi.dat
     open(mytmp, file='solver.schi.dat', form='formatted', status='unknown')

! write it
     do j=1,nband
         write(mytmp,'(a,i6)') '# flvr:', j
         do i=1,ntime
             write(mytmp,'(2f12.6)') tmesh(i), sschi(i,j)
         enddo ! over i={1,ntime} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over j={1,nband} loop

     write(mytmp,'(a,i6)') '# flvr:', 8888
     do i=1,ntime
         write(mytmp,'(2f12.6)') tmesh(i), schi(i) / real(nband)
     enddo ! over i={1,ntime} loop
     write(mytmp,*) ! write empty lines
     write(mytmp,*)

     write(mytmp,'(a,i6)') '# flvr:', 9999
     do i=1,ntime
         write(mytmp,'(2f12.6)') tmesh(i), sum( sschi(i,:) ) / real(nband)
     enddo ! over i={1,ntime} loop
     write(mytmp,*) ! write empty lines
     write(mytmp,*)

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_schi

!!>>> ctqmc_dump_ochi: write out the orbital-orbital correlation function
  subroutine ctqmc_dump_ochi(ochi, oochi)
     use constants, only : dp, mytmp

     use control, only : isvrt
     use control, only : norbs
     use control, only : ntime
     use context, only : tmesh

     implicit none

! external arguments
! orbital-orbital correlation function data, < N(0) N(\tau) >, totally-averaged
     real(dp), intent(in) :: ochi(ntime)

! orbital-orbital correlation function data, < N(0) N(\tau) >, orbital-resolved
     real(dp), intent(in) :: oochi(ntime,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! check if we need to dump orbital-orbital correlation function data
     if ( isvrt /= 3 ) RETURN

! open data file: solver.ochi.dat
     open(mytmp, file='solver.ochi.dat', form='formatted', status='unknown')

! write it
     do j=1,norbs
         write(mytmp,'(a,i5)') '# flvr:', j
         do i=1,ntime
             write(mytmp,'(2f12.6)') tmesh(i), oochi(i,j)
         enddo ! over i={1,ntime} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over j={1,norbs} loop

     write(mytmp,'(a,i5)') '# flvr:', 8888
     do i=1,ntime
         write(mytmp,'(2f12.6)') tmesh(i), ochi(i) / real(norbs)
     enddo ! over i={1,ntime} loop
     write(mytmp,*) ! write empty lines
     write(mytmp,*)

     write(mytmp,'(a,i5)') '# flvr:', 9999
     do i=1,ntime
         write(mytmp,'(2f12.6)') tmesh(i), sum( oochi(i,:) ) / real(norbs)
     enddo ! over i={1,ntime} loop
     write(mytmp,*) ! write empty lines
     write(mytmp,*)

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_ochi

!!>>> ctqmc_dump_twop: write out the two-particle green's function and
!!>>> vertex function
  subroutine ctqmc_dump_twop(g2_re, g2_im)
     use constants, only : dp, czero, mytmp

     use control, only : isvrt
     use control, only : norbs
     use control, only : nffrq, nbfrq
     use control, only : beta
     use context, only : grnf

     implicit none

! external arguments
! used to calculate two-particle green's function, real part
     real(dp), intent(in) :: g2_re(norbs,norbs,nffrq,nffrq,nbfrq)

! used to calculate two-particle green's function, imaginary part
     real(dp), intent(in) :: g2_im(norbs,norbs,nffrq,nffrq,nbfrq)

! local variables
! loop index for frequencies
     integer :: i
     integer :: j
     integer :: k
     integer :: p
     integer :: q

! loop index for orbitals
     integer :: m
     integer :: n

! dummy integer variables
     integer :: it
     integer :: jt

! dummy complex(dp) variables, used to store the correct green's function
     complex(dp) :: g1
     complex(dp) :: g2
     complex(dp) :: g3
     complex(dp) :: g4

! two-particle green's function, full record
     complex(dp) :: chit

! two-particle green's function, disconnected part
     complex(dp) :: chi0

! two-particle green's function, connected part
     complex(dp) :: chii

! check if we need to dump two-particle green's function and vertex
! function data to solver.twop.dat
     if ( isvrt /= 4 ) RETURN

! open data file: solver.twop.dat
     open(mytmp, file='solver.twop.dat', form='formatted', status='unknown')

! write it
     do m=1,norbs
         do n=1,norbs
             do k=1,nbfrq
                 write(mytmp,'(a,i5)') '# flvr1:', m
                 write(mytmp,'(a,i5)') '# flvr2:', n
                 write(mytmp,'(a,i5)') '# nbfrq:', k
                 do j=1,nffrq

! evaluate g2 and g1
                     if ( j <= nffrq/2 ) then
                         g2 = dconjg( grnf(nffrq/2-j+1,m,m) )
                     else
                         g2 = grnf(j-nffrq/2,m,m)
                     endif ! back if ( j <= nffrq/2 ) block
                     p = j + k - 1
                     if ( p <= nffrq/2 ) then
                         g1 = dconjg( grnf(nffrq/2-p+1,m,m) )
                     else
                         g1 = grnf(p-nffrq/2,m,m)
                     endif ! back if ( p <= nffrq/2 ) block

                     do i=1,nffrq

! evaluate g3 and g4
                         if ( i <= nffrq/2 ) then
                             g3 = dconjg( grnf(nffrq/2-i+1,n,n) )
                         else
                             g3 = grnf(i-nffrq/2,n,n)
                         endif ! back if ( i <= nffrq/2 ) block
                         q = i + k - 1
                         if ( q <= nffrq/2 ) then
                             g4 = dconjg( grnf(nffrq/2-q+1,m,m))
                         else
                             g4 = grnf(q-nffrq/2,m,m)
                         endif ! back if ( q <= nffrq/2 ) block

! evaluate chit
                         chit = dcmplx( g2_re(m,n,j,i,k), g2_im(m,n,j,i,k) )

! evaluate chi0
                         chi0 = czero
                         if ( k == 1 ) chi0 = chi0 + beta * g1 * g3
                         if ( i == j .and. m == n ) chi0 = chi0 - beta * g1 * g3

! evaluate chii, straightforward but less accurate
                         chii = chit - chi0

! jt: \omega
! it: \omega'
! chit: \chi_{tot}(\omega, \omega', \nu)
! chi0: \chi_{0}(\omega, \omega', \nu)
! chii: \chi_{irr}(\omega, \omega', \nu)
! chii/(g1*g2*g3*g4) : \gamma(\omega, \omega', \nu)
                         it = 2*i - nffrq - 1; jt = 2*j - nffrq - 1
                         write(mytmp,'(2i5,8f16.8)') jt, it, chit, chi0, chii, chii/(g1*g2*g3*g4)
                     enddo ! over i={1,nffrq} loop
                 enddo ! over j={1,nffrq} loop
                 write(mytmp,*) ! write empty lines
                 write(mytmp,*)
             enddo ! over k={1,nbfrq} loop
         enddo ! over n={1,norbs} loop
     enddo ! over m={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_twop

!!>>> ctqmc_dump_vrtx: write out the vertex function and two-particle
!!>>> green's function
  subroutine ctqmc_dump_vrtx(h2_re, h2_im)
     use constants, only : dp, czero, mytmp

     use control, only : isvrt
     use control, only : norbs
     use control, only : mfreq
     use control, only : nffrq, nbfrq
     use control, only : beta
     use context, only : g2_re, g2_im
     use context, only : grnf, frnf
     use context, only : sig2

     implicit none

! external arguments
! used to calculate vertex function, real part
     real(dp), intent(in) :: h2_re(norbs,norbs,nffrq,nffrq,nbfrq)

! used to calculate vertex function, imaginary part
     real(dp), intent(in) :: h2_im(norbs,norbs,nffrq,nffrq,nbfrq)

! local variables
! loop index for frequencies
     integer :: i
     integer :: j
     integer :: k
     integer :: p
     integer :: q

! loop index for orbitals
     integer :: m
     integer :: n

! dummy integer variables
     integer :: it
     integer :: jt

! dummy complex(dp) variables, used to store the correct green's function
     complex(dp) :: fw
     complex(dp) :: g1
     complex(dp) :: g2
     complex(dp) :: g3
     complex(dp) :: g4

! two-particle green's function, full record
     complex(dp) :: chit
     complex(dp) :: chih

! two-particle green's function, disconnected part
     complex(dp) :: chi0

! two-particle green's function, connected part
     complex(dp) :: chii

! check if we need to dump two-particle green's function and vertex
! function data to solver.vrtx.dat
     if ( isvrt /= 5 ) RETURN

! build frnf at first: F = G \Sigma
! in principle, F should be measured during the Monte Carlo procedure
     do m=1,norbs
         do k=1,mfreq
             frnf(k,m,m) = grnf(k,m,m) * sig2(k,m,m)
         enddo ! over k={1,mfreq} loop
     enddo ! over m={1,norbs} loop

! open data file: solver.vrtx.dat
     open(mytmp, file='solver.vrtx.dat', form='formatted', status='unknown')

! write it
     do m=1,norbs
         do n=1,norbs
             do k=1,nbfrq
                 write(mytmp,'(a,i5)') '# flvr1:', m
                 write(mytmp,'(a,i5)') '# flvr2:', n
                 write(mytmp,'(a,i5)') '# nbfrq:', k
                 do j=1,nffrq

! evaluate g2 and g1
                     if ( j <= nffrq/2 ) then
                         g2 = dconjg( grnf(nffrq/2-j+1,m,m) )
                     else
                         g2 = grnf(j-nffrq/2,m,m)
                     endif ! back if ( j <= nffrq/2 ) block
                     p = j + k - 1
                     if ( p <= nffrq/2 ) then
                         g1 = dconjg( grnf(nffrq/2-p+1,m,m) )
                     else
                         g1 = grnf(p-nffrq/2,m,m)
                     endif ! back if ( p <= nffrq/2 ) block

! evaluate fw
                     if ( p <= nffrq/2 ) then
                         fw = dconjg( frnf(nffrq/2-p+1,m,m) )
                     else
                         fw = frnf(p-nffrq/2,m,m)
                     endif ! back if ( p <= nffrq/2 ) block

                     do i=1,nffrq

! evaluate g3 and g4
                         if ( i <= nffrq/2 ) then
                             g3 = dconjg( grnf(nffrq/2-i+1,n,n) )
                         else
                             g3 = grnf(i-nffrq/2,n,n)
                         endif ! back if ( i <= nffrq/2 ) block
                         q = i + k - 1
                         if ( q <= nffrq/2 ) then
                             g4 = dconjg( grnf(nffrq/2-q+1,m,m))
                         else
                             g4 = grnf(q-nffrq/2,m,m)
                         endif ! back if ( q <= nffrq/2 ) block

! evaluate chih
                         chih = dcmplx( h2_re(m,n,j,i,k), h2_im(m,n,j,i,k) )

! evaluate chit
                         chit = dcmplx( g2_re(m,n,j,i,k), g2_im(m,n,j,i,k) )

! evaluate chi0
                         chi0 = czero
                         if ( k == 1 ) chi0 = chi0 + beta * g1 * g3
                         if ( i == j .and. m == n ) chi0 = chi0 - beta * g1 * g3

! evaluate chii, more accurate than that in ctqmc_dump_twop() subroutine
                         chii = g1 * chih - fw * chit

! jt: \omega
! it: \omega'
! chit: \chi_{tot}(\omega, \omega', \nu)
! chi0: \chi_{0}(\omega, \omega', \nu)
! chii: \chi_{irr}(\omega, \omega', \nu)
! chii/(g1*g2*g3*g4) : \gamma(\omega, \omega', \nu)
                         it = 2*i - nffrq - 1; jt = 2*j - nffrq - 1
                         write(mytmp,'(2i5,8f16.8)') jt, it, chit, chi0, chii, chii/(g1*g2*g3*g4)
                     enddo ! over i={1,nffrq} loop
                 enddo ! over j={1,nffrq} loop
                 write(mytmp,*) ! write empty lines
                 write(mytmp,*)
             enddo ! over k={1,nbfrq} loop
         enddo ! over n={1,norbs} loop
     enddo ! over m={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_vrtx
