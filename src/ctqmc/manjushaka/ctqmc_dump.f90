!!!-----------------------------------------------------------------------
!!! project : manjushaka
!!! program : ctqmc_dump_hist
!!!           ctqmc_dump_prob
!!!           ctqmc_dump_paux
!!!           ctqmc_dump_nmat <<<---
!!!           ctqmc_dump_gtau
!!!           ctqmc_dump_ftau
!!!           ctqmc_dump_htau
!!!           ctqmc_dump_wtau <<<---
!!!           ctqmc_dump_grnf
!!!           ctqmc_dump_frnf
!!!           ctqmc_dump_hybf
!!!           ctqmc_dump_wssf
!!!           ctqmc_dump_sigf <<<---
!!!           ctqmc_dump_kmat
!!!           ctqmc_dump_lrmm
!!!           ctqmc_dump_szpw <<<---
!!!           ctqmc_dump_twop
!!!           ctqmc_dump_pair <<<---
!!!           ctqmc_dump_diag <<<---
!!! source  : ctqmc_dump.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           05/16/2017 by li huang (last modified)
!!! purpose : dump key observables produced by the hybridization expansion
!!!           version continuous time quantum Monte Carlo (CTQMC) quantum
!!!           impurity solver and dynamical mean field theory (DMFT) self
!!!           -consistent engine to external files.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> dump data of physical observables 1                              <<<
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
     real(dp) :: hint(mkink)
     real(dp) :: haux(mkink)
     real(dp) :: htmp(mkink)

! evaluate hint, haux, and htmp at first, and then transform them
     hint = hist
     haux = hist / sum(hist)
     htmp = herr / sum(hist)

     hint = cshift(hint, -1)
     haux = cshift(haux, -1)
     htmp = cshift(htmp, -1)

! open data file: solver.hist.dat
     open(mytmp, file='solver.hist.dat', form='formatted', status='unknown')

! write it
     write(mytmp,'(a)') '# histogram: order | count | percent'
     do i=1,mkink
         write(mytmp,'(i6,i12,2f12.6)') i-1, int( hint(i) ), haux(i), htmp(i)
     enddo ! over i={1,mkink} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_hist

!!>>> ctqmc_dump_prob: write out the probability of eigenstates of local
!!>>> hamiltonian matrix
  subroutine ctqmc_dump_prob(prob, naux, saux, perr)
     use constants, only : dp, zero, eps6, mytmp

     use control, only : norbs, ncfgs

     use m_sect, only : nsect
     use m_sect, only : sectors

     implicit none

! external arguments
! probability data of eigenstates
     real(dp), intent(in) :: prob(ncfgs)
     real(dp), intent(in) :: perr(ncfgs)

! occupation number of eigenstates
     real(dp), intent(in) :: naux(ncfgs)

! net spin of eigenstates
     real(dp), intent(in) :: saux(ncfgs)

! local variables
! loop index
     integer  :: i
     integer  :: j

! number of individual spin values
     integer  :: ns

! start index of sectors
     integer  :: indx

! probability of occupation number distribution
     real(dp) :: oprob(0:norbs)

! dummy arrays, used to store spin of eigenstates
     real(dp) :: stmp1(ncfgs)
     real(dp) :: stmp2(ncfgs)

! probability of net spin distribution
     real(dp) :: sprob(ncfgs)

! probability of sectors
     real(dp) :: psect(nsect)
     real(dp) :: pserr(nsect)

! evaluate psect
     psect = zero
     pserr = zero
     do i=1,nsect
         indx = sectors(i)%istart
         do j=1,sectors(i)%ndim
             psect(i) = psect(i) + prob(indx+j-1)
             pserr(i) = pserr(i) + perr(indx+j-1)
         enddo ! over j={1,sectors(i)%ndim} loop
     enddo ! over i={1,nsect} loop

! evaluate oprob
     oprob = zero
     do i=1,ncfgs
         j = int( naux(i) )
         oprob(j) = oprob(j) + prob(i)
     enddo ! over i={1,ncfgs} loop

! sort all the spin values
     stmp1 = saux
     call s_sorter(ncfgs, stmp1)

! find out individual spin values, and store them into stmp2
     ns = 1
     stmp2 = zero
     stmp2(1) = stmp1(1)
     do i=2,ncfgs
         if ( stmp2(ns) < stmp1(i) ) then
             ns = ns + 1
             stmp2(ns) = stmp1(i)
         endif ! back if ( stmp2(ns) < stmp1(i) ) block
     enddo ! over i={2,ncfgs} loop

! evaluate sprob
     sprob = zero
     do i=1,ncfgs
         do j=1,ns
             if ( abs( stmp2(j) - saux(i) ) < eps6 ) then
                 sprob(j) = sprob(j) + prob(i); EXIT
             endif ! back if ( abs( stmp2(j) - saux(i) ) < eps6 ) block
         enddo ! over j={1,ns} loop
     enddo ! over i={1,ncfgs} loop

! open data file: solver.prob.dat
     open(mytmp, file='solver.prob.dat', form='formatted', status='unknown')

! write it
     write(mytmp,'(a)') '# state probability: index | prob | occupy | spin'
     do i=1,ncfgs
         write(mytmp,'(i6,4f12.6)') i, prob(i), naux(i), saux(i), perr(i)
     enddo ! over i={1,ncfgs} loop

     write(mytmp,'(a)') '# sector probability: index | occupy | prob'
     do i=1,nsect
         write(mytmp,'(i6,3f12.6)') i, real( sectors(i)%nele ), psect(i), pserr(i)
     enddo ! over i={1,nsect} loop
     write(mytmp,'(a6,12X,2f12.6)') 'sum', sum(psect), sum(pserr)

     write(mytmp,'(a)') '# orbital probability: index | occupy | prob'
     do i=0,norbs
         write(mytmp,'(i6,2f12.6)') i + 1, real(i), oprob(i)
     enddo ! over i={0,norbs} loop
     write(mytmp,'(a6,12X,f12.6)') 'sum', sum(oprob)

     write(mytmp,'(a)') '# spin probability: index | spin | prob'
     do i=1,ns
         write(mytmp,'(i6,2f12.6)') i, stmp2(i), sprob(i)
     enddo ! over i={1,ns} loop
     write(mytmp,'(a6,12X,f12.6)') 'sum', sum(sprob)

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_prob


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


!!>>> ctqmc_dump_nmat: write out the occupation matrix and double
!!>>> occupation matrix
  subroutine ctqmc_dump_nmat(nmat, nnmat, nerr, nnerr)
     use constants, only : dp, mytmp

     use control, only : nband, norbs

     implicit none

! external arguments
! occupation matrix data
     real(dp), intent(in) :: nmat(norbs)
     real(dp), intent(in) :: nerr(norbs)

! double occupation matrix data
     real(dp), intent(in) :: nnmat(norbs,norbs)
     real(dp), intent(in) :: nnerr(norbs,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: solver.nmat.dat
     open(mytmp, file='solver.nmat.dat', form='formatted', status='unknown')

! write it
     write(mytmp,'(a)') '#   < n_i >   data:'
     do i=1,norbs
         write(mytmp,'(i6,2f12.6)') i, nmat(i), nerr(i)
     enddo ! over i={1,norbs} loop
     write(mytmp,'(a6,2f12.6)') 'sup', sum( nmat(1:nband) ), sum( nerr(1:nband) )
     write(mytmp,'(a6,2f12.6)') 'sdn', sum( nmat(nband+1:norbs) ), sum( nerr(nband+1:norbs) )
     write(mytmp,'(a6,2f12.6)') 'sum', sum( nmat(1:norbs) ), sum( nerr(1:norbs) )

     write(mytmp,'(a)') '# < n_i n_j > data:'
     do i=1,norbs
         do j=1,norbs
             write(mytmp,'(2i6,2f12.6)') i, j, nnmat(i,j), nnerr(i,j)
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_nmat

!!>>> ctqmc_dump_kmat: write out the < k > and < k^2 >
  subroutine ctqmc_dump_kmat(kmat, kkmat, kerr, kkerr)
     use constants, only : dp, one, two, mytmp

     use control, only : issus
     use control, only : norbs

     implicit none

! external arguments
! number of operators, < k >
     real(dp), intent(in) :: kmat(norbs)
     real(dp), intent(in) :: kerr(norbs)

! square of number of operators, < k^2 >
     real(dp), intent(in) :: kkmat(norbs,norbs)
     real(dp), intent(in) :: kkerr(norbs,norbs)

! local variables
! loop index
     integer  :: i
     integer  :: j

! final value and corresponding error
     real(dp) :: f_val
     real(dp) :: f_err

! calculate f_val and f_err
     f_val = sum( kkmat ) - sum( kmat ) * ( one * sum( kmat ) + one )
     f_err = sum( kkerr ) - sum( kerr ) * ( two * sum( kmat ) + one )

! check if we need to dump the < k > and < k^2 > data
! to solver.kmat.dat
     if ( .not. btest(issus, 5) ) RETURN

! open data file: solver.kmat.dat
     open(mytmp, file='solver.kmat.dat', form='formatted', status='unknown')

! write it
     write(mytmp,'(a)') '# <  k  > data:'
     do i=1,norbs
         write(mytmp,'(i6,2f12.6)') i, kmat(i), kerr(i)
     enddo ! over i={1,norbs} loop
     write(mytmp,'(a6,2f12.6)') 'k_sum', sum( kmat ), sum( kerr )

     write(mytmp,'(a)') '# < k^2 > data:'
     do i=1,norbs
         do j=1,norbs
             write(mytmp,'(2i6,2f12.6)') i, j, kkmat(i,j), kkerr(i,j)
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop
     write(mytmp,'(a6,2f12.6)') 'kksum', sum( kkmat ), sum( kkerr )
     write(mytmp,'(a6,2f12.6)') 'final', f_val, f_err

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_kmat

!!>>> ctqmc_dump_lmat: write out the fidelity susceptibility
  subroutine ctqmc_dump_lmat(lmat, rmat, lrmat, lerr, rerr, lrerr)
     use constants, only : dp, mytmp

     use control, only : issus
     use control, only : norbs

     implicit none

! external arguments
! number of operators at left half axis, < k_l >
     real(dp), intent(in) :: lmat(norbs)
     real(dp), intent(in) :: lerr(norbs)

! number of operators at right half axis, < k_r >
     real(dp), intent(in) :: rmat(norbs)
     real(dp), intent(in) :: rerr(norbs)

! used to evaluate fidelity susceptibility, < k_l k_r >
     real(dp), intent(in) :: lrmat(norbs,norbs)
     real(dp), intent(in) :: lrerr(norbs,norbs)

! local variables
! loop index
     integer  :: i
     integer  :: j

! final value and corresponding error
     real(dp) :: f_val
     real(dp) :: f_err

! calculate f_val and f_err
     f_val = sum( lrmat ) - sum( lmat ) * sum( rmat )
     f_err = sum( lrerr ) - sum( rmat ) * sum( lerr ) - sum( lmat ) * sum( rerr )

! check if we need to dump the fidelity susceptibility data
! to solver.lmat.dat
     if ( .not. btest(issus, 6) ) RETURN

! open data file: solver.lmat.dat
     open(mytmp, file='solver.lmat.dat', form='formatted', status='unknown')

! write it
     write(mytmp,'(a)') '# < k_l > < k_r > data:'
     do i=1,norbs
         write(mytmp,'(i6,4f12.6)') i, lmat(i), rmat(i), lerr(i), rerr(i)
     enddo ! over i={1,norbs} loop
     write(mytmp,'(a6,2f12.6)') 'l_sum', sum( lmat ), sum( lerr )
     write(mytmp,'(a6,2f12.6)') 'r_sum', sum( rmat ), sum( rerr )

     write(mytmp,'(a)') '# < k_l k_r > data:'
     do i=1,norbs
         do j=1,norbs
             write(mytmp,'(2i6,2f12.6)') i, j, lrmat(i,j), lrerr(i,j)
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop
     write(mytmp,'(a6,2f12.6)') 'lrsum', sum( lrmat ), sum( lrerr )
     write(mytmp,'(a6,2f12.6)') 'fidel', f_val, f_err

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_lmat

!!>>> ctqmc_dump_twop: write out the two-particle green's function and
!!>>> full (reducible) vertex function
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
     real(dp), intent(in) :: g2_re(nffrq,nffrq,nbfrq,norbs,norbs)

! used to calculate two-particle green's function, imaginary part
     real(dp), intent(in) :: g2_im(nffrq,nffrq,nbfrq,norbs,norbs)

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

! check if we need to dump the two-particle green's function and vertex
! function data to solver.twop.dat
     if ( .not. btest(isvrt, 1) ) RETURN

! open data file: solver.twop.dat
     open(mytmp, file='solver.twop.dat', form='formatted', status='unknown')

! write it
     do m=1,norbs
         do n=1,m
             do k=1,nbfrq
                 write(mytmp,'(a,i6)') '# flvr1:', m
                 write(mytmp,'(a,i6)') '# flvr2:', n
                 write(mytmp,'(a,i6)') '# nbfrq:', k
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
                             g4 = dconjg( grnf(nffrq/2-q+1,n,n))
                         else
                             g4 = grnf(q-nffrq/2,n,n)
                         endif ! back if ( q <= nffrq/2 ) block

! evaluate chit
                         chit = dcmplx( g2_re(i,j,k,n,m), g2_im(i,j,k,n,m) )

! evaluate chi0
                         chi0 = czero
                         if ( k == 1 ) chi0 = chi0 + beta * g1 * g3
                         if ( i == j .and. m == n ) chi0 = chi0 - beta * g1 * g3

! evaluate chii, straightforward but less accurate
                         chii = chit - chi0

! jt: \omega, unit is \pi/\beta
! it: \omega', unit is \pi/\beta
! chit: \chi_{tot}(\omega, \omega', \nu), two-particle green's function
! chi0: \chi_{0}(\omega, \omega', \nu), bubble function
! chii: \chi_{irr}(\omega, \omega', \nu)
! chii/(g1*g2*g3*g4) : \gamma(\omega, \omega', \nu), full vertex function
                         it = 2*i - nffrq - 1; jt = 2*j - nffrq - 1
                         write(mytmp,'(2i6,8f16.8)') jt, it, chit, chi0, chii, chii/(g1*g2*g3*g4)
                     enddo ! over i={1,nffrq} loop
                 enddo ! over j={1,nffrq} loop
                 write(mytmp,*) ! write empty lines
                 write(mytmp,*)
             enddo ! over k={1,nbfrq} loop
         enddo ! over n={1,m} loop
     enddo ! over m={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_twop

!!>>> ctqmc_dump_pair: write out the particle-particle pair susceptibility
  subroutine ctqmc_dump_pair(ps_re, ps_im)
     use constants, only : dp, mytmp

     use control, only : isvrt
     use control, only : norbs
     use control, only : nffrq, nbfrq

     implicit none

! external arguments
! particle-particle pair susceptibility, real part
     real(dp), intent(in) :: ps_re(nffrq,nffrq,nbfrq,norbs,norbs)

! particle-particle pair susceptibility, imaginary part
     real(dp), intent(in) :: ps_im(nffrq,nffrq,nbfrq,norbs,norbs)

! local variables
! loop index for frequencies
     integer :: i
     integer :: j
     integer :: k

! loop index for orbitals
     integer :: m
     integer :: n

! dummy integer variables
     integer :: it
     integer :: jt

! check if we need to dump the particle-particle pair susceptibility
! to solver.pair.dat
     if ( .not. btest(isvrt, 3) ) RETURN

! open data file: solver.pair.dat
     open(mytmp, file='solver.pair.dat', form='formatted', status='unknown')

! write it
     do m=1,norbs
         do n=1,m
             do k=1,nbfrq
                 write(mytmp,'(a,i6)') '# flvr1:', m
                 write(mytmp,'(a,i6)') '# flvr2:', n
                 write(mytmp,'(a,i6)') '# nbfrq:', k
                 do j=1,nffrq
                     do i=1,nffrq
! jt: \omega, unit is \pi/\beta
! it: \omega', unit is \pi/\beta
                         it = 2*i - nffrq - 1; jt = 2*j - nffrq - 1
                         write(mytmp,'(2i6,2f16.8)') jt, it, ps_re(i,j,k,n,m), ps_im(i,j,k,n,m)
                     enddo ! over i={1,nffrq} loop
                 enddo ! over j={1,nffrq} loop
                 write(mytmp,*) ! write empty lines
                 write(mytmp,*)
             enddo ! over k={1,nbfrq} loop
         enddo ! over n={1,m} loop
     enddo ! over m={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_pair
