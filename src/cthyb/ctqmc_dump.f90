!!!-----------------------------------------------------------------------
!!! project : iqist @ manjushaka
!!! program : ctqmc_dump_ac_f
!!!           ctqmc_dump_hist
!!!           ctqmc_dump_prob
!!!           ctqmc_dump_paux
!!!           ctqmc_dump_nmat
!!!           ctqmc_dump_gtau
!!!           ctqmc_dump_ftau
!!!           ctqmc_dump_htau
!!!           ctqmc_dump_wtau
!!!           ctqmc_dump_grnf
!!!           ctqmc_dump_frnf
!!!           ctqmc_dump_hybf
!!!           ctqmc_dump_wssf
!!!           ctqmc_dump_sig2
!!!           ctqmc_dump_kmat
!!!           ctqmc_dump_lrmm
!!!           ctqmc_dump_szpw
!!!           ctqmc_dump_sp_t
!!!           ctqmc_dump_sp_w
!!!           ctqmc_dump_ch_t
!!!           ctqmc_dump_ch_w
!!!           ctqmc_dump_g2ph
!!!           ctqmc_dump_g2pp
!!!           ctqmc_dump_diag
!!! source  : ctqmc_dump.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli@caep.cn)
!!! history : 09/16/2009 by li huang (created)
!!!           06/20/2024 by li huang (last modified)
!!! purpose : dump key observables produced by the hybridization expansion
!!!           version continuous time quantum Monte Carlo (CTQMC) quantum
!!!           impurity solver and dynamical mean field theory (DMFT) self
!!!           -consistent engine to external files.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> dump data of autocorrelation time                                <<<
!!========================================================================

!!
!! @sub ctqmc_dump_ac_f
!!
!! write out the autocorrelation function for the total occupation number
!!
  subroutine ctqmc_dump_ac_f(ac_f)
     use constants, only : dp
     use constants, only : mytmp

     use control, only : ntime

     implicit none

!! external arguments
     ! autocorrelation function
     real(dp), intent(in) :: ac_f(ntime)

!! local variables
     ! loop index
     integer :: i

!! [body

     ! open data file: solver.ac_f.dat
     open(mytmp, file='solver.ac_f.dat', form='formatted', status='unknown')

     ! write it
     do i=1,ntime
         write(mytmp,'(i6,f12.6)') i, ac_f(i)
     enddo ! over i={1,ntime} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine ctqmc_dump_ac_f

!!========================================================================
!!>>> dump data of physical observables 1                              <<<
!!========================================================================

!!
!! @sub ctqmc_dump_hist
!!
!! write out the histogram for perturbation expansion series
!!
  subroutine ctqmc_dump_hist(hist, herr)
     use constants, only : dp
     use constants, only : mytmp

     use control, only : mkink

     implicit none

!! external arguments
     ! histogram data and its error bar
     real(dp), intent(in) :: hist(mkink)
     real(dp), intent(in) :: herr(mkink)

!! local variables
     ! loop index
     integer  :: i

     ! scaled histogram data
     real(dp) :: hint(mkink)
     real(dp) :: haux(mkink)
     real(dp) :: htmp(mkink)

!! [body

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

!! body]

     return
  end subroutine ctqmc_dump_hist

!!
!! @sub ctqmc_dump_prob
!!
!! write out the probability of atomic eigenstates of local hamiltonian
!!
  subroutine ctqmc_dump_prob(prob, perr)
     use constants, only : dp, zero, eps6, mytmp

     use control, only : norbs, ncfgs

     use context, only : naux, saux

     use m_sect, only : nsect
     use m_sect, only : sectors

     implicit none

! external arguments
! probability data of eigenstates
     real(dp), intent(in) :: prob(ncfgs)
     real(dp), intent(in) :: perr(ncfgs)

! local variables
! loop index
     integer  :: i
     integer  :: j

! number of individual spin values
     integer  :: ns

! start index of sectors
     integer  :: indx

! dummy arrays, used to store spin of eigenstates
     real(dp) :: stmp1(ncfgs)
     real(dp) :: stmp2(ncfgs)

! probability of sectors
     real(dp) :: psect(nsect)
     real(dp) :: pserr(nsect)

! probability of occupation number distribution
     real(dp) :: oprob(0:norbs)
     real(dp) :: operr(0:norbs) ! error bar

! probability of net spin distribution
     real(dp) :: sprob(ncfgs)
     real(dp) :: sperr(ncfgs) ! error bar

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
     operr = zero
     do i=1,ncfgs
         j = int( naux(i) )
         oprob(j) = oprob(j) + prob(i)
         operr(j) = operr(j) + perr(i)
     enddo ! over i={1,ncfgs} loop

! evaluate sprob
     sprob = zero
     sperr = zero
     do i=1,ncfgs
         SCAN_CYCLE: do j=1,ns
             if ( abs( stmp2(j) - saux(i) ) < eps6 ) then
                 sprob(j) = sprob(j) + prob(i)
                 sperr(j) = sperr(j) + perr(i)
                 EXIT SCAN_CYCLE
             endif ! back if ( abs( stmp2(j) - saux(i) ) < eps6 ) block
         enddo SCAN_CYCLE ! over j={1,ns} loop
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
     write(mytmp,'(a6,12X,f12.6)') 'sum', sum(psect)

     write(mytmp,'(a)') '# orbital probability: index | occupy | prob'
     do i=0,norbs
         write(mytmp,'(i6,3f12.6)') i + 1, real(i), oprob(i), operr(i)
     enddo ! over i={0,norbs} loop
     write(mytmp,'(a6,12X,f12.6)') 'sum', sum(oprob)

     write(mytmp,'(a)') '# spin probability: index | spin | prob'
     do i=1,ns
         write(mytmp,'(i6,3f12.6)') i, stmp2(i), sprob(i), sperr(i)
     enddo ! over i={1,ns} loop
     write(mytmp,'(a6,12X,f12.6)') 'sum', sum(sprob)

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine ctqmc_dump_prob

!!
!! @sub ctqmc_dump_paux
!!
!! write out the auxiliary physical observables
!!
  subroutine ctqmc_dump_paux(paux, perr)
     use constants, only : dp
     use constants, only : mytmp

     implicit none

!! external arguments
     ! auxiliary physical observables and its error bar
     real(dp), intent(in) :: paux(9)
     real(dp), intent(in) :: perr(9)

!! [body

     ! open data file: solver.paux.dat
     open(mytmp, file='solver.paux.dat', form='formatted', status='unknown')

     ! write it
     write(mytmp,'(a,2f12.6)') 'etot:', paux(1), perr(1)
     write(mytmp,'(a,2f12.6)') 'epot:', paux(2), perr(2)
     write(mytmp,'(a,2f12.6)') 'ekin:', paux(3), perr(3)
     write(mytmp,'(a,2f12.6)') '<Sz>:', paux(4), perr(4)
     write(mytmp,'(a,2f12.6)') '<N1>:', paux(5), perr(5)
     write(mytmp,'(a,2f12.6)') '<N2>:', paux(6), perr(6)
     write(mytmp,'(a,2e12.4)') '<K2>:', paux(7), perr(7)
     write(mytmp,'(a,2e12.4)') '<K3>:', paux(8), perr(8)
     write(mytmp,'(a,2e12.4)') '<K4>:', paux(9), perr(9)

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine ctqmc_dump_paux

!!
!! @sub ctqmc_dump_nmat
!!
!! write out the occupation number and double occupation matrix
!!
  subroutine ctqmc_dump_nmat(nimp, nmat, nerr, nbar)
     use constants, only : dp, mytmp

     use control, only : nband, norbs

     implicit none

!! external arguments
     ! occupation number data
     real(dp), intent(in) :: nimp(norbs)
     real(dp), intent(in) :: nerr(norbs)

     ! double occupation matrix data
     real(dp), intent(in) :: nmat(norbs,norbs)
     real(dp), intent(in) :: nbar(norbs,norbs)

!! local variables
     ! loop index
     integer :: i
     integer :: j

!! [body

     ! open data file: solver.nmat.dat
     open(mytmp, file='solver.nmat.dat', form='formatted', status='unknown')

     ! write it
     write(mytmp,'(a)') '#   < n_i >   data:'
     do i=1,norbs
         write(mytmp,'(i6,2f12.6)') i, nimp(i), nerr(i)
     enddo ! over i={1,norbs} loop
     write(mytmp,'(a6,2f12.6)') 'sup', sum( nimp(1:nband) ), sum( nerr(1:nband) )
     write(mytmp,'(a6,2f12.6)') 'sdn', sum( nimp(nband+1:norbs) ), sum( nerr(nband+1:norbs) )
     write(mytmp,'(a6,2f12.6)') 'sum', sum( nimp(1:norbs) ), sum( nerr(1:norbs) )

     write(mytmp,'(a)') '# < n_i n_j > data:'
     do i=1,norbs
         do j=1,norbs
             write(mytmp,'(2i6,2f12.6)') i, j, nmat(i,j), nbar(i,j)
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine ctqmc_dump_nmat

!!========================================================================
!!>>> dump data of physical observables 2: over imaginary time         <<<
!!========================================================================

!!
!! @sub ctqmc_dump_gtau
!!
!! write out impurity green's function in imaginary time space
!!
  subroutine ctqmc_dump_gtau(gtau, gerr)
     use constants, only : dp
     use constants, only : mytmp

     use control, only : norbs
     use control, only : ntime

     use context, only : tmesh

     implicit none

!! external arguments
     ! impurity green's function
     real(dp), intent(in) :: gtau(ntime,norbs,norbs)
     real(dp), intent(in) :: gerr(ntime,norbs,norbs)

!! local variables
     ! loop index
     integer  :: i
     integer  :: j

!! [body

     ! open data file: solver.green.dat
     open(mytmp, file='solver.green.dat', form='formatted', status='unknown')

     ! write it
     do i=1,norbs
         do j=1,ntime
             write(mytmp,'(2i6,3f12.6)') i, j, tmesh(j), gtau(j,i,i), gerr(j,i,i)
         enddo ! over j={1,ntime} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

     ! close data file
     close(mytmp)

!! body]

     return
  end subroutine ctqmc_dump_gtau

!!
!! @sub ctqmc_dump_ftau
!!
!! write out auxiliary correlation function in imaginary time space
!!
  subroutine ctqmc_dump_ftau(ftau, ferr)
     use constants, only : dp, mytmp

     use control, only : norbs
     use control, only : ntime

     use context, only : tmesh

     implicit none

! external arguments
! auxiliary correlation function
     real(dp), intent(in) :: ftau(ntime,norbs,norbs)
     real(dp), intent(in) :: ferr(ntime,norbs,norbs)

! local variables
! loop index
     integer  :: i
     integer  :: j

! scaled auxiliary correlation function
     real(dp) :: faux(ntime,norbs,norbs)
     real(dp) :: fbar(ntime,norbs,norbs)

! evaluate faux and fbar at first
     call ctqmc_make_ftau(tmesh, ftau, faux)
     call ctqmc_make_ftau(tmesh, ferr, fbar)

! open data file: solver.fcorr.dat
     open(mytmp, file='solver.fcorr.dat', form='formatted', status='unknown')

! write it
     do i=1,norbs
         do j=1,ntime
             write(mytmp,'(2i6,3f12.6)') i, j, tmesh(j), faux(j,i,i), fbar(j,i,i)
         enddo ! over j={1,ntime} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_ftau

!!
!! @sub ctqmc_dump_htau
!!
!! write out hybridization function in imaginary time space
!!
  subroutine ctqmc_dump_htau(htau)
     use constants, only : dp, zero, mytmp

     use control, only : norbs
     use control, only : ntime

     use context, only : tmesh

     implicit none

! external arguments
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

!!
!! @sub ctqmc_dump_wtau
!!
!! write out bath weiss's function in imaginary time space
!!
  subroutine ctqmc_dump_wtau(wtau)
     use constants, only : dp, zero, mytmp

     use control, only : norbs
     use control, only : ntime

     use context, only : tmesh

     implicit none

! external arguments
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

!!========================================================================
!!>>> dump data of physical observables 2: over matsubara frequency    <<<
!!========================================================================

!!
!! @sub ctqmc_dump_grnf
!!
!! write out impurity green's function in matsubara frequency space
!!
  subroutine ctqmc_dump_grnf(grnf, gerr)
     use constants, only : dp, mytmp

     use control, only : norbs
     use control, only : mfreq

     use context, only : rmesh

     implicit none

! external arguments
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
             write(mytmp,'(i6,5f16.8)') i, rmesh(j), grnf(j,i,i), gerr(j,i,i)
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_grnf

!!
!! @sub ctqmc_dump_frnf
!!
!! write out auxiliary correlation function in matsubara frequency space
!!
  subroutine ctqmc_dump_frnf(frnf)
     use constants, only : dp, czero, mytmp

     use control, only : norbs
     use control, only : mfreq

     use context, only : rmesh

     implicit none

! external arguments
! auxiliary correlation function
     complex(dp), intent(in) :: frnf(mfreq,norbs,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: solver.frn.dat
     open(mytmp, file='solver.frn.dat', form='formatted', status='unknown')

! write it
     do i=1,norbs
         do j=1,mfreq
             write(mytmp,'(i6,5f16.8)') i, rmesh(j), frnf(j,i,i), czero
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_frnf

!!
!! @sub ctqmc_dump_hybf
!!
!! write out hybridization function in matsubara frequency space
!!
  subroutine ctqmc_dump_hybf(hybf)
     use constants, only : dp, czero, mytmp

     use control, only : norbs
     use control, only : mfreq

     use context, only : rmesh

     implicit none

! external arguments
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
             write(mytmp,'(i6,5f16.8)') i, rmesh(j), hybf(j,i,i), czero
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_hybf

!!
!! @sub ctqmc_dump_wssf
!!
!! write out bath weiss's function in matsubara frequency space
!!
  subroutine ctqmc_dump_wssf(wssf)
     use constants, only : dp, czero, mytmp

     use control, only : norbs
     use control, only : mfreq

     use context, only : rmesh

     implicit none

! external arguments
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
             write(mytmp,'(i6,5f16.8)') i, rmesh(j), wssf(j,i,i), czero
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_wssf

!!
!! @sub ctqmc_dump_sigf
!!
!! write out self-energy function in matsubara frequency space
!!
  subroutine ctqmc_dump_sigf(sigf)
     use constants, only : dp, czero, mytmp

     use control, only : norbs
     use control, only : mfreq

     use context, only : rmesh

     implicit none

! external arguments
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
             write(mytmp,'(i6,5f16.8)') i, rmesh(j), sigf(j,i,i), czero
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_sigf

!!========================================================================
!!>>> dump data of physical observables 3                              <<<
!!========================================================================

!!
!! @sub ctqmc_dump_kmat
!!
!! write out the kinetic energy fluctuation
!!
  subroutine ctqmc_dump_kmat(knop, kmat, kerr, kbar)
     use constants, only : dp, one, two, mytmp

     use control, only : isobs
     use control, only : norbs

     implicit none

! external arguments
! number of operators, < k >
     real(dp), intent(in) :: knop(norbs)
     real(dp), intent(in) :: kerr(norbs)

! crossing product of k_i and k_j, < k_i k_j >
     real(dp), intent(in) :: kmat(norbs,norbs)
     real(dp), intent(in) :: kbar(norbs,norbs)

! local variables
! loop index
     integer  :: i
     integer  :: j

! final value and corresponding error
     real(dp) :: f_val
     real(dp) :: f_err

! calculate f_val and f_err
     f_val = sum( kmat ) - sum( knop ) * ( one * sum( knop ) + one )
     f_err = sum( kbar ) - sum( kerr ) * ( two * sum( knop ) + one )

! check if we need to dump the < k > and < k^2 > data
! to solver.kmat.dat
     if ( .not. btest(isobs, 1) ) RETURN

! open data file: solver.kmat.dat
     open(mytmp, file='solver.kmat.dat', form='formatted', status='unknown')

! write it
     write(mytmp,'(a)') '# <  k  > data:'
     do i=1,norbs
         write(mytmp,'(i6,2f12.6)') i, knop(i), kerr(i)
     enddo ! over i={1,norbs} loop
     write(mytmp,'(a6,2f12.6)') 'k_sum', sum( knop ), sum( kerr )

     write(mytmp,'(a)') '# < k^2 > data:'
     do i=1,norbs
         do j=1,norbs
             write(mytmp,'(2i6,2f12.6)') i, j, kmat(i,j), kbar(i,j)
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop
     write(mytmp,'(a6,2f12.6)') 'kksum', sum( kmat ), sum( kbar )
     write(mytmp,'(a6,2f12.6)') 'final', f_val, f_err

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_kmat

!!
!! @sub ctqmc_dump_lrmm
!!
!! write out the fidelity susceptibility
!!
  subroutine ctqmc_dump_lrmm(lnop, rnop, lrmm, lerr, rerr, lree)
     use constants, only : dp, mytmp

     use control, only : isobs
     use control, only : norbs

     implicit none

! external arguments
! number of operators at left half axis, < k_l >
     real(dp), intent(in) :: lnop(norbs)
     real(dp), intent(in) :: lerr(norbs)

! number of operators at right half axis, < k_r >
     real(dp), intent(in) :: rnop(norbs)
     real(dp), intent(in) :: rerr(norbs)

! crossing product of k_l and k_r, < k_l k_r >
     real(dp), intent(in) :: lrmm(norbs,norbs)
     real(dp), intent(in) :: lree(norbs,norbs)

! local variables
! loop index
     integer  :: i
     integer  :: j

! final value and corresponding error
     real(dp) :: f_val
     real(dp) :: f_err

! calculate f_val and f_err
     f_val = sum( lrmm ) - sum( lnop ) * sum( rnop )
     f_err = sum( lree ) - sum( rnop ) * sum( lerr ) - sum( lnop ) * sum( rerr )

! check if we need to dump the fidelity susceptibility data
! to solver.lrmm.dat
     if ( .not. btest(isobs, 2) ) RETURN

! open data file: solver.lrmm.dat
     open(mytmp, file='solver.lrmm.dat', form='formatted', status='unknown')

! write it
     write(mytmp,'(a)') '# < k_l > < k_r > data:'
     do i=1,norbs
         write(mytmp,'(i6,4f12.6)') i, lnop(i), rnop(i), lerr(i), rerr(i)
     enddo ! over i={1,norbs} loop
     write(mytmp,'(a6,2f12.6)') 'l_sum', sum( lnop ), sum( lerr )
     write(mytmp,'(a6,2f12.6)') 'r_sum', sum( rnop ), sum( rerr )

     write(mytmp,'(a)') '# < k_l k_r > data:'
     do i=1,norbs
         do j=1,norbs
             write(mytmp,'(2i6,2f12.6)') i, j, lrmm(i,j), lree(i,j)
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop
     write(mytmp,'(a6,2f12.6)') 'lrsum', sum( lrmm ), sum( lree )
     write(mytmp,'(a6,2f12.6)') 'fidel', f_val, f_err

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_lrmm

!!
!! @sub ctqmc_dump_szpw
!!
!! write out the powers of local magnetization, which can be used to
!! calculate the binder cumulant
!!
  subroutine ctqmc_dump_szpw(szpw, serr)
     use constants, only : dp, mytmp

     use control, only : isobs
     use control, only : nband, norbs

     implicit none

! external arguments
! powers of local magnetization
     real(dp), intent(in) :: szpw(4,norbs)
     real(dp), intent(in) :: serr(4,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! check if we need to dump the powers of local magnetization data
! to solver.szpw.dat
     if ( .not. btest(isobs, 3) ) RETURN

! open data file: solver.szpw.dat
     open(mytmp, file='solver.szpw.dat', form='formatted', status='unknown')

! write it
     do j=1,nband
         write(mytmp,'(a,i6)') '# flvr:', j
         do i=1,4
             write(mytmp,'(i4,2f12.6)') i, szpw(i,j), serr(i,j)
         enddo ! over i={1,4} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over j={1,nband} loop
     write(mytmp,'(a,i6)') '# flvr:', 8888
     do i=1,4
         write(mytmp,'(i4,2f12.6)') i, szpw(i,nband+1), serr(i,nband+1)
     enddo ! over i={1,4} loop
     write(mytmp,*) ! write empty lines
     write(mytmp,*)

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_szpw

!!========================================================================
!!>>> dump data of physical observables 4                              <<<
!!========================================================================

!!========================================================================
!!>>> dump data of physical observables 5                              <<<
!!========================================================================

!!
!! @sub ctqmc_dump_twop
!!
!! write out the two-particle green's function and full (reducible) vertex
!! function, the improved estimator was used to improve the accuracy
!!
  subroutine ctqmc_dump_twop(g2pw, h2pw, gerr, herr)
     use constants, only : dp, czero, mytmp

     use control, only : isvrt
     use control, only : norbs
     use control, only : nffrq, nbfrq

     use context, only : grnf, frnf

     implicit none

! external arguments
! two-particle green's functions
     complex(dp), intent(in) :: g2pw(nffrq,nffrq,nbfrq,norbs,norbs)
     complex(dp), intent(in) :: gerr(nffrq,nffrq,nbfrq,norbs,norbs)

! irreducible vertex functions
     complex(dp), intent(in) :: h2pw(nffrq,nffrq,nbfrq,norbs,norbs)
     complex(dp), intent(in) :: herr(nffrq,nffrq,nbfrq,norbs,norbs)

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
! jt: \omega, unit is \pi/\beta
! it: \omega', unit is \pi/\beta
     integer :: it
     integer :: jt

! dummy complex(dp) variables
! they are used to store the impurity green's function
     complex(dp) :: fw
     complex(dp) :: g1
     complex(dp) :: g2
     complex(dp) :: g3
     complex(dp) :: g4

! two-particle green's function, connected part, \chi_{irr}
     complex(dp) :: chic

! full vertex function, \gamma
     complex(dp) :: chig

! check whether we need to dump the two-particle green's function and
! irreducible vertex function data to solver.g2pw.dat and solver.h2pw.dat
     if ( .not. btest(isvrt, 1) ) RETURN

! task 1: dump two-particle green's function
!-------------------------------------------------------------------------

! open data file: solver.g2pw.dat
     open(mytmp, file='solver.g2pw.dat', form='formatted', status='unknown')

! write it
     do m=1,norbs
         do n=1,m
             do k=1,nbfrq
                 write(mytmp,'(a,i6)') '# flvr1:', m
                 write(mytmp,'(a,i6)') '# flvr2:', n
                 write(mytmp,'(a,i6)') '# nbfrq:', k
                 do j=1,nffrq
                     do i=1,nffrq
                         it = 2*i - nffrq - 1; jt = 2*j - nffrq - 1
                         write(mytmp,'(2i6,4f16.8)') jt, it, g2pw(i,j,k,n,m), gerr(i,j,k,n,m)
                     enddo ! over i={1,nffrq} loop
                 enddo ! over j={1,nffrq} loop
                 write(mytmp,*) ! write empty lines
                 write(mytmp,*)
             enddo ! over k={1,nbfrq} loop
         enddo ! over n={1,m} loop
     enddo ! over m={1,norbs} loop

! close data file
     close(mytmp)

! task 2: dump irreducible vertex function
!-------------------------------------------------------------------------

! open data file: solver.h2pw.dat
     open(mytmp, file='solver.h2pw.dat', form='formatted', status='unknown')

! write it
     do m=1,norbs
         do n=1,m
             do k=1,nbfrq
                 write(mytmp,'(a,i6)') '# flvr1:', m
                 write(mytmp,'(a,i6)') '# flvr2:', n
                 write(mytmp,'(a,i6)') '# nbfrq:', k
                 do j=1,nffrq
                     do i=1,nffrq
                         it = 2*i - nffrq - 1; jt = 2*j - nffrq - 1
                         write(mytmp,'(2i6,4f16.8)') jt, it, h2pw(i,j,k,n,m), herr(i,j,k,n,m)
                     enddo ! over i={1,nffrq} loop
                 enddo ! over j={1,nffrq} loop
                 write(mytmp,*) ! write empty lines
                 write(mytmp,*)
             enddo ! over k={1,nbfrq} loop
         enddo ! over n={1,m} loop
     enddo ! over m={1,norbs} loop

! close data file
     close(mytmp)

! task 3: dump irreducible vertex function
!-------------------------------------------------------------------------

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

! evaluate g2
                     if ( j <= nffrq/2 ) then
                         g2 = dconjg( grnf(nffrq/2-j+1,m,m) )
                     else
                         g2 = grnf(j-nffrq/2,m,m)
                     endif ! back if ( j <= nffrq/2 ) block

! evaluate g1 and fw
                     p = j + k - 1
                     if ( p <= nffrq/2 ) then
                         g1 = dconjg( grnf(nffrq/2-p+1,m,m) )
                         fw = dconjg( frnf(nffrq/2-p+1,m,m) )
                     else
                         g1 = grnf(p-nffrq/2,m,m)
                         fw = frnf(p-nffrq/2,m,m)
                     endif ! back if ( p <= nffrq/2 ) block

                     do i=1,nffrq

! evaluate g3
                         if ( i <= nffrq/2 ) then
                             g3 = dconjg( grnf(nffrq/2-i+1,n,n) )
                         else
                             g3 = grnf(i-nffrq/2,n,n)
                         endif ! back if ( i <= nffrq/2 ) block

! evaluate g4
                         q = i + k - 1
                         if ( q <= nffrq/2 ) then
                             g4 = dconjg( grnf(nffrq/2-q+1,n,n))
                         else
                             g4 = grnf(q-nffrq/2,n,n)
                         endif ! back if ( q <= nffrq/2 ) block

! evaluate chic
                         chic = g1 * h2pw(i,j,k,n,m) - fw * g2pw(i,j,k,n,m)

! evaluate chig
                         chig = chic / (g1 * g2 * g3 * g4)

                         it = 2*i - nffrq - 1; jt = 2*j - nffrq - 1
                         write(mytmp,'(2i6,4f16.8)') jt, it, chic, chig
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

!!
!! @sub ctqmc_dump_pair
!!
!! write out the particle-particle pairing susceptibility
!!
  subroutine ctqmc_dump_pair(p2pw, perr)
     use constants, only : dp, mytmp

     use control, only : isvrt
     use control, only : norbs
     use control, only : nffrq, nbfrq

     implicit none

! external arguments
! particle-particle pairing susceptibility
     complex(dp), intent(in) :: p2pw(nffrq,nffrq,nbfrq,norbs,norbs)
     complex(dp), intent(in) :: perr(nffrq,nffrq,nbfrq,norbs,norbs)

! local variables
! loop index for frequencies
     integer :: i
     integer :: j
     integer :: k

! loop index for orbitals
     integer :: m
     integer :: n

! dummy integer variables
! jt: \omega, unit is \pi/\beta
! it: \omega', unit is \pi/\beta
     integer :: it
     integer :: jt

! check if we need to dump the particle-particle pairing susceptibility
! to solver.pair.dat
     if ( .not. btest(isvrt, 2) ) RETURN

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
                         it = 2*i - nffrq - 1; jt = 2*j - nffrq - 1
                         write(mytmp,'(2i6,4f16.8)') jt, it, p2pw(i,j,k,n,m), perr(i,j,k,n,m)
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

!!========================================================================
!!>>> dump data of diagrammatic configuration                          <<<
!!========================================================================

!!
!! @sub ctqmc_dump_diag
!!
!! write out a snapshot for the current diagram configuration, the results
!! can be used to make a dynamical video.
!!
  subroutine ctqmc_dump_diag(iter, cstep)
     use constants, only : mystd, mytmp

     use control, only : norbs
     use control, only : niter
     use control, only : nwrite, nsweep

     use context, only : index_s, index_e
     use context, only : time_s, time_e
     use context, only : rank

     implicit none

! external arguments
! current self-consistent iteration number
     integer, intent(in) :: iter

! current QMC sweeping steps
     integer, intent(in) :: cstep

! local variables
! loop index for the flavor
     integer :: i

! loop index for the operator pair
     integer :: j

! setup the internal criterion
     if ( nsweep / nwrite < 100 ) RETURN

! write the snapshot
! open data file: solver.diag.dat
     open(mytmp, file='solver.diag.dat', form='formatted', status='unknown', position='append')

! write diagram info
     write(mytmp,'(2(a,i4))') '>> cur_iter:', iter, ' tot_iter:', niter
     write(mytmp,'(2(a,i4))') '>> cur_diag:', cstep/nwrite, ' tot_diag:', nsweep/nwrite

! write the position of operators
     do i=1,norbs
         write(mytmp,'(2(a,i4))') '# flvr:', i, ' rank:', rank(i)
         do j=1,rank(i)
             write(mytmp,'(i4,2f16.8)') i, time_s( index_s(j, i), i ), time_e( index_e(j, i), i )
         enddo ! over j={1,rank(i)} loop
     enddo ! over i={1,norbs} loop

! write two blank lines
     write(mytmp,*)
     write(mytmp,*)

! close data file
     close(mytmp)

! write the message to the terminal
     write(mystd,'(4X,a)') '>>> quantum impurity solver config: saving'

     return
  end subroutine ctqmc_dump_diag
