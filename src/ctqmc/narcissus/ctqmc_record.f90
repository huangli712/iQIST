!!!-----------------------------------------------------------------------
!!! project : narcissus
!!! program : ctqmc_record_gtau
!!!           ctqmc_record_ftau
!!!           ctqmc_record_grnf
!!!           ctqmc_record_hist
!!!           ctqmc_record_prob
!!!           ctqmc_record_nmat
!!!           ctqmc_record_kmat
!!!           ctqmc_record_lmat
!!!           ctqmc_record_szpw
!!!           ctqmc_record_schi
!!!           ctqmc_record_sfom
!!!           ctqmc_record_ochi
!!!           ctqmc_record_ofom
!!!           ctqmc_record_twop
!!!           ctqmc_record_vrtx
!!!           ctqmc_record_pair <<<---
!!!           ctqmc_reduce_gtau
!!!           ctqmc_reduce_ftau
!!!           ctqmc_reduce_grnf
!!!           ctqmc_reduce_hist
!!!           ctqmc_reduce_prob
!!!           ctqmc_reduce_nmat
!!!           ctqmc_reduce_kmat
!!!           ctqmc_reduce_lmat
!!!           ctqmc_reduce_szpw
!!!           ctqmc_reduce_schi
!!!           ctqmc_reduce_sfom
!!!           ctqmc_reduce_ochi
!!!           ctqmc_reduce_ofom
!!!           ctqmc_reduce_twop
!!!           ctqmc_reduce_vrtx
!!!           ctqmc_reduce_pair <<<---
!!!           ctqmc_symm_nmat
!!!           ctqmc_symm_gtau
!!!           ctqmc_symm_grnf
!!!           ctqmc_smth_sigf   <<<---
!!!           ctqmc_make_gtau
!!!           ctqmc_make_ftau   <<<---
!!!           ctqmc_make_iret
!!!           ctqmc_make_pref   <<<---
!!!           ctqmc_make_prod   <<<---
!!!           ctqmc_make_hub1
!!!           ctqmc_make_hub2   <<<---
!!! source  : ctqmc_record.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           08/17/2015 by li huang (last modified)
!!! purpose : measure, record, and postprocess the important observables
!!!           produced by the hybridization expansion version continuous
!!!           time quantum Monte Carlo (CTQMC) quantum impurity solver
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> measure physical observables                                     <<<
!!========================================================================

!!>>> ctqmc_record_gtau: record the impurity green's function in imaginary
!!>>> time axis
  subroutine ctqmc_record_gtau()
     use constants, only : dp, zero, one, two, pi

     use control, only : isort
     use control, only : norbs
     use control, only : lemax, legrd, chmax, chgrd
     use control, only : ntime
     use control, only : beta
     use context, only : index_s, index_e, time_s, time_e
     use context, only : ppleg, qqche
     use context, only : rank
     use context, only : mmat
     use context, only : gtau

     implicit none

! local variables
! loop indices for start and end points
     integer  :: is
     integer  :: ie

! loop index for flavor channel
     integer  :: flvr

! loop index for legendre polynomial
     integer  :: fleg

! loop index for chebyshev polynomial
     integer  :: fche

! index for imaginary time \tau
     integer  :: curr

! used to store the element of mmat matrix
     real(dp) :: maux

! imaginary time for start and end points
     real(dp) :: taus
     real(dp) :: taue

! length betweem taus and taue
     real(dp) :: dtau
     real(dp) :: daux

! interval for imaginary time slice
     real(dp) :: step

! select measurement method
     select case ( isort )

         case (1, 4)
             call cat_record_gtau1()

         case (2, 5)
             call cat_record_gtau2()

         case (3, 6)
             call cat_record_gtau3()

     end select

     return

  contains

!!>>> cat_record_gtau1: record impurity green's function using normal
!!>>> representation
  subroutine cat_record_gtau1()
     implicit none

! evaluate step at first
     step = real(ntime - 1) / beta

     CTQMC_FLAVOR_LOOP: do flvr=1,norbs

! get imaginary time value for segments
         do is=1,rank(flvr)
             taus = time_s( index_s(is, flvr), flvr )

             do ie=1,rank(flvr)
                 taue = time_e( index_e(ie, flvr), flvr )

! evaluate dtau
                 dtau = taue - taus

! get matrix element from mmat, pay special attention to the sign of dtau
                 maux = mmat(ie, is, flvr) * sign(one, dtau)

! adjust dtau, keep it stay in (zero, beta)
                 if ( dtau < zero ) then
                     dtau = dtau + beta
                 endif ! back if ( dtau < zero ) block

! determine index for imaginary time
                 curr = nint( dtau * step ) + 1

! special tricks for the first point and the last point
                 if ( curr == 1 .or. curr == ntime ) then
                     maux = two * maux
                 endif ! back if ( curr == 1 .or. curr == ntime ) block

! record gtau, we normalize gtau in ctqmc_make_gtau() subroutine
                 gtau(curr, flvr, flvr) = gtau(curr, flvr, flvr) - maux

             enddo ! over ie={1,rank(flvr)} loop
         enddo ! over is={1,rank(flvr)} loop

     enddo CTQMC_FLAVOR_LOOP ! over flvr={1,norbs} loop

     return
  end subroutine cat_record_gtau1

!!>>> cat_record_gtau2: record impurity green's function using legendre
!!>>> polynomial representation
  subroutine cat_record_gtau2()
     implicit none

! evaluate step at first
     step = real(legrd - 1) / two

     CTQMC_FLAVOR_LOOP: do flvr=1,norbs

! get imaginary time value for segments
         do is=1,rank(flvr)
             taus = time_s( index_s(is, flvr), flvr )

             do ie=1,rank(flvr)
                 taue = time_e( index_e(ie, flvr), flvr )

! evaluate dtau
                 dtau = taue - taus

! get matrix element from mmat, pay special attention to the sign of dtau
                 maux = mmat(ie, is, flvr) * sign(one, dtau)

! adjust dtau, keep it stay in (zero, beta)
                 if ( dtau < zero ) then
                     dtau = dtau + beta
                 endif ! back if ( dtau < zero ) block

! convert dtau in [0,\beta] to daux in [0,2]
                 daux = two * dtau / beta

! determine index for legendre polynomial interval
                 curr = nint( daux * step ) + 1

! record gtau, we normalize gtau in ctqmc_make_gtau() subroutine
                 CTQMC_FLALEG_LOOP: do fleg=1,lemax
                     dtau = sqrt(two * fleg - 1) * ppleg(curr,fleg)
                     gtau(fleg, flvr, flvr) = gtau(fleg, flvr, flvr) - maux * dtau
                 enddo CTQMC_FLALEG_LOOP ! over fleg={1,lemax} loop

             enddo ! over ie={1,rank(flvr)} loop
         enddo ! over is={1,rank(flvr)} loop

     enddo CTQMC_FLAVOR_LOOP ! over flvr={1,norbs} loop

     return
  end subroutine cat_record_gtau2

!!>>> cat_record_gtau3: record impurity green's function using chebyshev
!!>>> polynomial representation
  subroutine cat_record_gtau3()
     implicit none

! evaluate step at first
     step = real(chgrd - 1) / two

     CTQMC_FLAVOR_LOOP: do flvr=1,norbs

! get imaginary time value for segments
         do is=1,rank(flvr)
             taus = time_s( index_s(is, flvr), flvr )

             do ie=1,rank(flvr)
                 taue = time_e( index_e(ie, flvr), flvr )

! evaluate dtau
                 dtau = taue - taus

! get matrix element from mmat, pay special attention to the sign of dtau
                 maux = mmat(ie, is, flvr) * sign(one, dtau)

! adjust dtau, keep it stay in (zero, beta)
                 if ( dtau < zero ) then
                     dtau = dtau + beta
                 endif ! back if ( dtau < zero ) block

! convert dtau in [0,\beta] to daux in [0,2]
                 daux = two * dtau / beta

! determine index for chebyshev polynomial interval
                 curr = nint( daux * step ) + 1

! record gtau, we normalize gtau in ctqmc_make_gtau() subroutine
                 CTQMC_FLACHE_LOOP: do fche=1,chmax
                     dtau = (two / pi) * sqrt(one - (daux - one)**2) * qqche(curr,fche)
                     gtau(fche, flvr, flvr) = gtau(fche, flvr, flvr) - maux * dtau
                 enddo CTQMC_FLACHE_LOOP ! over fche={1,chmax} loop

             enddo ! over ie={1,rank(flvr)} loop
         enddo ! over is={1,rank(flvr)} loop

     enddo CTQMC_FLAVOR_LOOP ! over flvr={1,norbs} loop

     return
  end subroutine cat_record_gtau3
  end subroutine ctqmc_record_gtau

!!>>> ctqmc_record_ftau: record the auxiliary correlation function in
!!>>> imaginary time axis, F(\tau)
  subroutine ctqmc_record_ftau()
     use constants, only : dp, zero, one, two, pi

     use control, only : isort
     use control, only : norbs
     use control, only : lemax, legrd, chmax, chgrd
     use control, only : ntime
     use control, only : beta
     use context, only : index_s, index_e, time_s, time_e
     use context, only : ppleg, qqche
     use context, only : rank, pref
     use context, only : mmat
     use context, only : ftau

     implicit none

! local variables
! loop indices for start and end points
     integer  :: is
     integer  :: ie

! loop index for flavor channel
     integer  :: flvr

! loop index for legendre polynomial
     integer  :: fleg

! loop index for chebyshev polynomial
     integer  :: fche

! index for imaginary time \tau
     integer  :: curr

! used to store the element of mmat matrix
     real(dp) :: maux

! imaginary time for start and end points
     real(dp) :: taus
     real(dp) :: taue

! length betweem taus and taue
     real(dp) :: dtau
     real(dp) :: daux

! interval for imaginary time slice
     real(dp) :: step

! select measurement method
     select case ( isort )

         case (4)
             call cat_record_ftau1()

         case (5)
             call cat_record_ftau2()

         case (6)
             call cat_record_ftau3()

     end select

     return

  contains

!!>>> cat_record_ftau1: record auxiliary correlation function using normal
!!>>> representation
  subroutine cat_record_ftau1()
     implicit none

! calculate prefactor: pref
     call ctqmc_make_pref()

! evaluate step at first
     step = real(ntime - 1) / beta

     CTQMC_FLAVOR_LOOP: do flvr=1,norbs

! get imaginary time value for segments
         do is=1,rank(flvr)
             taus = time_s( index_s(is, flvr), flvr )

! get imaginary time value for segments
             do ie=1,rank(flvr)
                 taue = time_e( index_e(ie, flvr), flvr )

! evaluate dtau
                 dtau = taue - taus

! get matrix element from mmat, pay special attention to the sign of dtau
                 maux = mmat(ie, is, flvr) * sign(one, dtau) * pref(ie,flvr)

! adjust dtau, keep it stay in (zero, beta)
                 if ( dtau < zero ) then
                     dtau = dtau + beta
                 endif ! back if ( dtau < zero ) block

! determine index for imaginary time
                 curr = nint( dtau * step ) + 1

! special tricks for the first point and the last point
                 if ( curr == 1 .or. curr == ntime ) then
                     maux = two * maux
                 endif ! back if ( curr == 1 .or. curr == ntime ) block

! record ftau, we normalize ftau in ctqmc_make_ftau() subroutine
                 ftau(curr, flvr, flvr) = ftau(curr, flvr, flvr) - maux

             enddo ! over ie={1,rank(flvr)} loop
         enddo ! over is={1,rank(flvr)} loop

     enddo CTQMC_FLAVOR_LOOP ! over flvr={1,norbs} loop

     return
  end subroutine cat_record_ftau1

!!>>> cat_record_ftau2: record auxiliary correlation function using
!!>>> legendre polynomial representation
  subroutine cat_record_ftau2()
     implicit none

! calculate prefactor: pref
     call ctqmc_make_pref()

! evaluate step at first
     step = real(legrd - 1) / two

     CTQMC_FLAVOR_LOOP: do flvr=1,norbs

! get imaginary time value for segments
         do is=1,rank(flvr)
             taus = time_s( index_s(is, flvr), flvr )

             do ie=1,rank(flvr)
                 taue = time_e( index_e(ie, flvr), flvr )

! evaluate dtau
                 dtau = taue - taus

! get matrix element from mmat, pay special attention to the sign of dtau
                 maux = mmat(ie, is, flvr) * sign(one, dtau) * pref(ie,flvr)

! adjust dtau, keep it stay in (zero, beta)
                 if ( dtau < zero ) then
                     dtau = dtau + beta
                 endif ! back if ( dtau < zero ) block

! convert dtau in [0,\beta] to daux in [0,2]
                 daux = two * dtau / beta

! determine index for legendre polynomial interval
                 curr = nint( daux * step ) + 1

! record ftau, we normalize ftau in ctqmc_make_ftau() subroutine
                 CTQMC_FLALEG_LOOP: do fleg=1,lemax
                     dtau = sqrt(two * fleg - 1) * ppleg(curr,fleg)
                     ftau(fleg, flvr, flvr) = ftau(fleg, flvr, flvr) - maux * dtau
                 enddo CTQMC_FLALEG_LOOP ! over fleg={1,lemax} loop

             enddo ! over ie={1,rank(flvr)} loop
         enddo ! over is={1,rank(flvr)} loop

     enddo CTQMC_FLAVOR_LOOP ! over flvr={1,norbs} loop

     return
  end subroutine cat_record_ftau2

!!>>> cat_record_ftau3: record auxiliary correlation function using
!!>>> chebyshev polynomial representation
  subroutine cat_record_ftau3()
     implicit none

! calculate prefactor: pref
     call ctqmc_make_pref()

! evaluate step at first
     step = real(chgrd - 1) / two

     CTQMC_FLAVOR_LOOP: do flvr=1,norbs

! get imaginary time value for segments
         do is=1,rank(flvr)
             taus = time_s( index_s(is, flvr), flvr )

             do ie=1,rank(flvr)
                 taue = time_e( index_e(ie, flvr), flvr )

! evaluate dtau
                 dtau = taue - taus

! get matrix element from mmat, pay special attention to the sign of dtau
                 maux = mmat(ie, is, flvr) * sign(one, dtau) * pref(ie,flvr)

! adjust dtau, keep it stay in (zero, beta)
                 if ( dtau < zero ) then
                     dtau = dtau + beta
                 endif ! back if ( dtau < zero ) block

! convert dtau in [0,\beta] to daux in [0,2]
                 daux = two * dtau / beta

! determine index for chebyshev polynomial interval
                 curr = nint( daux * step ) + 1

! record ftau, we normalize ftau in ctqmc_make_ftau() subroutine
                 CTQMC_FLACHE_LOOP: do fche=1,chmax
                     dtau = (two / pi) * sqrt(one - (daux - one)**2) * qqche(curr,fche)
                     ftau(fche, flvr, flvr) = ftau(fche, flvr, flvr) - maux * dtau
                 enddo CTQMC_FLACHE_LOOP ! over fche={1,chmax} loop

             enddo ! over ie={1,rank(flvr)} loop
         enddo ! over is={1,rank(flvr)} loop

     enddo CTQMC_FLAVOR_LOOP ! over flvr={1,norbs} loop

     return
  end subroutine cat_record_ftau3
  end subroutine ctqmc_record_ftau

!!>>> ctqmc_record_grnf: record the impurity green's function in matsubara
!!>>> frequency space
  subroutine ctqmc_record_grnf()
     use control, only : norbs
     use control, only : nfreq
     use context, only : gmat
     use context, only : grnf

     implicit none

! local variables
! loop index over matsubara frequencies
     integer :: ifrq

! loop index for flavor channel
     integer :: flvr

! note: only the first nfreq points of grnf are modified
     do flvr=1,norbs
         do ifrq=1,nfreq
             grnf(ifrq, flvr, flvr) = grnf(ifrq, flvr, flvr) + gmat(ifrq, flvr, flvr)
         enddo ! over ifrq={1,nfreq} loop
     enddo ! over flvr={1,norbs} loop

     return
  end subroutine ctqmc_record_grnf

!!>>> ctqmc_record_hist: record the histogram of perturbation expansion series
  subroutine ctqmc_record_hist()
     use constants, only : one

     use control, only : mkink
     use context, only : ckink
     use context, only : hist

     implicit none

! note: if ckink == 0, we record its count in hist(mkink)
     if ( ckink > 0 ) then
         hist(ckink) = hist(ckink) + one
     else
         hist(mkink) = hist(mkink) + one
     endif ! back if ( ckink > 0 ) block

     return
  end subroutine ctqmc_record_hist

!!>>> ctqmc_record_prob: record the probability of atomic states
  subroutine ctqmc_record_prob()
     use constants, only : one

     use control, only : norbs
     use context, only : prob
     use context, only : stts

     implicit none

! local variables
! current flavor channel
     integer :: flvr

! atomic state index
     integer :: pstat

! current atomic state for segment representation
     integer :: state(norbs)

! generate current atomic state
     do flvr=1,norbs
         select case ( stts(flvr) )

             case (0:1)
                 state(flvr) = 0

             case (2:3)
                 state(flvr) = 1

         end select
     enddo ! over flvr={1,norbs} loop

! convert atomic state array to index
     call ctqmc_make_state(norbs, pstat, state)

! accumulate the data
     prob(pstat) = prob(pstat) + one

     return
  end subroutine ctqmc_record_prob

!!>>> ctqmc_record_nmat: record the occupation matrix, double occupation
!!>>> matrix, and auxiliary physical observables simulataneously
  subroutine ctqmc_record_nmat()
     use constants, only : dp, zero, two

     use control, only : nband, norbs
     use control, only : beta
     use context, only : ckink
     use context, only : index_s, index_e, time_s, time_e
     use context, only : paux, nmat, nnmat
     use context, only : rank, stts, uumat

     implicit none

! local variables
! loop index over segments
     integer  :: i

! loop index for flavor channel
     integer  :: flvr

! imaginary time for start and end points
     real(dp) :: ts
     real(dp) :: te

! total length of segments
     real(dp) :: sgmt(norbs)

! used to record overlaps between two segments
     real(dp) :: oaux(norbs)
     real(dp) :: ovlp(norbs,norbs)

! evaluate occupation matrix: < n_i >
!-------------------------------------------------------------------------
     do flvr=1,norbs

! case 1: null occupation
         if      ( stts(flvr) == 0 ) then
             sgmt(flvr) = zero

! case 2: partial occupation, segment scheme
         else if ( stts(flvr) == 1 ) then
             sgmt(flvr) = zero
             do i=1,rank(flvr)
                 ts = time_s(index_s(i, flvr), flvr)
                 te = time_e(index_e(i, flvr), flvr)
                 sgmt(flvr) = sgmt(flvr) + abs( te - ts )
             enddo ! over i={1,rank(flvr)} loop

! case 3: partial occupation, anti-segment scheme
         else if ( stts(flvr) == 2 ) then
             sgmt(flvr) = beta
             do i=1,rank(flvr)
                 ts = time_s(index_s(i, flvr), flvr)
                 te = time_e(index_e(i, flvr), flvr)
                 sgmt(flvr) = sgmt(flvr) - abs( ts - te )
             enddo ! over i={1,rank(flvr)} loop

! case 4: full occupation
         else if ( stts(flvr) == 3 ) then
             sgmt(flvr) = beta

         endif ! back if ( stts(flvr) == 0 ) block

         nmat(flvr) = nmat(flvr) + sgmt(flvr) / beta
     enddo ! over flvr={1,norbs} loop
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate double occupation matrix: < n_i n_j >
!-------------------------------------------------------------------------
     do flvr=1,norbs

! case 1: null occupation
         if      ( stts(flvr) == 0 ) then
             ovlp(flvr,:) = zero

! case 2: partial occupation, segment scheme
         else if ( stts(flvr) == 1 ) then
             ovlp(flvr,:) = zero
             do i=1,rank(flvr)
                 ts = time_s(index_s(i, flvr), flvr)
                 te = time_e(index_e(i, flvr), flvr)
                 call cat_ovlp_segments(flvr, ts, te, oaux)
                 ovlp(flvr,:) = ovlp(flvr,:) + oaux
             enddo ! over i={1,rank(flvr)} loop

! case 3: partial occupation, anti-segment scheme
! pay special attention to the head and tail parts
         else if ( stts(flvr) == 2 ) then
             ovlp(flvr,:) = zero
             do i=1,rank(flvr)-1
                 ts = time_s(index_s(i,   flvr), flvr)
                 te = time_e(index_e(i+1, flvr), flvr)
                 call cat_ovlp_segments(flvr, ts, te, oaux)
                 ovlp(flvr,:) = ovlp(flvr,:) + oaux
             enddo ! over i={1,rank(flvr)-1} loop

             te = time_e(index_e(1, flvr), flvr)
             call cat_ovlp_segments(flvr, zero, te, oaux)
             ovlp(flvr,:) = ovlp(flvr,:) + oaux

             ts = time_s(index_s(rank(flvr), flvr), flvr)
             call cat_ovlp_segments(flvr, ts, beta, oaux)
             ovlp(flvr,:) = ovlp(flvr,:) + oaux

! case 4: full occupation
         else if ( stts(flvr) == 3 ) then
             call cat_ovlp_segments(flvr, zero, beta, oaux)
             ovlp(flvr,:) = oaux

         endif ! back if ( stts(flvr) == 0 ) block

         nnmat(flvr,:) = nnmat(flvr,:) + ovlp(flvr,:) / beta
     enddo ! over flvr={1,norbs} loop
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate <K^4>
!-------------------------------------------------------------------------
     paux(9) = paux(9) + ( ckink * two )**4
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate <K^3>
!-------------------------------------------------------------------------
     paux(8) = paux(8) + ( ckink * two )**3
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate <K^2>
!-------------------------------------------------------------------------
     paux(7) = paux(7) + ( ckink * two )**2
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate <N^2>
!-------------------------------------------------------------------------
     paux(6) = paux(6) + ( sum(sgmt) / beta )**2
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate <N^1>
!-------------------------------------------------------------------------
     paux(5) = paux(5) + sum(sgmt) / beta
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate spin magnetization: < Sz >
!-------------------------------------------------------------------------
     do flvr=1,nband
         paux(4) = paux(4) + ( sgmt(flvr) - sgmt(flvr+nband) ) / beta
     enddo ! over flvr={1,nband} loop
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate kinetic energy: ekin
! equation : -T < k >
!-------------------------------------------------------------------------
     paux(3) = paux(3) - real(ckink * norbs) / beta
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate potential energy: epot
!-------------------------------------------------------------------------
     do flvr=1,norbs
         do i=1,flvr
             paux(2) = paux(2) + uumat(flvr,i) * ovlp(flvr,i) / beta
         enddo ! over i={1,flvr} loop
     enddo ! over flvr={1,norbs} loop
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate total energy: etot
!-------------------------------------------------------------------------
     paux(1) = paux(2) + paux(3)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     return
  end subroutine ctqmc_record_nmat

!!>>> ctqmc_record_kmat: record the < k^2 > - < k >^2
  subroutine ctqmc_record_kmat()
     use constants, only : dp

     use control, only : issus
     use control, only : norbs
     use context, only : kmat, kkmat
     use context, only : rank

     implicit none

! local variables
! loop index for flavor channel
     integer :: i
     integer :: j

! check whether there is conflict
     call s_assert( btest(issus, 5) )

! since rank means the number of operator pairs,
! so we have to multiply it with two
     do i=1,norbs
         kmat(i) = kmat(i) + rank(i) * 2.0_dp
     enddo ! over i={1,norbs} loop

     do i=1,norbs
         do j=1,norbs
             kkmat(i,j) = kkmat(i,j) + rank(i) * rank(j) * 4.0_dp
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

     return
  end subroutine ctqmc_record_kmat

!!>>> ctqmc_record_lmat: record the fidelity susceptibility
  subroutine ctqmc_record_lmat()
     use constants, only : dp, zero, one, two

     use control, only : issus
     use control, only : norbs
     use control, only : beta
     use context, only : index_s, index_e, time_s, time_e
     use context, only : lmat, rmat, lrmat
     use context, only : rank

     implicit none

! local variables
! loop index over segments
     integer  :: i

! loop index for flavor channel
     integer  :: flvr

! imaginary time for start and end points
     real(dp) :: ts
     real(dp) :: te

! number of operators at left half axis for the current configuration
     real(dp) :: kl(norbs)

! number of operators at right half axis for the current configuration
     real(dp) :: kr(norbs)

! check whether there is conflict
     call s_assert( btest(issus, 6) )

! init k_l and k_r
     kl = zero
     kr = zero

! loop over flavors and segments to calculate k_l and k_r
     do flvr=1,norbs
         do i=1,rank(flvr)
             ts = time_s(index_s(i, flvr), flvr)
             if ( ts < beta / two ) then
                 kl(flvr) = kl(flvr) + one
             else
                 kr(flvr) = kr(flvr) + one
             endif ! back if ( ts < beta / two ) block

             te = time_e(index_e(i, flvr), flvr)
             if ( te < beta / two ) then
                 kl(flvr) = kl(flvr) + one
             else
                 kr(flvr) = kr(flvr) + one
             endif ! back if ( te < beta / two ) block
         enddo ! over i={1,rank(flvr)} loop
     enddo ! over flvr={1,norbs} loop

! add contribution to < k_l > and < k_r >
     lmat = lmat + kl
     rmat = rmat + kr

! add contribution to < k_l k_r >
     do flvr=1,norbs
         do i=1,norbs
             lrmat(i,flvr) = lrmat(i,flvr) + kl(i) * kr(flvr)
         enddo ! over i={1,norbs} loop
     enddo ! over flvr={1,norbs} loop

     return
  end subroutine ctqmc_record_lmat

!!>>> ctqmc_record_szpw: record the powers of local magnetization
  subroutine ctqmc_record_szpw()
     use constants, only : dp, zero

     use control, only : issus
     use control, only : nband, norbs
     use control, only : ntime
     use control, only : beta
     use context, only : tmesh
     use context, only : szpow

     implicit none

! local variables
! loop index over times
     integer  :: i

! loop index for flavor channel
     integer  :: f1
     integer  :: f2

! \delta \tau
     real(dp) :: step

! integral of Sz(\tau), in order words:
!     sint = 1/\beta \int^{\beta}_0 Sz(\tau) d\tau
     real(dp) :: sint

! used to record occupations for current flavor channel and time
     real(dp) :: oaux(ntime,norbs)

! orbital-resolved Sz(\tau)
     real(dp) :: saux(ntime,nband)

! check whether there is conflict
     call s_assert( btest(issus, 7) )

! calculate oaux, obtain occupation status
! calculate saux, obtain Sz(\tau)
     oaux = zero
     saux = zero
     TIME_LOOP: do i=1,ntime
         do f1=1,norbs
             call cat_occupy_status(f1, tmesh(i), oaux(i,f1))
         enddo ! over f1={1,norbs} loop

         do f2=1,nband
             saux(i,f2) = oaux(i,f2) - oaux(i,f2+nband)
         enddo ! over f2={1,nband} loop
     enddo TIME_LOOP ! over i={1,ntime} loop

! accumulate szpow(1:4,1:nband)
! calculate \delta \tau
     step = ( tmesh(2) - tmesh(1) ) / 2.0
     BAND_LOOP: do f2=1,nband
! calculate sint using trapezoid algorithm 
         sint = zero
         do i=1,ntime-1
             sint = sint + ( saux(i,f2) + saux(i+1,f2) ) * step
         enddo ! over i={1,ntime-1} loop
         sint = sint / beta
! record the data
         szpow(1,f2) = szpow(1,f2) + sint**1.0
         szpow(2,f2) = szpow(2,f2) + sint**2.0
         szpow(3,f2) = szpow(3,f2) + sint**3.0
         szpow(4,f2) = szpow(4,f2) + sint**4.0
     enddo BAND_LOOP ! over f2={1,nband} loop

! accumulate szpow(1:4,nband+1)
! here we consider the contribution from all flavors
     sint = zero
     do i=1,ntime-1
         sint = sint + ( sum( saux(i,:) ) + sum( saux(i+1,:) ) ) * step
     enddo ! over i={1,ntime-1} loop
     sint = sint / beta
! record the data
     szpow(1,nband+1) = szpow(1,nband+1) + sint**1.0
     szpow(2,nband+1) = szpow(2,nband+1) + sint**2.0
     szpow(3,nband+1) = szpow(3,nband+1) + sint**3.0
     szpow(4,nband+1) = szpow(4,nband+1) + sint**4.0

     return
  end subroutine ctqmc_record_szpw

!!>>> ctqmc_record_schi: record the spin-spin correlation function
!!>>> imaginary-time version
  subroutine ctqmc_record_schi()
     use constants, only : dp, zero
     use spring, only : spring_sfmt_stream

     use control, only : issus
     use control, only : nband, norbs
     use control, only : ntime
     use context, only : tmesh
     use context, only : schi, sschi

     implicit none

! local parameters
! number of internal loop
! if you want to obtain more accurate results, please increase it
     integer, parameter :: num_try = 16

! local variables
! loop index over times
     integer  :: i
     integer  :: m
     integer  :: n

! loop index for flavor channel
     integer  :: f1

! used to record occupations for current flavor channel and time
     real(dp) :: oaux(ntime,norbs)

! check whether there is conflict
     call s_assert( btest(issus, 1) )

! calculate oaux, obtain occupation status
     oaux = zero
     TIME_LOOP: do i=1,ntime
         do f1=1,norbs
             call cat_occupy_status(f1, tmesh(i), oaux(i,f1))
         enddo ! over f1={1,norbs} loop
     enddo TIME_LOOP ! over i={1,ntime} loop
     oaux = oaux / real(num_try)

! calculate schi and sschi
     do f1=1,nband
         do i=1,num_try
             m = ceiling( spring_sfmt_stream() * ntime )
             if ( oaux(m,f1) > zero ) then
! n - m + ntime \in [ntime - m + 1, ntime]
                 do n=1,m
                     schi(n-m+ntime) = schi(n-m+ntime) + oaux(n,f1)
                     schi(n-m+ntime) = schi(n-m+ntime) - oaux(n,f1+nband)
                     sschi(n-m+ntime,f1) = sschi(n-m+ntime,f1) + oaux(n,f1)
                     sschi(n-m+ntime,f1) = sschi(n-m+ntime,f1) - oaux(n,f1+nband)
                 enddo ! over n={1,m} loop
! n - m \in [1, ntime - m]
                 do n=m+1,ntime
                     schi(n-m) = schi(n-m) + oaux(n,f1)
                     schi(n-m) = schi(n-m) - oaux(n,f1+nband)
                     sschi(n-m,f1) = sschi(n-m,f1) + oaux(n,f1)
                     sschi(n-m,f1) = sschi(n-m,f1) - oaux(n,f1+nband)
                 enddo ! over n={m+1,ntime} loop
             endif ! back if ( oaux(m,f1) > zero ) block

             if ( oaux(m,f1+nband) > zero ) then ! oaux(m,f1+nband) = one
! n - m + ntime \in [ntime - m + 1, ntime]
                 do n=1,m
                     schi(n-m+ntime) = schi(n-m+ntime) + oaux(n,f1+nband)
                     schi(n-m+ntime) = schi(n-m+ntime) - oaux(n,f1)
                     sschi(n-m+ntime,f1) = sschi(n-m+ntime,f1) + oaux(n,f1+nband)
                     sschi(n-m+ntime,f1) = sschi(n-m+ntime,f1) - oaux(n,f1)
                 enddo ! over n={1,m} loop
! n - m \in [1, ntime - m]
                 do n=m+1,ntime
                     schi(n-m) = schi(n-m) + oaux(n,f1+nband)
                     schi(n-m) = schi(n-m) - oaux(n,f1)
                     sschi(n-m,f1) = sschi(n-m,f1) + oaux(n,f1+nband)
                     sschi(n-m,f1) = sschi(n-m,f1) - oaux(n,f1)
                 enddo ! over n={m+1,ntime} loop
             endif ! back if ( oaux(m,f1+nband) > zero ) block
         enddo ! over i={1,num_try} loop
     enddo ! over f1={1,nband} loop

     return
  end subroutine ctqmc_record_schi

!!>>> ctqmc_record_sfom: record the spin-spin correlation function
!!>>> matsubara frequency version
  subroutine ctqmc_record_sfom()
     use constants, only : dp, zero, one, two, pi, czi

     use control, only : issus
     use control, only : nband, norbs
     use control, only : nbfrq
     use control, only : beta
     use context, only : index_s, index_e, time_s, time_e
     use context, only : ssfom
     use context, only : rank

     implicit none

! local variables
! loop index for flavor channel
     integer  :: f1
     integer  :: f2

! loop index for operators
     integer  :: it

! imaginary time for start and end points
     real(dp) :: taus
     real(dp) :: taue

! the first bosonic frequency
     complex(dp) :: dw

! used to record occupations for current flavor channel at \tau = 0
     real(dp) :: oaux(norbs)

! bosonic frequency mesh
     complex(dp) :: mesh(nbfrq)

! matsubara frequency exponents for create operators
     complex(dp) :: exps(nbfrq)

! matsubara frequency exponents for destroy operators
     complex(dp) :: expe(nbfrq)

! check whether there is conflict
     call s_assert( btest(issus, 3) )

! build bosonic frequency mesh
     dw = czi * two * pi / beta
     mesh = dw
     call s_cumsum_z(nbfrq, mesh, mesh)

! calculate oaux, obtain occupation status
     do f1=1,norbs
         call cat_occupy_status(f1, zero, oaux(f1))
     enddo ! over i={1,norbs} loop

! calculate ssfom, it must be real
! < Sz(t)Sz(0) > = < ( nu(t) - nd(t) ) * ( nu(0) - nd(0) ) >
     do f1=1,nband
         f2 = f1 + nband
! the contribution from oaux(f1) = one and oaux(f2) = one is zero
! the contribution from oaux(f1) = zero and oaux(f2) = zero is also zero
! here oaux(f1) = one; oaux(f2) = zero
         if ( oaux(f1) > zero .and. oaux(f2) < one ) then
! + nu(t)nu(0) term
             do it=1,rank(f1)
                 taus = time_s( index_s(it, f1), f1 )
                 taue = time_e( index_e(it, f1), f1 )
                 exps = exp( dw * taus )
                 expe = exp( dw * taue )
                 call s_cumprod_z(nbfrq, exps, exps)
                 call s_cumprod_z(nbfrq, expe, expe)
                 ssfom(:,f1) = ssfom(:,f1) + real( ( expe - exps ) / mesh )
             enddo ! over do it={1,rank(f1)} loop
! - nd(t)nu(0) term
             do it=1,rank(f2)
                 taus = time_s( index_s(it, f2), f2 )
                 taue = time_e( index_e(it, f2), f2 )
                 exps = exp( dw * taus )
                 expe = exp( dw * taue )
                 call s_cumprod_z(nbfrq, exps, exps)
                 call s_cumprod_z(nbfrq, expe, expe)
                 ssfom(:,f1) = ssfom(:,f1) - real( ( expe - exps ) / mesh )
             enddo ! over do it={1,rank(f2)} loop
         endif ! back if ( oaux(f1) > zero .and. oaux(f2) < one ) block

! here oaux(f2) = one; oaux(f1) = zero
         if ( oaux(f2) > zero .and. oaux(f1) < one ) then
! - nu(t)nd(0) term
             do it=1,rank(f1)
                 taus = time_s( index_s(it, f1), f1 )
                 taue = time_e( index_e(it, f1), f1 )
                 exps = exp( dw * taus )
                 expe = exp( dw * taue )
                 call s_cumprod_z(nbfrq, exps, exps)
                 call s_cumprod_z(nbfrq, expe, expe)
                 ssfom(:,f1) = ssfom(:,f1) - real( ( expe - exps ) / mesh )
             enddo ! over do it={1,rank(f1)} loop
! + nd(t)nd(0) term
             do it=1,rank(f2)
                 taus = time_s( index_s(it, f2), f2 )
                 taue = time_e( index_e(it, f2), f2 )
                 exps = exp( dw * taus )
                 expe = exp( dw * taue )
                 call s_cumprod_z(nbfrq, exps, exps)
                 call s_cumprod_z(nbfrq, expe, expe)
                 ssfom(:,f1) = ssfom(:,f1) + real( ( expe - exps ) / mesh )
             enddo ! over do it={1,rank(f2)} loop
         endif ! back if ( oaux(f2) > zero .and. oaux(f1) < one ) block
     enddo ! over f1={1,nband} loop

     return
  end subroutine ctqmc_record_sfom

!!>>> ctqmc_record_ochi: record the orbital-orbital correlation function
!!>>> imaginary-time version
  subroutine ctqmc_record_ochi()
     use constants, only : dp, zero
     use spring, only : spring_sfmt_stream

     use control, only : issus
     use control, only : norbs
     use control, only : ntime
     use context, only : tmesh
     use context, only : ochi, oochi

     implicit none

! local parameters
! number of internal loop
! if you want to obtain more accurate results, please increase it
     integer, parameter :: num_try = 16

! local variables
! loop index over times
     integer  :: i
     integer  :: m
     integer  :: n

! loop index for flavor channel
     integer  :: f1
     integer  :: f2

! used to record occupations for current flavor channel and time
     real(dp) :: oaux(ntime,norbs)

! check whether there is conflict
     call s_assert( btest(issus, 2) )

! calculate oaux, obtain occupation status
     oaux = zero
     TIME_LOOP: do i=1,ntime
         do f1=1,norbs
             call cat_occupy_status(f1, tmesh(i), oaux(i,f1))
         enddo ! over f1={1,norbs} loop
     enddo TIME_LOOP ! over i={1,ntime} loop
     oaux = oaux / real(num_try)

! calculate ochi and oochi
     do f1=1,norbs
         do f2=1,norbs
             do i=1,num_try
                 m = ceiling( spring_sfmt_stream() * ntime )
                 if ( oaux(m,f2) > zero ) then
! n - m + ntime \in [ntime - m + 1, ntime]
                     do n=1,m
                         ochi(n-m+ntime) = ochi(n-m+ntime) + oaux(n,f1)
                         oochi(n-m+ntime,f2,f1) = oochi(n-m+ntime,f2,f1) + oaux(n,f1)
                     enddo ! over n={1,m} loop
! n - m \in [1, ntime - m]
                     do n=m+1,ntime
                         ochi(n-m) = ochi(n-m) + oaux(n,f1)
                         oochi(n-m,f2,f1) = oochi(n-m,f2,f1) + oaux(n,f1)
                     enddo ! over n={m+1,ntime} loop
                 endif ! back if ( oaux(m,f2) > zero ) block
             enddo ! over i={1,num_try} loop
         enddo ! over f2={1,norbs} loop
     enddo ! over f1={1,norbs} loop

     return
  end subroutine ctqmc_record_ochi

!!>>> ctqmc_record_ofom: record the orbital-orbital correlation function
!!>>> matsubara frequency version
  subroutine ctqmc_record_ofom()
     use constants, only : dp, zero, two, pi, czi

     use control, only : issus
     use control, only : norbs
     use control, only : nbfrq
     use control, only : beta
     use context, only : index_s, index_e, time_s, time_e
     use context, only : oofom
     use context, only : rank

     implicit none

! local variables
! loop index for flavor channel
     integer  :: f1
     integer  :: f2

! loop index for operators
     integer  :: it

! imaginary time for start and end points
     real(dp) :: taus
     real(dp) :: taue

! the first bosonic frequency
     complex(dp) :: dw

! used to record occupations for current flavor channel at \tau = 0
     real(dp) :: oaux(norbs)

! bosonic frequency mesh
     complex(dp) :: mesh(nbfrq)

! matsubara frequency exponents for create operators
     complex(dp) :: exps(nbfrq)

! matsubara frequency exponents for destroy operators
     complex(dp) :: expe(nbfrq)

! check whether there is conflict
     call s_assert( btest(issus, 4) )

! build bosonic frequency mesh
     dw = czi * two * pi / beta
     mesh = dw
     call s_cumsum_z(nbfrq, mesh, mesh)

! calculate oaux, obtain occupation status
     do f1=1,norbs
         call cat_occupy_status(f1, zero, oaux(f1))
     enddo ! over i={1,norbs} loop

! calculate oofom, it must be real
     do f1=1,norbs
         do f2=1,f1
             if ( oaux(f2) > zero ) then
                 do it=1,rank(f1)
                     taus = time_s( index_s(it, f1), f1 )
                     taue = time_e( index_e(it, f1), f1 )
                     exps = exp( dw * taus )
                     expe = exp( dw * taue )
                     call s_cumprod_z(nbfrq, exps, exps)
                     call s_cumprod_z(nbfrq, expe, expe)
                     oofom(:,f2,f1) = oofom(:,f2,f1) + real( ( expe - exps ) / mesh )
                 enddo ! over do it={1,rank(f1)} loop
             endif ! back if ( oaux(f2) > zero ) block
             if ( f1 /= f2 ) then ! consider the symmetry
                 oofom(:,f1,f2) = oofom(:,f2,f1)
             endif ! back if ( f1 /= f2 ) block
         enddo ! over f2={1,f1} loop
     enddo ! over f1={1,norbs} loop

     return
  end subroutine ctqmc_record_ofom

!!>>> ctqmc_record_twop: record the two-particle green's function
  subroutine ctqmc_record_twop()
     use constants, only : dp, czero

     use control, only : isvrt
     use control, only : norbs
     use control, only : nffrq, nbfrq
     use control, only : beta
     use context, only : g2_re, g2_im
     use context, only : rank
     use context, only : mmat

     implicit none

! local variables
! loop indices for start and end points
     integer  :: is
     integer  :: ie

! loop index for flavor channel
     integer  :: f1
     integer  :: f2
     integer  :: flvr

! loop index for frequency
     integer  :: nfaux
     integer  :: wbn
     integer  :: w1n
     integer  :: w2n
     integer  :: w3n
     integer  :: w4n

! used to store the element of mmat matrix
     real(dp) :: maux

! dummy complex(dp) variables, used to calculate the g2_re and g2_im
     complex(dp) :: cmeas

! dummy complex(dp) arrays, used to store the intermediate results
     complex(dp), allocatable :: g2aux(:,:,:)
     complex(dp), allocatable :: caux1(:,:)
     complex(dp), allocatable :: caux2(:,:)

! check whether there is conflict
     call s_assert( btest(isvrt, 1) .and. .not. btest(isvrt, 2) )

! evaluate nfaux, determine the size of g2aux
     nfaux = nffrq + nbfrq - 1

! allocate memory for g2aux and then initialize it
     allocate( g2aux(nfaux, nfaux, norbs) ); g2aux = czero

! allocate memory for caux1 and caux2, and then initialize them
     allocate( caux1(nfaux, maxval(rank)) ); caux1 = czero
     allocate( caux2(nfaux, maxval(rank)) ); caux2 = czero

! calculate g2aux: see Eq. (52) in Phys. Rev. B 89, 235128 (2014)
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE (flvr, is, ie, maux, w2n, w1n, caux1, caux2)
     CTQMC_FLAVOR_LOOP: do flvr=1,norbs
         call ctqmc_make_prod(flvr, nfaux, maxval(rank), caux1, caux2)

         do is=1,rank(flvr)
             do ie=1,rank(flvr)

                 maux = mmat(ie, is, flvr)
                 do w2n=1,nfaux
                     do w1n=1,nfaux
                         g2aux(w1n,w2n,flvr) = g2aux(w1n,w2n,flvr) + maux * caux1(w2n,is) * caux2(w1n,ie)
                     enddo ! over w1n={1,nfaux} loop
                 enddo ! over w2n={1,nfaux} loop

             enddo ! over ie={1,rank(flvr)} loop
         enddo ! over is={1,rank(flvr)} loop

     enddo CTQMC_FLAVOR_LOOP ! over flvr={1,norbs} loop
!$OMP END DO

! calculate g2_re and g2_im
!$OMP DO PRIVATE (f1, f2, cmeas, wbn, w4n, w3n, w2n, w1n)
     CTQMC_ORBIT1_LOOP: do f1=1,norbs
         CTQMC_ORBIT2_LOOP: do f2=1,f1

             CTQMC_BOSONF_LOOP: do wbn=1,nbfrq

                 CTQMC_FERMI1_LOOP: do w2n=1,nffrq
                     CTQMC_FERMI2_LOOP: do w3n=1,nffrq
                         w1n = w2n + wbn - 1; w4n = w3n + wbn - 1

                         cmeas = g2aux(w1n,w2n,f1) * g2aux(w3n,w4n,f2)
                         if ( f1 == f2 ) then
                             cmeas = cmeas - g2aux(w1n,w4n,f1) * g2aux(w3n,w2n,f1)
                         endif ! back if ( f1 == f2 ) block
                         g2_re(w3n,w2n,wbn,f2,f1) = g2_re(w3n,w2n,wbn,f2,f1) +  real(cmeas) / beta
                         g2_im(w3n,w2n,wbn,f2,f1) = g2_im(w3n,w2n,wbn,f2,f1) + aimag(cmeas) / beta
                     enddo CTQMC_FERMI2_LOOP ! over w3n={1,nffrq} loop
                 enddo CTQMC_FERMI1_LOOP ! over w2n={1,nffrq} loop

             enddo CTQMC_BOSONF_LOOP ! over wbn={1,nbfrq} loop

         enddo CTQMC_ORBIT2_LOOP ! over f2={1,f1} loop
     enddo CTQMC_ORBIT1_LOOP ! over f1={1,norbs} loop
!$OMP END DO
!$OMP END PARALLEL

! deallocate memory
     deallocate( g2aux )
     deallocate( caux1 )
     deallocate( caux2 )

     return
  end subroutine ctqmc_record_twop

!!>>> ctqmc_record_vrtx: record the two-particle green's function
!!>>> improved estimator is used to improve the accuracy
  subroutine ctqmc_record_vrtx()
     use constants, only : dp, czero

     use control, only : isvrt
     use control, only : norbs
     use control, only : nffrq, nbfrq
     use control, only : beta
     use context, only : g2_re, g2_im, h2_re, h2_im
     use context, only : rank, pref
     use context, only : mmat

     implicit none

! local variables
! loop indices for start and end points
     integer  :: is
     integer  :: ie

! loop index for flavor channel
     integer  :: f1
     integer  :: f2
     integer  :: flvr

! loop index for frequency
     integer  :: nfaux
     integer  :: wbn
     integer  :: w1n
     integer  :: w2n
     integer  :: w3n
     integer  :: w4n

! used to store the element of mmat matrix
     real(dp) :: maux
     real(dp) :: naux

! dummy complex(dp) variables, used to calculate the h2_re and h2_im
     complex(dp) :: cmeas

! dummy complex(dp) arrays, used to store the intermediate results
     complex(dp), allocatable :: g2aux(:,:,:)
     complex(dp), allocatable :: h2aux(:,:,:)
     complex(dp), allocatable :: caux1(:,:)
     complex(dp), allocatable :: caux2(:,:)

! check whether there is conflict
     call s_assert( btest(isvrt, 2) .and. .not. btest(isvrt, 1) )

! evaluate nfaux, determine the size of g2aux and h2aux
     nfaux = nffrq + nbfrq - 1

! allocate memory for g2aux and then initialize it
     allocate( g2aux(nfaux, nfaux, norbs) ); g2aux = czero

! allocate memory for h2aux and then initialize it
     allocate( h2aux(nfaux, nfaux, norbs) ); h2aux = czero

! allocate memory for caux1 and caux2, and then initialize them
     allocate( caux1(nfaux, maxval(rank)) ); caux1 = czero
     allocate( caux2(nfaux, maxval(rank)) ); caux2 = czero

! calculate prefactor: pref
     call ctqmc_make_pref()

! calculate g2aux and h2aux: see Eq. (52) in Phys. Rev. B 89, 235128 (2014)
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE (flvr, is, ie, maux, naux, w2n, w1n, caux1, caux2)
     CTQMC_FLAVOR_LOOP: do flvr=1,norbs
         call ctqmc_make_prod(flvr, nfaux, maxval(rank), caux1, caux2)

         do is=1,rank(flvr)
             do ie=1,rank(flvr)

                 maux = mmat(ie, is, flvr)
                 naux = mmat(ie, is, flvr) * pref(ie,flvr)
                 do w2n=1,nfaux
                     do w1n=1,nfaux
                         g2aux(w1n,w2n,flvr) = g2aux(w1n,w2n,flvr) + maux * caux1(w2n,is) * caux2(w1n,ie)
                         h2aux(w1n,w2n,flvr) = h2aux(w1n,w2n,flvr) + naux * caux1(w2n,is) * caux2(w1n,ie)
                     enddo ! over w1n={1,nfaux} loop
                 enddo ! over w2n={1,nfaux} loop

             enddo ! over ie={1,rank(flvr)} loop
         enddo ! over is={1,rank(flvr)} loop

     enddo CTQMC_FLAVOR_LOOP ! over flvr={1,norbs} loop
!$OMP END DO

! calculate g2_re and g2_im, h2_re and h2_im
!$OMP DO PRIVATE (f1, f2, cmeas, wbn, w4n, w3n, w2n, w1n)
     CTQMC_ORBIT1_LOOP: do f1=1,norbs
         CTQMC_ORBIT2_LOOP: do f2=1,f1

             CTQMC_BOSONF_LOOP: do wbn=1,nbfrq

                 CTQMC_FERMI1_LOOP: do w2n=1,nffrq
                     CTQMC_FERMI2_LOOP: do w3n=1,nffrq
                         w1n = w2n + wbn - 1; w4n = w3n + wbn - 1

                         cmeas = g2aux(w1n,w2n,f1) * g2aux(w3n,w4n,f2)
                         if ( f1 == f2 ) then
                             cmeas = cmeas - g2aux(w1n,w4n,f1) * g2aux(w3n,w2n,f1)
                         endif ! back if ( f1 == f2 ) block
                         g2_re(w3n,w2n,wbn,f2,f1) = g2_re(w3n,w2n,wbn,f2,f1) +  real(cmeas) / beta
                         g2_im(w3n,w2n,wbn,f2,f1) = g2_im(w3n,w2n,wbn,f2,f1) + aimag(cmeas) / beta

                         cmeas = h2aux(w1n,w2n,f1) * g2aux(w3n,w4n,f2)
                         if ( f1 == f2 ) then
                             cmeas = cmeas - h2aux(w1n,w4n,f1) * g2aux(w3n,w2n,f1)
                         endif ! back if ( f1 == f2 ) block
                         h2_re(w3n,w2n,wbn,f2,f1) = h2_re(w3n,w2n,wbn,f2,f1) +  real(cmeas) / beta
                         h2_im(w3n,w2n,wbn,f2,f1) = h2_im(w3n,w2n,wbn,f2,f1) + aimag(cmeas) / beta
                     enddo CTQMC_FERMI2_LOOP ! over w3n={1,nffrq} loop
                 enddo CTQMC_FERMI1_LOOP ! over w2n={1,nffrq} loop

             enddo CTQMC_BOSONF_LOOP ! over wbn={1,nbfrq} loop

         enddo CTQMC_ORBIT2_LOOP ! over f2={1,f1} loop
     enddo CTQMC_ORBIT1_LOOP ! over f1={1,norbs} loop
!$OMP END DO
!$OMP END PARALLEL

! deallocate memory
     deallocate( g2aux )
     deallocate( h2aux )
     deallocate( caux1 )
     deallocate( caux2 )

     return
  end subroutine ctqmc_record_vrtx

!!>>> ctqmc_record_pair: record the particle-particle pair susceptibility
  subroutine ctqmc_record_pair()
     use constants, only : dp, czero

     use control, only : isvrt
     use control, only : norbs
     use control, only : nffrq, nbfrq
     use control, only : beta
     use context, only : ps_re, ps_im
     use context, only : rank
     use context, only : mmat

     implicit none

! local variables
! loop indices for start and end points
     integer  :: is
     integer  :: ie

! loop index for flavor channel
     integer  :: f1
     integer  :: f2
     integer  :: flvr

! loop index for frequency
     integer  :: nfaux
     integer  :: wbn
     integer  :: w1n
     integer  :: w2n
     integer  :: w3n
     integer  :: w4n

! used to store the element of mmat matrix
     real(dp) :: maux

! dummy complex(dp) variables, used to calculate the ps_re and ps_im
     complex(dp) :: cmeas

! dummy complex(dp) arrays, used to store the intermediate results
     complex(dp), allocatable :: g2aux(:,:,:)
     complex(dp), allocatable :: caux1(:,:)
     complex(dp), allocatable :: caux2(:,:)

! check whether there is conflict
     call s_assert( btest(isvrt, 3) )

! evaluate nfaux, determine the size of g2aux
     nfaux = nffrq + nbfrq - 1

! allocate memory for g2aux and then initialize it
     allocate( g2aux(nfaux, nfaux, norbs) ); g2aux = czero

! allocate memory for caux1 and caux2, and then initialize them
     allocate( caux1(nfaux, maxval(rank)) ); caux1 = czero
     allocate( caux2(nfaux, maxval(rank)) ); caux2 = czero

! calculate g2aux: see Eq. (52) in Phys. Rev. B 89, 235128 (2014)
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE (flvr, is, ie, maux, w2n, w1n, caux1, caux2)
     CTQMC_FLAVOR_LOOP: do flvr=1,norbs
         call ctqmc_make_prod(flvr, nfaux, maxval(rank), caux1, caux2)

         do is=1,rank(flvr)
             do ie=1,rank(flvr)

                 maux = mmat(ie, is, flvr)
                 do w2n=1,nfaux
                     do w1n=1,nfaux
                         g2aux(w1n,w2n,flvr) = g2aux(w1n,w2n,flvr) + maux * caux1(w2n,is) * caux2(w1n,ie)
                     enddo ! over w1n={1,nfaux} loop
                 enddo ! over w2n={1,nfaux} loop

             enddo ! over ie={1,rank(flvr)} loop
         enddo ! over is={1,rank(flvr)} loop

     enddo CTQMC_FLAVOR_LOOP ! over flvr={1,norbs} loop
!$OMP END DO

! calculate ps_re and ps_im
!$OMP DO PRIVATE (f1, f2, cmeas, wbn, w4n, w3n, w2n, w1n)
     CTQMC_ORBIT1_LOOP: do f1=1,norbs
         CTQMC_ORBIT2_LOOP: do f2=1,f1

             CTQMC_BOSONF_LOOP: do wbn=1,nbfrq

                 CTQMC_FERMI1_LOOP: do w2n=1,nffrq
                     CTQMC_FERMI2_LOOP: do w3n=1,nffrq
                         w1n = w2n + wbn - 1; w4n = w3n + wbn - 1

                         cmeas = czero
                         if ( f1 /= f2 ) then
                             cmeas = cmeas + g2aux(w1n,w4n,f1) * g2aux(nffrq-w2n+1,nffrq-w3n+1,f2)
                         endif ! back if ( f1 == f2 ) block
                         ps_re(w3n,w2n,wbn,f2,f1) = ps_re(w3n,w2n,wbn,f2,f1) +  real(cmeas) / beta
                         ps_im(w3n,w2n,wbn,f2,f1) = ps_im(w3n,w2n,wbn,f2,f1) + aimag(cmeas) / beta
                     enddo CTQMC_FERMI2_LOOP ! over w3n={1,nffrq} loop
                 enddo CTQMC_FERMI1_LOOP ! over w2n={1,nffrq} loop

             enddo CTQMC_BOSONF_LOOP ! over wbn={1,nbfrq} loop

         enddo CTQMC_ORBIT2_LOOP ! over f2={1,f1} loop
     enddo CTQMC_ORBIT1_LOOP ! over f1={1,norbs} loop
!$OMP END DO
!$OMP END PARALLEL

! deallocate memory
     deallocate( g2aux )
     deallocate( caux1 )
     deallocate( caux2 )

     return
  end subroutine ctqmc_record_pair

!!========================================================================
!!>>> reduce physical observables                                      <<<
!!========================================================================

!!>>> ctqmc_reduce_gtau: reduce the gtau from all children processes
  subroutine ctqmc_reduce_gtau(gtau_mpi, gtau_err)
     use constants, only : dp, zero
     use mmpi, only : mp_allreduce, mp_barrier

     use control, only : norbs
     use control, only : ntime
     use control, only : nprocs
     use context, only : gtau

     implicit none

! external arguments
! impurity green's function
     real(dp), intent(out) :: gtau_mpi(ntime,norbs,norbs)
     real(dp), intent(out) :: gtau_err(ntime,norbs,norbs)

! initialize gtau_mpi and gtau_err
     gtau_mpi = zero
     gtau_err = zero

! build gtau_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(gtau, gtau_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     gtau_mpi = gtau

# endif /* MPI */

! calculate the average
     gtau_mpi = gtau_mpi / real(nprocs)

! build gtau_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce((gtau - gtau_mpi)**2, gtau_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         gtau_err = sqrt( gtau_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

     return
  end subroutine ctqmc_reduce_gtau

!!>>> ctqmc_reduce_ftau: reduce the ftau from all children processes
  subroutine ctqmc_reduce_ftau(ftau_mpi, ftau_err)
     use constants, only : dp, zero
     use mmpi, only : mp_allreduce, mp_barrier

     use control, only : norbs
     use control, only : ntime
     use control, only : nprocs
     use context, only : ftau

     implicit none

! external arguments
! auxiliary correlation function, F(\tau)
     real(dp), intent(out) :: ftau_mpi(ntime,norbs,norbs)
     real(dp), intent(out) :: ftau_err(ntime,norbs,norbs)

! initialize ftau_mpi and ftau_err
     ftau_mpi = zero
     ftau_err = zero

! build ftau_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(ftau, ftau_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     ftau_mpi = ftau

# endif /* MPI */

! calculate the average
     ftau_mpi = ftau_mpi / real(nprocs)

! build ftau_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce((ftau - ftau_mpi)**2, ftau_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         ftau_err = sqrt( ftau_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

     return
  end subroutine ctqmc_reduce_ftau

!!>>> ctqmc_reduce_grnf: reduce the grnf from all children processes
  subroutine ctqmc_reduce_grnf(grnf_mpi, grnf_err)
     use constants, only : dp, zero, czero, czi
     use mmpi, only : mp_allreduce, mp_barrier

     use control, only : norbs
     use control, only : mfreq
     use control, only : nprocs
     use context, only : grnf

     implicit none

! external arguments
! impurity green's function
     complex(dp), intent(out) :: grnf_mpi(mfreq,norbs,norbs)
     complex(dp), intent(out) :: grnf_err(mfreq,norbs,norbs)

! local variables
! used to store the real and imaginary parts of impurity green's function
     real(dp), allocatable :: re_err(:,:,:)
     real(dp), allocatable :: im_err(:,:,:)

! allocate memory
     allocate(re_err(mfreq,norbs,norbs))
     allocate(im_err(mfreq,norbs,norbs))

! initialize re_err and im_err
     re_err = zero
     im_err = zero

! initialize grnf_mpi and grnf_err
     grnf_mpi = czero
     grnf_err = czero

! build grnf_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(grnf, grnf_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     grnf_mpi = grnf

# endif /* MPI */

! calculate the average
     grnf_mpi = grnf_mpi / real(nprocs)

! build grnf_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(( real(grnf - grnf_mpi))**2, re_err)
     call mp_allreduce((aimag(grnf - grnf_mpi))**2, im_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         re_err = sqrt( re_err / real( nprocs * ( nprocs - 1 ) ) )
         im_err = sqrt( im_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

! construct the final grnf_err
     grnf_err = re_err + im_err * czi

! deallocate memory
     deallocate(re_err)
     deallocate(im_err)

     return
  end subroutine ctqmc_reduce_grnf

!!>>> ctqmc_reduce_hist: reduce the hist from all children processes
  subroutine ctqmc_reduce_hist(hist_mpi, hist_err)
     use constants, only : dp, zero
     use mmpi, only : mp_allreduce, mp_barrier

     use control, only : mkink
     use control, only : nprocs
     use context, only : hist

     implicit none

! external arguments
! histogram for perturbation expansion series
     real(dp), intent(out) :: hist_mpi(mkink)
     real(dp), intent(out) :: hist_err(mkink)

! initialize hist_mpi and hist_err
     hist_mpi = zero
     hist_err = zero

! build hist_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(hist, hist_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     hist_mpi = hist

# endif /* MPI */

! calculate the average
     hist_mpi = hist_mpi / real(nprocs)

! build hist_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce((hist - hist_mpi)**2, hist_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         hist_err = sqrt( hist_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

     return
  end subroutine ctqmc_reduce_hist

!!>>> ctqmc_reduce_prob: reduce the prob from all children processes
  subroutine ctqmc_reduce_prob(prob_mpi, prob_err)
     use constants, only : dp, zero
     use mmpi, only : mp_allreduce, mp_barrier

     use control, only : ncfgs
     use control, only : nprocs
     use context, only : prob

     implicit none

! external arguments
! probability of atomic states
     real(dp), intent(out) :: prob_mpi(ncfgs)
     real(dp), intent(out) :: prob_err(ncfgs)

! initialize prob_mpi and prob_err
     prob_mpi = zero
     prob_err = zero

! build prob_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(prob, prob_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     prob_mpi = prob

# endif /* MPI */

! calculate the average
     prob_mpi = prob_mpi / real(nprocs)

! build prob_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce((prob - prob_mpi)**2, prob_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         prob_err = sqrt( prob_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

     return
  end subroutine ctqmc_reduce_prob

!!>>> ctqmc_reduce_nmat: reduce the nmat and nnmat from all children processes
  subroutine ctqmc_reduce_nmat(nmat_mpi, nnmat_mpi, nmat_err, nnmat_err)
     use constants, only : dp, zero
     use mmpi, only : mp_allreduce, mp_barrier

     use control, only : norbs
     use control, only : nprocs
     use context, only : nmat, nnmat

     implicit none

! external arguments
! occupation number matrix
     real(dp), intent(out) :: nmat_mpi(norbs)
     real(dp), intent(out) :: nmat_err(norbs)

! double occupation number matrix
     real(dp), intent(out) :: nnmat_mpi(norbs,norbs)
     real(dp), intent(out) :: nnmat_err(norbs,norbs)

! initialize nmat_mpi and nnmat_mpi, nmat_err and nnmat_err
     nmat_mpi = zero
     nnmat_mpi = zero

     nmat_err = zero
     nnmat_err = zero

! build nmat_mpi and nnmat_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(nmat, nmat_mpi)
     call mp_allreduce(nnmat, nnmat_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     nmat_mpi = nmat
     nnmat_mpi = nnmat

# endif /* MPI */

! calculate the average
     nmat_mpi = nmat_mpi / real(nprocs)
     nnmat_mpi = nnmat_mpi / real(nprocs)

! build nmat_err and nnmat_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce((nmat - nmat_mpi)**2, nmat_err)
     call mp_allreduce((nnmat - nnmat_mpi)**2, nnmat_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         nmat_err = sqrt( nmat_err / real( nprocs * ( nprocs - 1 ) ) )
         nnmat_err = sqrt( nnmat_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

     return
  end subroutine ctqmc_reduce_nmat

!!>>> ctqmc_reduce_kmat: reduce the kmat and kkmat from all children processes
  subroutine ctqmc_reduce_kmat(kmat_mpi, kkmat_mpi, kmat_err, kkmat_err)
     use constants, only : dp, zero
     use mmpi, only : mp_allreduce, mp_barrier

     use control, only : norbs
     use control, only : nprocs
     use context, only : kmat, kkmat

     implicit none

! external arguments
! number of operators
     real(dp), intent(out) :: kmat_mpi(norbs)
     real(dp), intent(out) :: kmat_err(norbs)

! square of number of operators
     real(dp), intent(out) :: kkmat_mpi(norbs,norbs)
     real(dp), intent(out) :: kkmat_err(norbs,norbs)

! initialize kmat_mpi and kkmat_mpi, kmat_err and kkmat_err
     kmat_mpi = zero
     kkmat_mpi = zero

     kmat_err = zero
     kkmat_err = zero

! build kmat_mpi and kkmat_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(kmat, kmat_mpi)
     call mp_allreduce(kkmat, kkmat_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     kmat_mpi = kmat
     kkmat_mpi = kkmat

# endif /* MPI */

! calculate the average
     kmat_mpi = kmat_mpi / real(nprocs)
     kkmat_mpi = kkmat_mpi / real(nprocs)

! build kmat_err and kkmat_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce((kmat - kmat_mpi)**2, kmat_err)
     call mp_allreduce((kkmat - kkmat_mpi)**2, kkmat_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         kmat_err = sqrt( kmat_err / real( nprocs * ( nprocs - 1 ) ) )
         kkmat_err = sqrt( kkmat_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

     return
  end subroutine ctqmc_reduce_kmat

!!>>> ctqmc_reduce_lmat: reduce the lmat, rmat, and lrmat from all children processes
  subroutine ctqmc_reduce_lmat(lmat_mpi, rmat_mpi, lrmat_mpi, lmat_err, rmat_err, lrmat_err)
     use constants, only : dp, zero
     use mmpi, only : mp_allreduce, mp_barrier

     use control, only : norbs
     use control, only : nprocs
     use context, only : lmat, rmat, lrmat

     implicit none

! external arguments
! number of operators at left half axis
     real(dp), intent(out) :: lmat_mpi(norbs)
     real(dp), intent(out) :: lmat_err(norbs)

! number of operators at right half axis
     real(dp), intent(out) :: rmat_mpi(norbs)
     real(dp), intent(out) :: rmat_err(norbs)

! used to evaluate fidelity susceptibility
     real(dp), intent(out) :: lrmat_mpi(norbs,norbs)
     real(dp), intent(out) :: lrmat_err(norbs,norbs)

! initialize lmat_mpi, rmat_mpi, and lrmat_mpi
! initialize lmat_err, rmat_err, and lrmat_err
     lmat_mpi = zero
     rmat_mpi = zero
     lrmat_mpi = zero

     lmat_err = zero
     rmat_err = zero
     lrmat_err = zero

! build lmat_mpi, rmat_mpi, and lrmat_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(lmat, lmat_mpi)
     call mp_allreduce(rmat, rmat_mpi)
     call mp_allreduce(lrmat, lrmat_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     lmat_mpi = lmat
     rmat_mpi = rmat
     lrmat_mpi = lrmat

# endif /* MPI */

! calculate the average
     lmat_mpi = lmat_mpi / real(nprocs)
     rmat_mpi = rmat_mpi / real(nprocs)
     lrmat_mpi = lrmat_mpi / real(nprocs)

! build lmat_err, rmat_err, and lrmat_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce((lmat - lmat_mpi)**2, lmat_err)
     call mp_allreduce((rmat - rmat_mpi)**2, rmat_err)
     call mp_allreduce((lrmat - lrmat_mpi)**2, lrmat_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         lmat_err = sqrt( lmat_err / real( nprocs * ( nprocs - 1 ) ) )
         rmat_err = sqrt( rmat_err / real( nprocs * ( nprocs - 1 ) ) )
         lrmat_err = sqrt( lrmat_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

     return
  end subroutine ctqmc_reduce_lmat

!!>>> ctqmc_reduce_szpw: reduce the szpow from all children processes
  subroutine ctqmc_reduce_szpw(szpow_mpi, szpow_err)
     use constants, only : dp, zero
     use mmpi, only : mp_allreduce, mp_barrier

     use control, only : norbs
     use control, only : nprocs
     use context, only : szpow

     implicit none

! external arguments
! powers of local magnetization, orbital-resolved
     real(dp), intent(out) :: szpow_mpi(4,norbs)
     real(dp), intent(out) :: szpow_err(4,norbs)

! initialize szpow_mpi and szpow_err
     szpow_mpi = zero
     szpow_err = zero

! build szpow_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(szpow, szpow_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     szpow_mpi = szpow

# endif /* MPI */

! calculate the average
     szpow_mpi = szpow_mpi / real(nprocs)

! build szpow_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce((szpow - szpow_mpi)**2, szpow_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         szpow_err = sqrt( szpow_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

     return
  end subroutine ctqmc_reduce_szpw

!!>>> ctqmc_reduce_schi: reduce the schi and sschi from all children processes
  subroutine ctqmc_reduce_schi(schi_mpi, sschi_mpi, schi_err, sschi_err)
     use constants, only : dp, zero
     use mmpi, only : mp_allreduce, mp_barrier

     use control, only : nband
     use control, only : ntime
     use control, only : nprocs
     use context, only : schi, sschi

     implicit none

! external arguments
! spin-spin correlation function, totally-averaged
     real(dp), intent(out) :: schi_mpi(ntime)
     real(dp), intent(out) :: schi_err(ntime)

! spin-spin correlation function, orbital-resolved
     real(dp), intent(out) :: sschi_mpi(ntime,nband)
     real(dp), intent(out) :: sschi_err(ntime,nband)

! initialize schi_mpi and sschi_mpi, schi_err and sschi_err
     schi_mpi = zero
     sschi_mpi = zero

     schi_err = zero
     sschi_err = zero

! build schi_mpi and sschi_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(schi, schi_mpi)
     call mp_allreduce(sschi, sschi_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     schi_mpi = schi
     sschi_mpi = sschi

# endif /* MPI */

! calculate the average
     schi_mpi = schi_mpi / real(nprocs)
     sschi_mpi = sschi_mpi / real(nprocs)

! build schi_err and sschi_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce((schi - schi_mpi)**2, schi_err)
     call mp_allreduce((sschi - sschi_mpi)**2, sschi_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         schi_err = sqrt( schi_err / real( nprocs * ( nprocs - 1 ) ) )
         sschi_err = sqrt( sschi_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

     return
  end subroutine ctqmc_reduce_schi

!!>>> ctqmc_reduce_sfom: reduce the ssfom from all children processes
  subroutine ctqmc_reduce_sfom(ssfom_mpi, ssfom_err)
     use constants, only : dp, zero
     use mmpi, only : mp_allreduce, mp_barrier

     use control, only : nband
     use control, only : nbfrq
     use control, only : nprocs
     use context, only : ssfom

     implicit none

! external arguments
! spin-spin correlation function, orbital-resolved
     real(dp), intent(out) :: ssfom_mpi(nbfrq,nband)
     real(dp), intent(out) :: ssfom_err(nbfrq,nband)

! initialize ssfom_mpi and ssfom_err
     ssfom_mpi = zero
     ssfom_err = zero

! build ssfom_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(ssfom, ssfom_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     ssfom_mpi = ssfom

# endif /* MPI */

! calculate the average
     ssfom_mpi = ssfom_mpi / real(nprocs)

! build ssfom_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce((ssfom - ssfom_mpi)**2, ssfom_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         ssfom_err = sqrt( ssfom_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

     return
  end subroutine ctqmc_reduce_sfom

!!>>> ctqmc_reduce_ochi: reduce the ochi and oochi from all children processes
  subroutine ctqmc_reduce_ochi(ochi_mpi, oochi_mpi, ochi_err, oochi_err)
     use constants, only : dp, zero
     use mmpi, only : mp_allreduce, mp_barrier

     use control, only : norbs
     use control, only : ntime
     use control, only : nprocs
     use context, only : ochi, oochi

     implicit none

! external arguments
! orbital-orbital correlation function, totally-averaged
     real(dp), intent(out) :: ochi_mpi(ntime)
     real(dp), intent(out) :: ochi_err(ntime)

! orbital-orbital correlation function, orbital-resolved
     real(dp), intent(out) :: oochi_mpi(ntime,norbs,norbs)
     real(dp), intent(out) :: oochi_err(ntime,norbs,norbs)

! initialize ochi_mpi and oochi_mpi, ochi_err and oochi_err
     ochi_mpi = zero
     oochi_mpi = zero

     ochi_err = zero
     oochi_err = zero

! build ochi_mpi and oochi_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(ochi, ochi_mpi)
     call mp_allreduce(oochi, oochi_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     ochi_mpi = ochi
     oochi_mpi = oochi

# endif /* MPI */

! calculate the average
     ochi_mpi = ochi_mpi / real(nprocs)
     oochi_mpi = oochi_mpi / real(nprocs)

! build ochi_err and oochi_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce((ochi - ochi_mpi)**2, ochi_err)
     call mp_allreduce((oochi - oochi_mpi)**2, oochi_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         ochi_err = sqrt( ochi_err / real( nprocs * ( nprocs - 1 ) ) )
         oochi_err = sqrt( oochi_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

     return
  end subroutine ctqmc_reduce_ochi

!!>>> ctqmc_reduce_ofom: reduce the oofom from all children processes
  subroutine ctqmc_reduce_ofom(oofom_mpi, oofom_err)
     use constants, only : dp, zero
     use mmpi, only : mp_allreduce, mp_barrier

     use control, only : norbs
     use control, only : nbfrq
     use control, only : nprocs
     use context, only : oofom

     implicit none

! external arguments
! orbital-orbital correlation function, orbital-resolved
     real(dp), intent(out) :: oofom_mpi(nbfrq,norbs,norbs)
     real(dp), intent(out) :: oofom_err(nbfrq,norbs,norbs)

! initialize oofom_mpi and oofom_err
     oofom_mpi = zero
     oofom_err = zero

! build oofom_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(oofom, oofom_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     oofom_mpi = oofom

# endif /* MPI */

! calculate the average
     oofom_mpi = oofom_mpi / real(nprocs)

! build oofom_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce((oofom - oofom_mpi)**2, oofom_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         oofom_err = sqrt( oofom_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

     return
  end subroutine ctqmc_reduce_ofom

!!>>> ctqmc_reduce_twop: reduce the g2_re_mpi and g2_im_mpi from all
!!>>> children processes
  subroutine ctqmc_reduce_twop(g2_re_mpi, g2_im_mpi)
     use constants, only : dp, zero
     use mmpi, only : mp_allreduce, mp_barrier

     use control, only : norbs
     use control, only : nffrq, nbfrq
     use control, only : nprocs
     use context, only : g2_re, g2_im

     implicit none

! external arguments
! two-particle green's function, real part
     real(dp), intent(out) :: g2_re_mpi(nffrq,nffrq,nbfrq,norbs,norbs)

! two-particle green's function, imaginary part
     real(dp), intent(out) :: g2_im_mpi(nffrq,nffrq,nbfrq,norbs,norbs)

! initialize g2_re_mpi and g2_im_mpi
     g2_re_mpi = zero
     g2_im_mpi = zero

! build g2_re_mpi and g2_im_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(g2_re, g2_re_mpi)
     call mp_allreduce(g2_im, g2_im_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     g2_re_mpi = g2_re
     g2_im_mpi = g2_im

# endif /* MPI */

! calculate the average
     g2_re_mpi = g2_re_mpi / real(nprocs)
     g2_im_mpi = g2_im_mpi / real(nprocs)

     return
  end subroutine ctqmc_reduce_twop

!!>>> ctqmc_reduce_vrtx: reduce the h2_re_mpi and h2_im_mpi from all
!!>>> children processes
  subroutine ctqmc_reduce_vrtx(h2_re_mpi, h2_im_mpi)
     use constants, only : dp, zero
     use mmpi, only : mp_allreduce, mp_barrier

     use control, only : norbs
     use control, only : nffrq, nbfrq
     use control, only : nprocs
     use context, only : h2_re, h2_im

     implicit none

! external arguments
! two-particle green's function, real part
     real(dp), intent(out) :: h2_re_mpi(nffrq,nffrq,nbfrq,norbs,norbs)

! two-particle green's function, imaginary part
     real(dp), intent(out) :: h2_im_mpi(nffrq,nffrq,nbfrq,norbs,norbs)

! initialize h2_re_mpi and h2_im_mpi
     h2_re_mpi = zero
     h2_im_mpi = zero

! build h2_re_mpi and h2_im_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(h2_re, h2_re_mpi)
     call mp_allreduce(h2_im, h2_im_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     h2_re_mpi = h2_re
     h2_im_mpi = h2_im

# endif /* MPI */

! calculate the average
     h2_re_mpi = h2_re_mpi / real(nprocs)
     h2_im_mpi = h2_im_mpi / real(nprocs)

     return
  end subroutine ctqmc_reduce_vrtx

!!>>> ctqmc_reduce_pair: reduce the ps_re_mpi and ps_im_mpi from all
!!>>> children processes
  subroutine ctqmc_reduce_pair(ps_re_mpi, ps_im_mpi)
     use constants, only : dp, zero
     use mmpi, only : mp_allreduce, mp_barrier

     use control, only : norbs
     use control, only : nffrq, nbfrq
     use control, only : nprocs
     use context, only : ps_re, ps_im

     implicit none

! external arguments
! particle-particle pair susceptibility, real part
     real(dp), intent(out) :: ps_re_mpi(nffrq,nffrq,nbfrq,norbs,norbs)

! particle-particle pair susceptibility, imaginary part
     real(dp), intent(out) :: ps_im_mpi(nffrq,nffrq,nbfrq,norbs,norbs)

! initialize ps_re_mpi and ps_im_mpi
     ps_re_mpi = zero
     ps_im_mpi = zero

! build ps_re_mpi and ps_im_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(ps_re, ps_re_mpi)
     call mp_allreduce(ps_im, ps_im_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     ps_re_mpi = ps_re
     ps_im_mpi = ps_im

# endif /* MPI */

! calculate the average
     ps_re_mpi = ps_re_mpi / real(nprocs)
     ps_im_mpi = ps_im_mpi / real(nprocs)

     return
  end subroutine ctqmc_reduce_pair

!!========================================================================
!!>>> symmetrize physical observables                                  <<<
!!========================================================================

!!>>> ctqmc_symm_nmat: symmetrize the nmat according to symm vector
  subroutine ctqmc_symm_nmat(symm, nmat)
     use constants, only : dp, zero, two

     use control, only : issun, isspn
     use control, only : nband, norbs

     implicit none

! external arguments
! symmetry vector
     integer, intent(in) :: symm(norbs)

! occupation number
     real(dp), intent(inout) :: nmat(norbs)

! local variables
! loop index over bands
     integer  :: ibnd
     integer  :: jbnd

! dummy variables
     real(dp) :: raux

! histogram vector
! note: it is NOT the global one
     integer  :: hist(norbs)

! build histogram
     hist = 0
     do ibnd=1,norbs
         hist(symm(ibnd)) = hist(symm(ibnd)) + 1
     enddo ! over ibnd={1,norbs} loop

! perform symmetrization for those orbitals which symm index are identity
     if ( issun == 2 ) then
         do ibnd=1,norbs
             if ( hist(ibnd) > 0 ) then         ! need to enforce symmetry
                 raux = zero

                 do jbnd=1,norbs                ! gather the data
                     if ( symm(jbnd) == ibnd ) then
                         raux = raux + nmat(jbnd)
                     endif ! back if ( symm(jbnd) == ibnd ) block
                 enddo ! over jbnd={1,norbs} loop

                 raux = raux / real(hist(ibnd)) ! calculate average value

                 do jbnd=1,norbs                ! setup it
                     if ( symm(jbnd) == ibnd ) then
                         nmat(jbnd) = raux
                     endif ! back if ( symm(jbnd) == ibnd ) block
                 enddo ! over jbnd={1,norbs} loop
             endif ! back if ( hist(ibnd) > 0 ) block
         enddo ! over ibnd={1,norbs} loop
     endif ! back if ( issun == 2 ) block

! symmetrize nmat over spin
     if ( isspn == 1 ) then
         do jbnd=1,nband
             raux = ( nmat(jbnd) + nmat(jbnd+nband) ) / two
             nmat(jbnd) = raux
             nmat(jbnd+nband) = raux
         enddo ! over jbnd={1,nband} loop
     endif ! back if ( isspn == 1 ) block

     return
  end subroutine ctqmc_symm_nmat

!!>>> ctqmc_symm_gtau: symmetrize the gtau according to symm vector
!!>>> only the diagonal elements are taken into considerations
  subroutine ctqmc_symm_gtau(symm, gtau)
     use constants, only : dp, zero, two

     use control, only : issun, isspn
     use control, only : nband, norbs
     use control, only : ntime

     implicit none

! external arguments
! symmetry vector
     integer, intent(in) :: symm(norbs)

! impurity green's function
     real(dp), intent(inout) :: gtau(ntime,norbs,norbs)

! local variables
! loop index over bands
     integer  :: ibnd
     integer  :: jbnd

! loop index over imaginary-time points
     integer  :: ktau

! dummy variables
     real(dp) :: raux

! histogram vector
! note: it is NOT the global one
     integer  :: hist(norbs)

! build histogram
     hist = 0
     do ibnd=1,norbs
         hist(symm(ibnd)) = hist(symm(ibnd)) + 1
     enddo ! over ibnd={1,norbs} loop

! perform symmetrization for those orbitals which symm index are identity
     if ( issun == 2 ) then
         do ktau=1,ntime
             do ibnd=1,norbs
                 if ( hist(ibnd) > 0 ) then         ! need to enforce symmetry
                     raux = zero

                     do jbnd=1,norbs                ! gather the data
                         if ( symm(jbnd) == ibnd ) then
                             raux = raux + gtau(ktau,jbnd,jbnd)
                         endif ! back if ( symm(jbnd) == ibnd ) block
                     enddo ! over jbnd={1,norbs} loop

                     raux = raux / real(hist(ibnd)) ! calculate average value

                     do jbnd=1,norbs                ! setup it
                         if ( symm(jbnd) == ibnd ) then
                             gtau(ktau,jbnd,jbnd) = raux
                         endif ! back if ( symm(jbnd) == ibnd ) block
                     enddo ! over jbnd={1,norbs} loop
                 endif ! back if ( hist(ibnd) > 0 ) block
             enddo ! over ibnd={1,norbs} loop
         enddo ! over ktau={1,ntime} loop
     endif ! back if ( issun == 2 ) block

! symmetrize gtau over spin
     if ( isspn == 1 ) then
         do ktau=1,ntime
             do jbnd=1,nband
                 raux = ( gtau(ktau,jbnd,jbnd) + gtau(ktau,jbnd+nband,jbnd+nband) ) / two
                 gtau(ktau,jbnd,jbnd) = raux
                 gtau(ktau,jbnd+nband,jbnd+nband) = raux
             enddo ! over jbnd={1,nband} loop
         enddo ! over ktau={1,ntime} loop
     endif ! back if ( isspn == 1 ) block

     return
  end subroutine ctqmc_symm_gtau

!!>>> ctqmc_symm_grnf: symmetrize the grnf according to symm vector
!!>>> only the diagonal elements are taken into considerations
  subroutine ctqmc_symm_grnf(symm, grnf)
     use constants, only : dp, two, czero

     use control, only : issun, isspn
     use control, only : nband, norbs
     use control, only : mfreq

     implicit none

! external arguments
! symmetry vector
     integer, intent(in) :: symm(norbs)

! impurity green's function
     complex(dp), intent(inout) :: grnf(mfreq,norbs,norbs)

! local variables
! loop index over bands
     integer :: ibnd
     integer :: jbnd

! loop index over matsubara frequencies
     integer :: kfrq

! dummy variables
     complex(dp) :: caux

! histogram vector
! note: it is NOT the global one
     integer :: hist(norbs)

! build histogram
     hist = 0
     do ibnd=1,norbs
         hist(symm(ibnd)) = hist(symm(ibnd)) + 1
     enddo ! over ibnd={1,norbs} loop

! perform symmetrization for those orbitals which symm index are identity
     if ( issun == 2 ) then
         do kfrq=1,mfreq
             do ibnd=1,norbs
                 if ( hist(ibnd) > 0 ) then         ! need to enforce symmetry
                     caux = czero

                     do jbnd=1,norbs                ! gather the data
                         if ( symm(jbnd) == ibnd ) then
                             caux = caux + grnf(kfrq,jbnd,jbnd)
                         endif ! back if ( symm(jbnd) == ibnd ) block
                     enddo ! over jbnd={1,norbs} loop

                     caux = caux / real(hist(ibnd)) ! calculate average value

                     do jbnd=1,norbs                ! setup it
                         if ( symm(jbnd) == ibnd ) then
                             grnf(kfrq,jbnd,jbnd) = caux
                         endif ! back if ( symm(jbnd) == ibnd ) block
                     enddo ! over jbnd={1,norbs} loop
                 endif ! back if ( hist(ibnd) > 0 ) block
             enddo ! over ibnd={1,norbs} loop
         enddo ! over kfrq={1,mfreq} loop
     endif ! back if ( issun == 2 ) block

! symmetrize grnf over spin
     if ( isspn == 1 ) then
         do kfrq=1,mfreq
             do jbnd=1,nband
                 caux = ( grnf(kfrq,jbnd,jbnd) + grnf(kfrq,jbnd+nband,jbnd+nband) ) / two
                 grnf(kfrq,jbnd,jbnd) = caux
                 grnf(kfrq,jbnd+nband,jbnd+nband) = caux
             enddo ! over jbnd={1,nband} loop
         enddo ! over kfrq={1,mfreq} loop
     endif ! back if ( isspn == 1 ) block

     return
  end subroutine ctqmc_symm_grnf

!!>>> ctqmc_smth_sigf: smooth impurity self-energy function in low
!!>>> frequency region
  subroutine ctqmc_smth_sigf(sigf)
     use constants, only : dp, czero

     use control, only : nfreq

     implicit none

! external arguments
! impurity self-energy function to be smoothen
     complex(dp), intent(inout) :: sigf(nfreq)

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k

! smooth radius
     integer  :: lrad
     integer  :: srad

! imaginary part of self-energy function
     real(dp) :: ti

! real part of self-energy function
     real(dp) :: tr

! dummy variables for addition
     complex(dp) :: saux

! dummy self-energy function
     complex(dp) :: stmp(nfreq)

! determine smooth radius
     lrad = nfreq / 4  ! large radius
     srad = nfreq / 16 ! small radius

! |---------|---------|----------------------|
! 1         lrad      2*lrad                 nfreq
! deal with [1,lrad], head part
     do k=1,lrad
         stmp(k) = sigf(k)
     enddo ! over k={1,lrad} loop

! deal with [lrad+1,2*lrad], intermediate part
     do i=1,lrad
         k = lrad + i
         saux = czero
         do j=-srad,srad
             saux = saux + sigf(k+j)
         enddo ! over j={-srad,srad} loop
         stmp(k) = saux / real(2 * srad + 1)
         stmp(k) = ( (lrad - i) * sigf(k) + i * stmp(k) ) / real(lrad)
     enddo ! over i={1,lrad} loop

! deal with [nfreq-2*lrad+1,nfreq], tail part
     do k=nfreq-2*lrad+1,nfreq
         tr =  real( stmp(nfreq-2*lrad) )
         ti = aimag( stmp(nfreq-2*lrad) ) * real(nfreq - 2 * lrad) / real(k)
         stmp(k) = dcmplx(tr,ti)
     enddo ! over k={nfreq-2*lrad+1,nfreq} loop

! copy stmp to sigf
     do k=1,nfreq
         sigf(k) = stmp(k)
     enddo ! over k={1,nfreq} loop

     return
  end subroutine ctqmc_smth_sigf

!!========================================================================
!!>>> postprocess physical observables                                 <<<
!!========================================================================

!!>>> ctqmc_make_gtau: build imaginary green's function using orthogonal
!!>>> polynomial representation
  subroutine ctqmc_make_gtau(tmesh, gtau, gaux)
     use constants, only : dp, zero, one, two, pi

     use control, only : isort
     use control, only : norbs
     use control, only : lemax, legrd, chmax, chgrd
     use control, only : ntime
     use control, only : beta
     use context, only : ppleg, qqche

     implicit none

! external arguments
! imaginary time mesh
     real(dp), intent(in)  :: tmesh(ntime)

! impurity green's function/orthogonal polynomial coefficients
     real(dp), intent(in)  :: gtau(ntime,norbs,norbs)

! calculated impurity green's function
     real(dp), intent(out) :: gaux(ntime,norbs,norbs)

! local parameters
! scheme of integral kernel used to damp the Gibbs oscillation
! damp = 0, Dirichlet   mode
! damp = 1, Jackson     mode, preferred
! damp = 2, Lorentz     mode
! damp = 3, Fejer       mode
! damp = 4, Wang-Zunger mode
     integer, parameter :: damp = 0

! local variables
! loop index
     integer  :: i
     integer  :: j

! loop index for legendre polynomial
     integer  :: fleg

! loop index for chebyshev polynomial
     integer  :: fche

! index for imaginary time \tau
     integer  :: curr

! interval for imaginary time slice
     real(dp) :: step

! dummy variables
     real(dp) :: raux

! initialize gaux
     gaux = zero

! select calculation method
     select case ( isort )

         case (1, 4)
             call cat_make_gtau1()

         case (2, 5)
             call cat_make_gtau2()

         case (3, 6)
             call cat_make_gtau3()

     end select

     return

  contains

!!>>> cat_make_kpm: build the integral kernel function
  subroutine cat_make_kpm(kdim, kern)
     implicit none

! external arguments
! dimension of integral kernel function
     integer, intent(in)   :: kdim

! integral kernel function
     real(dp), intent(out) :: kern(kdim)

! local variables
! loop index
     integer :: kcur

     kern = zero
     do kcur=1,kdim
         select case ( damp )

! Dirichlet mode
             case (0)
                 kern(kcur) = one

! Jackson mode
             case (1)
                 i = kcur - 1
                 curr = kdim + 1
                 raux = pi * i / curr
                 kern(kcur) = ( (curr-i) * cos(raux) + sin(raux) / tan(pi/curr) ) / curr

! Lorentz mode
             case (2)
                 kern(kcur) = sinh( one - (kcur - one) / real( kdim ) ) / sinh(one)

! Fejer mode
             case (3)
                 kern(kcur) = one - ( kcur - one ) / real( kdim )

! Wang-Zunger mode
             case (4)
                 kern(kcur) = exp( - ( (kcur - one) / real( kdim ) )**4 )

         end select
     enddo ! over kcur={1,kdim} loop

     return
  end subroutine cat_make_kpm

!!>>> cat_make_gtau1: build impurity green's function using normal
!!>>> representation
  subroutine cat_make_gtau1()
     implicit none

     raux = real(ntime) / (beta * beta)
     do i=1,norbs
         do j=1,ntime
             gaux(j,i,i) = gtau(j,i,i) * raux
         enddo ! over j={1,ntime} loop
     enddo ! over i={1,norbs} loop

     return
  end subroutine cat_make_gtau1

!!>>> cat_make_gtau2: build impurity green's function using legendre
!!>>> polynomial representation
  subroutine cat_make_gtau2()
     implicit none

! integral kernel
     real(dp) :: ker1(lemax)

! build kernel function at first
     ker1 = one; call cat_make_kpm(lemax, ker1)

! reconstruct green's function
     step = real(legrd - 1) / two
     do i=1,norbs
         do j=1,ntime
             raux = two * tmesh(j) / beta
             curr = nint(raux * step) + 1
             do fleg=1,lemax
                 raux = sqrt(two * fleg - 1) / (beta * beta) * ker1(fleg)
                 gaux(j,i,i) = gaux(j,i,i) + raux * gtau(fleg,i,i) * ppleg(curr,fleg)
             enddo ! over fleg={1,lemax} loop
         enddo ! over j={1,ntime} loop
     enddo ! over i={1,norbs} loop

     return
  end subroutine cat_make_gtau2

!!>>> cat_make_gtau3: build impurity green's function using chebyshev
!!>>> polynomial representation
  subroutine cat_make_gtau3()
     implicit none

! integral kernel
     real(dp) :: ker2(chmax)

! build kernel function at first
     ker2 = one; call cat_make_kpm(chmax, ker2)

! reconstruct green's function
     step = real(chgrd - 1) / two
     do i=1,norbs
         do j=1,ntime
             raux = two * tmesh(j) / beta
             curr = nint(raux * step) + 1
             do fche=1,chmax
                 raux = two / (beta * beta) * ker2(fche)
                 gaux(j,i,i) = gaux(j,i,i) + raux * gtau(fche,i,i) * qqche(curr,fche)
             enddo ! over fche={1,chmax} loop
         enddo ! over j={1,ntime} loop
     enddo ! over i={1,norbs} loop

     return
  end subroutine cat_make_gtau3
  end subroutine ctqmc_make_gtau

!!>>> ctqmc_make_ftau: build auxiliary correlation function using
!!>>> orthogonal polynomial representation, F(\tau)
  subroutine ctqmc_make_ftau(tmesh, ftau, faux)
     use constants, only : dp, zero, two

     use control, only : isort
     use control, only : norbs
     use control, only : lemax, legrd, chmax, chgrd
     use control, only : ntime
     use control, only : beta
     use context, only : ppleg, qqche

     implicit none

! external arguments
! imaginary time mesh
     real(dp), intent(in)  :: tmesh(ntime)

! auxiliary correlation function/orthogonal polynomial coefficients
     real(dp), intent(in)  :: ftau(ntime,norbs,norbs)

! calculated auxiliary correlation function
     real(dp), intent(out) :: faux(ntime,norbs,norbs)

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k

! loop index for legendre polynomial
     integer  :: fleg

! loop index for chebyshev polynomial
     integer  :: fche

! index for imaginary time \tau
     integer  :: curr

! interval for imaginary time slice
     real(dp) :: step

! dummy variables
     real(dp) :: raux

! initialize faux
     faux = zero

! select calculation method
     select case ( isort )

         case (4)
             call cat_make_ftau1()

         case (5)
             call cat_make_ftau2()

         case (6)
             call cat_make_ftau3()

     end select

     return

  contains

!!>>> cat_make_ftau1: build auxiliary correlation function using normal
!!>>> representation
  subroutine cat_make_ftau1()
     implicit none

     raux = real(ntime) / (beta * beta)
     do i=1,norbs
         do j=1,norbs
             do k=1,ntime
                 faux(k,j,i) = ftau(k,j,i) * raux
             enddo ! over k={1,ntime} loop
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

     return
  end subroutine cat_make_ftau1

!!>>> cat_make_ftau2: build auxiliary correlation function using legendre
!!>>> polynomial representation
  subroutine cat_make_ftau2()
     implicit none

     step = real(legrd - 1) / two
     do i=1,norbs
         do j=1,norbs
             do k=1,ntime
                 raux = two * tmesh(k) / beta
                 curr = nint(raux * step) + 1
                 do fleg=1,lemax
                     raux = sqrt(two * fleg - 1) / (beta * beta)
                     faux(k,j,i) = faux(k,j,i) + raux * ftau(fleg,j,i) * ppleg(curr,fleg)
                 enddo ! over fleg={1,lemax} loop
             enddo ! over k={1,ntime} loop
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

     return
  end subroutine cat_make_ftau2

!!>>> cat_make_ftau3: build auxiliary correlation function using chebyshev
!!>>> polynomial representation
  subroutine cat_make_ftau3()
     implicit none

     step = real(chgrd - 1) / two
     do i=1,norbs
         do j=1,norbs
             do k=1,ntime
                 raux = two * tmesh(k) / beta
                 curr = nint(raux * step) + 1
                 raux = two / (beta * beta)
                 do fche=1,chmax
                     faux(k,j,i) = faux(k,j,i) + raux * ftau(fche,j,i) * qqche(curr,fche)
                 enddo ! over fche={1,chmax} loop
             enddo ! over k={1,ntime} loop
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

     return
  end subroutine cat_make_ftau3
  end subroutine ctqmc_make_ftau

!!========================================================================
!!>>> build prefactor for improved estimator                           <<<
!!========================================================================

!!>>> ctqmc_make_iret: to calculate the integral I(\tau_end) which is very
!!>>> important when retarded interaction is included
  subroutine ctqmc_make_iret(time, iret)
     use constants, only : dp, zero, two

     use control, only : norbs
     use context, only : index_s, index_e, time_s, time_e
     use context, only : rank

     implicit none

! external arguments
! imaginary time point, in principle, it is \tau_end
     real(dp), intent(in)  :: time

! integral value for I(\tau_end)
     real(dp), intent(out) :: iret

! local variables
! loop index for start and end points
     integer  :: it

! loop index for flavor channel
     integer  :: flvr

! length betweem two time points
     real(dp) :: dtau
     real(dp) :: daux

! calculate integral I(\tau_end), the equation is
! Eq. (39) in Phys. Rev. B 89, 235128 (2014)
     iret = zero

     do flvr=1,norbs
         do it=1,rank(flvr)
! contribution from create operators
             dtau = time_s( index_s(it, flvr), flvr ) - time
             if ( dtau >= zero ) then
                 call cat_weight_kernel(2, +dtau, daux)
                 iret = iret + daux
             else
                 call cat_weight_kernel(2, -dtau, daux)
                 iret = iret - daux
             endif ! back if ( dtau >= zero ) block

! contribution from destroy operators
             dtau = time_e( index_e(it, flvr), flvr ) - time
             if ( dtau >= zero ) then
                 call cat_weight_kernel(2, +dtau, daux)
                 iret = iret - daux
             else
                 call cat_weight_kernel(2, -dtau, daux)
                 iret = iret + daux
             endif ! back if ( dtau >= zero ) block
         enddo ! over it={1,rank(flvr)} loop
     enddo ! over flvr={1,norbs} loop

! add additional term
     call cat_weight_kernel(2, zero, daux)
     iret = -iret - two * daux

     return
  end subroutine ctqmc_make_iret

!!>>> ctqmc_make_pref: to calculate the prefactor used by the improved
!!>>> estimator for self-energy function and vertex function
  subroutine ctqmc_make_pref()
     use constants, only : dp, zero, half

     use control, only : isscr
     use control, only : norbs
     use context, only : index_e, time_e
     use context, only : rank, pref, uumat

     implicit none

! local variables
! loop index for start and end points
     integer  :: it

! loop index for flavor channel
     integer  :: flvr

! loop index for colour channel
     integer  :: clur

! occupation number at \tau_end
     real(dp) :: occu

! integral value for I(\tau_end)
     real(dp) :: iret

! if it is holstein-hubbard model, the improved estimator algorithm
! can not be used now
     call s_assert2(isscr /= 2,'sorry, isscr = 2 is not compatible with improved estimator')

     do flvr=1,norbs
         do it=1,rank(flvr)

! reset the prefactor
             pref(it,flvr) = zero

! calculate normal contribution
             do clur=1,norbs
                 call cat_occupy_status(clur, time_e( index_e(it, flvr), flvr ), occu)
                 pref(it,flvr) = pref(it,flvr) + half * ( uumat(flvr,clur) + uumat(clur,flvr) ) * occu
             enddo ! over clur={1,norbs} loop

! if retarded interaction is considered, we have to include the
! contribution from I(\tau_end)
             if ( isscr > 1 ) then
                 call ctqmc_make_iret(time_e( index_e(it, flvr), flvr ), iret)
                 pref(it,flvr) = pref(it,flvr) + iret
             endif ! back if ( isscr > 1 ) block
         enddo ! over it={1,rank(flvr)} loop
     enddo ! over flvr={1,norbs} loop

     return
  end subroutine ctqmc_make_pref

!!========================================================================
!!>>> build auxiliary two-particle related variables                   <<<
!!========================================================================

!!>>> ctqmc_make_prod: try to calculate the product of matsubara
!!>>> frequency exponents exp(i \omega_n \tau)
  subroutine ctqmc_make_prod(flvr, nfaux, mrank, caux1, caux2)
     use constants, only : dp, two, pi, czi

     use control, only : nffrq
     use control, only : beta
     use context, only : index_s, index_e, time_s, time_e
     use context, only : rank

     implicit none

! external arguments
! current flavor channel
     integer, intent(in) :: flvr

! combination of nffrq and nbfrq
     integer, intent(in) :: nfaux

! maximum number of operators in different flavor channels
     integer, intent(in) :: mrank

! matsubara frequency exponents for create operators
     complex(dp), intent(out) :: caux1(nfaux,mrank)

! matsubara frequency exponents for destroy operators
     complex(dp), intent(out) :: caux2(nfaux,mrank)

! local variables
! loop indices for start and end points
     integer  :: is
     integer  :: ie

! imaginary time for start and end points
     real(dp) :: taus
     real(dp) :: taue

! for create operators
     do is=1,rank(flvr)
         taus = time_s( index_s(is, flvr), flvr )
         caux1(:,is) = exp(-two * czi * pi * taus / beta)
         call s_cumprod_z(nfaux, caux1(:,is), caux1(:,is))
         caux1(:,is) = caux1(:,is) * exp(+(nffrq + 1) * czi * pi * taus / beta)
     enddo ! over is={1,rank(flvr)} loop

! for destroy operators
     do ie=1,rank(flvr)
         taue = time_e( index_e(ie, flvr), flvr )
         caux2(:,ie) = exp(+two * czi * pi * taue / beta)
         call s_cumprod_z(nfaux, caux2(:,ie), caux2(:,ie))
         caux2(:,ie) = caux2(:,ie) * exp(-(nffrq + 1) * czi * pi * taue / beta)
     enddo ! over ie={1,rank(flvr)} loop

     return
  end subroutine ctqmc_make_prod

!!========================================================================
!!>>> build self-energy function                                       <<<
!!========================================================================

!!>>> ctqmc_make_hub1: build atomic green's function and self-energy
!!>>> function using improved Hubbard-I approximation, and then make
!!>>> interpolation for self-energy function between low frequency QMC
!!>>> data and high frequency Hubbard-I approximation data, the full
!!>>> impurity green's function can be obtained by using dyson's equation
!!>>> finally
  subroutine ctqmc_make_hub1()
     use constants, only : dp, zero, one, two, czi, czero

     use control, only : norbs, ncfgs
     use control, only : mfreq
     use control, only : nfreq
     use control, only : mune
     use control, only : myid, master
     use context, only : rmesh
     use context, only : prob, nmat
     use context, only : eimp, uumat
     use context, only : grnf
     use context, only : hybf
     use context, only : sig2

     implicit none

! local parameters
! maximum allowable number of non-zero elements in F matrix
     integer, parameter :: nzero = 1024

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k
     integer  :: m
     integer  :: n

! dummy integer variables, used to build F matrix
     integer  :: start
     integer  :: value
     integer  :: permute

! dummy real variables, used to interpolate self-energy function
     real(dp) :: ob, oe
     real(dp) :: d0, d1
     real(dp) :: shift

! dummy complex variables, used to interpolate self-energy function
     complex(dp) :: cb, ce
     complex(dp) :: sinf

! dummy atomic states: alpha, beta, gamma
     integer  :: sa(norbs)
     integer  :: sb(norbs)
     integer  :: sc(norbs)

! atomic basis sets
     integer  :: basis(ncfgs,norbs)

! F matrix, < alpha | f_{n} | beta >
     integer  :: fcounter(norbs)
     integer  :: fa(nzero,norbs)
     integer  :: fb(nzero,norbs)
     integer  :: fv(nzero,norbs)

! eigenvalues for local hmailtonian
     real(dp) :: eaux(ncfgs)

! dummy imurity green's function: G^{-1}
     complex(dp) :: gaux(norbs,norbs)

! atomic green's function and self-energy function in Hubbard-I approximation
     complex(dp) :: ghub(mfreq,norbs)
     complex(dp) :: shub(mfreq,norbs)

! evaluate the shift for the Coulomb interaction and chemical potential
! if the retarded interaction is used
     call ctqmc_eval_shift(shift)

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

! evaluate atomic eigenvalues directly
     eaux = zero
     do i=1,ncfgs
         do j=1,norbs
             eaux(i) = eaux(i) + ( eimp(j) - mune ) * basis(i,j)
         enddo ! over j={1,norbs} loop
         do j=1,norbs-1
             do k=j+1,norbs
                 if ( basis(i,j) == 1 .and. basis(i,k) == 1 ) then
                     eaux(i) = eaux(i) + uumat(j,k)
                 endif ! back if ( basis(i,j) == 1 .and. basis(i,k) == 1 ) block
             enddo ! over k={j+1,norbs} loop
         enddo ! over j={1,norbs-1} loop
     enddo ! over i={1,ncfgs} loop

! build F matrix < alpha | f_{n} | beta >
! note 1: to save the memory and accelerate the computation, we only store
! the non-zero element of F matrix
! note 2: it is crucial to check whether the number of non-zero elements
! exceed limit (nzero)
     fcounter = 0
     alpha_loop: do i=1,ncfgs
         sa = basis(i,:)
         beta_loop: do j=1,ncfgs
             sb = basis(j,:)

             orbital_loop: do m=1,norbs
                 sc = sb

                 if ( sc(m) == 1 ) then
                     permute = 1
                     do n=1,m-1
                         if ( sc(n) == 1 ) permute = -permute
                     enddo ! over n={1,m-1} loop
                     sc(m) = 0

                     value = 1
                     do n=1,norbs
                         if ( sa(n) /= sc(n) ) value = 0
                     enddo ! over n={1,norbs} loop
                     value = value * permute
                 else
                     value = 0
                 endif ! back if ( sc(m) == 1 ) block

                 if ( value /= 0 ) then
                     fcounter(m) = fcounter(m) + 1
                     if ( fcounter(m) > nzero ) then
                         call s_print_error('ctqmc_make_hub1','non-zero elements exceed limit')
                     endif ! back if ( fcounter(m) > nzero ) block
                     fa(fcounter(m),m) = i
                     fb(fcounter(m),m) = j
                     fv(fcounter(m),m) = value
                 endif ! back if ( value /= 0 ) block
             enddo orbital_loop ! over m={1,norbs} loop

         enddo beta_loop ! over j={1,ncfgs} loop
     enddo alpha_loop ! over i={1,ncfgs} loop

! calculate atomic green's function using Hubbard-I approximation
     do i=1,norbs
         do k=1,mfreq
             ghub(k,i) = czero
             do m=1,fcounter(i)
                 ob = fv(m,i) * fv(m,i) * ( prob(fa(m,i)) + prob(fb(m,i)) )
                 cb = czi * rmesh(k) + eaux(fa(m,i)) - eaux(fb(m,i))
                 ghub(k,i) = ghub(k,i) + ob / cb
             enddo ! over m={1,fcounter(i)} loop
         enddo ! over k={1,mfreq} loop
     enddo ! over i={1,norbs} loop

! calculate atomic self-energy function using dyson's equation
     do i=1,norbs
         do k=1,mfreq
             shub(k,i) = czi * rmesh(k) + mune - eimp(i) - one / ghub(k,i)
             shub(k,i) = shub(k,i) - nmat(i) * shift
         enddo ! over k={1,mfreq} loop
     enddo ! over i={1,norbs} loop

! dump the ghub and shub, only for reference, only the master node can do it
     if ( myid == master ) then
         call ctqmc_dump_hub1(rmesh, ghub, shub)
     endif ! back if ( myid == master ) block

! build self-energy function at low frequency region
!-------------------------------------------------------------------------
! filter grnf to suppress the fluctuation of its real part
!-------------------------------------------------------------------------
!<     do k=1,nfreq
!<         do i=1,norbs
!<             ob =  real( grnf(k,i,i) )
!<             oe = aimag( grnf(k,i,i) )
!<             grnf(k,i,i) = dcmplx( zero, oe )
!<         enddo ! over i={1,norbs} loop
!<     enddo ! over k={1,nfreq} loop
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     do k=1,nfreq
         gaux = grnf(k,:,:)
         call s_inv_z(norbs, gaux)
         do i=1,norbs
             sig2(k,i,i) = czi * rmesh(k) + mune - eimp(i) - gaux(i,i) - hybf(k,i,i)
             sig2(k,i,i) = sig2(k,i,i) - shift / two
         enddo ! over i={1,norbs} loop
     enddo ! over k={1,nfreq} loop
!-------------------------------------------------------------------------
! filter sig2 to suppress the fluctuation of its imaginary part
!-------------------------------------------------------------------------
     do k=1,16
         do i=1,norbs
             call ctqmc_smth_sigf( sig2(1:nfreq,i,i) ) ! smooth it 16 times
         enddo ! over i={1,norbs} loop
     enddo ! over k={1,16} loop
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! interpolates self-energy function between low energy QMC data and high
! energy Hubbard-I approximation
     do i=1,norbs

! determine the base point, its value is calculated by using five points
         cb = czero
         do k=nfreq-4,nfreq
             cb = cb + sig2(k,i,i)
         enddo ! over k={nfreq-4,nfreq} loop

         cb = cb / real(5)
         ob = rmesh(nfreq-2)

! step A: for the imaginary part
! determine the intermediate region [nfreq+1,start] at first
         start = 0
         do k=nfreq+1,mfreq
             start = k
             d0 = aimag( shub(k,i) - cb ) / ( rmesh(k) - ob )
             d1 = aimag( shub(k,i) - shub(k-1,i) ) / ( rmesh(k) - rmesh(k-1) )
             if ( abs( d0 - d1 ) < 0.02_dp ) EXIT
         enddo ! over k={nfreq+1,mfreq} loop

! we just constrain start \in [nfreq + 32, nfreq + 128]
         if ( start - nfreq <  32 ) start = nfreq +  32
         if ( start - nfreq > 128 ) start = nfreq + 128

         ce = shub(start,i)
         oe = rmesh(start)

! deal with the intermediate region, using linear interpolation
         do k=nfreq+1,start
             sig2(k,i,i) = dcmplx( zero, aimag(cb) + aimag( ce - cb ) * ( rmesh(k) - ob ) / ( oe - ob ) )
         enddo ! over k={nfreq+1,start} loop

! deal with the tail region, using atomic self-energy function directly
         do k=start+1,mfreq
             sig2(k,i,i) = dcmplx( zero, aimag( shub(k,i) ) )
         enddo ! over k={start+1,mfreq} loop

! step B: for the real part
         sinf = shub(mfreq,i)
         do k=nfreq+1,mfreq
             sig2(k,i,i) = sig2(k,i,i) + real(sinf) + ( ob / rmesh(k) )**2 * real( cb - sinf )
         enddo ! over k={nfreq+1,mfreq} loop

     enddo ! over i={1,norbs} loop

! calculate final impurity green's function using dyson's equation
     do k=1,mfreq
         gaux = czero
         do i=1,norbs
             gaux(i,i) = czi * rmesh(k) + mune - eimp(i) - sig2(k,i,i) - hybf(k,i,i)
             gaux(i,i) = gaux(i,i) - shift / two
         enddo ! over i={1,norbs} loop
         call s_inv_z(norbs, gaux)
         grnf(k,:,:) = gaux
     enddo ! over k={1,mfreq} loop

     return
  end subroutine ctqmc_make_hub1

!!>>> ctqmc_make_hub2: build atomic green's function and self-energy
!!>>> function using improved Hubbard-I approximation, and then make
!!>>> forward fourier transformation for impurity green's function and
!!>>> auxiliary correlation function. then the final self-energy function
!!>>> is obtained by analytical formula.
  subroutine ctqmc_make_hub2()
     use constants, only : dp, zero, one, two, pi, czi, czero

     use control, only : isort
     use control, only : norbs, ncfgs
     use control, only : lemax
     use control, only : mfreq
     use control, only : ntime
     use control, only : mune, beta
     use control, only : myid, master
     use context, only : tmesh, rmesh
     use context, only : prob, nmat
     use context, only : eimp, uumat
     use context, only : gtau, ftau, grnf, frnf
     use context, only : sig2

     implicit none

! local parameters
! maximum allowable number of non-zero elements in F matrix
     integer, parameter :: nzero = 1024

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k
     integer  :: m
     integer  :: n

! dummy integer variables, used to build F matrix
     integer  :: value
     integer  :: permute

! dummy real variables, used to build atomic green's function
     real(dp) :: ob
     real(dp) :: shift

! dummy complex variables, used to build atomic green's function
     complex(dp) :: cb

! dummy atomic states: alpha, beta, gamma
     integer  :: sa(norbs)
     integer  :: sb(norbs)
     integer  :: sc(norbs)

! atomic basis sets
     integer  :: basis(ncfgs,norbs)

! F matrix, < alpha | f_{n} | beta >
     integer  :: fcounter(norbs)
     integer  :: fa(nzero,norbs)
     integer  :: fb(nzero,norbs)
     integer  :: fv(nzero,norbs)

! eigenvalues for local hmailtonian
     real(dp) :: eaux(ncfgs)

! spherical Bessel functions
     real(dp) :: jaux(mfreq,lemax)

! imaginary time green's function
     real(dp) :: gaux(ntime,norbs,norbs)

! auxiliary correlation function on imaginary time axis
     real(dp) :: faux(ntime,norbs,norbs)

! unitary transformation matrix for legendre polynomial
     complex(dp) :: taux(mfreq,lemax)

! atomic green's function and self-energy function in Hubbard-I approximation
     complex(dp) :: ghub(mfreq,norbs)
     complex(dp) :: shub(mfreq,norbs)

! evaluate the shift for the Coulomb interaction and chemical potential
! if the retarded interaction is used
     call ctqmc_eval_shift(shift)

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

! evaluate atomic eigenvalues directly
     eaux = zero
     do i=1,ncfgs
         do j=1,norbs
             eaux(i) = eaux(i) + ( eimp(j) - mune ) * basis(i,j)
         enddo ! over j={1,norbs} loop
         do j=1,norbs-1
             do k=j+1,norbs
                 if ( basis(i,j) == 1 .and. basis(i,k) == 1 ) then
                     eaux(i) = eaux(i) + uumat(j,k)
                 endif ! back if ( basis(i,j) == 1 .and. basis(i,k) == 1 ) block
             enddo ! over k={j+1,norbs} loop
         enddo ! over j={1,norbs-1} loop
     enddo ! over i={1,ncfgs} loop

! build F matrix < alpha | f_{n} | beta >
! note 1: to save the memory and accelerate the computation, we only store
! the non-zero element of F matrix
! note 2: it is crucial to check whether the number of non-zero elements
! exceed limit (nzero)
     fcounter = 0
     alpha_loop: do i=1,ncfgs
         sa = basis(i,:)
         beta_loop: do j=1,ncfgs
             sb = basis(j,:)

             orbital_loop: do m=1,norbs
                 sc = sb

                 if ( sc(m) == 1 ) then
                     permute = 1
                     do n=1,m-1
                         if ( sc(n) == 1 ) permute = -permute
                     enddo ! over n={1,m-1} loop
                     sc(m) = 0

                     value = 1
                     do n=1,norbs
                         if ( sa(n) /= sc(n) ) value = 0
                     enddo ! over n={1,norbs} loop
                     value = value * permute
                 else
                     value = 0
                 endif ! back if ( sc(m) == 1 ) block

                 if ( value /= 0 ) then
                     fcounter(m) = fcounter(m) + 1
                     if ( fcounter(m) > nzero ) then
                         call s_print_error('ctqmc_make_hub2','non-zero elements exceed limit')
                     endif ! back if ( fcounter(m) > nzero ) block
                     fa(fcounter(m),m) = i
                     fb(fcounter(m),m) = j
                     fv(fcounter(m),m) = value
                 endif ! back if ( value /= 0 ) block
             enddo orbital_loop ! over m={1,norbs} loop

         enddo beta_loop ! over j={1,ncfgs} loop
     enddo alpha_loop ! over i={1,ncfgs} loop

! calculate atomic green's function using Hubbard-I approximation
     do i=1,norbs
         do k=1,mfreq
             ghub(k,i) = czero
             do m=1,fcounter(i)
                 ob = fv(m,i) * fv(m,i) * ( prob(fa(m,i)) + prob(fb(m,i)) )
                 cb = czi * rmesh(k) + eaux(fa(m,i)) - eaux(fb(m,i))
                 ghub(k,i) = ghub(k,i) + ob / cb
             enddo ! over m={1,fcounter(i)} loop
         enddo ! over k={1,mfreq} loop
     enddo ! over i={1,norbs} loop

! calculate atomic self-energy function using dyson's equation
     do i=1,norbs
         do k=1,mfreq
             shub(k,i) = czi * rmesh(k) + mune - eimp(i) - one / ghub(k,i)
             shub(k,i) = shub(k,i) - nmat(i) * shift
         enddo ! over k={1,mfreq} loop
     enddo ! over i={1,norbs} loop

! dump the ghub and shub, only for reference, only the master node can do it
     if ( myid == master ) then
         call ctqmc_dump_hub1(rmesh, ghub, shub)
     endif ! back if ( myid == master ) block

! build final impurity green's function and then transform them into
! matsubara frequency axis
! note: only for isort == 4 .or. isort == 6 cases
     if ( isort /= 5 ) then
         call ctqmc_make_gtau(tmesh, gtau, gaux)
         call ctqmc_four_htau(gaux, grnf)
     endif ! back if ( isort /= 5 ) block

! build final auxiliary correlation function and then transform them into
! matsubara frequency axis
! note: only for isort == 4 .or. isort == 6 cases
     if ( isort /= 5 ) then
         call ctqmc_make_ftau(tmesh, ftau, faux)
         call ctqmc_four_htau(faux, frnf)
     endif ! back if ( isort /= 5 ) block

! special consideration must be taken for legendre representation, we can
! calculate grnf and frnf directly by using legendre coefficients, instead
! of performing fourier transformation
     if ( isort == 5 ) then
! build spherical Bessel functions: jaux
         jaux = zero
         do k=1,mfreq
             ob = (two * k - one) * pi / two
             call s_sbessel(lemax-1, ob, jaux(k,:))
         enddo ! over k={1,mfreq} loop

! build unitary transformation matrix: taux
         taux = czero
         do i=1,lemax
             do k=1,mfreq
                 ob = (-one)**(k - 1) * sqrt(two * i - one)
                 cb = czi**i
                 taux(k,i) = jaux(k,i) * ob * cb
             enddo ! over k={1,mfreq} loop
         enddo ! over i={1,lemax} loop

! rebuild impurity green's function on matsubara frequency (grnf) using
! orthogonal polynomial representation, G(i\omega)
! rebuild auxiliary correlation function on matsubara frequency (frnf)
! using orthogonal polynomial representation, F(i\omega)
         grnf = czero
         frnf = czero
         do i=1,norbs
             do j=1,lemax
                 do k=1,mfreq
                     grnf(k,i,i) = grnf(k,i,i) + taux(k,j) * gtau(j,i,i) / beta
                     frnf(k,i,i) = frnf(k,i,i) + taux(k,j) * ftau(j,i,i) / beta
                 enddo ! over k={1,mfreq} loop
             enddo ! over j={1,lemax} loop
         enddo ! over i={1,norbs} loop
     endif ! back if ( isort == 5 ) block

! build full self-energy function by using frnf and grnf
     do i=1,norbs
         do k=1,mfreq
             sig2(k,i,i) = frnf(k,i,i) / grnf(k,i,i)
         enddo ! over k={1,nfreq} loop
     enddo ! over i={1,norbs} loop

     return
  end subroutine ctqmc_make_hub2
