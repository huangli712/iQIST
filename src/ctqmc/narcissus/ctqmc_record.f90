!!!-----------------------------------------------------------------------
!!! project : narcissus
!!! program : ctqmc_record_gtau
!!!           ctqmc_record_ftau
!!!           ctqmc_record_grnf
!!!           ctqmc_record_hist
!!!           ctqmc_record_prob
!!!           ctqmc_record_nmat
!!!           ctqmc_record_schi
!!!           ctqmc_record_ochi
!!!           ctqmc_record_twop
!!!           ctqmc_record_vrtx
!!!           ctqmc_record_pair <<<---
!!!           ctqmc_reduce_gtau
!!!           ctqmc_reduce_ftau
!!!           ctqmc_reduce_grnf
!!!           ctqmc_reduce_hist
!!!           ctqmc_reduce_prob
!!!           ctqmc_reduce_nmat
!!!           ctqmc_reduce_schi
!!!           ctqmc_reduce_ochi
!!!           ctqmc_reduce_twop
!!!           ctqmc_reduce_vrtx
!!!           ctqmc_reduce_pair <<<---
!!!           ctqmc_symm_nmat
!!!           ctqmc_symm_gtau
!!!           ctqmc_symm_grnf
!!!           ctqmc_smth_sigf   <<<---
!!!           ctqmc_make_gtau
!!!           ctqmc_make_ftau   <<<---
!!!           ctqmc_make_hub1
!!!           ctqmc_make_hub2   <<<---
!!! source  : ctqmc_record.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 09/16/2009 by li huang
!!!           09/29/2010 by li huang
!!!           09/18/2014 by li huang
!!!           10/13/2014 by li huang
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
                 endif

! determine index for imaginary time
                 curr = nint( dtau * step ) + 1

! special tricks for the first point and the last point
                 if ( curr == 1 .or. curr == ntime ) then
                     maux = two * maux
                 endif

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
                 endif

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
                 endif

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
!!>>> imaginary time axis
  subroutine ctqmc_record_ftau()
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
     use context, only : ftau

     implicit none

! local variables
! loop indices for start and end points
     integer  :: is
     integer  :: ie

! loop index for flavor channel
     integer  :: flvr

! loop index for colour channel
     integer  :: clur

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

! occupation number at taus
     real(dp) :: occu

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

! evaluate step at first
     step = real(ntime - 1) / beta

     CTQMC_FLAVOR_LOOP: do flvr=1,norbs
     CTQMC_COLOUR_LOOP: do clur=1,norbs

! skip diagonal term
         if ( flvr == clur ) CYCLE

! get imaginary time value for segments
         do is=1,rank(flvr)
             taus = time_s( index_s(is, flvr), flvr )

! evaluate occu, and then check it
             call ctqmc_spin_counter(clur, taus, occu); if ( occu < one ) CYCLE

! get imaginary time value for segments
             do ie=1,rank(flvr)
                 taue = time_e( index_e(ie, flvr), flvr )

! evaluate dtau
                 dtau = taue - taus

! get matrix element from mmat, pay special attention to the sign of dtau
                 maux = mmat(ie, is, flvr) * sign(one, dtau) * occu

! adjust dtau, keep it stay in (zero, beta)
                 if ( dtau < zero ) then
                     dtau = dtau + beta
                 endif

! determine index for imaginary time
                 curr = nint( dtau * step ) + 1

! special tricks for the first point and the last point
                 if ( curr == 1 .or. curr == ntime ) then
                     maux = two * maux
                 endif

! record ftau, we normalize ftau in ctqmc_make_ftau() subroutine
                 ftau(curr, clur, flvr) = ftau(curr, clur, flvr) - maux

             enddo ! over ie={1,rank(flvr)} loop
         enddo ! over is={1,rank(flvr)} loop

     enddo CTQMC_COLOUR_LOOP ! over clur={1,norbs} loop
     enddo CTQMC_FLAVOR_LOOP ! over flvr={1,norbs} loop

     return
  end subroutine cat_record_ftau1

!!>>> cat_record_ftau2: record auxiliary correlation function using
!!>>> legendre polynomial representation
  subroutine cat_record_ftau2()
     implicit none

! evaluate step at first
     step = real(legrd - 1) / two

     CTQMC_FLAVOR_LOOP: do flvr=1,norbs
     CTQMC_COLOUR_LOOP: do clur=1,norbs

! skip diagonal term
         if ( flvr == clur ) CYCLE

! get imaginary time value for segments
         do is=1,rank(flvr)
             taus = time_s( index_s(is, flvr), flvr )

! evaluate occu, and then check it
             call ctqmc_spin_counter(clur, taus, occu); if ( occu < one ) CYCLE

             do ie=1,rank(flvr)
                 taue = time_e( index_e(ie, flvr), flvr )

! evaluate dtau
                 dtau = taue - taus

! get matrix element from mmat, pay special attention to the sign of dtau
                 maux = mmat(ie, is, flvr) * sign(one, dtau) * occu

! adjust dtau, keep it stay in (zero, beta)
                 if ( dtau < zero ) then
                     dtau = dtau + beta
                 endif

! convert dtau in [0,\beta] to daux in [0,2]
                 daux = two * dtau / beta

! determine index for legendre polynomial interval
                 curr = nint( daux * step ) + 1

! record ftau, we normalize ftau in ctqmc_make_ftau() subroutine
                 CTQMC_FLALEG_LOOP: do fleg=1,lemax
                     dtau = sqrt(two * fleg - 1) * ppleg(curr,fleg)
                     ftau(fleg, clur, flvr) = ftau(fleg, clur, flvr) - maux * dtau
                 enddo CTQMC_FLALEG_LOOP ! over fleg={1,lemax} loop

             enddo ! over ie={1,rank(flvr)} loop
         enddo ! over is={1,rank(flvr)} loop

     enddo CTQMC_COLOUR_LOOP ! over clur={1,norbs} loop
     enddo CTQMC_FLAVOR_LOOP ! over flvr={1,norbs} loop

     return
  end subroutine cat_record_ftau2

!!>>> cat_record_ftau3: record auxiliary correlation function using
!!>>> chebyshev polynomial representation
  subroutine cat_record_ftau3()
     implicit none

! evaluate step at first
     step = real(chgrd - 1) / two

     CTQMC_FLAVOR_LOOP: do flvr=1,norbs
     CTQMC_COLOUR_LOOP: do clur=1,norbs

! skip diagonal term
         if ( flvr == clur ) CYCLE

! get imaginary time value for segments
         do is=1,rank(flvr)
             taus = time_s( index_s(is, flvr), flvr )

! evaluate occu, and then check it
             call ctqmc_spin_counter(clur, taus, occu); if ( occu < one ) CYCLE

             do ie=1,rank(flvr)
                 taue = time_e( index_e(ie, flvr), flvr )

! evaluate dtau
                 dtau = taue - taus

! get matrix element from mmat, pay special attention to the sign of dtau
                 maux = mmat(ie, is, flvr) * sign(one, dtau) * occu

! adjust dtau, keep it stay in (zero, beta)
                 if ( dtau < zero ) then
                     dtau = dtau + beta
                 endif

! convert dtau in [0,\beta] to daux in [0,2]
                 daux = two * dtau / beta

! determine index for chebyshev polynomial interval
                 curr = nint( daux * step ) + 1

! record ftau, we normalize ftau in ctqmc_make_ftau() subroutine
                 CTQMC_FLACHE_LOOP: do fche=1,chmax
                     dtau = (two / pi) * sqrt(one - (daux - one)**2) * qqche(curr,fche)
                     ftau(fche, clur, flvr) = ftau(fche, clur, flvr) - maux * dtau
                 enddo CTQMC_FLACHE_LOOP ! over fche={1,chmax} loop

             enddo ! over ie={1,rank(flvr)} loop
         enddo ! over is={1,rank(flvr)} loop

     enddo CTQMC_COLOUR_LOOP ! over clur={1,norbs} loop
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
     use control, only : mkink
     use context, only : ckink
     use context, only : hist

     implicit none

! note: if ckink == 0, we record its count in hist(mkink)
     if ( ckink > 0 ) then
         hist(ckink) = hist(ckink) + 1
     else
         hist(mkink) = hist(mkink) + 1
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
     use constants, only : dp, zero

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
                 call ctqmc_make_overlap(flvr, ts, te, oaux)
                 ovlp(flvr,:) = ovlp(flvr,:) + oaux
             enddo ! over i={1,rank(flvr)} loop

! case 3: partial occupation, anti-segment scheme
! pay special attention to the head and tail parts
         else if ( stts(flvr) == 2 ) then
             ovlp(flvr,:) = zero
             do i=1,rank(flvr)-1
                 ts = time_s(index_s(i,   flvr), flvr)
                 te = time_e(index_e(i+1, flvr), flvr)
                 call ctqmc_make_overlap(flvr, ts, te, oaux)
                 ovlp(flvr,:) = ovlp(flvr,:) + oaux
             enddo ! over i={1,rank(flvr)-1} loop

             te = time_e(index_e(1, flvr), flvr)
             call ctqmc_make_overlap(flvr, zero, te, oaux)
             ovlp(flvr,:) = ovlp(flvr,:) + oaux

             ts = time_s(index_s(rank(flvr), flvr), flvr)
             call ctqmc_make_overlap(flvr, ts, beta, oaux)
             ovlp(flvr,:) = ovlp(flvr,:) + oaux

! case 4: full occupation
         else if ( stts(flvr) == 3 ) then
             call ctqmc_make_overlap(flvr, zero, beta, oaux)
             ovlp(flvr,:) = oaux

         endif ! back if ( stts(flvr) == 0 ) block

         nnmat(flvr,:) = nnmat(flvr,:) + ovlp(flvr,:) / beta
     enddo ! over flvr={1,norbs} loop
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

!!>>> ctqmc_record_schi: record the spin-spin correlation function
  subroutine ctqmc_record_schi()
     use constants, only : dp, zero

     use control, only : nband, norbs
     use control, only : ntime
     use context, only : tmesh
     use context, only : schi, sschi

     implicit none

! local variables
! loop index over segments
     integer  :: i

! loop index for flavor channel
     integer  :: flvr

! Sz(0) and Sz(\tau)
     real(dp) :: sz1_s
     real(dp) :: sz2_s
     real(dp) :: sz1_i(nband)
     real(dp) :: sz2_i(nband)

! used to record occupations for current flavor channel and time
     real(dp) :: oaux(norbs)

     TIME_LOOP: do i=1,ntime

! obtain occupation status
         oaux = zero
         do flvr=1,norbs
             call ctqmc_spin_counter(flvr, tmesh(i), oaux(flvr))
         enddo ! over flvr={1,norbs} loop

! calculate schi
! evaluate Sz(\tau)
         sz2_s = zero
         do flvr=1,nband
             sz2_s = sz2_s + oaux(flvr) - oaux(flvr+nband)
         enddo ! over flvr={1,nband} loop

! evaluate Sz(0)
         if ( i == 1 ) then
             sz1_s = sz2_s
         endif ! back if ( i == 1 ) block

! sum up the contribution to schi
         schi(i) = schi(i) + sz1_s * sz2_s

! calculate sschi
         BAND_LOOP: do flvr=1,nband

! evaluate Sz(\tau)
             sz2_i(flvr) = oaux(flvr) - oaux(flvr+nband)

! evaluate Sz(0)
             if ( i == 1 ) then
                 sz1_i(flvr) = sz2_i(flvr)
             endif ! back if ( i == 1 ) block

! sum up the contribution to sschi
             sschi(i,flvr) = sschi(i,flvr) + sz1_i(flvr) * sz2_i(flvr)

         enddo BAND_LOOP ! over flvr={1,nband} loop

     enddo TIME_LOOP ! over i={1,ntime} loop

     return
  end subroutine ctqmc_record_schi

!!>>> ctqmc_record_ochi: record the orbital-orbital correlation function
  subroutine ctqmc_record_ochi()
     use constants, only : dp, zero

     use control, only : norbs
     use control, only : ntime
     use context, only : tmesh
     use context, only : ochi, oochi

     implicit none

! local variables
! loop index over segments
     integer  :: i

! loop index for flavor channel
     integer  :: flvr

! N(0) and N(\tau)
     real(dp) :: nt_s
     real(dp) :: nz_s
     real(dp) :: nt_i(norbs)
     real(dp) :: nz_i(norbs)

! used to record occupations for current flavor channel and time
     real(dp) :: oaux(norbs)

     TIME_LOOP: do i=1,ntime

! obtain occupation status
         oaux = zero
         do flvr=1,norbs
             call ctqmc_spin_counter(flvr, tmesh(i), oaux(flvr))
         enddo ! over flvr={1,norbs} loop

! calculate ochi
! evaluate N(\tau)
         nt_s = sum( oaux )

! evaluate N(0)
         if ( i == 1 ) then
             nz_s = nt_s
         endif ! back if ( i == 1 ) block

! sum up the contribution to ochi
         ochi(i) = ochi(i) + nz_s * nt_s

! calculate oochi
         BAND_LOOP: do flvr=1,norbs

! evaluate N_{\alpha}(\tau)
             nt_i(flvr) = oaux(flvr)

! evaluate N_{\alpha}(0)
             if ( i == 1 ) then
                 nz_i(flvr) = nt_i(flvr)
             endif ! back if ( i == 1 ) block

! sum up the contribution to oochi
             oochi(i,flvr) = oochi(i,flvr) + nz_i(flvr) * nt_i(flvr)

         enddo BAND_LOOP ! over flvr={1,norbs} loop

     enddo TIME_LOOP ! over i={1,ntime} loop

     return
  end subroutine ctqmc_record_ochi

!!>>> ctqmc_record_twop: record the two-particle green's function
  subroutine ctqmc_record_twop()
     use constants, only : dp, two, pi, czi, czero

     use control, only : norbs
     use control, only : nffrq, nbfrq
     use control, only : beta
     use context, only : index_s, index_e, time_s, time_e
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

! imaginary time for start and end points
     real(dp) :: taus
     real(dp) :: taue

! dummy complex(dp) variables, used to calculate the exponential function
     complex(dp) :: cmeas
     complex(dp) :: dexp1
     complex(dp) :: dexp2
     complex(dp) :: iexp1
     complex(dp) :: iexp2
     complex(dp) :: cexp1
     complex(dp) :: cexp2

! dummy complex(dp) arrays, used to store the intermediate results
     complex(dp), allocatable :: g2aux(:,:,:)

! evaluate nfaux, determine the size of g2aux
     nfaux = nffrq + nbfrq - 1

! allocate memory for g2aux and then initialize it
     allocate( g2aux(norbs, nfaux, nfaux) ); g2aux = czero

     CTQMC_FLAVOR_LOOP: do flvr=1,norbs

! get imaginary time value for segments
         do is=1,rank(flvr)
             taus = time_s( index_s(is, flvr), flvr )

             do ie=1,rank(flvr)
                 taue = time_e( index_e(ie, flvr), flvr )

! get matrix element from mmat
                 maux = mmat(ie, is, flvr)

! calculate g2aux
                 dexp1 = exp(+    two     * czi * pi * taue / beta)
                 dexp2 = exp(-    two     * czi * pi * taus / beta)
                 iexp1 = exp(-(nffrq + 1) * czi * pi * taue / beta)
                 iexp2 = exp(+(nffrq + 1) * czi * pi * taus / beta)

                 cexp1 = iexp1
                 do w1n=1,nfaux
                     cexp1 = cexp1 * dexp1

                     cexp2 = iexp2
                     do w2n=1,nfaux
                         cexp2 = cexp2 * dexp2

                         g2aux(flvr, w1n, w2n) = g2aux(flvr, w1n, w2n) + maux * cexp1 * cexp2
                     enddo ! over w2n={1,nfaux} loop
                 enddo ! over w1n={1,nfaux} loop

             enddo ! over ie={1,rank(flvr)} loop
         enddo ! over is={1,rank(flvr)} loop

     enddo CTQMC_FLAVOR_LOOP ! over flvr={1,norbs} loop

! calculate g2_re and g2_im
     CTQMC_ORBIT1_LOOP: do f1=1,norbs
         CTQMC_ORBIT2_LOOP: do f2=1,norbs

             CTQMC_FERMI1_LOOP: do w2n=1,nffrq
                 CTQMC_FERMI2_LOOP: do w3n=1,nffrq

                     CTQMC_BOSONF_LOOP: do wbn=1,nbfrq
                         w1n = w2n + wbn - 1; w4n = w3n + wbn - 1
                         cmeas = g2aux(f1,w1n,w2n) * g2aux(f2,w3n,w4n)
                         if ( f1 == f2 ) then
                             cmeas = cmeas - g2aux(f1,w1n,w4n) * g2aux(f1,w3n,w2n)
                         endif ! back if ( f1 == f2 ) block
                         g2_re(f1,f2,w2n,w3n,wbn) = g2_re(f1,f2,w2n,w3n,wbn) +  real(cmeas) / beta
                         g2_im(f1,f2,w2n,w3n,wbn) = g2_im(f1,f2,w2n,w3n,wbn) + aimag(cmeas) / beta
                     enddo CTQMC_BOSONF_LOOP ! over wbn={1,nbfrq} loop

                 enddo CTQMC_FERMI2_LOOP ! over w3n={1,nffrq} loop
             enddo CTQMC_FERMI1_LOOP ! over w2n={1,nffrq} loop

         enddo CTQMC_ORBIT2_LOOP ! over f2={1,norbs} loop
     enddo CTQMC_ORBIT1_LOOP ! over f1={1,norbs} loop

! deallocate memory
     deallocate( g2aux )

     return
  end subroutine ctqmc_record_twop

!!>>> ctqmc_record_vrtx: record the fake vertex function
  subroutine ctqmc_record_vrtx()
     use constants, only : dp, zero, one, two, half, pi, czi, czero

     use control, only : norbs
     use control, only : nffrq, nbfrq
     use control, only : beta
     use context, only : index_s, index_e, time_s, time_e
     use context, only : g2_re, g2_im, h2_re, h2_im
     use context, only : rank, uumat
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

! imaginary time for start and end points
     real(dp) :: taus
     real(dp) :: taue

! occupation number at taus
     real(dp) :: occu
     real(dp) :: oaux

! dummy complex(dp) variables, used to calculate the exponential function
     complex(dp) :: cmeas
     complex(dp) :: dexp1
     complex(dp) :: dexp2
     complex(dp) :: iexp1
     complex(dp) :: iexp2
     complex(dp) :: cexp1
     complex(dp) :: cexp2

! dummy complex(dp) arrays, used to store the intermediate results
     complex(dp), allocatable :: g2aux(:,:,:)
     complex(dp), allocatable :: h2aux(:,:,:)

! evaluate nfaux, determine the size of g2aux and h2aux
     nfaux = nffrq + nbfrq - 1

! allocate memory for g2aux and then initialize it
     allocate( g2aux(norbs, nfaux, nfaux) ); g2aux = czero

! allocate memory for h2aux and then initialize it
     allocate( h2aux(norbs, nfaux, nfaux) ); h2aux = czero

     CTQMC_FLAVOR_LOOP: do flvr=1,norbs

! get imaginary time value for segments
         do ie=1,rank(flvr)
             taue = time_e( index_e(ie, flvr), flvr )

! evaluate occu, and then check it
             oaux = zero
             do f1=1,norbs
                 call ctqmc_spin_counter(f1, taue, occu); if ( occu < one ) CYCLE
                 oaux = oaux + half * ( uumat(f1,flvr) + uumat(flvr,f1) ) * occu
             enddo ! over f1={1,norbs} loop

             do is=1,rank(flvr)
                 taus = time_s( index_s(is, flvr), flvr )

! get matrix element from mmat
                 maux = mmat(ie, is, flvr)

! calculate g2aux and h2aux
                 dexp1 = exp(+    two     * czi * pi * taue / beta)
                 dexp2 = exp(-    two     * czi * pi * taus / beta)
                 iexp1 = exp(-(nffrq + 1) * czi * pi * taue / beta)
                 iexp2 = exp(+(nffrq + 1) * czi * pi * taus / beta)

                 cexp1 = iexp1
                 do w1n=1,nfaux
                     cexp1 = cexp1 * dexp1

                     cexp2 = iexp2
                     do w2n=1,nfaux
                         cexp2 = cexp2 * dexp2

                         g2aux(flvr, w1n, w2n) = g2aux(flvr, w1n, w2n) + maux * cexp1 * cexp2
                         h2aux(flvr, w1n, w2n) = h2aux(flvr, w1n, w2n) + maux * cexp1 * cexp2 * oaux
                     enddo ! over w2n={1,nfaux} loop
                 enddo ! over w1n={1,nfaux} loop

             enddo ! over is={1,rank(flvr)} loop
         enddo ! over ie={1,rank(flvr)} loop

     enddo CTQMC_FLAVOR_LOOP ! over flvr={1,norbs} loop

! calculate g2_re and g2_im, h2_re and h2_im
     CTQMC_ORBIT1_LOOP: do f1=1,norbs
         CTQMC_ORBIT2_LOOP: do f2=1,norbs

             CTQMC_FERMI1_LOOP: do w2n=1,nffrq
                 CTQMC_FERMI2_LOOP: do w3n=1,nffrq

                     CTQMC_BOSONF_LOOP: do wbn=1,nbfrq
                         w1n = w2n + wbn - 1; w4n = w3n + wbn - 1

                         cmeas = g2aux(f1,w1n,w2n) * g2aux(f2,w3n,w4n)
                         if ( f1 == f2 ) then
                             cmeas = cmeas - g2aux(f1,w1n,w4n) * g2aux(f1,w3n,w2n)
                         endif ! back if ( f1 == f2 ) block
                         g2_re(f1,f2,w2n,w3n,wbn) = g2_re(f1,f2,w2n,w3n,wbn) +  real(cmeas) / beta
                         g2_im(f1,f2,w2n,w3n,wbn) = g2_im(f1,f2,w2n,w3n,wbn) + aimag(cmeas) / beta

                         cmeas = h2aux(f1,w1n,w2n) * g2aux(f2,w3n,w4n)
                         if ( f1 == f2 ) then
                             cmeas = cmeas - h2aux(f1,w1n,w4n) * g2aux(f1,w3n,w2n)
                         endif ! back if ( f1 == f2 ) block
                         h2_re(f1,f2,w2n,w3n,wbn) = h2_re(f1,f2,w2n,w3n,wbn) +  real(cmeas) / beta
                         h2_im(f1,f2,w2n,w3n,wbn) = h2_im(f1,f2,w2n,w3n,wbn) + aimag(cmeas) / beta

                     enddo CTQMC_BOSONF_LOOP ! over wbn={1,nbfrq} loop

                 enddo CTQMC_FERMI2_LOOP ! over w3n={1,nffrq} loop
             enddo CTQMC_FERMI1_LOOP ! over w2n={1,nffrq} loop

         enddo CTQMC_ORBIT2_LOOP ! over f2={1,norbs} loop
     enddo CTQMC_ORBIT1_LOOP ! over f1={1,norbs} loop

! deallocate memory
     deallocate( g2aux )
     deallocate( h2aux )

     return
  end subroutine ctqmc_record_vrtx

!!========================================================================
!!>>> reduce physical observables                                      <<<
!!========================================================================

!!>>> ctqmc_reduce_gtau: reduce the gtau from all children processes
  subroutine ctqmc_reduce_gtau(gtau_mpi)
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

! initialize gtau_mpi
     gtau_mpi = zero

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

     return
  end subroutine ctqmc_reduce_gtau

!!>>> ctqmc_reduce_ftau: reduce the ftau from all children processes
  subroutine ctqmc_reduce_ftau(ftau_mpi)
     use constants, only : dp, zero
     use mmpi, only : mp_allreduce, mp_barrier

     use control, only : norbs
     use control, only : ntime
     use control, only : nprocs
     use context, only : ftau

     implicit none

! external arguments
! auxiliary correlation function
     real(dp), intent(out) :: ftau_mpi(ntime,norbs,norbs)

! initialize ftau_mpi
     ftau_mpi = zero

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

     return
  end subroutine ctqmc_reduce_ftau

!!>>> ctqmc_reduce_grnf: reduce the grnf from all children processes
  subroutine ctqmc_reduce_grnf(grnf_mpi)
     use constants, only : dp, czero
     use mmpi, only : mp_allreduce, mp_barrier

     use control, only : norbs
     use control, only : mfreq
     use control, only : nprocs
     use context, only : grnf

     implicit none

! external arguments
! impurity green's function
     complex(dp), intent(out) :: grnf_mpi(mfreq,norbs,norbs)

! initialize grnf_mpi
     grnf_mpi = czero

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

     return
  end subroutine ctqmc_reduce_grnf

!!>>> ctqmc_reduce_hist: reduce the hist from all children processes
!!>>> note: since hist_mpi and hist are integer (kind=4) type, it is
!!>>> important to avoid data overflow in them
  subroutine ctqmc_reduce_hist(hist_mpi)
     use mmpi, only : mp_allreduce, mp_barrier

     use control, only : mkink
     use control, only : nprocs
     use context, only : hist

     implicit none

! external arguments
! histogram for perturbation expansion series
     integer, intent(out) :: hist_mpi(mkink)

! initialize hist_mpi
     hist_mpi = 0

! build hist_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(hist / nprocs, hist_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     hist_mpi = hist

# endif /* MPI */

! calculate the average
     hist_mpi = hist_mpi / 1

     return
  end subroutine ctqmc_reduce_hist

!!>>> ctqmc_reduce_prob: reduce the prob from all children processes
  subroutine ctqmc_reduce_prob(prob_mpi)
     use constants, only : dp, zero
     use mmpi, only : mp_allreduce, mp_barrier

     use control, only : ncfgs
     use control, only : nprocs
     use context, only : prob

     implicit none

! external arguments
! probability of atomic states
     real(dp), intent(out) :: prob_mpi(ncfgs)

! initialize prob_mpi
     prob_mpi = zero

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

     return
  end subroutine ctqmc_reduce_prob

!!>>> ctqmc_reduce_nmat: reduce the nmat and nnmat from all children processes
  subroutine ctqmc_reduce_nmat(nmat_mpi, nnmat_mpi)
     use constants, only : dp, zero
     use mmpi, only : mp_allreduce, mp_barrier

     use control, only : norbs
     use control, only : nprocs
     use context, only : nmat, nnmat

     implicit none

! external arguments
! occupation number matrix
     real(dp), intent(out) :: nmat_mpi(norbs)

! double occupation number matrix
     real(dp), intent(out) :: nnmat_mpi(norbs,norbs)

! initialize nmat_mpi and nnmat_mpi
     nmat_mpi = zero
     nnmat_mpi = zero

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

     return
  end subroutine ctqmc_reduce_nmat

!!>>> ctqmc_reduce_schi: reduce the schi and sschi from all children processes
  subroutine ctqmc_reduce_schi(schi_mpi, sschi_mpi)
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

! spin-spin correlation function, orbital-resolved
     real(dp), intent(out) :: sschi_mpi(ntime,nband)

! initialize schi_mpi and sschi_mpi
     schi_mpi = zero
     sschi_mpi = zero

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

     return
  end subroutine ctqmc_reduce_schi

!!>>> ctqmc_reduce_ochi: reduce the ochi and oochi from all children processes
  subroutine ctqmc_reduce_ochi(ochi_mpi, oochi_mpi)
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

! orbital-orbital correlation function, orbital-resolved
     real(dp), intent(out) :: oochi_mpi(ntime,norbs)

! initialize ochi_mpi and oochi_mpi
     ochi_mpi = zero
     oochi_mpi = zero

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

     return
  end subroutine ctqmc_reduce_ochi

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
     real(dp), intent(out) :: g2_re_mpi(norbs,norbs,nffrq,nffrq,nbfrq)

! two-particle green's function, imaginary part
     real(dp), intent(out) :: g2_im_mpi(norbs,norbs,nffrq,nffrq,nbfrq)

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
! vertex function, real part
     real(dp), intent(out) :: h2_re_mpi(norbs,norbs,nffrq,nffrq,nbfrq)

! vertex function, imaginary part
     real(dp), intent(out) :: h2_im_mpi(norbs,norbs,nffrq,nffrq,nbfrq)

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
                     endif
                 enddo ! over jbnd={1,norbs} loop

                 raux = raux / real(hist(ibnd)) ! calculate average value

                 do jbnd=1,norbs                ! setup it
                     if ( symm(jbnd) == ibnd ) then
                         nmat(jbnd) = raux
                     endif
                 enddo ! over jbnd={1,norbs} loop
             endif
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
                         endif
                     enddo ! over jbnd={1,norbs} loop

                     raux = raux / real(hist(ibnd)) ! calculate average value

                     do jbnd=1,norbs                ! setup it
                         if ( symm(jbnd) == ibnd ) then
                             gtau(ktau,jbnd,jbnd) = raux
                         endif
                     enddo ! over jbnd={1,norbs} loop
                 endif
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
                         endif
                     enddo ! over jbnd={1,norbs} loop

                     caux = caux / real(hist(ibnd)) ! calculate average value

                     do jbnd=1,norbs                ! setup it
                         if ( symm(jbnd) == ibnd ) then
                             grnf(kfrq,jbnd,jbnd) = caux
                         endif
                     enddo ! over jbnd={1,norbs} loop
                 endif
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

!!>>> cat_make_kernel: build the integral kernel function
  subroutine cat_make_kernel(kdim, kern)
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
                 kern(kcur) = ( (curr-i) * cos(raux) + sin(raux) * cotan(pi/curr) ) / curr

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
  end subroutine cat_make_kernel

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
     ker1 = one; call cat_make_kernel(lemax, ker1)

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
     ker2 = one; call cat_make_kernel(chmax, ker2)

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
!!>>> orthogonal polynomial representation
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

!>>> build auxiliary correlation function using normal representation
  subroutine cat_make_ftau1()
     implicit none

     raux = real(ntime) / (beta * beta)
     do i=1,norbs
         do j=1,norbs
             if ( i == j ) CYCLE

             do k=1,ntime
                 faux(k,j,i) = ftau(k,j,i) * raux
             enddo ! over k={1,ntime} loop
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

     return
  end subroutine cat_make_ftau1

!>>> build auxiliary correlation function using legendre polynomial representation
  subroutine cat_make_ftau2()
     implicit none

     step = real(legrd - 1) / two
     do i=1,norbs
         do j=1,norbs
             if ( i == j ) CYCLE

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

!>>> build auxiliary correlation function using chebyshev polynomial representation
  subroutine cat_make_ftau3()
     implicit none

     step = real(chgrd - 1) / two
     do i=1,norbs
         do j=1,norbs
             if ( i == j ) CYCLE

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
!!>>> build self-energy function                                       <<<
!!========================================================================

!!>>> ctqmc_make_hub1: build atomic green's function and self-energy
!!>>> function using improved Hubbard-I approximation, and then make
!!>>> interpolation for self-energy function between low frequency QMC
!!>>> data and high frequency Hubbard-I approximation data, the full
!!>>> impurity green's function can be obtained by using dyson's equation
!!>>> finally
  subroutine ctqmc_make_hub1()
     use constants, only : dp, zero, one, czi, czero

     use control, only : norbs, ncfgs
     use control, only : mfreq
     use control, only : nfreq
     use control, only : mune
     use control, only : myid, master
     use context, only : rmesh
     use context, only : prob
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

! dummy complex variables, used to interpolate self-energy function
     complex(dp) :: cb, ce
     complex(dp) :: sinf
     complex(dp) :: caux

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

! build atomic basis set, we do not order them according to their
! occupation numbers
     do i=1,ncfgs
         do j=1,norbs
             if ( btest(i-1,j-1) .eqv. .true. ) then
                 basis(i,j) = 1
             else
                 basis(i,j) = 0
             endif
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
                 endif
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
                     endif
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
             caux = czero
             do m=1,fcounter(i)
                 ob = fv(m,i) * fv(m,i) * ( prob(fa(m,i)) + prob(fb(m,i)) )
                 cb = czi * rmesh(k) + eaux(fa(m,i)) - eaux(fb(m,i))
                 caux = caux +  ob / cb
             enddo ! over m={1,fcounter(i)} loop
             ghub(k,i) = caux
         enddo ! over k={1,mfreq} loop
     enddo ! over i={1,norbs} loop

! calculate atomic self-energy function using dyson's equation
     do i=1,norbs
         do k=1,mfreq
             shub(k,i) = czi * rmesh(k) + mune - eimp(i) - one / ghub(k,i)
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
     use constants, only : dp, zero, one, two, half, pi, czi, czero

     use control, only : isort
     use control, only : norbs, ncfgs
     use control, only : lemax
     use control, only : mfreq
     use control, only : ntime
     use control, only : mune, beta
     use control, only : myid, master
     use context, only : tmesh, rmesh
     use context, only : prob
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

! build atomic basis set, we do not order them according to their
! occupation numbers
     do i=1,ncfgs
         do j=1,norbs
             if ( btest(i-1,j-1) .eqv. .true. ) then
                 basis(i,j) = 1
             else
                 basis(i,j) = 0
             endif
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
                 endif
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
                     endif
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
                 ghub(k,i) = ghub(k,i) +  ob / cb
             enddo ! over m={1,fcounter(i)} loop
         enddo ! over k={1,mfreq} loop
     enddo ! over i={1,norbs} loop

! calculate atomic self-energy function using dyson's equation
     do i=1,norbs
         do k=1,mfreq
             shub(k,i) = czi * rmesh(k) + mune - eimp(i) - one / ghub(k,i)
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
! orthogonal polynomial representation
         grnf = czero
         do i=1,norbs
             do k=1,mfreq
                 do j=1,lemax
                     grnf(k,i,i) = grnf(k,i,i) + taux(k,j) * gtau(j,i,i) / beta
                 enddo ! over j={1,lemax} loop
             enddo ! over k={1,mfreq} loop
         enddo ! over i={1,norbs} loop

! rebuild auxiliary correlation function on matsubara frequency (frnf)
! using orthogonal polynomial representation
         frnf = czero
         do i=1,norbs
             do j=1,norbs
                 if ( i == j ) CYCLE
                 do k=1,mfreq
                     do m=1,lemax
                         frnf(k,j,i) = frnf(k,j,i) + taux(k,m) * ftau(m,j,i) / beta
                     enddo ! over m={1,lemax} loop
                 enddo ! over k={1,mfreq} loop
             enddo ! over j={1,norbs} loop
         enddo ! over i={1,norbs} loop
     endif ! back if ( isort == 5 ) block

! build full self-energy function by using frnf and grnf
     do i=1,norbs
         do k=1,mfreq
             sig2(k,i,i) = czero
             do j=1,norbs
                 sig2(k,i,i) = sig2(k,i,i) + ( uumat(j,i) + uumat(i,j) ) * frnf(k,j,i)
             enddo ! over j={1,norbs} loop
             sig2(k,i,i) = half * sig2(k,i,i) / grnf(k,i,i)
         enddo ! over k={1,nfreq} loop
     enddo ! over i={1,norbs} loop

     return
  end subroutine ctqmc_make_hub2
