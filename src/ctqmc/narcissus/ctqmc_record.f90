!!!-----------------------------------------------------------------------
!!! project : narcissus
!!! program : ctqmc_record_hist
!!!           ctqmc_record_prob
!!!           ctqmc_record_paux
!!!           ctqmc_record_nmat <<<---
!!!           ctqmc_record_gtau
!!!           ctqmc_record_ftau
!!!           ctqmc_record_grnf <<<---
!!!           ctqmc_record_kmat
!!!           ctqmc_record_lmat
!!!           ctqmc_record_szpw <<<---
!!!           ctqmc_record_schi
!!!           ctqmc_record_sfom
!!!           ctqmc_record_ochi
!!!           ctqmc_record_ofom <<<---
!!!           ctqmc_record_twop
!!!           ctqmc_record_pair <<<---
!!!           ctqmc_reduce_hist
!!!           ctqmc_reduce_prob
!!!           ctqmc_reduce_paux
!!!           ctqmc_reduce_nmat <<<---
!!!           ctqmc_reduce_gtau
!!!           ctqmc_reduce_ftau
!!!           ctqmc_reduce_grnf <<<---
!!!           ctqmc_reduce_kmat
!!!           ctqmc_reduce_lmat
!!!           ctqmc_reduce_szpw <<<---
!!!           ctqmc_reduce_schi
!!!           ctqmc_reduce_sfom
!!!           ctqmc_reduce_ochi
!!!           ctqmc_reduce_ofom <<<---
!!!           ctqmc_reduce_twop
!!!           ctqmc_reduce_pair <<<---
!!! source  : ctqmc_record.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           04/26/2017 by li huang (last modified)
!!! purpose : measure, record, and postprocess the important observables
!!!           produced by the hybridization expansion version continuous
!!!           time quantum Monte Carlo (CTQMC) quantum impurity solver
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> measure physical observables 1                                   <<<
!!========================================================================

!!
!! @sub ctqmc_record_hist
!!
!! record the histogram of perturbation expansion series
!!
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

!!
!! @sub ctqmc_record_prob
!!
!! record the probability of atomic eigenstates
!!
  subroutine ctqmc_record_prob()
     use constants, only : one

     use control, only : norbs
     use context, only : prob
     use context, only : stts

     implicit none

! local variables
! current flavor channel
     integer :: flvr

! atomic eigenstate index
     integer :: pstat

! current atomic eigenstate for segment representation
     integer :: state(norbs)

! generate current atomic eigenstate
     do flvr=1,norbs
         select case ( stts(flvr) )

             case (0:1)
                 state(flvr) = 0

             case (2:3)
                 state(flvr) = 1

         end select
     enddo ! over flvr={1,norbs} loop

! convert atomic eigenstate array to index
     call ctqmc_make_state(norbs, pstat, state)

! accumulate the data
     prob(pstat) = prob(pstat) + one

     return
  end subroutine ctqmc_record_prob

!!
!! @sub ctqmc_record_paux
!!
!! record some auxiliary physical observables. the occupation matrix and
!! double occupation matrix are measured at the same time in order to
!! improve the computational efficiency
!!
  subroutine ctqmc_record_paux()
     use constants, only : dp, zero, two

     use control, only : nband, norbs
     use control, only : beta
     use context, only : ckink
     use context, only : index_s, index_e
     use context, only : time_s, time_e
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

!-------------------------------------------------------------------------
! prepare sgmt array
!-------------------------------------------------------------------------
     SGMT_BLOCK: do flvr=1,norbs

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

     enddo SGMT_BLOCK ! over flvr={1,norbs} loop
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-------------------------------------------------------------------------
! prepare ovlp matrix
!-------------------------------------------------------------------------
     OVLP_BLOCK: do flvr=1,norbs

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

     enddo OVLP_BLOCK ! over flvr={1,norbs} loop
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate <K^4>
     paux(9) = paux(9) + ( ckink * two )**4

! evaluate <K^3>
     paux(8) = paux(8) + ( ckink * two )**3

! evaluate <K^2>
     paux(7) = paux(7) + ( ckink * two )**2

! evaluate <N^2>
     paux(6) = paux(6) + ( sum(sgmt) / beta )**2

! evaluate <N^1>
     paux(5) = paux(5) + sum(sgmt) / beta

! evaluate spin magnetization: < Sz >
     do flvr=1,nband
         paux(4) = paux(4) + ( sgmt(flvr) - sgmt(flvr+nband) ) / beta
     enddo ! over flvr={1,nband} loop

! evaluate kinetic energy: ekin
     paux(3) = paux(3) - real(ckink * norbs) / beta

! evaluate potential energy: epot
     do flvr=1,norbs
         do i=1,flvr
             paux(2) = paux(2) + uumat(flvr,i) * ovlp(flvr,i) / beta
         enddo ! over i={1,flvr} loop
     enddo ! over flvr={1,norbs} loop

! evaluate total energy: etot
     paux(1) = paux(2) + paux(3)

! evaluate occupation matrix: < n_i >
     do flvr=1,norbs
         nmat(flvr) = nmat(flvr) + sgmt(flvr) / beta
     enddo ! over flvr={1,norbs} loop

! evaluate double occupation matrix: < n_i n_j >
     do flvr=1,norbs
         nnmat(flvr,:) = nnmat(flvr,:) + ovlp(flvr,:) / beta
     enddo ! over flvr={1,norbs} loop

     return
  end subroutine ctqmc_record_paux

!!
!! @sub ctqmc_record_nmat
!!
!! record the occupation matrix, double occupation matrix. note that it is
!! a empty subroutine. all features are implemented in ctqmc_record_paux()
!!
  subroutine ctqmc_record_nmat()
     implicit none

     return
  end subroutine ctqmc_record_nmat

!!========================================================================
!!>>> measure physical observables 2                                   <<<
!!========================================================================

!!
!! @sub ctqmc_record_gtau
!!
!! record the impurity green's function in imaginary time axis
!!
  subroutine ctqmc_record_gtau()
     use constants, only : dp, zero, one, two

     use control, only : isort
     use control, only : norbs
     use control, only : lemax, legrd
     use control, only : ntime
     use control, only : beta
     use context, only : index_s, index_e
     use context, only : time_s, time_e
     use context, only : rep_l
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

! evaluate step at first
     if ( isort == 1 ) step = real(ntime - 1) / beta
     if ( isort == 2 ) step = real(legrd - 1) / two

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

!-------------------------------------------------------------------------
! using normal representation
!-------------------------------------------------------------------------
                 STD_BLOCK: if ( isort == 1 ) then

! determine index for imaginary time
                     curr = nint( dtau * step ) + 1

! special tricks for the first point and the last point
                     if ( curr == 1 .or. curr == ntime ) then
                         maux = two * maux
                     endif ! back if ( curr == 1 .or. curr == ntime ) block

! record gtau, we normalize gtau in ctqmc_make_gtau() subroutine
                     gtau(curr, flvr, flvr) = gtau(curr, flvr, flvr) - maux

                 endif STD_BLOCK ! back if ( isort == 1 ) block
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-------------------------------------------------------------------------
! using legendre polynomial representation
!-------------------------------------------------------------------------
                 LEG_BLOCK: if ( isort == 2 ) then

! convert dtau in [0,\beta] to daux in [0,2]
                     daux = two * dtau / beta

! determine index for legendre polynomial interval
                     curr = nint( daux * step ) + 1

! record gtau, we normalize gtau in ctqmc_make_gtau() subroutine
                     CTQMC_FLALEG_LOOP: do fleg=1,lemax
                         dtau = sqrt(two * fleg - 1) * rep_l(curr,fleg)
                         gtau(fleg, flvr, flvr) = gtau(fleg, flvr, flvr) - maux * dtau
                     enddo CTQMC_FLALEG_LOOP ! over fleg={1,lemax} loop

                 endif LEG_BLOCK ! back if ( isort == 2 ) block
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

             enddo ! over ie={1,rank(flvr)} loop
         enddo ! over is={1,rank(flvr)} loop

     enddo CTQMC_FLAVOR_LOOP ! over flvr={1,norbs} loop

     return
  end subroutine ctqmc_record_gtau

!!
!! @sub ctqmc_record_ftau
!!
!! record the auxiliary correlation function in imaginary time axis,
!! F(\tau). latter, we will use it to compute the self-energy function
!!
  subroutine ctqmc_record_ftau()
     use constants, only : dp, zero, one, two

     use control, only : isort
     use control, only : norbs
     use control, only : lemax, legrd
     use control, only : ntime
     use control, only : beta
     use context, only : index_s, index_e
     use context, only : time_s, time_e
     use context, only : rep_l
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

!-------------------------------------------------------------------------
! using normal representation
!-------------------------------------------------------------------------

! calculate prefactor: pref
     call ctqmc_make_pref()

! evaluate step at first
     step = real(ntime - 1) / beta
     step = real(legrd - 1) / two

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

!-------------------------------------------------------------------------
! using legendre polynomial representation
!-------------------------------------------------------------------------
! convert dtau in [0,\beta] to daux in [0,2]
                 daux = two * dtau / beta

! determine index for legendre polynomial interval
                 curr = nint( daux * step ) + 1

! record ftau, we normalize ftau in ctqmc_make_ftau() subroutine
                 CTQMC_FLALEG_LOOP: do fleg=1,lemax
                     dtau = sqrt(two * fleg - 1) * rep_l(curr,fleg)
                     ftau(fleg, flvr, flvr) = ftau(fleg, flvr, flvr) - maux * dtau
                 enddo CTQMC_FLALEG_LOOP ! over fleg={1,lemax} loop


             enddo ! over ie={1,rank(flvr)} loop
         enddo ! over is={1,rank(flvr)} loop

     enddo CTQMC_FLAVOR_LOOP ! over flvr={1,norbs} loop

     return
  end subroutine ctqmc_record_ftau

!!
!! @sub ctqmc_record_grnf
!!
!! record the impurity green's function in matsubara frequency space
!!
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

!!========================================================================
!!>>> measure physical observables 3                                   <<<
!!========================================================================

!!
!! @sub ctqmc_record_kmat
!!
!! record the < k^2 > - < k >^2
!!
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

!!
!! @sub ctqmc_record_lmat
!!
!! record the fidelity susceptibility
!!
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

!!
!! @sub ctqmc_record_szpw
!!
!! record the powers of local magnetization
!!
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

!!========================================================================
!!>>> measure physical observables 4                                   <<<
!!========================================================================

!!
!! @sub ctqmc_record_schi
!!
!! record the spin-spin correlation function imaginary-time version
!!
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

!!
!! @sub ctqmc_record_sfom
!!
!! record the spin-spin correlation function matsubara frequency version
!!
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

!!
!! @sub ctqmc_record_ochi
!!
!! record the orbital-orbital correlation function imaginary-time version
!!
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

!!
!! @sub ctqmc_record_ofom
!!
!! record the orbital-orbital correlation function matsubara frequency version
!!
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

!!========================================================================
!!>>> measure physical observables 5                                   <<<
!!========================================================================

!!
!! @sub ctqmc_record_twop
!!
!! record the two-particle green's function
!!
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

!!
!! @sub ctqmc_record_vrtx
!!
!! record the two-particle green's function improved estimator is used to
!! improve the accuracy
!!
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

!!
!! @sub ctqmc_record_pair
!!
!! record the particle-particle pair susceptibility
!!
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

!!
!! @sub ctqmc_reduce_hist
!!
!! reduce the hist from all children processes
!!
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

!!
!! @sub ctqmc_reduce_prob
!!
!! reduce the prob from all children processes
!!
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

!!
!! @sub ctqmc_reduce_nmat
!!
!! reduce the nmat and nnmat from all children processes
!!
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

!!
!! @sub ctqmc_reduce_gtau
!!
!! reduce the gtau from all children processes
!!
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

!!
!! @sub ctqmc_reduce_ftau
!!
!! reduce the ftau from all children processes
!!
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

!!
!! @sub ctqmc_reduce_grnf
!!
!! reduce the grnf from all children processes
!!
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

!!
!! @sub ctqmc_reduce_kmat
!!
!! reduce the kmat and kkmat from all children processes
!!
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

!!
!! @sub ctqmc_reduce_lmat
!!
!! reduce the lmat, rmat, and lrmat from all children processes
!!
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

!!
!! @sub ctqmc_reduce_szpw
!!
!! reduce the szpow from all children processes
!!
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

!!
!! @sub ctqmc_reduce_schi
!!
!! reduce the schi and sschi from all children processes
!!
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

!!
!! @sub ctqmc_reduce_sfom
!!
!! reduce the ssfom from all children processes
!!
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

!!
!! @sub ctqmc_reduce_ochi
!! 
!! reduce the ochi and oochi from all children processes
!!
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

!!
!! @sub ctqmc_reduce_ofom
!!
!! reduce the oofom from all children processes
!!
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

!!
!! @sub ctqmc_reduce_twop
!!
!! reduce the g2_re_mpi and g2_im_mpi from all children processes
!!
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

!!
!! @sub ctqmc_reduce_vrtx
!!
!! reduce the h2_re_mpi and h2_im_mpi from all children processes
!!
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

!!
!! @sub ctqmc_reduce_pair
!!
!! reduce the ps_re_mpi and ps_im_mpi from all children processes
!!
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
