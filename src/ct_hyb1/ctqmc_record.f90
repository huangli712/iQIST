!!!-----------------------------------------------------------------------
!!! project : narcissus
!!! program : ctqmc_record_ac_f <<<---
!!!           ctqmc_record_hist
!!!           ctqmc_record_prob
!!!           ctqmc_record_paux
!!!           ctqmc_record_nmat <<<---
!!!           ctqmc_record_gtau
!!!           ctqmc_record_ftau
!!!           ctqmc_record_grnf
!!!           ctqmc_record_frnf
!!!           ctqmc_record_sig2 <<<---
!!!           ctqmc_record_kmat
!!!           ctqmc_record_lrmm
!!!           ctqmc_record_szpw <<<---
!!!           ctqmc_record_sp_t
!!!           ctqmc_record_sp_w
!!!           ctqmc_record_ch_t
!!!           ctqmc_record_ch_w <<<---
!!!           ctqmc_record_g2ph
!!!           ctqmc_record_g2pp <<<---
!!!           ctqmc_reduce_ac_f <<<---
!!!           ctqmc_reduce_hist
!!!           ctqmc_reduce_prob
!!!           ctqmc_reduce_paux
!!!           ctqmc_reduce_nmat <<<---
!!!           ctqmc_reduce_gtau
!!!           ctqmc_reduce_ftau
!!!           ctqmc_reduce_grnf
!!!           ctqmc_reduce_frnf
!!!           ctqmc_reduce_sig2 <<<---
!!!           ctqmc_reduce_kmat
!!!           ctqmc_reduce_lrmm
!!!           ctqmc_reduce_szpw <<<---
!!!           ctqmc_reduce_sp_t
!!!           ctqmc_reduce_sp_w
!!!           ctqmc_reduce_ch_t
!!!           ctqmc_reduce_ch_w <<<---
!!!           ctqmc_reduce_g2ph
!!!           ctqmc_reduce_g2pp <<<---
!!! source  : ctqmc_record.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           07/21/2017 by li huang (last modified)
!!! purpose : measure and collect physical observables produced by the
!!!           hybridization expansion version continuous time quantum
!!!           Monte Carlo (CTQMC) quantum impurity solver.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> measure autocorrelation function                                 <<<
!!========================================================================

!!
!! @sub ctqmc_record_ac_f
!!
!! record autocorrelation function for the total occupation number
!!
  subroutine ctqmc_record_ac_f()
     use constants, only : dp

     use control, only : norbs
     use control, only : ntime
     use control, only : beta

     use context, only : ac_v, ac_f

     implicit none

! local variables
! used to record how many times this subroutine were called
     integer, save :: starter = 0

! loop index
     integer  :: i
     integer  :: j

! memory position for current observable, which is used to measure the
! autocorrelation function
     integer  :: p

! total length of segments for all flavor channels
     real(dp) :: sgmt(norbs)

! calculate occupation status
     call cat_occupy_single(sgmt)

! record autocorrelation function: <A_{n} A_{n+k}>
! s1: increase the counter
     starter = starter + 1

! s2: determine memory location used to store the observable
     p = mod(starter, ntime)
     if ( p == 0 ) p = ntime

! s3: measure the autocorrelation function, starter is used to ensure that
! the vector ac_v is filled completely
     if ( starter > ntime ) then
         i = 0
         do j=p,ntime
             i = i + 1
             ac_f(i) = ac_f(i) + ac_v(p) * ac_v(j)
         enddo ! over j={p,ntime} loop
         do j=1,p-1
             i = i + 1
             ac_f(i) = ac_f(i) + ac_v(p) * ac_v(j)
         enddo ! over j={1,p-1} loop
         call s_assert( i == ntime )
     endif ! back if ( starter > ntime ) block

! s4: save the specified observable (the total occupation number) in ac_v
     ac_v(p) = sum(sgmt) / beta

! s5: here, ac_f(ntime + 1) is used to store the mean value for the total
! occupation number, it is very important
     ac_f(ntime + 1) = ac_f(ntime + 1) + sum(sgmt) / beta
     ac_f(ntime + 2) = ac_f(ntime + 2) + ( sum(sgmt) / beta )**2

     return
  end subroutine ctqmc_record_ac_f

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

! if ckink == 0, we record its count in hist(mkink)
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

! index of atomic eigenstate
     integer :: pstat

! current atomic eigenstate from segment representation
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
     call ctqmc_make_fock(norbs, pstat, state)

! accumulate the data
     prob(pstat) = prob(pstat) + one

     return
  end subroutine ctqmc_record_prob

!!
!! @sub ctqmc_record_paux
!!
!! record some auxiliary physical observables. the occupation number
!! and double occupation matrix are measured at the same time in order
!! to save the computational time
!!
  subroutine ctqmc_record_paux()
     use constants, only : dp
     use constants, only : two

     use control, only : nband, norbs
     use control, only : beta

     use context, only : ckink
     use context, only : paux
     use context, only : nimp, nmat
     use context, only : umat

     implicit none

! local variables
! loop index for flavor channel
     integer  :: flvr
     integer  :: i

! total length of segments for all flavor channels
     real(dp) :: sgmt(norbs)

! overlap length of segments between different flavors
     real(dp) :: ovlp(norbs,norbs)

! prepare sgmt array
     call cat_occupy_single(sgmt)

! prepare ovlp array
     call cat_occupy_double(ovlp)

!-------------------------------------------------------------------------

     CALC_PAUX: BLOCK

! evaluate < K^4 >
         paux(9) = paux(9) + ( ckink * two )**4

! evaluate < K^3 >
         paux(8) = paux(8) + ( ckink * two )**3

! evaluate < K^2 >
         paux(7) = paux(7) + ( ckink * two )**2

! evaluate < N^2 >
         paux(6) = paux(6) + ( sum(sgmt) / beta )**2

! evaluate < N^1 >
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
                 paux(2) = paux(2) + umat(flvr,i) * ovlp(flvr,i) / beta
             enddo ! over i={1,flvr} loop
         enddo ! over flvr={1,norbs} loop

! evaluate total energy: etot
         paux(1) = paux(2) + paux(3)

     END BLOCK CALC_PAUX

!-------------------------------------------------------------------------

     CALC_NMAT: BLOCK

! evaluate occupation matrix: < n_i >
         nimp = nimp + sgmt / beta

! evaluate double occupation matrix: < n_i n_j >
         nmat = nmat + ovlp / beta

     END BLOCK CALC_NMAT

     return
  end subroutine ctqmc_record_paux

!!
!! @sub ctqmc_record_nmat
!!
!! record the occupation number, double occupation matrix. actually it is
!! an empty subroutine. this required feature is already implemented in
!! the ctqmc_record_paux() subroutine
!!
  subroutine ctqmc_record_nmat()
     implicit none

     CONTINUE

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
     use constants, only : dp
     use constants, only : zero, one, two

     use control, only : isort
     use control, only : norbs
     use control, only : lemax, legrd
     use control, only : svmax, svgrd
     use control, only : ntime
     use control, only : beta

     use context, only : index_s, index_e
     use context, only : time_s, time_e
     use context, only : rep_l, rep_s
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

! loop index for orthogonal polynomial
     integer  :: fleg
     integer  :: fsvd

! index for imaginary time \tau
     integer  :: curr

! used to store the element of mmat matrix
     real(dp) :: maux

! imaginary time for start and end points
     real(dp) :: taus
     real(dp) :: taue

! distance betweem taus and taue
     real(dp) :: dtau
     real(dp) :: daux

! interval for imaginary time slice
     real(dp) :: step

! evaluate step at first
     select case ( isort )

         case (1)
             step = real(ntime - 1) / beta

         case (2)
             step = real(legrd - 1) / two

         case (3)
             step = real(svgrd - 1) / two

     end select

     FLVR_CYCLE: do flvr=1,norbs

! get imaginary time value for segments
         do is=1,rank(flvr)
             taus = time_s( index_s(is, flvr), flvr )

! get imaginary time value for segments
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
! using standard representation
!-------------------------------------------------------------------------
                 STD_BLOCK: if ( isort == 1 ) then

! determine index for imaginary time
                     curr = nint( dtau * step ) + 1

! special tricks for the first point and the last point
                     if ( curr == 1 .or. curr == ntime ) then
                         maux = two * maux
                     endif ! back if ( curr == 1 .or. curr == ntime ) block

! record gtau, we normalize gtau in ctqmc_tran_gtau() subroutine
                     gtau(curr, flvr, flvr) = gtau(curr, flvr, flvr) - maux

                 endif STD_BLOCK ! back if ( isort == 1 ) block
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-------------------------------------------------------------------------
! using legendre orthogonal polynomial representation
!-------------------------------------------------------------------------
                 LEG_BLOCK: if ( isort == 2 ) then

! convert dtau in [0,\beta] to daux in [0,2]
                     daux = two * dtau / beta

! determine index for legendre orthogonal polynomial interval
                     curr = nint( daux * step ) + 1

! special tricks for the first point and the last point
                     if ( curr == 1 .or. curr == legrd ) then
                         maux = two * maux
                     endif ! back if ( curr == 1 .or. curr == legrd ) block

! record gtau, we normalize gtau in ctqmc_tran_gtau() subroutine
                     LEG_CYCLE: do fleg=1,lemax
                         dtau = sqrt(two * fleg - 1) * rep_l(curr,fleg)
                         gtau(fleg, flvr, flvr) = gtau(fleg, flvr, flvr) - maux * dtau
                     enddo LEG_CYCLE ! over fleg={1,lemax} loop

                 endif LEG_BLOCK ! back if ( isort == 2 ) block
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-------------------------------------------------------------------------
! using svd orthogonal polynomial representation
!-------------------------------------------------------------------------
                 SVD_BLOCK: if ( isort == 3 ) then

! convert dtau in [0,\beta] to daux in [-1,1]
                     daux = two * dtau / beta - one

! determine index for svd orthogonal polynomial interval
                     call s_svd_point(daux, step, curr)

! record gtau, we normalize gtau in ctqmc_tran_gtau() subroutine
                     SVD_CYCLE: do fsvd=1,svmax
                         dtau = rep_s(curr,fsvd)
                         gtau(fsvd, flvr, flvr) = gtau(fsvd, flvr, flvr) - maux * dtau
                     enddo SVD_CYCLE ! over fsvd={1,svmax} loop

                 endif SVD_BLOCK ! back if ( isort == 3 ) block
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

             enddo ! over ie={1,rank(flvr)} loop
         enddo ! over is={1,rank(flvr)} loop

     enddo FLVR_CYCLE ! over flvr={1,norbs} loop

     return
  end subroutine ctqmc_record_gtau

!!
!! @sub ctqmc_record_ftau
!!
!! record the auxiliary correlation function in imaginary time axis.
!! latter, we will use it to compute the self-energy function
!!
  subroutine ctqmc_record_ftau()
     use constants, only : dp
     use constants, only : zero, one, two

     use control, only : isort
     use control, only : norbs
     use control, only : lemax, legrd
     use control, only : svmax, svgrd
     use control, only : ntime
     use control, only : beta

     use context, only : index_s, index_e
     use context, only : time_s, time_e
     use context, only : rep_l, rep_s
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

! loop index for orthogonal polynomial
     integer  :: fleg
     integer  :: fsvd

! index for imaginary time \tau
     integer  :: curr

! used to store the element of mmat matrix
     real(dp) :: maux

! imaginary time for start and end points
     real(dp) :: taus
     real(dp) :: taue

! distance betweem taus and taue
     real(dp) :: dtau
     real(dp) :: daux

! interval for imaginary time slice
     real(dp) :: step

! evaluate step at first
     select case ( isort )

         case (1)
             step = real(ntime - 1) / beta

         case (2)
             step = real(legrd - 1) / two

         case (3)
             step = real(svgrd - 1) / two

     end select

! calculate prefactor: pref
     call ctqmc_make_pref()

     FLVR_CYCLE: do flvr=1,norbs

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

!-------------------------------------------------------------------------
! using standard representation
!-------------------------------------------------------------------------
                 STD_BLOCK: if ( isort == 1 ) then

! determine index for imaginary time
                     curr = nint( dtau * step ) + 1

! special tricks for the first point and the last point
                     if ( curr == 1 .or. curr == ntime ) then
                         maux = two * maux
                     endif ! back if ( curr == 1 .or. curr == ntime ) block

! record ftau, we normalize ftau in ctqmc_tran_gtau() subroutine
                     ftau(curr, flvr, flvr) = ftau(curr, flvr, flvr) - maux

                 endif STD_BLOCK ! back if ( isort == 1 ) block
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-------------------------------------------------------------------------
! using legendre orthogonal polynomial representation
!-------------------------------------------------------------------------
                 LEG_BLOCK: if ( isort == 2 ) then

! convert dtau in [0,\beta] to daux in [0,2]
                     daux = two * dtau / beta

! determine index for legendre orthogonal polynomial interval
                     curr = nint( daux * step ) + 1

! special tricks for the first point and the last point
                     if ( curr == 1 .or. curr == legrd ) then
                         maux = two * maux
                     endif ! back if ( curr == 1 .or. curr == legrd ) block

! record ftau, we normalize ftau in ctqmc_tran_gtau() subroutine
                     LEG_CYCLE: do fleg=1,lemax
                         dtau = sqrt(two * fleg - 1) * rep_l(curr,fleg)
                         ftau(fleg, flvr, flvr) = ftau(fleg, flvr, flvr) - maux * dtau
                     enddo LEG_CYCLE ! over fleg={1,lemax} loop

                 endif LEG_BLOCK ! back if ( isort == 2 ) block
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-------------------------------------------------------------------------
! using svd orthogonal polynomial representation
!-------------------------------------------------------------------------
                 SVD_BLOCK: if ( isort == 3 ) then

! convert dtau in [0,\beta] to daux in [-1,1]
                     daux = two * dtau / beta - one

! determine index for svd orthogonal polynomial interval
                     call s_svd_point(daux, step, curr)

! record ftau, we normalize ftau in ctqmc_tran_gtau() subroutine
                     SVD_CYCLE: do fsvd=1,svmax
                         dtau = rep_s(curr,fsvd)
                         ftau(fsvd, flvr, flvr) = ftau(fsvd, flvr, flvr) - maux * dtau
                     enddo SVD_CYCLE ! over fsvd={1,svmax} loop

                 endif SVD_BLOCK ! back if ( isort == 3 ) block
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

             enddo ! over ie={1,rank(flvr)} loop
         enddo ! over is={1,rank(flvr)} loop

     enddo FLVR_CYCLE ! over flvr={1,norbs} loop

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

! only the first nfreq points of grnf are modified
     do flvr=1,norbs
         do ifrq=1,nfreq
             grnf(ifrq, flvr, flvr) = grnf(ifrq, flvr, flvr) + gmat(ifrq, flvr, flvr)
         enddo ! over ifrq={1,nfreq} loop
     enddo ! over flvr={1,norbs} loop

     return
  end subroutine ctqmc_record_grnf

!!
!! @sub ctqmc_record_frnf
!!
!! record the auxiliary correlation function in matsubara frequency space.
!! the required feature is implemented in the ctqmc_make_hub2() subroutine
!!
  subroutine ctqmc_record_frnf()
     implicit none

     CONTINUE

     return
  end subroutine ctqmc_record_frnf

!!
!! @sub ctqmc_record_sig2
!!
!! record the self-energy function in matsubara frequency space.
!! the required feature is implemented in the ctqmc_make_hub2() subroutine
!!
  subroutine ctqmc_record_sig2()
     implicit none

     CONTINUE

     return
  end subroutine ctqmc_record_sig2

!!========================================================================
!!>>> measure physical observables 3                                   <<<
!!========================================================================

!!
!! @sub ctqmc_record_kmat
!!
!! record the kinetic energy fluctuation < k^2 > - < k >^2
!!
  subroutine ctqmc_record_kmat()
     use constants, only : dp

     use control, only : isobs
     use control, only : norbs

     use context, only : knop, kmat
     use context, only : rank

     implicit none

! local variables
! loop index for flavor channel
     integer :: i
     integer :: j

! check whether there is conflict
     call s_assert2( btest(isobs, 1), 'in ctqmc_record_kmat' )

! since rank means the number of operator pairs, so we have to multiply
! it with two
     do i=1,norbs
         knop(i) = knop(i) + rank(i) * 2.0_dp
     enddo ! over i={1,norbs} loop

     do j=1,norbs
         do i=1,norbs
             kmat(i,j) = kmat(i,j) + rank(i) * rank(j) * 4.0_dp
         enddo ! over i={1,norbs} loop
     enddo ! over j={1,norbs} loop

     return
  end subroutine ctqmc_record_kmat

!!
!! @sub ctqmc_record_lrmm
!!
!! record the fidelity susceptibility
!!
  subroutine ctqmc_record_lrmm()
     use constants, only : dp
     use constants, only : zero, one, two

     use control, only : isobs
     use control, only : norbs
     use control, only : beta

     use context, only : index_s, index_e
     use context, only : time_s, time_e
     use context, only : lnop, rnop, lrmm
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
     call s_assert2( btest(isobs, 2), 'in ctqmc_record_lrmm' )

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
     lnop = lnop + kl
     rnop = rnop + kr

! add contribution to < k_l k_r >
     do flvr=1,norbs
         do i=1,norbs
             lrmm(i,flvr) = lrmm(i,flvr) + kl(i) * kr(flvr)
         enddo ! over i={1,norbs} loop
     enddo ! over flvr={1,norbs} loop

     return
  end subroutine ctqmc_record_lrmm

!!
!! @sub ctqmc_record_szpw
!!
!! record the powers of local magnetization, which will be used to compute
!! the binder cumulant
!!
  subroutine ctqmc_record_szpw()
     use constants, only : dp
     use constants, only : zero, two

     use control, only : isobs
     use control, only : nband, norbs
     use control, only : ntime
     use control, only : beta

     use context, only : tmesh
     use context, only : szpw

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
!     sint = 1/\beta \int^{\beta}_{0} Sz(\tau) d\tau
     real(dp) :: sint

! used to record occupations for current flavor channel and time
     real(dp) :: oaux(ntime,norbs)

! orbital-resolved Sz(\tau)
     real(dp) :: saux(ntime,nband)

! check whether there is conflict
     call s_assert2( btest(isobs, 3), 'in ctqmc_record_szpw' )

! calculate oaux, obtain occupation status
! calculate saux, obtain Sz(\tau)
     oaux = zero
     saux = zero
     TIME_CYCLE: do i=1,ntime
         do f1=1,norbs
             call cat_occupy_status(f1, tmesh(i), oaux(i,f1))
         enddo ! over f1={1,norbs} loop

         do f2=1,nband
             saux(i,f2) = oaux(i,f2) - oaux(i,f2+nband)
         enddo ! over f2={1,nband} loop
     enddo TIME_CYCLE ! over i={1,ntime} loop

! accumulate szpw(1:4,1:nband)
! calculate \delta \tau
     step = ( tmesh(2) - tmesh(1) ) / two
     FLVR_CYCLE: do f2=1,nband
! calculate sint using trapezoid algorithm
         sint = zero
         do i=1,ntime-1
             sint = sint + ( saux(i,f2) + saux(i+1,f2) ) * step
         enddo ! over i={1,ntime-1} loop
         sint = sint / beta
! record the data
         szpw(1,f2) = szpw(1,f2) + sint**1.0
         szpw(2,f2) = szpw(2,f2) + sint**2.0
         szpw(3,f2) = szpw(3,f2) + sint**3.0
         szpw(4,f2) = szpw(4,f2) + sint**4.0
     enddo FLVR_CYCLE ! over f2={1,nband} loop

! accumulate szpw(1:4,nband+1)
! here we consider the contribution from all flavors
     sint = zero
     do i=1,ntime-1
         sint = sint + ( sum( saux(i,:) ) + sum( saux(i+1,:) ) ) * step
     enddo ! over i={1,ntime-1} loop
     sint = sint / beta
! record the data
     szpw(1,nband+1) = szpw(1,nband+1) + sint**1.0
     szpw(2,nband+1) = szpw(2,nband+1) + sint**2.0
     szpw(3,nband+1) = szpw(3,nband+1) + sint**3.0
     szpw(4,nband+1) = szpw(4,nband+1) + sint**4.0

     return
  end subroutine ctqmc_record_szpw

!!========================================================================
!!>>> measure physical observables 4                                   <<<
!!========================================================================

!!
!! @sub ctqmc_record_sp_t
!!
!! record the spin-spin correlation function in imaginary time axis
!!
  subroutine ctqmc_record_sp_t()
     use constants, only : dp
     use constants, only : zero

     use spring, only : spring_sfmt_stream

     use control, only : issus
     use control, only : nband, norbs
     use control, only : ntime

     use context, only : tmesh
     use context, only : schi, sp_t

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

! occupation status for current flavor channel and time
     real(dp) :: oaux(ntime,norbs)

! check whether there is conflict
     call s_assert2( btest(issus, 1), 'in ctqmc_record_sp_t' )

! calculate oaux, obtain occupation status
     oaux = zero
     TIME_CYCLE: do i=1,ntime
         do f1=1,norbs
             call cat_occupy_status(f1, tmesh(i), oaux(i,f1))
         enddo ! over f1={1,norbs} loop
     enddo TIME_CYCLE ! over i={1,ntime} loop
     oaux = oaux / real(num_try)

! calculate schi and sp_t
     FLVR_CYCLE: do f1=1,nband
         do i=1,num_try
             m = ceiling( spring_sfmt_stream() * ntime )
             if ( oaux(m,f1) > zero ) then
! when n - m + ntime \in [ntime - m + 1, ntime]
                 do n=1,m
                     associate ( val => oaux(n,f1) - oaux(n,f1+nband) )
                         schi(n-m+ntime) = schi(n-m+ntime) + val
                         sp_t(n-m+ntime,f1) = sp_t(n-m+ntime,f1) + val
                     end associate
                 enddo ! over n={1,m} loop
! when n - m \in [1, ntime - m]
                 do n=m+1,ntime
                     associate ( val => oaux(n,f1) - oaux(n,f1+nband) )
                         schi(n-m) = schi(n-m) + val
                         sp_t(n-m,f1) = sp_t(n-m,f1) + val
                     end associate
                 enddo ! over n={m+1,ntime} loop
             endif ! back if ( oaux(m,f1) > zero ) block

             if ( oaux(m,f1+nband) > zero ) then
! when n - m + ntime \in [ntime - m + 1, ntime]
                 do n=1,m
                     associate ( val => oaux(n,f1+nband) - oaux(n,f1) )
                         schi(n-m+ntime) = schi(n-m+ntime) + val
                         sp_t(n-m+ntime,f1) = sp_t(n-m+ntime,f1) + val
                     end associate
                 enddo ! over n={1,m} loop
! when n - m \in [1, ntime - m]
                 do n=m+1,ntime
                     associate ( val => oaux(n,f1+nband) - oaux(n,f1) )
                         schi(n-m) = schi(n-m) + val
                         sp_t(n-m,f1) = sp_t(n-m,f1) + val
                     end associate
                 enddo ! over n={m+1,ntime} loop
             endif ! back if ( oaux(m,f1+nband) > zero ) block
         enddo ! over i={1,num_try} loop
     enddo FLVR_CYCLE ! over f1={1,nband} loop

     return
  end subroutine ctqmc_record_sp_t

!!
!! @sub ctqmc_record_sp_w
!!
!! record the spin-spin correlation function in matsubara frequency axis
!!
  subroutine ctqmc_record_sp_w()
     use constants, only : dp
     use constants, only : pi, zero, one, two, czi

     use control, only : issus
     use control, only : nband, norbs
     use control, only : nbfrq
     use control, only : beta

     use context, only : index_s, index_e
     use context, only : exp_s, exp_e
     use context, only : sp_w
     use context, only : rank

     implicit none

! local variables
! loop index for flavor channel
     integer  :: f1
     integer  :: f2

! loop index for operators
     integer  :: it

! the first bosonic frequency
     complex(dp) :: dw

! occupation status for current flavor channel at \tau = 0
     real(dp) :: oaux(norbs)

! total length of segments
     real(dp) :: sgmt(norbs)

! bosonic frequency mesh
     complex(dp) :: mesh(2:nbfrq)

! matsubara frequency exponents for creation operators
     complex(dp) :: bexp_s(2:nbfrq)

! matsubara frequency exponents for annihilation operators
     complex(dp) :: bexp_e(2:nbfrq)

! check whether there is conflict
     call s_assert2( btest(issus, 3), 'in ctqmc_record_sp_w' )

! build bosonic frequency mesh, zero frequency is not included
     dw = czi * two * pi / beta
     mesh = dw
     call s_cumsum_z(nbfrq - 1, mesh, mesh)

! calculate oaux, obtain occupation status
     do f1=1,norbs
         call cat_occupy_status(f1, zero, oaux(f1))
     enddo ! over i={1,norbs} loop

! calculate sgmt, obtain total length of segments
     call cat_occupy_single(sgmt)

! calculate sp_w, it must be real
! < Sz(t)Sz(0) > = < ( nu(t) - nd(t) ) * ( nu(0) - nd(0) ) >
     FLVR_CYCLE: do f1=1,nband
         f2 = f1 + nband

!
! note:
!
! the contribution from oaux(f1) = one and oaux(f2) = one is zero, and
! the contribution from oaux(f1) = zero and oaux(f2) = zero is also zero
!
! here oaux(f1) = one; oaux(f2) = zero
         if ( oaux(f1) > zero .and. oaux(f2) < one ) then
! + nu(t)nu(0) term
!-------------------------------------------------------------------------
             sp_w(1,f1) = sp_w(1,f1) + sgmt(f1)
! loop over the segments
! here we do not calculate bexp_s and bexp_e directly. we just utilize
! the available data in exp_s and exp_e. so that the required exponent
! functions are obtained via only a few complex number multiplications
             do it=1,rank(f1)
                 dw = exp_s(1,index_s(it, f1), f1)
                 bexp_s = dw * exp_s(1:nbfrq-1,index_s(it, f1), f1)
                 dw = exp_e(1,index_e(it, f1), f1)
                 bexp_e = dw * exp_e(1:nbfrq-1,index_e(it, f1), f1)
                 sp_w(2:,f1) = sp_w(2:,f1) + real( ( bexp_e - bexp_s ) / mesh )
             enddo ! over do it={1,rank(f1)} loop

! - nd(t)nu(0) term
!-------------------------------------------------------------------------
             sp_w(1,f1) = sp_w(1,f1) - sgmt(f2)
! loop over the segments
! here we do not calculate bexp_s and bexp_e directly. we just utilize
! the available data in exp_s and exp_e. so that the required exponent
! functions are obtained via only a few complex number multiplications
             do it=1,rank(f2)
                 dw = exp_s(1,index_s(it, f2), f2)
                 bexp_s = dw * exp_s(1:nbfrq-1,index_s(it, f2), f2)
                 dw = exp_e(1,index_e(it, f2), f2)
                 bexp_e = dw * exp_e(1:nbfrq-1,index_e(it, f2), f2)
                 sp_w(2:,f1) = sp_w(2:,f1) - real( ( bexp_e - bexp_s ) / mesh )
             enddo ! over do it={1,rank(f2)} loop
         endif ! back if ( oaux(f1) > zero .and. oaux(f2) < one ) block

! here oaux(f2) = one; oaux(f1) = zero
         if ( oaux(f2) > zero .and. oaux(f1) < one ) then
! - nu(t)nd(0) term
!-------------------------------------------------------------------------
             sp_w(1,f1) = sp_w(1,f1) - sgmt(f1)
! loop over the segments
! here we do not calculate bexp_s and bexp_e directly. we just utilize
! the available data in exp_s and exp_e. so that the required exponent
! functions are obtained via only a few complex number multiplications
             do it=1,rank(f1)
                 dw = exp_s(1,index_s(it, f1), f1)
                 bexp_s = dw * exp_s(1:nbfrq-1,index_s(it, f1), f1)
                 dw = exp_e(1,index_e(it, f1), f1)
                 bexp_e = dw * exp_e(1:nbfrq-1,index_e(it, f1), f1)
                 sp_w(2:,f1) = sp_w(2:,f1) - real( ( bexp_e - bexp_s ) / mesh )
             enddo ! over do it={1,rank(f1)} loop

! + nd(t)nd(0) term
!-------------------------------------------------------------------------
             sp_w(1,f1) = sp_w(1,f1) + sgmt(f2)
! loop over the segments
! here we do not calculate bexp_s and bexp_e directly. we just utilize
! the available data in exp_s and exp_e. so that the required exponent
! functions are obtained via only a few complex number multiplications
             do it=1,rank(f2)
                 dw = exp_s(1,index_s(it, f2), f2)
                 bexp_s = dw * exp_s(1:nbfrq-1,index_s(it, f2), f2)
                 dw = exp_e(1,index_e(it, f2), f2)
                 bexp_e = dw * exp_e(1:nbfrq-1,index_e(it, f2), f2)
                 sp_w(2:,f1) = sp_w(2:,f1) + real( ( bexp_e - bexp_s ) / mesh )
             enddo ! over do it={1,rank(f2)} loop
         endif ! back if ( oaux(f2) > zero .and. oaux(f1) < one ) block

     enddo FLVR_CYCLE ! over f1={1,nband} loop

     return
  end subroutine ctqmc_record_sp_w

!!
!! @sub ctqmc_record_ch_t
!!
!! record the charge-charge correlation function in imaginary time axis
!!
  subroutine ctqmc_record_ch_t()
     use constants, only : dp
     use constants, only : zero, one, two

     use spring, only : spring_sfmt_stream

     use control, only : issus
     use control, only : norbs
     use control, only : ntime

     use context, only : tmesh
     use context, only : cchi, ch_t

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

! factor for orbital symmetry
     real(dp) :: fa

! occupation status for current flavor channel and time
     real(dp) :: oaux(ntime,norbs)

! check whether there is conflict
     call s_assert2( btest(issus, 2), 'in ctqmc_record_ch_t' )

! calculate oaux, obtain occupation status
     oaux = zero
     TIME_CYCLE: do i=1,ntime
         do f1=1,norbs
             call cat_occupy_status(f1, tmesh(i), oaux(i,f1))
         enddo ! over f1={1,norbs} loop
     enddo TIME_CYCLE ! over i={1,ntime} loop
     oaux = oaux / real(num_try)

! calculate cchi and ch_t
     FLVR_CYCLE: do f1=1,norbs
         do f2=1,f1
! evaluate the symmetry factor
             if ( f1 /= f2 ) then
                 fa = two
             else
                 fa = one
             endif ! back if ( f1 /= f2 ) block
             do i=1,num_try
                 m = ceiling( spring_sfmt_stream() * ntime )
                 if ( oaux(m,f2) > zero ) then
! when n - m + ntime \in [ntime - m + 1, ntime]
                     do n=1,m
                         cchi(n-m+ntime) = cchi(n-m+ntime) + fa * oaux(n,f1)
                         ch_t(n-m+ntime,f2,f1) = ch_t(n-m+ntime,f2,f1) + oaux(n,f1)
                     enddo ! over n={1,m} loop
! when n - m \in [1, ntime - m]
                     do n=m+1,ntime
                         cchi(n-m) = cchi(n-m) + fa * oaux(n,f1)
                         ch_t(n-m,f2,f1) = ch_t(n-m,f2,f1) + oaux(n,f1)
                     enddo ! over n={m+1,ntime} loop
                 endif ! back if ( oaux(m,f2) > zero ) block
             enddo ! over i={1,num_try} loop
             if ( f1 /= f2 ) then ! consider the symmetry
                 ch_t(:,f1,f2) = ch_t(:,f2,f1)
             endif ! back if ( f1 /= f2 ) block
         enddo ! over f2={1,f1} loop
     enddo FLVR_CYCLE ! over f1={1,norbs} loop

     return
  end subroutine ctqmc_record_ch_t

!!
!! @sub ctqmc_record_ch_w
!!
!! record the charge-charge correlation function in matsubara frequency axis
!!
  subroutine ctqmc_record_ch_w()
     use constants, only : dp
     use constants, only : pi, zero, two, czi

     use control, only : issus
     use control, only : norbs
     use control, only : nbfrq
     use control, only : beta

     use context, only : index_s, index_e
     use context, only : exp_s, exp_e
     use context, only : ch_w
     use context, only : rank

     implicit none

! local variables
! loop index for flavor channel
     integer  :: f1
     integer  :: f2

! loop index for operators
     integer  :: it

! the first bosonic frequency
     complex(dp) :: dw

! occupation status for current flavor channel at \tau = 0
     real(dp) :: oaux(norbs)

! total length of segments
     real(dp) :: sgmt(norbs)

! bosonic frequency mesh
     complex(dp) :: mesh(2:nbfrq)

! matsubara frequency exponents for creation operators
     complex(dp) :: bexp_s(2:nbfrq)

! matsubara frequency exponents for annihilation operators
     complex(dp) :: bexp_e(2:nbfrq)

! check whether there is conflict
     call s_assert2( btest(issus, 4), 'in ctqmc_record_ch_w' )

! build bosonic frequency mesh, zero frequency is not included
     dw = czi * two * pi / beta
     mesh = dw
     call s_cumsum_z(nbfrq - 1, mesh, mesh)

! calculate oaux, obtain occupation status
     do f1=1,norbs
         call cat_occupy_status(f1, zero, oaux(f1))
     enddo ! over i={1,norbs} loop

! calculate sgmt, obtain total length of segments
     call cat_occupy_single(sgmt)

! calculate ch_w, it must be real
     FLVR_CYCLE: do f1=1,norbs
         do f2=1,f1
             if ( oaux(f2) > zero ) then
! special treatment for the first frequency
                 ch_w(1,f2,f1) = ch_w(1,f2,f1) + sgmt(f1)
! loop over the segments
! here we do not calculate bexp_s and bexp_e directly. we just utilize
! the available data in exp_s and exp_e. so that the required exponent
! functions are obtained via only a few complex number multiplications
                 do it=1,rank(f1)
                     dw = exp_s(1,index_s(it, f1), f1)
                     bexp_s = dw * exp_s(1:nbfrq-1,index_s(it, f1), f1)
                     dw = exp_e(1,index_e(it, f1), f1)
                     bexp_e = dw * exp_e(1:nbfrq-1,index_e(it, f1), f1)
                     ch_w(2:,f2,f1) = ch_w(2:,f2,f1) + real( ( bexp_e - bexp_s ) / mesh )
                 enddo ! over do it={1,rank(f1)} loop
             endif ! back if ( oaux(f2) > zero ) block
             if ( f1 /= f2 ) then ! consider the symmetry
                 ch_w(:,f1,f2) = ch_w(:,f2,f1)
             endif ! back if ( f1 /= f2 ) block
         enddo ! over f2={1,f1} loop
     enddo FLVR_CYCLE ! over f1={1,norbs} loop

     return
  end subroutine ctqmc_record_ch_w

!!========================================================================
!!>>> measure physical observables 5                                   <<<
!!========================================================================

!!
!! @sub ctqmc_record_g2ph
!!
!! record the two-particle green's and vertex functions in the ph channel.
!! here improved estimator is used to improve the accuracy
!!
!! notation:
!!
!!     G^{(2)}_{\alpha\beta\gamma\delta} (\tau_1, \tau_2, \tau_3, \tau_4)
!!         = \langle T_\tau
!!               c_{\alpha} (\tau_1) c^{\dagger}_{\beta} (\tau_2)
!!               c_{\gamma} (\tau_3) c^{\dagger}_{\delta} (\tau_4)
!!           \rangle
!!
!!     G^{(2)}_{\alpha\beta\gamma\delta,ph} (\nu, \nu', \omega)
!!         = \langle
!!               c_{\alpha} (\nu + \omega) c^{*}_{\beta} (\nu)
!!               c_{\gamma} (\nu') c^{*}_{\delta} (\nu' + \omega)
!!           \rangle
!!
!!     \alpha, \beta, \gamma, \delta: orbital index
!!     \nu and \nu': fermionic matsubara frequency
!!     \omega: bosonic matsubara frequency
!!
!!       in             out
!!        \              /
!!     v+w \            / v'+w
!!          \          /
!!         i \--------/ l
!!           |        |
!!           |        |
!!           |        |
!!         j /--------\ k
!!          /          \
!!       v /            \ v'
!!        /              \
!!       out             in
!!
  subroutine ctqmc_record_g2ph()
     use control, only : isort
     use control, only : isvrt

     implicit none

! check whether there is conflict
! this subroutine is only designed for the particle-hole channel
     call s_assert2( btest(isvrt, 1) .or. btest(isvrt, 2), 'in ctqmc_record_g2ph' )

! you can not calculate the AABB and ABBA components at the same time
     call s_assert2( .not. ( btest(isvrt, 1) .and. btest(isvrt, 2) ), 'in ctqmc_record_g2ph' )

     select case ( isort )

         case (1)
             call cat_record_g2ph_std()

         case (2)
             call cat_record_g2ph_leg()

         case (3)
             call cat_record_g2ph_svd()

         case default
             call s_print_error('ctqmc_record_g2ph','this feature is not implemented')

     end select

     return
  end subroutine ctqmc_record_g2ph

!!
!! @sub cat_record_g2ph_std
!!
!! record the two-particle green's and vertex functions in the ph channel.
!! here improved estimator is used to improve the accuracy
!!
!! note:
!!
!!     we try to measure the two-particle green's and vertex functions in
!!     the particle-hole channel and matsubara frequency representation
!!     in this subroutine. in order to simplify the calculations, we just
!!     consider the block structure of G^{(2)}
!!
!!     G^{(2)}_{abcd,AABB,ph} (\nu, \nu', \omega) = \frac{1}{\beta}
!!         \langle
!!             \sum^{K_A}_{ij=1} \sum^{K_B}_{kl=1}
!!             ( M^{A}_{ij} M^{B}_{kl} - \delta_{AB} M^{A}_{il} M^{B}_{kj} )
!!             exp [ i (\nu + \omega) \tau'_i ]
!!             exp [ -i \nu \tau_j ]
!!             exp [ i \nu' \tau'_k ]
!!             exp [ -i (\nu' + \omega) \tau_l ]
!!             \delta_{a,i} \delta_{b,j} \delta_{c,k} \delta_{d,l}
!!         \rangle
!!
!!     G^{(2)}_{abcd,ABBA,ph} (\nu, \nu', \omega) = \frac{1}{\beta}
!!         \langle
!!             \sum^{K_A}_{ij=1} \sum^{K_B}_{kl=1}
!!             ( \delta_{AB} M^{A}_{ij} M^{B}_{kl} - M^{A}_{il} M^{B}_{kj} )
!!             exp [ i (\nu + \omega) \tau'_i ]
!!             exp [ -i \nu \tau_j ]
!!             exp [ i \nu' \tau'_k ]
!!             exp [ -i (\nu' + \omega) \tau_l ]
!!             \delta_{a,i} \delta_{b,j} \delta_{c,k} \delta_{d,l}
!!         \rangle
!!
!!     \tau'_i and \tau'_k: imaginary time for annihilation operators
!!     \tau_j and \tau_l: imaginary time for creation operators
!!     \nu and \nu': fermionic matsubara frequency
!!     \omega: bosonic matsubara frequency
!!
  subroutine cat_record_g2ph_std()
     use constants, only : dp
     use constants, only : czero

     use control, only : isvrt
     use control, only : norbs
     use control, only : nffrq, nbfrq
     use control, only : beta

     use context, only : g2ph
     use context, only : h2ph
     use context, only : rank, pref
     use context, only : mmat

     implicit none

! local variables
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

! loop indices for start and end points
     integer  :: is
     integer  :: ie

! used to store the element of mmat matrix
     real(dp) :: maux
     real(dp) :: naux

! dummy complex(dp) variables, used to calculate the g2ph and h2ph
     complex(dp) :: zg
     complex(dp) :: zh

! exp [i \omega_n \tau_s] and exp [i \omega_n \tau_e]
     complex(dp), allocatable :: caux1(:,:)
     complex(dp), allocatable :: caux2(:,:)

! \sum_{ij=1} exp [i \omega_m \tau'_i ] M_{ij} exp [ i \omega_n \tau_j ]
! where m and n are the first two frequency indices for g2aux and h2aux
     complex(dp), allocatable :: g2aux(:,:,:)
     complex(dp), allocatable :: h2aux(:,:,:)

! evaluate nfaux, determine the size of g2aux and h2aux
     nfaux = nffrq + nbfrq - 1

! allocate memory for caux1 and caux2, and then initialize them
     allocate( caux1(nfaux, maxval(rank)) ); caux1 = czero
     allocate( caux2(nfaux, maxval(rank)) ); caux2 = czero

! allocate memory for g2aux and h2aux, and then initialize them
     allocate( g2aux(nfaux, nfaux, norbs) ); g2aux = czero
     allocate( h2aux(nfaux, nfaux, norbs) ); h2aux = czero

! calculate prefactor: pref
     call ctqmc_make_pref()

! calculate g2aux and h2aux
! see Eq. (52) in Phys. Rev. B 89, 235128 (2014)
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE (flvr, is, ie, w2n, w1n, maux, naux, caux1, caux2)
     FLVR_CYCLE: do flvr=1,norbs
         call ctqmc_make_fexp(flvr, nfaux, maxval(rank), caux1, caux2)

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

     enddo FLVR_CYCLE ! over flvr={1,norbs} loop
!$OMP END DO

! calculate g2ph and h2ph
!
! note (for G2_PH_AABB component):
!
!     g2aux(w1n,w2n,f1) ->
!         exp [ i (\nu + \omega) \tau'_i ] exp [ -i \nu \tau_j ]
!
!     g2aux(w3n,w4n,f2) ->
!         exp [ i \nu' \tau'_k ] exp [ -i (\nu' + \omega) \tau_l ]
!
!     g2aux(w1n,w4n,f1) ->
!         exp [ i (\nu + \omega) \tau'_i ] exp [ -i (\nu' + \omega) \tau_l ]
!
!     g2aux(w3n,w2n,f1) ->
!         exp [ i \nu' \tau'_k ] exp [ -i \nu \tau_j ]
!
! note (for G2_PH_ABBA component):
!
!     g2aux(w1n,w2n,f1) ->
!         exp [ i (\nu + \omega) \tau'_i ] exp [ -i \nu \tau_j ]
!
!     g2aux(w3n,w4n,f1) ->
!         exp [ i \nu' \tau'_k ] exp [ -i (\nu' + \omega) \tau_l ]
!
!     g2aux(w1n,w4n,f1) ->
!         exp [ i (\nu + \omega) \tau'_i ] exp [ -i (\nu' + \omega) \tau_l ]
!
!     g2aux(w3n,w2n,f2) ->
!         exp [ i \nu' \tau'_k ] exp [ -i \nu \tau_j ]
!
!$OMP DO PRIVATE (f1, f2, wbn, w4n, w3n, w2n, w1n, zg, zh)
     ORB1_CYCLE: do f1=1,norbs                 ! block index: A
         ORB2_CYCLE: do f2=1,f1                ! block index: B
                                               !
             WB_CYCLE: do wbn=1,nbfrq          ! bosonic matsubara frequency: w
                                               !
                 WF1_CYCLE: do w2n=1,nffrq     ! fermionic matsubara frequency: v
                     WF2_CYCLE: do w3n=1,nffrq ! fermionic matsubara frequency: v'
                         w1n = w2n + wbn - 1 ! think it carefully
                         w4n = w3n + wbn - 1

                         zg = czero; zh = czero

! G2_PH_AABB component
!-------------------------------------------------------------------------
                     CALC_G2_PH_AABB: BLOCK

                         if ( btest(isvrt,1) ) then
                             zg = zg + g2aux(w1n,w2n,f1) * g2aux(w3n,w4n,f2)
                             zh = zh + h2aux(w1n,w2n,f1) * g2aux(w3n,w4n,f2)

                             if ( f1 == f2 ) then
                                 zg = zg - g2aux(w1n,w4n,f1) * g2aux(w3n,w2n,f1)
                                 zh = zh - h2aux(w1n,w4n,f1) * g2aux(w3n,w2n,f1)
                             endif ! back if ( f1 == f2 ) block
                         endif ! back if ( btest(isvrt,1) ) block

                     END BLOCK CALC_G2_PH_AABB

! G2_PH_ABBA component
!-------------------------------------------------------------------------
                     CALC_G2_PH_ABBA: BLOCK

                         if ( btest(isvrt,2) ) then
                             zg = zg - g2aux(w1n,w4n,f1) * g2aux(w3n,w2n,f2)
                             zh = zh - h2aux(w1n,w4n,f1) * g2aux(w3n,w2n,f2)

                             if ( f1 == f2 ) then
                                 zg = zg + g2aux(w1n,w2n,f1) * g2aux(w3n,w4n,f1)
                                 zh = zh + h2aux(w1n,w2n,f1) * g2aux(w3n,w4n,f1)
                             endif ! back if ( f1 == f2 ) block
                         endif ! back if ( btest(isvrt,2) ) block

                     END BLOCK CALC_G2_PH_ABBA

                         g2ph(w3n,w2n,wbn,f2,f1) = g2ph(w3n,w2n,wbn,f2,f1) + zg / beta
                         h2ph(w3n,w2n,wbn,f2,f1) = h2ph(w3n,w2n,wbn,f2,f1) + zh / beta
                     enddo WF2_CYCLE ! over w3n={1,nffrq} loop
                 enddo WF1_CYCLE ! over w2n={1,nffrq} loop

             enddo WB_CYCLE ! over wbn={1,nbfrq} loop

         enddo ORB2_CYCLE ! over f2={1,f1} loop
     enddo ORB1_CYCLE ! over f1={1,norbs} loop
!$OMP END DO
!$OMP END PARALLEL

! deallocate memory
     deallocate( caux1 )
     deallocate( caux2 )
     deallocate( g2aux )
     deallocate( h2aux )

     return
  end subroutine cat_record_g2ph_std

!!
!! @sub cat_record_g2ph_leg
!!
!! record the two-particle green's and vertex functions in the ph channel.
!! here improved estimator is used to improve the accuracy
!!
!! note:
!!
!!     we try to measure the two-particle green's and vertex functions in
!!     the particle-hole channel and legendre/matsubara representation
!!     in this subroutine. in order to simplify the calculations, we just
!!     consider the block structure of G^{(2)}
!!
!!     G^{(2)}_{abcd,AABB,ph} (l, l', \omega) =  (-1)^l'
!!         \frac{ \sqrt{2l - 1} \sqrt{2l' - 1} }{ \beta }
!!         \langle
!!             \sum^{K_A}_{ij} \sum^{K_B}_{kl}
!!             ( M^{A}_{ij} M^{B}_{kl} - \delta_{AB} M^{A}_{il} M^{B}_{kj} )
!!             p_l( x(\tau'_i - \tau_j) ) p_l'( x(\tau'_k - \tau_l) )
!!             exp [ i \omega (\tau'_i - \tau_l) ]
!!             \delta_{a,i} \delta_{b,j} \delta_{c,k} \delta_{d,l}
!!         \rangle
!!
!!     G^{(2)}_{abcd,ABBA,ph} (l, l', \omega) =  (-1)^l'
!!         \frac{ \sqrt{2l - 1} \sqrt{2l' - 1} }{ \beta }
!!         \langle
!!             \sum^{K_A}_{il} \sum^{K_B}_{kj}
!!             ( \delta_{AB} M^{A}_{ij} M^{B}_{kl} - M^{A}_{il} M^{B}_{kj} )
!!             p_l( x(\tau'_i - \tau_j) ) p_l'( x(\tau'_k - \tau_l) )
!!             exp [ i \omega (\tau'_i - \tau_l) ]
!!             \delta_{a,i} \delta_{b,j} \delta_{c,k} \delta_{d,l}
!!         \rangle
!!
!!     \tau'_i and \tau'_k: imaginary time for annihilation operators
!!     \tau_j and \tau_l: imaginary time for creation operators
!!     p_l and p_l': legendre polynomial
!!     \omega: bosonic matsubara frequency
!!
  subroutine cat_record_g2ph_leg()
     use constants, only : dp
     use constants, only : zero, one, two, czero

     use control, only : isvrt
     use control, only : norbs
     use control, only : lemax, legrd
     use control, only : nbfrq
     use control, only : beta

     use context, only : index_s, index_e
     use context, only : time_s, time_e
     use context, only : rep_l
     use context, only : g2ph, h2ph
     use context, only : rank, pref
     use context, only : mmat

     implicit none

! local variables
! loop index for flavor channel
     integer  :: f1
     integer  :: f2

! loop index for frequency
     integer  :: wbn
     integer  :: l1
     integer  :: l2

! loop indices for start and end points
     integer  :: is1
     integer  :: is2
     integer  :: ie1
     integer  :: ie2

! index for imaginary time \tau
     integer  :: curr

! interval for imaginary time slice
     real(dp) :: step

! distance betweem \tau_s and \tau_e
     real(dp) :: dt

! sign for p_l(x(\tau))
     real(dp) :: ms

! real(dp) dummy variables
     real(dp) :: mm
     real(dp) :: pp

! complex(dp) dummy variables
     complex(dp) :: ee

! sqrt(2l+1) sqrt(2l'+1) (-1)^{(l'+1)}
     real(dp), allocatable :: lfun(:,:)

! p_l(x(\tau_e - \tau_s))
     real(dp), allocatable :: pfun(:,:,:,:,:)

! exp [i \omega_n \tau_s] and exp [i \omega_n \tau_e]
! note here \omega_n is bosonic
     complex(dp), allocatable :: caux1(:,:,:)
     complex(dp), allocatable :: caux2(:,:,:)

! allocate memory
     allocate( lfun(lemax,lemax) ); lfun = zero
     allocate( pfun(lemax, maxval(rank), maxval(rank), norbs, norbs)); pfun = zero

     allocate( caux1(nbfrq, maxval(rank), norbs) ); caux1 = czero
     allocate( caux2(nbfrq, maxval(rank), norbs) ); caux2 = czero

! calculate prefactor: pref
     call ctqmc_make_pref()

! prepare some important arrays: lfun
     do l1=1,lemax     ! legendre polynomial index: l
         do l2=1,lemax ! legendre polynomial index: l'
             lfun(l1,l2) = sqrt(two * l1 - one) * sqrt(two * l2 - one) * ( (-one)**l2 )
         enddo ! over l2={1,lemax} loop
     enddo ! over l1={1,lemax} loop

! prepare some important arrays: pfun
     step = real(legrd - 1) / two
     do f1=1,norbs
         do is1=1,rank(f1)
             do f2=1,norbs
                 do ie2=1,rank(f2)
! determine dt (distance) and ms (sign)
                     dt = time_e( index_e(ie2, f2), f2 ) - time_s( index_s(is1, f1), f1 )
                     ms = sign(one, dt)

! adjust dt, keep it stay in (zero, beta)
                     if ( dt < zero ) then
                         dt = dt + beta
                     endif ! back if ( dt < zero ) block

! determine index for imaginary time
                     curr = nint( ( two * dt / beta ) * step ) + 1

! special tricks for the first point and the last point
                     if ( curr == 1 .or. curr == legrd ) then
                         ms = two * ms
                     endif ! back if ( curr == 1 .or. curr == legrd ) block

! fill pfun
                     do l1=1,lemax
                         pfun(l1,ie2,is1,f2,f1) = ms * rep_l(curr,l1)
                     enddo ! over l1={1,lemax} loop
                 enddo ! over ie2={1,rank(f2)} loop
             enddo ! over f2={1,norbs} loop
         enddo ! over is1={1,rank(f1)} loop
     enddo ! over f1={1,norbs} loop

! prepare some important arrays: caux1 and caux2
     do f1=1,norbs
         call ctqmc_make_bexp(f1, nbfrq, maxval(rank), caux1(:,:,f1), caux2(:,:,f1))
     enddo ! over f1={1,norbs} loop

! calculate g2ph and h2ph
!
! G2_PH_AABB component
!-------------------------------------------------------------------------
     CALC_G2_PH_AABB: BLOCK

         if ( btest(isvrt,1) ) then

             do f1=1,norbs                         ! block index: A
                 do f2=1,f1                        ! block index: B
                     do is1=1,rank(f1)             ! \beta  -> j: creation operator
                         do ie1=1,rank(f1)         ! \alpha -> i: annihilation operator
                             do is2=1,rank(f2)     ! \delta -> l: creation operator
                                 do ie2=1,rank(f2) ! \gamma -> k: annihilation operator
             !-------------------!
             do wbn=1,nbfrq                        ! bosonic matsubara frequency: w
                 do l1=1,lemax                     ! legendre polynomial index: l
                     do l2=1,lemax                 ! legendre polynomial index: l'
                         ee = caux2(wbn,ie1,f1) * caux1(wbn,is2,f2)
                         pp = pfun(l1,ie1,is1,f1,f1) * pfun(l2,ie2,is2,f2,f2) * lfun(l1,l2)
                         mm = mmat(ie1, is1, f1) * mmat(ie2, is2, f2)
                         if ( f1 == f2 ) then
                             mm = mm - mmat(ie1, is2, f1) * mmat(ie2, is1, f1)
                         endif ! back if ( f1 == f2 ) block

                         g2ph(l2,l1,wbn,f2,f1) = g2ph(l2,l1,wbn,f2,f1) + mm * pp * ee / beta
                         h2ph(l2,l1,wbn,f2,f1) = h2ph(l2,l1,wbn,f2,f1) + mm * pp * ee / beta * pref(ie1,f1)
                     enddo ! over l2={1,lemax} loop
                 enddo ! over l1={1,lemax} loop
             enddo ! over wbn={1,nbfrq} loop
             !-------------------!
                                 enddo ! over ie2={1,rank(f2)} loop
                             enddo ! over is2={1,rank(f2)} loop
                         enddo ! over ie1={1,rank(f1)} loop
                     enddo ! is1={1,rank(f1)} loop
                 enddo ! over f2={1,f1} loop
             enddo ! over f1={1,norbs} loop

         endif ! back if ( btest(isvrt,1) ) block

     END BLOCK CALC_G2_PH_AABB

! G2_PH_ABBA component
!-------------------------------------------------------------------------
     CALC_G2_PH_ABBA: BLOCK

         if ( btest(isvrt,2) ) then

             do f1=1,norbs                         ! block index: A
                 do f2=1,f1                        ! block index: B
                     do is1=1,rank(f1)             ! \delta -> l: creation operator
                         do ie1=1,rank(f1)         ! \alpha -> i: annihilation operator
                             do is2=1,rank(f2)     ! \beta  -> j: creation operator
                                 do ie2=1,rank(f2) ! \gamma -> k: annihilation operator
             !-------------------!
             do wbn=1,nbfrq                        ! bosonic matsubara frequency: w
                 do l1=1,lemax                     ! legendre polynomial index: l
                     do l2=1,lemax                 ! legendre polynomial index: l'
                         ee = caux2(wbn,ie1,f1) * caux1(wbn,is1,f1)
                         pp = pfun(l1,ie1,is2,f1,f2) * pfun(l2,ie2,is1,f2,f1) * lfun(l1,l2)
                         mm = mmat(ie1, is1, f1) * mmat(ie2, is2, f2)
                         if ( f1 == f2 ) then
                             mm = mm - mmat(ie1, is2, f1) * mmat(ie2, is1, f1)
                         endif ! back if ( f1 == f2 ) block

                         g2ph(l2,l1,wbn,f2,f1) = g2ph(l2,l1,wbn,f2,f1) - mm * pp * ee / beta
                         h2ph(l2,l1,wbn,f2,f1) = h2ph(l2,l1,wbn,f2,f1) - mm * pp * ee / beta * pref(ie1,f1)
                     enddo ! over l2={1,lemax} loop
                 enddo ! over l1={1,lemax} loop
             enddo ! over wbn={1,nbfrq} loop
             !-------------------!
                                 enddo ! over ie2={1,rank(f2)} loop
                             enddo ! over is2={1,rank(f2)} loop
                         enddo ! over ie1={1,rank(f1)} loop
                     enddo ! is1={1,rank(f1)} loop
                 enddo ! over f2={1,f1} loop
             enddo ! over f1={1,norbs} loop

         endif ! back if ( btest(isvrt,2) ) block

     END BLOCK CALC_G2_PH_ABBA

! deallocate memory
     deallocate( lfun  )
     deallocate( pfun  )
     deallocate( caux1 )
     deallocate( caux2 )

     return
  end subroutine cat_record_g2ph_leg

!!
!! @sub cat_record_g2ph_svd
!!
!! record the two-particle green's and vertex functions in the ph channel.
!! here improved estimator is used to improve the accuracy
!!
!! note:
!!
!!     we try to measure the two-particle green's and vertex functions in
!!     the particle-hole channel and intermediate/matsubara representation
!!     in this subroutine. in order to simplify the calculations, we just
!!     consider the block structure of G^{(2)}
!!
!!     G^{(2)}_{abcd,AABB,ph} (l, l', \omega) =  (-1)^l'
!!         \frac{ 1 }{ \beta }
!!         \langle
!!             \sum^{K_A}_{ij} \sum^{K_B}_{kl}
!!             ( M^{A}_{ij} M^{B}_{kl} - \delta_{AB} M^{A}_{il} M^{B}_{kj} )
!!             u_l( x(\tau'_i - \tau_j) ) u_l'( x(\tau'_k - \tau_l) )
!!             exp [ i \omega (\tau'_i - \tau_l) ]
!!             \delta_{a,i} \delta_{b,j} \delta_{c,k} \delta_{d,l}
!!         \rangle
!!
!!     G^{(2)}_{abcd,ABBA,ph} (l, l', \omega) =  (-1)^l'
!!         \frac{ \sqrt{2l - 1} \sqrt{2l' - 1} }{ \beta }
!!         \langle
!!             \sum^{K_A}_{il} \sum^{K_B}_{kj}
!!             ( \delta_{AB} M^{A}_{ij} M^{B}_{kl} - M^{A}_{il} M^{B}_{kj} )
!!             u_l( x(\tau'_i - \tau_j) ) u_l'( x(\tau'_k - \tau_l) )
!!             exp [ i \omega (\tau'_i - \tau_l) ]
!!             \delta_{a,i} \delta_{b,j} \delta_{c,k} \delta_{d,l}
!!         \rangle
!!
!!     \tau'_i and \tau'_k: imaginary time for annihilation operators
!!     \tau_j and \tau_l: imaginary time for creation operators
!!     u_l and u_l': svd polynomial
!!     \omega: bosonic matsubara frequency
!!
  subroutine cat_record_g2ph_svd()
     use constants, only : dp
     use constants, only : zero, one, two, czero

     use control, only : isvrt
     use control, only : norbs
     use control, only : svmax, svgrd
     use control, only : nbfrq
     use control, only : beta

     use context, only : index_s, index_e
     use context, only : time_s, time_e
     use context, only : rep_s
     use context, only : g2ph, h2ph
     use context, only : rank, pref
     use context, only : mmat

     implicit none

! local variables
! loop index for flavor channel
     integer  :: f1
     integer  :: f2

! loop index for frequency
     integer  :: wbn
     integer  :: l1
     integer  :: l2

! loop indices for start and end points
     integer  :: is1
     integer  :: is2
     integer  :: ie1
     integer  :: ie2

! index for imaginary time \tau
     integer  :: curr

! interval for imaginary time slice
     real(dp) :: step

! distance betweem \tau_s and \tau_e
     real(dp) :: dt

! sign for u_l(x(\tau))
     real(dp) :: ms

! real(dp) dummy variables
     real(dp) :: mm
     real(dp) :: uu

! complex(dp) dummy variables
     complex(dp) :: ee

! u_l(x(\tau_e - \tau_s))
     real(dp), allocatable :: ufun(:,:,:,:,:)

! exp [i \omega_n \tau_s] and exp [i \omega_n \tau_e]
! note here \omega_n is bosonic
     complex(dp), allocatable :: caux1(:,:,:)
     complex(dp), allocatable :: caux2(:,:,:)

! allocate memory
     allocate( ufun(svmax, maxval(rank), maxval(rank), norbs, norbs)); ufun = zero

     allocate( caux1(nbfrq, maxval(rank), norbs) ); caux1 = czero
     allocate( caux2(nbfrq, maxval(rank), norbs) ); caux2 = czero

! calculate prefactor: pref
     call ctqmc_make_pref()

! prepare some important arrays: ufun
     step = real(svgrd - 1) / two
     do f1=1,norbs
         do is1=1,rank(f1)
             do f2=1,norbs
                 do ie2=1,rank(f2)
! determine dt (distance) and ms (sign)
                     dt = time_e( index_e(ie2, f2), f2 ) - time_s( index_s(is1, f1), f1 )
                     ms = sign(one, dt)

! adjust dt, keep it stay in (zero, beta)
                     if ( dt < zero ) then
                         dt = dt + beta
                     endif ! back if ( dt < zero ) block

! determine index for imaginary time
                     call s_svd_point(two * dt / beta - one, step, curr)

! fill ufun
                     do l1=1,svmax
                         ufun(l1,ie2,is1,f2,f1) = ms * rep_s(curr,l1)
                     enddo ! over l1={1,svmax} loop
                 enddo ! over ie2={1,rank(f2)} loop
             enddo ! over f2={1,norbs} loop
         enddo ! over is1={1,rank(f1)} loop
     enddo ! over f1={1,norbs} loop

! prepare some important arrays: caux1 and caux2
     do f1=1,norbs
         call ctqmc_make_bexp(f1, nbfrq, maxval(rank), caux1(:,:,f1), caux2(:,:,f1))
     enddo ! over f1={1,norbs} loop

! calculate g2ph and h2ph
!
! G2_PH_AABB component
!-------------------------------------------------------------------------
     CALC_G2_PH_AABB: BLOCK

         if ( btest(isvrt,1) ) then

             do f1=1,norbs                         ! block index: A
                 do f2=1,f1                        ! block index: B
                     do is1=1,rank(f1)             ! \beta  -> j: creation operator
                         do ie1=1,rank(f1)         ! \alpha -> i: annihilation operator
                             do is2=1,rank(f2)     ! \delta -> l: creation operator
                                 do ie2=1,rank(f2) ! \gamma -> k: annihilation operator
             !-------------------!
             do wbn=1,nbfrq                        ! bosonic matsubara frequency: w
                 do l1=1,svmax                     ! svd polynomial index: l
                     do l2=1,svmax                 ! svd polynomial index: l'
                         ee = caux2(wbn,ie1,f1) * caux1(wbn,is2,f2)
                         uu = ufun(l1,ie1,is1,f1,f1) * ufun(l2,ie2,is2,f2,f2)
                         mm = mmat(ie1, is1, f1) * mmat(ie2, is2, f2)
                         if ( f1 == f2 ) then
                             mm = mm - mmat(ie1, is2, f1) * mmat(ie2, is1, f1)
                         endif ! back if ( f1 == f2 ) block

                         g2ph(l2,l1,wbn,f2,f1) = g2ph(l2,l1,wbn,f2,f1) + mm * uu * ee / beta
                         h2ph(l2,l1,wbn,f2,f1) = h2ph(l2,l1,wbn,f2,f1) + mm * uu * ee / beta * pref(ie1,f1)
                     enddo ! over l2={1,svmax} loop
                 enddo ! over l1={1,svmax} loop
             enddo ! over wbn={1,nbfrq} loop
             !-------------------!
                                 enddo ! over ie2={1,rank(f2)} loop
                             enddo ! over is2={1,rank(f2)} loop
                         enddo ! over ie1={1,rank(f1)} loop
                     enddo ! is1={1,rank(f1)} loop
                 enddo ! over f2={1,f1} loop
             enddo ! over f1={1,norbs} loop

         endif ! back if ( btest(isvrt,1) ) block

     END BLOCK CALC_G2_PH_AABB

! G2_PH_ABBA component
!-------------------------------------------------------------------------

     CALC_G2_PH_ABBA: BLOCK

         if ( btest(isvrt,2) ) then

             do f1=1,norbs                         ! block index: A
                 do f2=1,f1                        ! block index: B
                     do is1=1,rank(f1)             ! \delta -> l: creation operator
                         do ie1=1,rank(f1)         ! \alpha -> i: annihilation operator
                             do is2=1,rank(f2)     ! \beta  -> j: creation operator
                                 do ie2=1,rank(f2) ! \gamma -> k: annihilation operator
             !-------------------!
             do wbn=1,nbfrq                        ! bosonic matsubara frequency: w
                 do l1=1,svmax                     ! svd polynomial index: l
                     do l2=1,svmax                 ! svd polynomial index: l'
                         ee = caux2(wbn,ie1,f1) * caux1(wbn,is1,f1)
                         uu = ufun(l1,ie1,is2,f1,f2) * ufun(l2,ie2,is1,f2,f1)
                         mm = mmat(ie1, is1, f1) * mmat(ie2, is2, f2)
                         if ( f1 == f2 ) then
                             mm = mm - mmat(ie1, is2, f1) * mmat(ie2, is1, f1)
                         endif ! back if ( f1 == f2 ) block

                         g2ph(l2,l1,wbn,f2,f1) = g2ph(l2,l1,wbn,f2,f1) - mm * uu * ee / beta
                         h2ph(l2,l1,wbn,f2,f1) = h2ph(l2,l1,wbn,f2,f1) - mm * uu * ee / beta * pref(ie1,f1)
                     enddo ! over l2={1,svmax} loop
                 enddo ! over l1={1,svmax} loop
             enddo ! over wbn={1,nbfrq} loop
             !-------------------!
                                 enddo ! over ie2={1,rank(f2)} loop
                             enddo ! over is2={1,rank(f2)} loop
                         enddo ! over ie1={1,rank(f1)} loop
                     enddo ! is1={1,rank(f1)} loop
                 enddo ! over f2={1,f1} loop
             enddo ! over f1={1,norbs} loop

         endif ! back if ( btest(isvrt,2) ) block

     END BLOCK CALC_G2_PH_ABBA

! deallocate memory
     deallocate( ufun  )
     deallocate( caux1 )
     deallocate( caux2 )

     return
  end subroutine cat_record_g2ph_svd

!!
!! @sub ctqmc_record_g2pp
!!
!! record the two-particle green's and vertex functions in the pp channel.
!! here improved estimator is used to improve the accuracy
!!
!! notation:
!!
!!     G^{(2)}_{\alpha\beta\gamma\delta} (\tau_1, \tau_2, \tau_3, \tau_4)
!!         = \langle T_\tau
!!               c_{\alpha} (\tau_1) c^{\dagger}_{\beta} (\tau_2)
!!               c_{\gamma} (\tau_3) c^{\dagger}_{\delta} (\tau_4)
!!           \rangle
!!
!!     G^{(2)}_{\alpha\beta\gamma\delta,pp} (\nu, \nu', \omega)
!!         = \langle
!!               c_{\alpha} (\omega - \nu') c^{*}_{\beta} (\nu)
!!               c_{\gamma} (\nu') c^{*}_{\delta} (\omega - nu)
!!           \rangle
!!
!!     \nu and \nu': fermionic matsubara frequency
!!     \omega: bosonic matsubara frequency
!!
!!        in             out
!!         \              /
!!     w-v' \            / w-v
!!           \          /
!!          i \--------/ l
!!            |        |
!!            |        |
!!            |        |
!!          j /--------\ k
!!           /          \
!!        v /            \ v'
!!         /              \
!!        out             in
!!
  subroutine ctqmc_record_g2pp()
     use control, only : isort
     use control, only : isvrt

     implicit none

! check whether there is conflict
! this subroutine is only designed for the particle-particle channel
     call s_assert2( btest(isvrt, 3) .or. btest(isvrt, 4), 'in ctqmc_record_g2pp' )

! you can not calculate the AABB and ABBA components at the same time
     call s_assert2( .not. ( btest(isvrt, 3) .and. btest(isvrt, 4) ), 'in ctqmc_record_g2pp' )

     select case ( isort )

         case (1)
             call cat_record_g2pp_std()

         case (2)
             call cat_record_g2pp_leg()

         case (3)
             call cat_record_g2pp_svd()

         case default
             call s_print_error('ctqmc_record_g2pp','this feature is not implemented')

     end select

     return
  end subroutine ctqmc_record_g2pp

!!
!! @sub cat_record_g2pp_std
!!
!! record the two-particle green's and vertex functions in the pp channel.
!! here improved estimator is used to improve the accuracy
!!
!! note:
!!
!!     we try to measure the two-particle green's and vertex functions in
!!     the particle-particle channel and matsubara frequency representation
!!     in this subroutine. in order to simplify the calculations, we just
!!     consider the block structure of G^{(2)}
!!
!!     G^{(2)}_{abcd,AABB,pp} (\nu, \nu', \omega) = \frac{1}{\beta}
!!         \langle
!!             \sum^{K_A}_{ij=1} \sum^{K_B}_{kl=1}
!!             ( M^{A}_{ij} M^{B}_{kl} - \delta_{AB} M^{A}_{il} M^{B}_{kj} )
!!             exp [ i (\omega - \nu') \tau'_i ]
!!             exp [ -i \nu \tau_j ]
!!             exp [ i \nu' \tau'_k ]
!!             exp [ -i (\omega - \nu) \tau_l ]
!!             \delta_{a,i} \delta_{b,j} \delta_{c,k} \delta_{d,l}
!!         \rangle
!!
!!     G^{(2)}_{abcd,ABBA,pp} (\nu, \nu', \omega) = \frac{1}{\beta}
!!         \langle
!!             \sum^{K_A}_{ij=1} \sum^{K_B}_{kl=1}
!!             ( \delta_{AB} M^{A}_{ij} M^{B}_{kl} - M^{A}_{il} M^{B}_{kj} )
!!             exp [ i (\omega - \nu') \tau'_i ]
!!             exp [ -i \nu \tau_j ]
!!             exp [ i \nu' \tau'_k ]
!!             exp [ -i (\omega - \nu) \tau_l ]
!!             \delta_{a,i} \delta_{b,j} \delta_{c,k} \delta_{d,l}
!!         \rangle
!!
!!     \tau'_i and \tau'_k: imaginary time for annihilation operator
!!     \tau_j and \tau_l: imaginary time for creation operator
!!     \nu and \nu': fermionic matsubara frequency
!!     \omega: bosonic matsubara frequency
!!
  subroutine cat_record_g2pp_std()
     use constants, only : dp
     use constants, only : czero

     use control, only : isvrt
     use control, only : norbs
     use control, only : nffrq, nbfrq
     use control, only : beta

     use context, only : g2pp
     use context, only : h2pp
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

! dummy complex(dp) variables, used to calculate the g2pp and h2pp
     complex(dp) :: zg
     complex(dp) :: zh

! exp [i \omega_n \tau_s] and exp [i \omega_n \tau_e]
     complex(dp), allocatable :: caux1(:,:)
     complex(dp), allocatable :: caux2(:,:)

! \sum_{ij=1} exp [i \omega_m \tau'_i ] M_{ij} exp [ i \omega_n \tau_j ]
! where m and n are the first two frequency indices for g2aux and h2aux
     complex(dp), allocatable :: g2aux(:,:,:)
     complex(dp), allocatable :: h2aux(:,:,:)

! evaluate nfaux, determine the size of g2aux and h2aux
     nfaux = nffrq + nbfrq - 1

! allocate memory for caux1 and caux2, and then initialize them
     allocate( caux1(nfaux, maxval(rank)) ); caux1 = czero
     allocate( caux2(nfaux, maxval(rank)) ); caux2 = czero

! allocate memory for g2aux and h2aux, and then initialize them
     allocate( g2aux(nfaux, nfaux, norbs) ); g2aux = czero
     allocate( h2aux(nfaux, nfaux, norbs) ); h2aux = czero

! calculate prefactor: pref
     call ctqmc_make_pref()

! calculate g2aux and h2aux
! see Eq. (52) in Phys. Rev. B 89, 235128 (2014)
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO PRIVATE (flvr, is, ie, w2n, w1n, maux, naux, caux1, caux2)
     FLVR_CYCLE: do flvr=1,norbs
         call ctqmc_make_fexp(flvr, nfaux, maxval(rank), caux1, caux2)

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

     enddo FLVR_CYCLE ! over flvr={1,norbs} loop
!$OMP END DO

! calculate g2pp and h2pp
!
! note (for G2_PP_AABB component):
!
!     g2aux(w1n,w2n,f1) ->
!         exp [ i (\omega - \nu') \tau'_i ] exp [ -i \nu \tau_j ]
!
!     g2aux(w3n,w4n,f2) ->
!         exp [ i \nu' \tau'_k ] exp [ -i (\omega - \nu) \tau_l ]
!
!     g2aux(w1n,w4n,f1) ->
!         exp [ i (\omega - \nu') \tau'_i ] exp [ -i (\omega - nu) \tau_l ]
!
!     g2aux(w3n,w2n,f1) ->
!         exp [ i \nu' \tau'_k ] exp [ -i \nu \tau_j ]
!
! note (for G2_PP_ABBA component):
!
!     g2aux(w1n,w2n,f1) ->
!         exp [ i (\omega - \nu') \tau'_i ] exp [ -i \nu \tau_j ]
!
!     g2aux(w3n,w4n,f1) ->
!         exp [ i \nu' \tau'_k ] exp [ -i (\omega - \nu) \tau_l ]
!
!     g2aux(w1n,w4n,f1) ->
!         exp [ i (\omega - \nu') \tau'_i ] exp [ -i (\omega - nu) \tau_l ]
!
!     g2aux(w3n,w2n,f2) ->
!         exp [ i \nu' \tau'_k ] exp [ -i \nu \tau_j ]
!
!$OMP DO PRIVATE (f1, f2, wbn, w4n, w3n, w2n, w1n, zg, zh)
     ORB1_CYCLE: do f1=1,norbs                 ! block index: A
         ORB2_CYCLE: do f2=1,f1                ! block index: B
                                               !
             WB_CYCLE: do wbn=1,nbfrq          ! bosonic matsubara frequency: w
                                               !
                 WF1_CYCLE: do w2n=1,nffrq     ! fermionic matsubara frequency: v
                     WF2_CYCLE: do w3n=1,nffrq ! fermionic matsubara frequency: v'
                         w1n = wbn - w3n + nffrq ! think it carefully
                         w4n = wbn - w2n + nffrq

                         zg = czero; zh = czero

! G2_PP_AABB component
!-------------------------------------------------------------------------
                     CALC_G2_PP_AABB: BLOCK

                         if ( btest(isvrt,3) ) then
                             zg = zg + g2aux(w1n,w2n,f1) * g2aux(w3n,w4n,f2)
                             zh = zh + h2aux(w1n,w2n,f1) * g2aux(w3n,w4n,f2)

                             if ( f1 == f2 ) then
                                 zg = zg - g2aux(w1n,w4n,f1) * g2aux(w3n,w2n,f1)
                                 zh = zh - h2aux(w1n,w4n,f1) * g2aux(w3n,w2n,f1)
                             endif ! back if ( f1 == f2 ) block
                         endif ! back if ( btest(isvrt,3) ) block

                     END BLOCK CALC_G2_PP_AABB

! G2_PP_ABBA component
!-------------------------------------------------------------------------
                     CALC_G2_PP_ABBA: BLOCK

                         if ( btest(isvrt,4) ) then
                             zg = zg - g2aux(w1n,w4n,f1) * g2aux(w3n,w2n,f2)
                             zh = zh - h2aux(w1n,w4n,f1) * g2aux(w3n,w2n,f2)

                             if ( f1 == f2 ) then
                                 zg = zg + g2aux(w1n,w2n,f1) * g2aux(w3n,w4n,f1)
                                 zh = zh + h2aux(w1n,w2n,f1) * g2aux(w3n,w4n,f1)
                             endif ! back if ( f1 == f2 ) block
                         endif ! back if ( btest(isvrt,4) ) block

                     END BLOCK CALC_G2_PP_ABBA

                         g2pp(w3n,w2n,wbn,f2,f1) = g2pp(w3n,w2n,wbn,f2,f1) + zg / beta
                         h2pp(w3n,w2n,wbn,f2,f1) = h2pp(w3n,w2n,wbn,f2,f1) + zh / beta
                     enddo WF2_CYCLE ! over w3n={1,nffrq} loop
                 enddo WF1_CYCLE ! over w2n={1,nffrq} loop

             enddo WB_CYCLE ! over wbn={1,nbfrq} loop

         enddo ORB2_CYCLE ! over f2={1,f1} loop
     enddo ORB1_CYCLE ! over f1={1,norbs} loop
!$OMP END DO
!$OMP END PARALLEL

! deallocate memory
     deallocate( caux1 )
     deallocate( caux2 )
     deallocate( g2aux )
     deallocate( h2aux )

     return
  end subroutine cat_record_g2pp_std

!!
!! @sub cat_record_g2pp_leg
!!
!! record the two-particle green's and vertex functions in the pp channel.
!! here improved estimator is used to improve the accuracy
!!
!! note:
!!
!!     we try to measure the two-particle green's and vertex functions in
!!     the particle-particle channel and legendre/matsubara representation
!!     in this subroutine. in order to simplify the calculations, we just
!!     consider the block structure of G^{(2)}
!!
!!     G^{(2)}_{abcd,AABB,pp} (l, l', \omega) =  (-1)^l'
!!         \frac{ \sqrt{2l - 1} \sqrt{2l' - 1} }{ \beta }
!!         \langle
!!             \sum^{K_A}_{ij} \sum^{K_B}_{kl}
!!             ( M^{A}_{ij} M^{B}_{kl} - \delta_{AB} M^{A}_{il} M^{B}_{kj} )
!!             p_l( x(\tau_l - \tau_j) ) p_l'( x(\tau'_k - \tau'_i) )
!!             exp [ i \omega (\tau'_i - \tau_l) ]
!!             \delta_{a,i} \delta_{b,j} \delta_{c,k} \delta_{d,l}
!!         \rangle
!!
!!     G^{(2)}_{abcd,ABBA,pp} (l, l', \omega) =  (-1)^l'
!!         \frac{ \sqrt{2l - 1} \sqrt{2l' - 1} }{ \beta }
!!         \langle
!!             \sum^{K_A}_{il} \sum^{K_B}_{kj}
!!             ( \delta_{AB} M^{A}_{ij} M^{B}_{kl} - M^{A}_{il} M^{B}_{kj} )
!!             p_l( x(\tau_l - \tau_j) ) p_l'( x(\tau'_k - \tau'_i) )
!!             exp [ i \omega (\tau'_i - \tau_l) ]
!!             \delta_{a,i} \delta_{b,j} \delta_{c,k} \delta_{d,l}
!!         \rangle
!!
!!     \tau'_i and \tau'_k: imaginary time for annihilation operators
!!     \tau_j and \tau_l: imaginary time for creation operators
!!     p_l and p_l': legendre polynomial
!!     \omega: bosonic matsubara frequency
!!
  subroutine cat_record_g2pp_leg()
     use constants, only : dp
     use constants, only : zero, one, two, czero

     use control, only : isvrt
     use control, only : norbs
     use control, only : lemax, legrd
     use control, only : nbfrq
     use control, only : beta

     use context, only : index_s, index_e
     use context, only : time_s, time_e
     use context, only : rep_l
     use context, only : g2pp, h2pp
     use context, only : rank, pref
     use context, only : mmat

     implicit none

! local variables
! loop index for flavor channel
     integer  :: f1
     integer  :: f2

! loop index for frequency
     integer  :: wbn
     integer  :: l1
     integer  :: l2

! loop indices for start and end points
     integer  :: is1
     integer  :: is2
     integer  :: ie1
     integer  :: ie2

! index for imaginary time \tau
     integer  :: curr

! interval for imaginary time slice
     real(dp) :: step

! distance betweem \tau_s and \tau_e
     real(dp) :: dt

! sign for p_l(x(\tau))
     real(dp) :: ms

! real(dp) dummy variables
     real(dp) :: mm
     real(dp) :: pp

! complex(dp) dummy variables
     complex(dp) :: ee

! sqrt(2l+1) sqrt(2l'+1) (-1)^{(l'+1)}
     real(dp), allocatable :: lfun(:,:)

! p_l(x(\tau_s2 - \tau_s1))
     real(dp), allocatable :: pl_s(:,:,:,:,:)

! p_l(x(\tau_e2 - \tau_e1))
     real(dp), allocatable :: pl_e(:,:,:,:,:)

! exp [i \omega_n \tau_s] and exp [i \omega_n \tau_e]
! note here \omega_n is bosonic
     complex(dp), allocatable :: caux1(:,:,:)
     complex(dp), allocatable :: caux2(:,:,:)

! allocate memory
     allocate( lfun(lemax,lemax) ); lfun = zero
     allocate( pl_s(lemax, maxval(rank), maxval(rank), norbs, norbs)); pl_s = zero
     allocate( pl_e(lemax, maxval(rank), maxval(rank), norbs, norbs)); pl_e = zero

     allocate( caux1(nbfrq, maxval(rank), norbs) ); caux1 = czero
     allocate( caux2(nbfrq, maxval(rank), norbs) ); caux2 = czero

! calculate prefactor: pref
     call ctqmc_make_pref()

! prepare some important arrays: lfun
     do l1=1,lemax     ! legendre polynomial index: l
         do l2=1,lemax ! legendre polynomial index: l'
             lfun(l1,l2) = sqrt(two * l1 - one) * sqrt(two * l2 - one) * ( (-one)**l2 )
         enddo ! over l2={1,lemax} loop
     enddo ! over l1={1,lemax} loop

! prepare some important arrays: pl_s
     step = real(legrd - 1) / two
     do f1=1,norbs
         do is1=1,rank(f1)
             do f2=1,norbs
                 do is2=1,rank(f2)
! determine dt (distance) and ms (sign)
                     dt = time_s( index_s(is2, f2), f2 ) - time_s( index_s(is1, f1), f1 )
                     ms = sign(one, dt)

! adjust dt, keep it stay in (zero, beta)
                     if ( dt < zero ) then
                         dt = dt + beta
                     endif ! back if ( dt < zero ) block

! determine index for imaginary time
                     curr = nint( ( two * dt / beta ) * step ) + 1

! special tricks for the first point and the last point
                     if ( curr == 1 .or. curr == legrd ) then
                         ms = two * ms
                     endif ! back if ( curr == 1 .or. curr == legrd ) block

! fill pl_s
                     do l1=1,lemax
                         pl_s(l1,is2,is1,f2,f1) = ms * rep_l(curr,l1)
                     enddo ! over l1={1,lemax} loop
                 enddo ! over is2={1,rank(f2)} loop
             enddo ! over f2={1,norbs} loop
         enddo ! over is1={1,rank(f1)} loop
     enddo ! over f1={1,norbs} loop

! prepare some important arrays: pl_e
     step = real(legrd - 1) / two
     do f1=1,norbs
         do ie1=1,rank(f1)
             do f2=1,norbs
                 do ie2=1,rank(f2)
! determine dt (distance) and ms (sign)
                     dt = time_e( index_e(ie2, f2), f2 ) - time_e( index_e(ie1, f1), f1 )
                     ms = sign(one, dt)

! adjust dt, keep it stay in (zero, beta)
                     if ( dt < zero ) then
                         dt = dt + beta
                     endif ! back if ( dt < zero ) block

! determine index for imaginary time
                     curr = nint( ( two * dt / beta ) * step ) + 1

! special tricks for the first point and the last point
                     if ( curr == 1 .or. curr == legrd ) then
                         ms = two * ms
                     endif ! back if ( curr == 1 .or. curr == legrd ) block

! fill pl_e
                     do l1=1,lemax
                         pl_e(l1,ie2,ie1,f2,f1) = ms * rep_l(curr,l1)
                     enddo ! over l1={1,lemax} loop
                 enddo ! over ie2={1,rank(f2)} loop
             enddo ! over f2={1,norbs} loop
         enddo ! over ie1={1,rank(f1)} loop
     enddo ! over f1={1,norbs} loop

! prepare some important arrays: caux1 and caux2
     do f1=1,norbs
         call ctqmc_make_bexp(f1, nbfrq, maxval(rank), caux1(:,:,f1), caux2(:,:,f1))
     enddo ! over f1={1,norbs} loop

! calculate g2pp and h2pp
!
! G2_PP_AABB component
!-------------------------------------------------------------------------
     CALC_G2_PP_AABB: BLOCK

         if ( btest(isvrt,3) ) then

             do f1=1,norbs                         ! block index: A
                 do f2=1,f1                        ! block index: B
                     do is1=1,rank(f1)             ! \beta  -> j: creation operator
                         do ie1=1,rank(f1)         ! \alpha -> i: annihilation operator
                             do is2=1,rank(f2)     ! \delta -> l: creation operator
                                 do ie2=1,rank(f2) ! \gamma -> k: annihilation operator
             !-------------------!
             do wbn=1,nbfrq                        ! bosonic matsubara frequency: w
                 do l1=1,lemax                     ! legendre polynomial index: l
                     do l2=1,lemax                 ! legendre polynomial index: l'
                         ee = caux2(wbn,ie1,f1) * caux1(wbn,is2,f2)
                         pp = pl_s(l1,is2,is1,f2,f1) * pl_e(l2,ie2,ie1,f2,f1) * lfun(l1,l2)
                         mm = mmat(ie1, is1, f1) * mmat(ie2, is2, f2)
                         if ( f1 == f2 ) then
                             mm = mm - mmat(ie1, is2, f1) * mmat(ie2, is1, f1)
                         endif ! back if ( f1 == f2 ) block

                         g2pp(l2,l1,wbn,f2,f1) = g2pp(l2,l1,wbn,f2,f1) + mm * pp * ee / beta
                         h2pp(l2,l1,wbn,f2,f1) = h2pp(l2,l1,wbn,f2,f1) + mm * pp * ee / beta * pref(ie1,f1)
                     enddo ! over l2={1,lemax} loop
                 enddo ! over l1={1,lemax} loop
             enddo ! over wbn={1,nbfrq} loop
             !-------------------!
                                 enddo ! over ie2={1,rank(f2)} loop
                             enddo ! over is2={1,rank(f2)} loop
                         enddo ! over ie1={1,rank(f1)} loop
                     enddo ! is1={1,rank(f1)} loop
                 enddo ! over f2={1,f1} loop
             enddo ! over f1={1,norbs} loop

         endif ! back if ( btest(isvrt,3) ) block

     END BLOCK CALC_G2_PP_AABB

! G2_PP_ABBA component
!-------------------------------------------------------------------------
     CALC_G2_PP_ABBA: BLOCK

         if ( btest(isvrt,4) ) then

             do f1=1,norbs                         ! block index: A
                 do f2=1,f1                        ! block index: B
                     do is1=1,rank(f1)             ! \delta -> l: creation operator
                         do ie1=1,rank(f1)         ! \alpha -> i: annihilation operator
                             do is2=1,rank(f2)     ! \beta  -> j: creation operator
                                 do ie2=1,rank(f2) ! \gamma -> k: annihilation operator
             !-------------------!
             do wbn=1,nbfrq                        ! bosonic matsubara frequency: w
                 do l1=1,lemax                     ! legendre polynomial index: l
                     do l2=1,lemax                 ! legendre polynomial index: l'
                         ee = caux2(wbn,ie1,f1) * caux1(wbn,is1,f1)
                         pp = pl_s(l1,is1,is2,f1,f2) * pl_e(l2,ie2,ie1,f2,f1) * lfun(l1,l2)
                         mm = mmat(ie1, is1, f1) * mmat(ie2, is2, f2)
                         if ( f1 == f2 ) then
                             mm = mm - mmat(ie1, is2, f1) * mmat(ie2, is1, f1)
                         endif ! back if ( f1 == f2 ) block

                         g2pp(l2,l1,wbn,f2,f1) = g2pp(l2,l1,wbn,f2,f1) - mm * pp * ee / beta
                         h2pp(l2,l1,wbn,f2,f1) = h2pp(l2,l1,wbn,f2,f1) - mm * pp * ee / beta * pref(ie1,f1)
                     enddo ! over l2={1,lemax} loop
                 enddo ! over l1={1,lemax} loop
             enddo ! over wbn={1,nbfrq} loop
             !-------------------!
                                 enddo ! over ie2={1,rank(f2)} loop
                             enddo ! over is2={1,rank(f2)} loop
                         enddo ! over ie1={1,rank(f1)} loop
                     enddo ! is1={1,rank(f1)} loop
                 enddo ! over f2={1,f1} loop
             enddo ! over f1={1,norbs} loop

         endif ! back if ( btest(isvrt,4) ) block

     END BLOCK CALC_G2_PP_ABBA

! deallocate memory
     deallocate( lfun  )
     deallocate( pl_s  )
     deallocate( pl_e  )
     deallocate( caux1 )
     deallocate( caux2 )

     return
  end subroutine cat_record_g2pp_leg

!!
!! @sub cat_record_g2pp_svd
!!
!! record the two-particle green's and vertex functions in the pp channel.
!! here improved estimator is used to improve the accuracy
!!
!! note:
!!
!!     we try to measure the two-particle green's and vertex functions in
!!     the particle-particle channel and intermediate/matsubara representation
!!     in this subroutine. in order to simplify the calculations, we just
!!     consider the block structure of G^{(2)}
!!
!!     G^{(2)}_{abcd,AABB,pp} (l, l', \omega) =  (-1)^l'
!!         \frac{ 1 }{ \beta }
!!         \langle
!!             \sum^{K_A}_{ij} \sum^{K_B}_{kl}
!!             ( M^{A}_{ij} M^{B}_{kl} - \delta_{AB} M^{A}_{il} M^{B}_{kj} )
!!             u_l( x(\tau_l - \tau_j) ) u_l'( x(\tau'_k - \tau'_i) )
!!             exp [ i \omega (\tau'_i - \tau_l) ]
!!             \delta_{a,i} \delta_{b,j} \delta_{c,k} \delta_{d,l}
!!         \rangle
!!
!!     G^{(2)}_{abcd,ABBA,pp} (l, l', \omega) =  (-1)^l'
!!         \frac{ 1 }{ \beta }
!!         \langle
!!             \sum^{K_A}_{il} \sum^{K_B}_{kj}
!!             ( \delta_{AB} M^{A}_{ij} M^{B}_{kl} - M^{A}_{il} M^{B}_{kj} )
!!             u_l( x(\tau_l - \tau_j) ) u_l'( x(\tau'_k - \tau'_i) )
!!             exp [ i \omega (\tau'_i - \tau_l) ]
!!             \delta_{a,i} \delta_{b,j} \delta_{c,k} \delta_{d,l}
!!         \rangle
!!
!!     \tau'_i and \tau'_k: imaginary time for annihilation operators
!!     \tau_j and \tau_l: imaginary time for creation operators
!!     u_l and u_l': svd polynomial
!!     \omega: bosonic matsubara frequency
!!
  subroutine cat_record_g2pp_svd()
     use constants, only : dp
     use constants, only : zero, one, two, czero

     use control, only : isvrt
     use control, only : norbs
     use control, only : svmax, svgrd
     use control, only : nbfrq
     use control, only : beta

     use context, only : index_s, index_e
     use context, only : time_s, time_e
     use context, only : rep_s
     use context, only : g2pp, h2pp
     use context, only : rank, pref
     use context, only : mmat

     implicit none

! local variables
! loop index for flavor channel
     integer  :: f1
     integer  :: f2

! loop index for frequency
     integer  :: wbn
     integer  :: l1
     integer  :: l2

! loop indices for start and end points
     integer  :: is1
     integer  :: is2
     integer  :: ie1
     integer  :: ie2

! index for imaginary time \tau
     integer  :: curr

! interval for imaginary time slice
     real(dp) :: step

! distance betweem \tau_s and \tau_e
     real(dp) :: dt

! sign for u_l(x(\tau))
     real(dp) :: ms

! real(dp) dummy variables
     real(dp) :: mm
     real(dp) :: uu

! complex(dp) dummy variables
     complex(dp) :: ee

! u_l(x(\tau_s2 - \tau_s1))
     real(dp), allocatable :: ul_s(:,:,:,:,:)

! u_l(x(\tau_e2 - \tau_e1))
     real(dp), allocatable :: ul_e(:,:,:,:,:)

! exp [i \omega_n \tau_s] and exp [i \omega_n \tau_e]
! note here \omega_n is bosonic
     complex(dp), allocatable :: caux1(:,:,:)
     complex(dp), allocatable :: caux2(:,:,:)

! allocate memory
     allocate( ul_s(svmax, maxval(rank), maxval(rank), norbs, norbs)); ul_s = zero
     allocate( ul_e(svmax, maxval(rank), maxval(rank), norbs, norbs)); ul_e = zero

     allocate( caux1(nbfrq, maxval(rank), norbs) ); caux1 = czero
     allocate( caux2(nbfrq, maxval(rank), norbs) ); caux2 = czero

! calculate prefactor: pref
     call ctqmc_make_pref()

! prepare some important arrays: ul_s
     step = real(svgrd - 1) / two
     do f1=1,norbs
         do is1=1,rank(f1)
             do f2=1,norbs
                 do is2=1,rank(f2)
! determine dt (distance) and ms (sign)
                     dt = time_s( index_s(is2, f2), f2 ) - time_s( index_s(is1, f1), f1 )
                     ms = sign(one, dt)

! adjust dt, keep it stay in (zero, beta)
                     if ( dt < zero ) then
                         dt = dt + beta
                     endif ! back if ( dt < zero ) block

! determine index for imaginary time
                     call s_svd_point(two * dt / beta - one, step, curr)

! fill ul_s
                     do l1=1,svmax
                         ul_s(l1,is2,is1,f2,f1) = ms * rep_s(curr,l1)
                     enddo ! over l1={1,svmax} loop
                 enddo ! over is2={1,rank(f2)} loop
             enddo ! over f2={1,norbs} loop
         enddo ! over is1={1,rank(f1)} loop
     enddo ! over f1={1,norbs} loop

! prepare some important arrays: ul_e
     step = real(svgrd - 1) / two
     do f1=1,norbs
         do ie1=1,rank(f1)
             do f2=1,norbs
                 do ie2=1,rank(f2)
! determine dt (distance) and ms (sign)
                     dt = time_e( index_e(ie2, f2), f2 ) - time_e( index_e(ie1, f1), f1 )
                     ms = sign(one, dt)

! adjust dt, keep it stay in (zero, beta)
                     if ( dt < zero ) then
                         dt = dt + beta
                     endif ! back if ( dt < zero ) block

! determine index for imaginary time
                     call s_svd_point(two * dt / beta - one, step, curr)

! fill ul_e
                     do l1=1,svmax
                         ul_e(l1,ie2,ie1,f2,f1) = ms * rep_s(curr,l1)
                     enddo ! over l1={1,svmax} loop
                 enddo ! over ie2={1,rank(f2)} loop
             enddo ! over f2={1,norbs} loop
         enddo ! over ie1={1,rank(f1)} loop
     enddo ! over f1={1,norbs} loop

! prepare some important arrays: caux1 and caux2
     do f1=1,norbs
         call ctqmc_make_bexp(f1, nbfrq, maxval(rank), caux1(:,:,f1), caux2(:,:,f1))
     enddo ! over f1={1,norbs} loop

! calculate g2pp and h2pp
!
! G2_PP_AABB component
!-------------------------------------------------------------------------
     CALC_G2_PP_AABB: BLOCK

         if ( btest(isvrt,3) ) then

             do f1=1,norbs                         ! block index: A
                 do f2=1,f1                        ! block index: B
                     do is1=1,rank(f1)             ! \beta  -> j: creation operator
                         do ie1=1,rank(f1)         ! \alpha -> i: annihilation operator
                             do is2=1,rank(f2)     ! \delta -> l: creation operator
                                 do ie2=1,rank(f2) ! \gamma -> k: annihilation operator
             !-------------------!
             do wbn=1,nbfrq                        ! bosonic matsubara frequency: w
                 do l1=1,svmax                     ! svd polynomial index: l
                     do l2=1,svmax                 ! svd polynomial index: l'
                         ee = caux2(wbn,ie1,f1) * caux1(wbn,is2,f2)
                         uu = ul_s(l1,is2,is1,f2,f1) * ul_e(l2,ie2,ie1,f2,f1)
                         mm = mmat(ie1, is1, f1) * mmat(ie2, is2, f2)
                         if ( f1 == f2 ) then
                             mm = mm - mmat(ie1, is2, f1) * mmat(ie2, is1, f1)
                         endif ! back if ( f1 == f2 ) block

                         g2pp(l2,l1,wbn,f2,f1) = g2pp(l2,l1,wbn,f2,f1) + mm * uu * ee / beta
                         h2pp(l2,l1,wbn,f2,f1) = h2pp(l2,l1,wbn,f2,f1) + mm * uu * ee / beta * pref(ie1,f1)
                     enddo ! over l2={1,svmax} loop
                 enddo ! over l1={1,svmax} loop
             enddo ! over wbn={1,nbfrq} loop
             !-------------------!
                                 enddo ! over ie2={1,rank(f2)} loop
                             enddo ! over is2={1,rank(f2)} loop
                         enddo ! over ie1={1,rank(f1)} loop
                     enddo ! is1={1,rank(f1)} loop
                 enddo ! over f2={1,f1} loop
             enddo ! over f1={1,norbs} loop

         endif ! back if ( btest(isvrt,3) ) block

     END BLOCK CALC_G2_PP_AABB

! G2_PP_ABBA component
!-------------------------------------------------------------------------
     CALC_G2_PP_ABBA: BLOCK

         if ( btest(isvrt,4) ) then

             do f1=1,norbs                         ! block index: A
                 do f2=1,f1                        ! block index: B
                     do is1=1,rank(f1)             ! \delta -> l: creation operator
                         do ie1=1,rank(f1)         ! \alpha -> i: annihilation operator
                             do is2=1,rank(f2)     ! \beta  -> j: creation operator
                                 do ie2=1,rank(f2) ! \gamma -> k: annihilation operator
             !-------------------!
             do wbn=1,nbfrq                        ! bosonic matsubara frequency: w
                 do l1=1,svmax                     ! svd polynomial index: l
                     do l2=1,svmax                 ! svd polynomial index: l'
                         ee = caux2(wbn,ie1,f1) * caux1(wbn,is1,f1)
                         uu = ul_s(l1,is1,is2,f1,f2) * ul_e(l2,ie2,ie1,f2,f1)
                         mm = mmat(ie1, is1, f1) * mmat(ie2, is2, f2)
                         if ( f1 == f2 ) then
                             mm = mm - mmat(ie1, is2, f1) * mmat(ie2, is1, f1)
                         endif ! back if ( f1 == f2 ) block

                         g2pp(l2,l1,wbn,f2,f1) = g2pp(l2,l1,wbn,f2,f1) - mm * uu * ee / beta
                         h2pp(l2,l1,wbn,f2,f1) = h2pp(l2,l1,wbn,f2,f1) - mm * uu * ee / beta * pref(ie1,f1)
                     enddo ! over l2={1,svmax} loop
                 enddo ! over l1={1,svmax} loop
             enddo ! over wbn={1,nbfrq} loop
             !-------------------!
                                 enddo ! over ie2={1,rank(f2)} loop
                             enddo ! over is2={1,rank(f2)} loop
                         enddo ! over ie1={1,rank(f1)} loop
                     enddo ! is1={1,rank(f1)} loop
                 enddo ! over f2={1,f1} loop
             enddo ! over f1={1,norbs} loop

         endif ! back if ( btest(isvrt,4) ) block

     END BLOCK CALC_G2_PP_ABBA

! deallocate memory
     deallocate( ul_s  )
     deallocate( ul_e  )
     deallocate( caux1 )
     deallocate( caux2 )

     return
  end subroutine cat_record_g2pp_svd

!!========================================================================
!!>>> reduce autocorrelation function                                  <<<
!!========================================================================

!!
!! @sub ctqmc_reduce_ac_f
!!
!! reduce the ac_f from all children processes
!!
  subroutine ctqmc_reduce_ac_f(ac_f_mpi, ac_f_err)
     use constants, only : dp
     use constants, only : zero

     use mmpi, only : mp_allreduce
     use mmpi, only : mp_barrier

     use control, only : ntime
     use control, only : nprocs

     use context, only : ac_f

     implicit none

! external arguments
! autocorrelation function
     real(dp), intent(out) :: ac_f_mpi(ntime + 2)
     real(dp), intent(out) :: ac_f_err(ntime + 2)

! initialize ac_f_mpi and ac_f_err
     ac_f_mpi = zero
     ac_f_err = zero

! build ac_f_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(ac_f, ac_f_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     ac_f_mpi = ac_f

# endif /* MPI */

! calculate the average
     ac_f_mpi = ac_f_mpi / real(nprocs)

! build ac_f_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce((ac_f - ac_f_mpi)**2, ac_f_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         ac_f_err = sqrt( ac_f_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

     return
  end subroutine ctqmc_reduce_ac_f

!!========================================================================
!!>>> reduce physical observables 1                                    <<<
!!========================================================================

!!
!! @sub ctqmc_reduce_hist
!!
!! reduce the hist from all children processes
!!
  subroutine ctqmc_reduce_hist(hist_mpi, hist_err)
     use constants, only : dp
     use constants, only : zero

     use mmpi, only : mp_allreduce
     use mmpi, only : mp_barrier

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
     use constants, only : dp
     use constants, only : zero

     use mmpi, only : mp_allreduce
     use mmpi, only : mp_barrier

     use control, only : ncfgs
     use control, only : nprocs

     use context, only : prob

     implicit none

! external arguments
! probability of atomic eigenstates
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
!! @sub ctqmc_reduce_paux
!!
!! reduce the paux from all children processes
!!
  subroutine ctqmc_reduce_paux(paux_mpi, paux_err)
     use constants, only : dp
     use constants, only : zero

     use mmpi, only : mp_allreduce
     use mmpi, only : mp_barrier

     use control, only : nprocs

     use context, only : paux

     implicit none

! external arguments
! auxiliary physical observables
     real(dp), intent(out) :: paux_mpi(  9  )
     real(dp), intent(out) :: paux_err(  9  )

! initialize paux_mpi and paux_err
     paux_mpi = zero
     paux_err = zero

! build paux_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(paux, paux_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     paux_mpi = paux

# endif /* MPI */

! calculate the average
     paux_mpi = paux_mpi / real(nprocs)

! build paux_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce((paux - paux_mpi)**2, paux_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         paux_err = sqrt( paux_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

     return
  end subroutine ctqmc_reduce_paux

!!
!! @sub ctqmc_reduce_nmat
!!
!! reduce the nimp and nmat from all children processes
!!
  subroutine ctqmc_reduce_nmat(nimp_mpi, nmat_mpi, nimp_err, nmat_err)
     use constants, only : dp
     use constants, only : zero

     use mmpi, only : mp_allreduce
     use mmpi, only : mp_barrier

     use control, only : norbs
     use control, only : nprocs

     use context, only : nimp, nmat

     implicit none

! external arguments
! occupation number
     real(dp), intent(out) :: nimp_mpi(norbs)
     real(dp), intent(out) :: nimp_err(norbs)

! double occupation number matrix
     real(dp), intent(out) :: nmat_mpi(norbs,norbs)
     real(dp), intent(out) :: nmat_err(norbs,norbs)

! initialize nimp_mpi and nmat_mpi, nimp_err and nmat_err
     nimp_mpi = zero
     nmat_mpi = zero

     nimp_err = zero
     nmat_err = zero

! build nimp_mpi and nmat_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(nimp, nimp_mpi)
     call mp_allreduce(nmat, nmat_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     nimp_mpi = nimp
     nmat_mpi = nmat

# endif /* MPI */

! calculate the average
     nimp_mpi = nimp_mpi / real(nprocs)
     nmat_mpi = nmat_mpi / real(nprocs)

! build nimp_err and nmat_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce((nimp - nimp_mpi)**2, nimp_err)
     call mp_allreduce((nmat - nmat_mpi)**2, nmat_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         nimp_err = sqrt( nimp_err / real( nprocs * ( nprocs - 1 ) ) )
         nmat_err = sqrt( nmat_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

     return
  end subroutine ctqmc_reduce_nmat

!!========================================================================
!!>>> reduce physical observables 2                                    <<<
!!========================================================================

!!
!! @sub ctqmc_reduce_gtau
!!
!! reduce the gtau from all children processes
!!
  subroutine ctqmc_reduce_gtau(gtau_mpi, gtau_err)
     use constants, only : dp
     use constants, only : zero

     use mmpi, only : mp_allreduce
     use mmpi, only : mp_barrier

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
     use constants, only : dp
     use constants, only : zero

     use mmpi, only : mp_allreduce
     use mmpi, only : mp_barrier

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
     use constants, only : dp
     use constants, only : zero, czero, czi

     use mmpi, only : mp_allreduce
     use mmpi, only : mp_barrier

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
     real(dp), allocatable :: g_re_err(:,:,:)
     real(dp), allocatable :: g_im_err(:,:,:)

! allocate memory
     allocate(g_re_err(mfreq,norbs,norbs))
     allocate(g_im_err(mfreq,norbs,norbs))

! initialize g_re_err and g_im_err
     g_re_err = zero
     g_im_err = zero

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
     call mp_allreduce(( real(grnf - grnf_mpi))**2, g_re_err)
     call mp_allreduce((aimag(grnf - grnf_mpi))**2, g_im_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         g_re_err = sqrt( g_re_err / real( nprocs * ( nprocs - 1 ) ) )
         g_im_err = sqrt( g_im_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

! construct the final grnf_err
     grnf_err = g_re_err + g_im_err * czi

! deallocate memory
     deallocate( g_re_err )
     deallocate( g_im_err )

     return
  end subroutine ctqmc_reduce_grnf

!!
!! @sub ctqmc_reduce_frnf
!!
!! reduce the frnf from all children processes
!!
  subroutine ctqmc_reduce_frnf(frnf_mpi, frnf_err)
     use constants, only : dp
     use constants, only : zero, czero, czi

     use mmpi, only : mp_allreduce
     use mmpi, only : mp_barrier

     use control, only : norbs
     use control, only : mfreq
     use control, only : nprocs

     use context, only : frnf

     implicit none

! external arguments
! auxiliary correlation function
     complex(dp), intent(out) :: frnf_mpi(mfreq,norbs,norbs)
     complex(dp), intent(out) :: frnf_err(mfreq,norbs,norbs)

! local variables
! used to store the real and imaginary parts of auxiliary correlation function
     real(dp), allocatable :: f_re_err(:,:,:)
     real(dp), allocatable :: f_im_err(:,:,:)

! allocate memory
     allocate(f_re_err(mfreq,norbs,norbs))
     allocate(f_im_err(mfreq,norbs,norbs))

! initialize f_re_err and f_im_err
     f_re_err = zero
     f_im_err = zero

! initialize frnf_mpi and frnf_err
     frnf_mpi = czero
     frnf_err = czero

! build frnf_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(frnf, frnf_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     frnf_mpi = frnf

# endif /* MPI */

! calculate the average
     frnf_mpi = frnf_mpi / real(nprocs)

! build frnf_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(( real(frnf - frnf_mpi))**2, f_re_err)
     call mp_allreduce((aimag(frnf - frnf_mpi))**2, f_im_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         f_re_err = sqrt( f_re_err / real( nprocs * ( nprocs - 1 ) ) )
         f_im_err = sqrt( f_im_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

! construct the final frnf_err
     frnf_err = f_re_err + f_im_err * czi

! deallocate memory
     deallocate( f_re_err )
     deallocate( f_im_err )

     return
  end subroutine ctqmc_reduce_frnf

!!
!! @sub ctqmc_reduce_sig2
!!
!! reduce the sig2 from all children processes
!!
  subroutine ctqmc_reduce_sig2(sig2_mpi, sig2_err)
     use constants, only : dp
     use constants, only : zero, czero, czi

     use mmpi, only : mp_allreduce
     use mmpi, only : mp_barrier

     use control, only : norbs
     use control, only : mfreq
     use control, only : nprocs

     use context, only : sig2

     implicit none

! external arguments
! self-energy function
     complex(dp), intent(out) :: sig2_mpi(mfreq,norbs,norbs)
     complex(dp), intent(out) :: sig2_err(mfreq,norbs,norbs)

! local variables
! used to store the real and imaginary parts of self-energy function
     real(dp), allocatable :: s_re_err(:,:,:)
     real(dp), allocatable :: s_im_err(:,:,:)

! allocate memory
     allocate(s_re_err(mfreq,norbs,norbs))
     allocate(s_im_err(mfreq,norbs,norbs))

! initialize s_re_err and s_im_err
     s_re_err = zero
     s_im_err = zero

! initialize sig2_mpi and sig2_err
     sig2_mpi = czero
     sig2_err = czero

! build sig2_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(sig2, sig2_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     sig2_mpi = sig2

# endif /* MPI */

! calculate the average
     sig2_mpi = sig2_mpi / real(nprocs)

! build sig2_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(( real(sig2 - sig2_mpi))**2, s_re_err)
     call mp_allreduce((aimag(sig2 - sig2_mpi))**2, s_im_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         s_re_err = sqrt( s_re_err / real( nprocs * ( nprocs - 1 ) ) )
         s_im_err = sqrt( s_im_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

! construct the final sig2_err
     sig2_err = s_re_err + s_im_err * czi

! deallocate memory
     deallocate( s_re_err )
     deallocate( s_im_err )

     return
  end subroutine ctqmc_reduce_sig2

!!========================================================================
!!>>> reduce physical observables 3                                    <<<
!!========================================================================

!!
!! @sub ctqmc_reduce_kmat
!!
!! reduce the knop and kmat from all children processes
!!
  subroutine ctqmc_reduce_kmat(knop_mpi, kmat_mpi, knop_err, kmat_err)
     use constants, only : dp
     use constants, only : zero

     use mmpi, only : mp_allreduce
     use mmpi, only : mp_barrier

     use control, only : isobs
     use control, only : norbs
     use control, only : nprocs

     use context, only : knop, kmat

     implicit none

! external arguments
! number of operators
     real(dp), intent(out) :: knop_mpi(norbs)
     real(dp), intent(out) :: knop_err(norbs)

! crossing product of k_i and k_j
     real(dp), intent(out) :: kmat_mpi(norbs,norbs)
     real(dp), intent(out) :: kmat_err(norbs,norbs)

! check whether this observable has been measured
     if ( .not. btest(isobs, 1) ) RETURN

! initialize knop_mpi and kmat_mpi, knop_err and kmat_err
     knop_mpi = zero
     kmat_mpi = zero

     knop_err = zero
     kmat_err = zero

! build knop_mpi and kmat_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(knop, knop_mpi)
     call mp_allreduce(kmat, kmat_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     knop_mpi = knop
     kmat_mpi = kmat

# endif /* MPI */

! calculate the average
     knop_mpi = knop_mpi / real(nprocs)
     kmat_mpi = kmat_mpi / real(nprocs)

! build knop_err and kmat_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce((knop - knop_mpi)**2, knop_err)
     call mp_allreduce((kmat - kmat_mpi)**2, kmat_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         knop_err = sqrt( knop_err / real( nprocs * ( nprocs - 1 ) ) )
         kmat_err = sqrt( kmat_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

     return
  end subroutine ctqmc_reduce_kmat

!!
!! @sub ctqmc_reduce_lrmm
!!
!! reduce the lnop, rnop, and lrmm from all children processes
!!
  subroutine ctqmc_reduce_lrmm(lnop_mpi, rnop_mpi, lrmm_mpi, lnop_err, rnop_err, lrmm_err)
     use constants, only : dp
     use constants, only : zero

     use mmpi, only : mp_allreduce
     use mmpi, only : mp_barrier

     use control, only : isobs
     use control, only : norbs
     use control, only : nprocs

     use context, only : lnop, rnop, lrmm

     implicit none

! external arguments
! number of operators at left half axis
     real(dp), intent(out) :: lnop_mpi(norbs)
     real(dp), intent(out) :: lnop_err(norbs)

! number of operators at right half axis
     real(dp), intent(out) :: rnop_mpi(norbs)
     real(dp), intent(out) :: rnop_err(norbs)

! crossing product of k_l and k_r
     real(dp), intent(out) :: lrmm_mpi(norbs,norbs)
     real(dp), intent(out) :: lrmm_err(norbs,norbs)

! check whether this observable has been measured
     if ( .not. btest(isobs, 2) ) RETURN

! initialize lnop_mpi, rnop_mpi, and lrmm_mpi
! initialize lnop_err, rnop_err, and lrmm_err
     lnop_mpi = zero
     rnop_mpi = zero
     lrmm_mpi = zero

     lnop_err = zero
     rnop_err = zero
     lrmm_err = zero

! build lnop_mpi, rnop_mpi, and lrmm_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(lnop, lnop_mpi)
     call mp_allreduce(rnop, rnop_mpi)
     call mp_allreduce(lrmm, lrmm_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     lnop_mpi = lnop
     rnop_mpi = rnop
     lrmm_mpi = lrmm

# endif /* MPI */

! calculate the average
     lnop_mpi = lnop_mpi / real(nprocs)
     rnop_mpi = rnop_mpi / real(nprocs)
     lrmm_mpi = lrmm_mpi / real(nprocs)

! build lnop_err, rnop_err, and lrmm_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce((lnop - lnop_mpi)**2, lnop_err)
     call mp_allreduce((rnop - rnop_mpi)**2, rnop_err)
     call mp_allreduce((lrmm - lrmm_mpi)**2, lrmm_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         lnop_err = sqrt( lnop_err / real( nprocs * ( nprocs - 1 ) ) )
         rnop_err = sqrt( rnop_err / real( nprocs * ( nprocs - 1 ) ) )
         lrmm_err = sqrt( lrmm_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

     return
  end subroutine ctqmc_reduce_lrmm

!!
!! @sub ctqmc_reduce_szpw
!!
!! reduce the szpw from all children processes
!!
  subroutine ctqmc_reduce_szpw(szpw_mpi, szpw_err)
     use constants, only : dp
     use constants, only : zero

     use mmpi, only : mp_allreduce
     use mmpi, only : mp_barrier

     use control, only : isobs
     use control, only : norbs
     use control, only : nprocs

     use context, only : szpw

     implicit none

! external arguments
! powers of local magnetization, orbital-resolved
     real(dp), intent(out) :: szpw_mpi(4,norbs)
     real(dp), intent(out) :: szpw_err(4,norbs)

! check whether this observable has been measured
     if ( .not. btest(isobs, 3) ) RETURN

! initialize szpw_mpi and szpw_err
     szpw_mpi = zero
     szpw_err = zero

! build szpw_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(szpw, szpw_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     szpw_mpi = szpw

# endif /* MPI */

! calculate the average
     szpw_mpi = szpw_mpi / real(nprocs)

! build szpw_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce((szpw - szpw_mpi)**2, szpw_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         szpw_err = sqrt( szpw_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

     return
  end subroutine ctqmc_reduce_szpw

!!========================================================================
!!>>> reduce physical observables 4                                    <<<
!!========================================================================

!!
!! @sub ctqmc_reduce_sp_t
!!
!! reduce the schi and sp_t from all children processes
!!
  subroutine ctqmc_reduce_sp_t(schi_mpi, sp_t_mpi, schi_err, sp_t_err)
     use constants, only : dp
     use constants, only : zero

     use mmpi, only : mp_allreduce
     use mmpi, only : mp_barrier

     use control, only : issus
     use control, only : nband
     use control, only : ntime
     use control, only : nprocs

     use context, only : schi, sp_t

     implicit none

! external arguments
! spin-spin correlation function, totally-averaged
     real(dp), intent(out) :: schi_mpi(ntime)
     real(dp), intent(out) :: schi_err(ntime)

! spin-spin correlation function, orbital-resolved
     real(dp), intent(out) :: sp_t_mpi(ntime,nband)
     real(dp), intent(out) :: sp_t_err(ntime,nband)

! check whether this observable has been measured
     if ( .not. btest(issus, 1) ) RETURN

! initialize schi_mpi and sp_t_mpi, schi_err and sp_t_err
     schi_mpi = zero
     sp_t_mpi = zero

     schi_err = zero
     sp_t_err = zero

! build schi_mpi and sp_t_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(schi, schi_mpi)
     call mp_allreduce(sp_t, sp_t_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     schi_mpi = schi
     sp_t_mpi = sp_t

# endif /* MPI */

! calculate the average
     schi_mpi = schi_mpi / real(nprocs)
     sp_t_mpi = sp_t_mpi / real(nprocs)

! build schi_err and sp_t_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce((schi - schi_mpi)**2, schi_err)
     call mp_allreduce((sp_t - sp_t_mpi)**2, sp_t_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         schi_err = sqrt( schi_err / real( nprocs * ( nprocs - 1 ) ) )
         sp_t_err = sqrt( sp_t_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

     return
  end subroutine ctqmc_reduce_sp_t

!!
!! @sub ctqmc_reduce_sp_w
!!
!! reduce the sp_w from all children processes
!!
  subroutine ctqmc_reduce_sp_w(sp_w_mpi, sp_w_err)
     use constants, only : dp
     use constants, only : zero

     use mmpi, only : mp_allreduce
     use mmpi, only : mp_barrier

     use control, only : issus
     use control, only : nband
     use control, only : nbfrq
     use control, only : nprocs

     use context, only : sp_w

     implicit none

! external arguments
! spin-spin correlation function, orbital-resolved
     real(dp), intent(out) :: sp_w_mpi(nbfrq,nband)
     real(dp), intent(out) :: sp_w_err(nbfrq,nband)

! check whether this observable has been measured
     if ( .not. btest(issus, 3) ) RETURN

! initialize sp_w_mpi and sp_w_err
     sp_w_mpi = zero
     sp_w_err = zero

! build sp_w_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(sp_w, sp_w_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     sp_w_mpi = sp_w

# endif /* MPI */

! calculate the average
     sp_w_mpi = sp_w_mpi / real(nprocs)

! build sp_w_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce((sp_w - sp_w_mpi)**2, sp_w_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         sp_w_err = sqrt( sp_w_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

     return
  end subroutine ctqmc_reduce_sp_w

!!
!! @sub ctqmc_reduce_ch_t
!!
!! reduce the cchi and ch_t from all children processes
!!
  subroutine ctqmc_reduce_ch_t(cchi_mpi, ch_t_mpi, cchi_err, ch_t_err)
     use constants, only : dp
     use constants, only : zero

     use mmpi, only : mp_allreduce
     use mmpi, only : mp_barrier

     use control, only : issus
     use control, only : norbs
     use control, only : ntime
     use control, only : nprocs

     use context, only : cchi, ch_t

     implicit none

! external arguments
! charge-charge correlation function, totally-averaged
     real(dp), intent(out) :: cchi_mpi(ntime)
     real(dp), intent(out) :: cchi_err(ntime)

! charge-charge correlation function, orbital-resolved
     real(dp), intent(out) :: ch_t_mpi(ntime,norbs,norbs)
     real(dp), intent(out) :: ch_t_err(ntime,norbs,norbs)

! check whether this observable has been measured
     if ( .not. btest(issus, 2) ) RETURN

! initialize cchi_mpi and ch_t_mpi, cchi_err and ch_t_err
     cchi_mpi = zero
     ch_t_mpi = zero

     cchi_err = zero
     ch_t_err = zero

! build cchi_mpi and ch_t_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(cchi, cchi_mpi)
     call mp_allreduce(ch_t, ch_t_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     cchi_mpi = cchi
     ch_t_mpi = ch_t

# endif /* MPI */

! calculate the average
     cchi_mpi = cchi_mpi / real(nprocs)
     ch_t_mpi = ch_t_mpi / real(nprocs)

! build cchi_err and ch_t_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce((cchi - cchi_mpi)**2, cchi_err)
     call mp_allreduce((ch_t - ch_t_mpi)**2, ch_t_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         cchi_err = sqrt( cchi_err / real( nprocs * ( nprocs - 1 ) ) )
         ch_t_err = sqrt( ch_t_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

     return
  end subroutine ctqmc_reduce_ch_t

!!
!! @sub ctqmc_reduce_ch_w
!!
!! reduce the ch_w from all children processes
!!
  subroutine ctqmc_reduce_ch_w(ch_w_mpi, ch_w_err)
     use constants, only : dp
     use constants, only : zero

     use mmpi, only : mp_allreduce
     use mmpi, only : mp_barrier

     use control, only : issus
     use control, only : norbs
     use control, only : nbfrq
     use control, only : nprocs

     use context, only : ch_w

     implicit none

! external arguments
! charge-charge correlation function, orbital-resolved
     real(dp), intent(out) :: ch_w_mpi(nbfrq,norbs,norbs)
     real(dp), intent(out) :: ch_w_err(nbfrq,norbs,norbs)

! check whether this observable has been measured
     if ( .not. btest(issus, 4) ) RETURN

! initialize ch_w_mpi and ch_w_err
     ch_w_mpi = zero
     ch_w_err = zero

! build ch_w_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(ch_w, ch_w_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     ch_w_mpi = ch_w

# endif /* MPI */

! calculate the average
     ch_w_mpi = ch_w_mpi / real(nprocs)

! build ch_w_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce((ch_w - ch_w_mpi)**2, ch_w_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         ch_w_err = sqrt( ch_w_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

     return
  end subroutine ctqmc_reduce_ch_w

!!========================================================================
!!>>> reduce physical observables 5                                    <<<
!!========================================================================

!!
!! @sub ctqmc_reduce_g2ph
!!
!! reduce the g2ph and h2ph from all children processes
!!
  subroutine ctqmc_reduce_g2ph(g2ph_mpi, h2ph_mpi, g2ph_err, h2ph_err)
     use constants, only : dp
     use constants, only : zero, czero, czi

     use mmpi, only : mp_allreduce
     use mmpi, only : mp_barrier

     use control, only : isvrt
     use control, only : norbs
     use control, only : nffrq, nbfrq
     use control, only : nprocs

     use context, only : g2ph
     use context, only : h2ph

     implicit none

! external arguments
! two-particle green's function
     complex(dp), intent(out) :: g2ph_mpi(nffrq,nffrq,nbfrq,norbs,norbs)
     complex(dp), intent(out) :: g2ph_err(nffrq,nffrq,nbfrq,norbs,norbs)

! two-particle vertex function
     complex(dp), intent(out) :: h2ph_mpi(nffrq,nffrq,nbfrq,norbs,norbs)
     complex(dp), intent(out) :: h2ph_err(nffrq,nffrq,nbfrq,norbs,norbs)

! local variables
! used to store the real and imaginary parts of green's function
     real(dp), allocatable :: g_re_err(:,:,:,:,:)
     real(dp), allocatable :: g_im_err(:,:,:,:,:)

! used to store the real and imaginary parts of vertex function
     real(dp), allocatable :: h_re_err(:,:,:,:,:)
     real(dp), allocatable :: h_im_err(:,:,:,:,:)

! check whether this observable has been measured
     if ( .not. ( btest(isvrt, 1) .or. btest(isvrt, 2) ) ) RETURN

! allocate memory
     allocate(g_re_err(nffrq,nffrq,nbfrq,norbs,norbs))
     allocate(g_im_err(nffrq,nffrq,nbfrq,norbs,norbs))
     allocate(h_re_err(nffrq,nffrq,nbfrq,norbs,norbs))
     allocate(h_im_err(nffrq,nffrq,nbfrq,norbs,norbs))

! initialize g_re_err and g_im_err
     g_re_err = zero
     g_im_err = zero

! initialize h_re_err and h_im_err
     h_re_err = zero
     h_im_err = zero

! initialize g2ph_mpi and g2ph_err
     g2ph_mpi = czero
     g2ph_err = czero

! initialize h2ph_mpi and h2ph_err
     h2ph_mpi = czero
     h2ph_err = czero

! build g2ph_mpi and h2ph_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(g2ph, g2ph_mpi)
     call mp_allreduce(h2ph, h2ph_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     g2ph_mpi = g2ph
     h2ph_mpi = h2ph

# endif /* MPI */

! calculate the average
     g2ph_mpi = g2ph_mpi / real(nprocs)
     h2ph_mpi = h2ph_mpi / real(nprocs)

! build g2ph_err and h2ph_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(( real(g2ph - g2ph_mpi))**2, g_re_err)
     call mp_allreduce((aimag(g2ph - g2ph_mpi))**2, g_im_err)
     call mp_allreduce(( real(h2ph - h2ph_mpi))**2, h_re_err)
     call mp_allreduce((aimag(h2ph - h2ph_mpi))**2, h_im_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         g_re_err = sqrt( g_re_err / real( nprocs * ( nprocs - 1 ) ) )
         g_im_err = sqrt( g_im_err / real( nprocs * ( nprocs - 1 ) ) )
         h_re_err = sqrt( h_re_err / real( nprocs * ( nprocs - 1 ) ) )
         h_im_err = sqrt( h_im_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

! construct the final g2ph_err and h2ph_err
     g2ph_err = g_re_err + g_im_err * czi
     h2ph_err = h_re_err + h_im_err * czi

! deallocate memory
     deallocate(g_re_err)
     deallocate(g_im_err)
     deallocate(h_re_err)
     deallocate(h_im_err)

     return
  end subroutine ctqmc_reduce_g2ph

!!
!! @sub ctqmc_reduce_g2pp
!!
!! reduce the g2pp and h2pp from all children processes
!!
  subroutine ctqmc_reduce_g2pp(g2pp_mpi, h2pp_mpi, g2pp_err, h2pp_err)
     use constants, only : dp
     use constants, only : zero, czero, czi

     use mmpi, only : mp_allreduce
     use mmpi, only : mp_barrier

     use control, only : isvrt
     use control, only : norbs
     use control, only : nffrq, nbfrq
     use control, only : nprocs

     use context, only : g2pp
     use context, only : h2pp

     implicit none

! external arguments
! two-particle green's function
     complex(dp), intent(out) :: g2pp_mpi(nffrq,nffrq,nbfrq,norbs,norbs)
     complex(dp), intent(out) :: g2pp_err(nffrq,nffrq,nbfrq,norbs,norbs)

! two-particle vertex function
     complex(dp), intent(out) :: h2pp_mpi(nffrq,nffrq,nbfrq,norbs,norbs)
     complex(dp), intent(out) :: h2pp_err(nffrq,nffrq,nbfrq,norbs,norbs)

! local variables
! used to store the real and imaginary parts of green's function
     real(dp), allocatable :: g_re_err(:,:,:,:,:)
     real(dp), allocatable :: g_im_err(:,:,:,:,:)

! used to store the real and imaginary parts of vertex function
     real(dp), allocatable :: h_re_err(:,:,:,:,:)
     real(dp), allocatable :: h_im_err(:,:,:,:,:)

! check whether this observable has been measured
     if ( .not. ( btest(isvrt, 3) .or. btest(isvrt, 4) ) ) RETURN

! allocate memory
     allocate(g_re_err(nffrq,nffrq,nbfrq,norbs,norbs))
     allocate(g_im_err(nffrq,nffrq,nbfrq,norbs,norbs))
     allocate(h_re_err(nffrq,nffrq,nbfrq,norbs,norbs))
     allocate(h_im_err(nffrq,nffrq,nbfrq,norbs,norbs))

! initialize g_re_err and g_im_err
     g_re_err = zero
     g_im_err = zero

! initialize h_re_err and h_im_err
     h_re_err = zero
     h_im_err = zero

! initialize g2pp_mpi and g2pp_err
     g2pp_mpi = czero
     g2pp_err = czero

! initialize h2pp_mpi and h2pp_err
     h2pp_mpi = czero
     h2pp_err = czero

! build g2pp_mpi and h2pp_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(g2pp, g2pp_mpi)
     call mp_allreduce(h2pp, h2pp_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     g2pp_mpi = g2pp
     h2pp_mpi = h2pp

# endif /* MPI */

! calculate the average
     g2pp_mpi = g2pp_mpi / real(nprocs)
     h2pp_mpi = h2pp_mpi / real(nprocs)

! build g2pp_err and h2pp_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(( real(g2pp - g2pp_mpi))**2, g_re_err)
     call mp_allreduce((aimag(g2pp - g2pp_mpi))**2, g_im_err)
     call mp_allreduce(( real(h2pp - h2pp_mpi))**2, h_re_err)
     call mp_allreduce((aimag(h2pp - h2pp_mpi))**2, h_im_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         g_re_err = sqrt( g_re_err / real( nprocs * ( nprocs - 1 ) ) )
         g_im_err = sqrt( g_im_err / real( nprocs * ( nprocs - 1 ) ) )
         h_re_err = sqrt( h_re_err / real( nprocs * ( nprocs - 1 ) ) )
         h_im_err = sqrt( h_im_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

! construct the final g2pp_err and h2pp_err
     g2pp_err = g_re_err + g_im_err * czi
     h2pp_err = h_re_err + h_im_err * czi

! deallocate memory
     deallocate(g_re_err)
     deallocate(g_im_err)
     deallocate(h_re_err)
     deallocate(h_im_err)

     return
  end subroutine ctqmc_reduce_g2pp
