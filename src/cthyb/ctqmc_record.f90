!!!-----------------------------------------------------------------------
!!! project : manjushaka
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
!!!           08/15/2017 by li huang (last modified)
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
     implicit none

     CONTINUE

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

     use context, only : ckink, csign, caves
     use context, only : hist

     implicit none

! record current sign as a byproduct
     caves = caves + csign

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
     use control, only : ncfgs
     use context, only : csign, c_mtr
     use context, only : prob
     use context, only : diag

     implicit none

! local variables
! loop index
     integer :: i

     do i=1,ncfgs
         prob(i) = prob(i) + csign * diag(i,2) / c_mtr
     enddo ! over i={1,ncfgs} loop

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
     use constants, only : dp, zero, two

     use control, only : norbs, ncfgs
     use control, only : Uc, mune, beta
     use context, only : ckink, csign, c_mtr
     use context, only : paux, nimp, nmat
     use context, only : diag, eigs

     use m_sect, only : nsect
     use m_sect, only : sectors

     implicit none

! local variables
! loop index
     integer  :: i
     integer  :: j

! start index of sectors
     integer  :: indx

! current occupation number and Sz
     real(dp) :: nele
     real(dp) :: sz

! current probability for eigenstates
     real(dp) :: cprob(ncfgs)

! current probability for sectors
     real(dp) :: sprob(nsect)

! evaluate cprob at first, it is current atomic probability
     do i=1,ncfgs
         cprob(i) = diag(i,2) / c_mtr
     enddo ! over i={1,ncfgs} loop

! evaluate sprob, it is current sector prob
     sprob = zero
     do i=1,nsect
         indx = sectors(i)%istart
         do j=1,sectors(i)%ndim
             sprob(i) = sprob(i) + cprob(indx+j-1)
         enddo ! over j={1,sectors(i)%ndim} loop
     enddo ! over i={1,nsect} loop

! evaluate the total occupation number
! this algorithm is somewhat rough, not very accurate
     nele = zero
     do i=1,nsect
         nele = nele + sectors(i)%nele * sprob(i)
     enddo ! over i={1,nsect} loop

! evaluate the total Sz
! this algorithm is somewhat rough, and only useful when the Sz quantum
! number is used to generate the atom.cix
     sz = zero
     do i=1,nsect
         sz = sz + sectors(i)%sz * sprob(i)
     enddo ! over i={1,nsect} loop

! evaluate occupation matrix: < n_i >
! equation : Tr ( e^{- \beta H} c^{\dag}_i c_i ) / Tr ( e^{- \beta H} )
!-------------------------------------------------------------------------
     nimp = zero
! this feature will not be implemented for majushaka code
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate double occupation matrix: < n_i n_j >
! equation : Tr ( e^{- \beta H} c^{\dag}_i c_i c^{\dag}_j c_j ) / Tr ( e^{- \beta H} )
!-------------------------------------------------------------------------
     nmat = zero
! this feature will not be implemented for manjushaka code
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
     paux(6) = paux(6) + csign * nele ** 2
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate <N^1>
!-------------------------------------------------------------------------
     paux(5) = paux(5) + csign * nele
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate spin magnetization: < Sz >
!-------------------------------------------------------------------------
     paux(4) = paux(4) + csign * sz
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate kinetic energy: ekin
! equation : -T < k >
!-------------------------------------------------------------------------
     paux(3) = paux(3) - csign * real(ckink * norbs) / beta
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate potential energy: epot
! it is < H_{loc} > in fact, not equal to the definition in azalea project
! equation : \sum_m P_m E_m
! note: here U denotes as energy zero point
!-------------------------------------------------------------------------
     do i=1,ncfgs
         paux(2) = paux(2) + csign * cprob(i) * ( eigs(i) + Uc )
     enddo ! over i={1,ncfgs} loop
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate total energy: etot
! equation : E_{tot} = < H_{loc} > - T < k > + \mu N
!-------------------------------------------------------------------------
     paux(1) = paux(2) + paux(3) + mune * nele
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

     use context, only : c_sgn
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

! get imaginary time value for creation operators
         do is=1,rank(flvr)
             taus = time_s( index_s(is, flvr), flvr )

! get imaginary time value for annihilation operators
             do ie=1,rank(flvr)
                 taue = time_e( index_e(ie, flvr), flvr )

! evaluate dtau
                 dtau = taue - taus

! get matrix element from mmat, pay special attention to the sign of dtau
                 maux = mmat(ie, is, flvr) * sign(one, dtau) * c_sgn

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
     implicit none

     CONTINUE

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

     use context, only : c_sgn
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
             grnf(ifrq, flvr, flvr) = grnf(ifrq, flvr, flvr) + c_sgn * gmat(ifrq, flvr, flvr)
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

     use context, only : c_sgn
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
         knop(i) = knop(i) + rank(i) * 2.0_dp * c_sgn
     enddo ! over i={1,norbs} loop

     do j=1,norbs
         do i=1,norbs
             kmat(i,j) = kmat(i,j) + rank(i) * rank(j) * 4.0_dp * c_sgn
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

     use context, only : c_sgn
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
     lnop = lnop + kl * c_sgn
     rnop = rnop + kr * c_sgn

! add contribution to < k_l k_r >
     do flvr=1,norbs
         do i=1,norbs
             lrmm(i,flvr) = lrmm(i,flvr) + kl(i) * kr(flvr) * c_sgn
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
     implicit none

     CONTINUE

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
     implicit none

     CONTINUE

     return
  end subroutine ctqmc_record_sp_t

!!
!! @sub ctqmc_record_sp_w
!!
!! record the spin-spin correlation function in matsubara frequency axis
!!
  subroutine ctqmc_record_sp_w()
     implicit none

     CONTINUE

     return
  end subroutine ctqmc_record_sp_w

!!
!! @sub ctqmc_record_ch_t
!!
!! record the charge-charge correlation function in imaginary time axis
!!
  subroutine ctqmc_record_ch_t()
     implicit none

     CONTINUE

     return
  end subroutine ctqmc_record_ch_t

!!
!! @sub ctqmc_record_ch_w
!!
!! record the charge-charge correlation function in matsubara frequency axis
!!
  subroutine ctqmc_record_ch_w()
     implicit none

     CONTINUE

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
