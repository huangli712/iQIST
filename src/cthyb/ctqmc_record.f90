!!!-----------------------------------------------------------------------
!!! project : manjushaka
!!! program : ctqmc_record_hist
!!!           ctqmc_record_prob
!!!           ctqmc_record_paux
!!!           ctqmc_record_nmat <<<---
!!!           ctqmc_record_gtau
!!!           ctqmc_record_ftau
!!!           ctqmc_record_grnf <<<---
!!!           ctqmc_record_kmat
!!!           ctqmc_record_lrmm
!!!           ctqmc_record_szpw <<<---
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
!!!           ctqmc_reduce_lrmm
!!!           ctqmc_reduce_szpw <<<---
!!!           ctqmc_reduce_twop
!!!           ctqmc_reduce_pair <<<---
!!! source  : ctqmc_record.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           05/18/2017 by li huang (last modified)
!!! purpose : measure and collect physical observables produced by the
!!!           hybridization expansion version continuous time quantum
!!!           Monte Carlo (CTQMC) quantum impurity solver.
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

!!>>> ctqmc_record_prob: record the probability of atomic states
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

  subroutine ctqmc_record_paux()
  end subroutine ctqmc_record_paux

!!>>> ctqmc_record_nmat: record the occupation matrix, double occupation
!!>>> matrix, and auxiliary physical observables simulataneously
  subroutine ctqmc_record_nmat()
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

     use context, only : csign
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

! distance betweem taus and taue
     real(dp) :: dtau
     real(dp) :: daux

! interval for imaginary time slice
     real(dp) :: step

! evaluate step at first
     if ( isort == 1 ) step = real(ntime - 1) / beta
     if ( isort == 2 ) step = real(legrd - 1) / two

     FLVR_CYCLE: do flvr=1,norbs

! get imaginary time value for creation and annihilation operators
         do is=1,rank(flvr)
             taus = time_s( index_s(is, flvr), flvr )

             do ie=1,rank(flvr)
                 taue = time_e( index_e(ie, flvr), flvr )

! evaluate dtau
                 dtau = taue - taus

! get matrix element from mmat, pay special attention to the sign of dtau
                 maux = mmat(ie, is, flvr) * sign(one, dtau) * csign

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
                     LEG_CYCLE: do fleg=1,lemax
                         dtau = sqrt(two * fleg - 1) * rep_l(curr,fleg)
                         gtau(fleg, flvr, flvr) = gtau(fleg, flvr, flvr) - maux * dtau
                     enddo LEG_CYCLE ! over fleg={1,lemax} loop

                 endif LEG_BLOCK ! back if ( isort == 2 ) block
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
  end subroutine ctqmc_record_ftau

!!
!! @sub ctqmc_record_grnf
!!
!! record the impurity green's function in matsubara frequency space
!!
  subroutine ctqmc_record_grnf()
     use control, only : norbs
     use control, only : nfreq

     use context, only : csign
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
             grnf(ifrq, flvr, flvr) = grnf(ifrq, flvr, flvr) + csign * gmat(ifrq, flvr, flvr)
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
!! record the kinetic energy fluctuation < k^2 > - < k >^2
!!
  subroutine ctqmc_record_kmat()
     use constants, only : dp

     use control, only : isobs
     use control, only : norbs

     use context, only : csign
     use context, only : knop, kmat
     use context, only : rank

     implicit none

! local variables
! loop index for flavor channel
     integer :: i
     integer :: j

! check whether there is conflict
     call s_assert( btest(isobs, 1) )

! since rank means the number of operator pairs, so we have to multiply
! it with two
     do i=1,norbs
         knop(i) = knop(i) + rank(i) * 2.0_dp * csign
     enddo ! over i={1,norbs} loop

     do j=1,norbs
         do i=1,norbs
             kmat(i,j) = kmat(i,j) + rank(i) * rank(j) * 4.0_dp * csign
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
     use constants, only : dp, zero, one, two

     use control, only : isobs
     use control, only : norbs
     use control, only : beta

     use context, only : csign
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
     call s_assert( btest(isobs, 2) )

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
     lnop = lnop + kl * csign
     rnop = rnop + kr * csign

! add contribution to < k_l k_r >
     do flvr=1,norbs
         do i=1,norbs
             lrmm(i,flvr) = lrmm(i,flvr) + kl(i) * kr(flvr) * csign
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
  end subroutine ctqmc_record_szpw

!!========================================================================
!!>>> measure physical observables 4                                   <<<
!!========================================================================

!!========================================================================
!!>>> measure physical observables 5                                   <<<
!!========================================================================

!!
!! @sub ctqmc_record_twop
!!
!! record the two-particle green's function. here improved estimator is
!! used to improve the accuracy
!!
  subroutine ctqmc_record_twop()
     use constants, only : dp, czero

     use control, only : isvrt
     use control, only : norbs
     use control, only : nffrq, nbfrq
     use control, only : beta

     use context, only : csign
     use context, only : g2pw
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

! dummy complex(dp) variables, used to calculate the g2pw
     complex(dp) :: cmeas

! dummy complex(dp) arrays, used to store the intermediate results
     complex(dp), allocatable :: g2aux(:,:,:)
     complex(dp), allocatable :: caux1(:,:)
     complex(dp), allocatable :: caux2(:,:)

! check whether there is conflict
     call s_assert( btest(isvrt, 1) )

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
     FLVR_CYCLE: do flvr=1,norbs
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

     enddo FLVR_CYCLE ! over flvr={1,norbs} loop
!$OMP END DO

! calculate g2pw
!$OMP DO PRIVATE (f1, f2, cmeas, wbn, w4n, w3n, w2n, w1n)
     ORB1_CYCLE: do f1=1,norbs
         ORB2_CYCLE: do f2=1,f1

             WB_CYCLE: do wbn=1,nbfrq

                 WF1_CYCLE: do w2n=1,nffrq
                     WF2_CYCLE: do w3n=1,nffrq
                         w1n = w2n + wbn - 1; w4n = w3n + wbn - 1

                         cmeas = g2aux(w1n,w2n,f1) * g2aux(w3n,w4n,f2)
                         if ( f1 == f2 ) then
                             cmeas = cmeas - g2aux(w1n,w4n,f1) * g2aux(w3n,w2n,f1)
                         endif ! back if ( f1 == f2 ) block
                         g2pw(w3n,w2n,wbn,f2,f1) = g2pw(w3n,w2n,wbn,f2,f1) + cmeas * csign / beta
                     enddo WF2_CYCLE ! over w3n={1,nffrq} loop
                 enddo WF1_CYCLE ! over w2n={1,nffrq} loop

             enddo WB_CYCLE ! over wbn={1,nbfrq} loop

         enddo ORB2_CYCLE ! over f2={1,f1} loop
     enddo ORB1_CYCLE ! over f1={1,norbs} loop
!$OMP END DO
!$OMP END PARALLEL

! deallocate memory
     deallocate( g2aux )
     deallocate( caux1 )
     deallocate( caux2 )

     return
  end subroutine ctqmc_record_twop

!!
!! @sub ctqmc_record_pair
!!
!! record the particle-particle pairing susceptibility
!!
  subroutine ctqmc_record_pair()
     use constants, only : dp, czero

     use control, only : isvrt
     use control, only : norbs
     use control, only : nffrq, nbfrq
     use control, only : beta

     use context, only : csign
     use context, only : p2pw
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

! dummy complex(dp) variables, used to calculate the p2pw
     complex(dp) :: cmeas

! dummy complex(dp) arrays, used to store the intermediate results
     complex(dp), allocatable :: g2aux(:,:,:)
     complex(dp), allocatable :: caux1(:,:)
     complex(dp), allocatable :: caux2(:,:)

! check whether there is conflict
     call s_assert( btest(isvrt, 2) )

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
     FLVR_CYCLE: do flvr=1,norbs
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

     enddo FLVR_CYCLE ! over flvr={1,norbs} loop
!$OMP END DO

! calculate p2pw
!$OMP DO PRIVATE (f1, f2, cmeas, wbn, w4n, w3n, w2n, w1n)
     ORB1_CYCLE: do f1=1,norbs
         ORB2_CYCLE: do f2=1,f1

             WB_CYCLE: do wbn=1,nbfrq

                 WF1_CYCLE: do w2n=1,nffrq
                     WF2_CYCLE: do w3n=1,nffrq
                         w1n = w2n + wbn - 1; w4n = w3n + wbn - 1

                         cmeas = czero
                         if ( f1 /= f2 ) then
                             cmeas = cmeas + g2aux(w1n,w4n,f1) * g2aux(nffrq-w2n+1,nffrq-w3n+1,f2)
                         endif ! back if ( f1 == f2 ) block
                         p2pw(w3n,w2n,wbn,f2,f1) = p2pw(w3n,w2n,wbn,f2,f1) + cmeas * csign / beta
                     enddo WF2_CYCLE ! over w3n={1,nffrq} loop
                 enddo WF1_CYCLE ! over w2n={1,nffrq} loop

             enddo WB_CYCLE ! over wbn={1,nbfrq} loop

         enddo ORB2_CYCLE ! over f2={1,f1} loop
     enddo ORB1_CYCLE ! over f1={1,norbs} loop
!$OMP END DO
!$OMP END PARALLEL

! deallocate memory
     deallocate( g2aux )
     deallocate( caux1 )
     deallocate( caux2 )

     return
  end subroutine ctqmc_record_pair

!!========================================================================
!!>>> reduce physical observables 1                                    <<<
!!========================================================================

!!
!! @sub ctqmc_reduce_hist
!!
!! reduce the hist from all children processes
!!
  subroutine ctqmc_reduce_hist(hist_mpi, hist_err)
     use constants, only : dp, zero

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
     use constants, only : dp, zero

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
     use constants, only : dp, zero

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
     use constants, only : dp, zero

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
     use constants, only : dp, zero

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
     use constants, only : dp, zero

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
     use constants, only : dp, zero, czero, czi

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

!!========================================================================
!!>>> reduce physical observables 3                                    <<<
!!========================================================================

!!
!! @sub ctqmc_reduce_kmat
!!
!! reduce the knop and kmat from all children processes
!!
  subroutine ctqmc_reduce_kmat(knop_mpi, kmat_mpi, knop_err, kmat_err)
     use constants, only : dp, zero

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
     use constants, only : dp, zero

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
     use constants, only : dp, zero

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

!!========================================================================
!!>>> reduce physical observables 5                                    <<<
!!========================================================================

!!
!! @sub ctqmc_reduce_twop
!!
!! reduce the g2pw and h2pw from all children processes
!!
  subroutine ctqmc_reduce_twop(g2pw_mpi, h2pw_mpi, g2pw_err, h2pw_err)
     use constants, only : dp, zero, czero, czi

     use mmpi, only : mp_allreduce
     use mmpi, only : mp_barrier

     use control, only : isvrt
     use control, only : norbs
     use control, only : nffrq, nbfrq
     use control, only : nprocs

     use context, only : g2pw
     use context, only : h2pw

     implicit none

! external arguments
! two-particle green's function
     complex(dp), intent(out) :: g2pw_mpi(nffrq,nffrq,nbfrq,norbs,norbs)
     complex(dp), intent(out) :: g2pw_err(nffrq,nffrq,nbfrq,norbs,norbs)

! irreducible vertex function
     complex(dp), intent(out) :: h2pw_mpi(nffrq,nffrq,nbfrq,norbs,norbs)
     complex(dp), intent(out) :: h2pw_err(nffrq,nffrq,nbfrq,norbs,norbs)

! local variables
! used to store the real and imaginary parts of green's function
     real(dp), allocatable :: g_re_err(:,:,:,:,:)
     real(dp), allocatable :: g_im_err(:,:,:,:,:)

! used to store the real and imaginary parts of vertex function
     real(dp), allocatable :: h_re_err(:,:,:,:,:)
     real(dp), allocatable :: h_im_err(:,:,:,:,:)

! check whether this observable has been measured
     if ( .not. btest(isvrt, 1) ) RETURN

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

! initialize g2pw_mpi and g2pw_err
     g2pw_mpi = czero
     g2pw_err = czero

! initialize h2pw_mpi and h2pw_err
     h2pw_mpi = czero
     h2pw_err = czero

! build g2pw_mpi and h2pw_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(g2pw, g2pw_mpi)
     call mp_allreduce(h2pw, h2pw_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     g2pw_mpi = g2pw
     h2pw_mpi = h2pw

# endif /* MPI */

! calculate the average
     g2pw_mpi = g2pw_mpi / real(nprocs)
     h2pw_mpi = h2pw_mpi / real(nprocs)

! build g2pw_err and h2pw_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(( real(g2pw - g2pw_mpi))**2, g_re_err)
     call mp_allreduce((aimag(g2pw - g2pw_mpi))**2, g_im_err)
     call mp_allreduce(( real(h2pw - h2pw_mpi))**2, h_re_err)
     call mp_allreduce((aimag(h2pw - h2pw_mpi))**2, h_im_err)

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

! construct the final g2pw_err and h2pw_err
     g2pw_err = g_re_err + g_im_err * czi
     h2pw_err = h_re_err + h_im_err * czi

! deallocate memory
     deallocate(g_re_err)
     deallocate(g_im_err)
     deallocate(h_re_err)
     deallocate(h_im_err)

     return
  end subroutine ctqmc_reduce_twop

!!
!! @sub ctqmc_reduce_pair
!!
!! reduce the p2pw from all children processes
!!
  subroutine ctqmc_reduce_pair(p2pw_mpi, p2pw_err)
     use constants, only : dp, zero, czero, czi

     use mmpi, only : mp_allreduce
     use mmpi, only : mp_barrier

     use control, only : isvrt
     use control, only : norbs
     use control, only : nffrq, nbfrq
     use control, only : nprocs

     use context, only : p2pw

     implicit none

! external arguments
! particle-particle pairing susceptibility
     complex(dp), intent(out) :: p2pw_mpi(nffrq,nffrq,nbfrq,norbs,norbs)
     complex(dp), intent(out) :: p2pw_err(nffrq,nffrq,nbfrq,norbs,norbs)

! local variables
! used to store the real and imaginary parts of pairing susceptibility
     real(dp), allocatable :: p_re_err(:,:,:,:,:)
     real(dp), allocatable :: p_im_err(:,:,:,:,:)

! check whether this observable has been measured
     if ( .not. btest(isvrt, 2) ) RETURN

! allocate memory
     allocate(p_re_err(nffrq,nffrq,nbfrq,norbs,norbs))
     allocate(p_im_err(nffrq,nffrq,nbfrq,norbs,norbs))

! initialize p_re_err and p_im_err
     p_re_err = zero
     p_im_err = zero

! initialize p2pw_mpi and p2pw_err
     p2pw_mpi = czero
     p2pw_err = czero

! build p2pw_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(p2pw, p2pw_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     p2pw_mpi = p2pw

# endif /* MPI */

! calculate the average
     p2pw_mpi = p2pw_mpi / real(nprocs)

! build p2pw_err, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(( real(p2pw - p2pw_mpi))**2, p_re_err)
     call mp_allreduce((aimag(p2pw - p2pw_mpi))**2, p_im_err)

! block until all processes have reached here
     call mp_barrier()

# endif /* MPI */

! calculate standard deviation
     if ( nprocs > 1 ) then
         p_re_err = sqrt( p_re_err / real( nprocs * ( nprocs - 1 ) ) )
         p_im_err = sqrt( p_im_err / real( nprocs * ( nprocs - 1 ) ) )
     endif ! back if ( nprocs > 1 ) block

! construct the final p2pw_err
     p2pw_err = p_re_err + p_im_err * czi

! deallocate memory
     deallocate(p_re_err)
     deallocate(p_im_err)

     return
  end subroutine ctqmc_reduce_pair
