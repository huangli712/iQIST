!!!-----------------------------------------------------------------------
!!! project : begonia
!!! program : ctqmc_record_gtau
!!!           ctqmc_record_grnf
!!!           ctqmc_record_hist
!!!           ctqmc_record_prob
!!!           ctqmc_record_nmat <<<---
!!!           ctqmc_reduce_gtau
!!!           ctqmc_reduce_grnf
!!!           ctqmc_reduce_hist
!!!           ctqmc_reduce_prob
!!!           ctqmc_reduce_nmat <<<---
!!!           ctqmc_symm_nmat
!!!           ctqmc_symm_gtau
!!!           ctqmc_symm_grnf
!!!           ctqmc_smth_sigf   <<<---
!!!           ctqmc_make_gtau   <<<---
!!!           ctqmc_make_hub1   <<<---
!!! source  : ctqmc_record.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 09/16/2009 by li huang
!!!           09/29/2010 by li huang
!!!           11/08/2014 by li huang
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
     use constants, only : dp, zero, one, two

     use control, only : norbs
     use control, only : ntime
     use control, only : beta
     use context, only : csign
     use context, only : index_s, index_e, time_s, time_e
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

! index for imaginary time \tau
     integer  :: curr

! used to store the element of mmat matrix
     real(dp) :: maux

! imaginary time for start and end points
     real(dp) :: taus
     real(dp) :: taue

! length betweem taus and taue
     real(dp) :: dtau

! interval for imaginary time slice
     real(dp) :: step

! evaluate step at first
     step = real(ntime - 1) / beta

     CTQMC_FLAVOR_LOOP: do flvr=1,norbs

! get imaginary time value for create and destroy operators
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
  end subroutine ctqmc_record_gtau

!!>>> ctqmc_record_grnf: record the impurity green's function in matsubara
!!>>> frequency space
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

! note: only the first nfreq points of grnf are modified
     do flvr=1,norbs
         do ifrq=1,nfreq
             grnf(ifrq, flvr, flvr) = grnf(ifrq, flvr, flvr) + csign * gmat(ifrq, flvr, flvr)
         enddo ! over ifrq={1,nfreq} loop
     enddo ! over flvr={1,norbs} loop

     return
  end subroutine ctqmc_record_grnf

!!>>> ctqmc_record_hist: record the histogram of perturbation expansion series
  subroutine ctqmc_record_hist()
     use constants, only : one

     use control, only : mkink
     use context, only : ckink, csign, caves
     use context, only : hist

     implicit none

! record current sign as a byproduct
     caves = caves + csign

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
     use control, only : ncfgs
     use context, only : csign, matrix_ptrace
     use context, only : prob
     use context, only : diag

     implicit none

! local variables
! loop index
     integer :: i

     do i=1,ncfgs
         prob(i) = prob(i) + csign * diag(i,2) / matrix_ptrace
     enddo ! over i={1,ncfgs} loop

     return
  end subroutine ctqmc_record_prob

!!>>> ctqmc_record_nmat: record the occupation matrix, double occupation
!!>>> matrix, and auxiliary physical observables simulataneously
  subroutine ctqmc_record_nmat()
     use constants, only : dp, zero
     use sparse, only : sp_csr_cp_elm, sp_csr_mm_csr

     use control, only : nband, norbs, ncfgs, nzero
     use control, only : U, mune, beta
     use context, only : ckink, csign, matrix_ptrace
     use context, only : paux, nmat, nnmat
     use context, only : diag, eigs
     use context, only : T_spmat, ctqmc_new_spmat, ctqmc_del_spmat
     use context, only : spm_s, spm_n, spm_m

     implicit none

! local variables
! loop index
     integer  :: i
     integer  :: j

! loop index for flavor channel
     integer  :: flvr

! dummy variables
     real(dp) :: raux1
     real(dp) :: raux2

! dummy array, denote as current occupation number
     real(dp) :: nvec(norbs)

! current probability for eigenstates
     real(dp) :: cprob(ncfgs)

! dummy sparse matrix, used to calculate nmat and nnmat
     type (T_spmat) :: spm_t

! allocate memory for the dummy sparse matrix
     call ctqmc_new_spmat(spm_t)

! evaluate cprob at first, it is current atomic probability
     do i=1,ncfgs
         cprob(i) = diag(i,2) / matrix_ptrace
     enddo ! over i={1,ncfgs} loop

! evaluate raux2, it is Tr ( e^{- \beta H} )
! i think it is equal to matrix_ptrace, to be checked
     raux2 = zero
     do i=1,ncfgs
         raux2 = raux2 + sp_csr_cp_elm( i, i, &
                                ncfgs, nzero, &
                                 spm_s(2)%vv, &
                                 spm_s(2)%jv, &
                                 spm_s(2)%iv )
     enddo ! over i={1,ncfgs} loop

! check validity of raux2
!<     if ( abs(raux2) < epss ) then
!<         call s_print_exception('ctqmc_record_nmat()','Z trace is too small')
!<     endif ! back if ( abs(raux2) < epss ) block

! evaluate occupation matrix: < n_i >
! equation : Tr ( e^{- \beta H} c^{\dag}_i c_i ) / Tr ( e^{- \beta H} )
!-------------------------------------------------------------------------
     do flvr=1,norbs
         raux1 = zero
         call sp_csr_mm_csr(        ncfgs, ncfgs, ncfgs, nzero, &
                         spm_s(2)%vv, spm_s(2)%jv, spm_s(2)%iv, &
                spm_n(flvr)%vv, spm_n(flvr)%jv, spm_n(flvr)%iv, &
                                  spm_t%vv, spm_t%jv, spm_t%iv )
         do i=1,ncfgs
             raux1 = raux1 + sp_csr_cp_elm( i, i, ncfgs, nzero, &
                                  spm_t%vv, spm_t%jv, spm_t%iv )
         enddo ! over i={1,ncfgs} loop
         nvec(flvr) = raux1 / raux2
     enddo ! over flvr={1,norbs} loop

! update nmat
     nmat = nmat + csign * nvec
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate double occupation matrix: < n_i n_j >
! equation : Tr ( e^{- \beta H} c^{\dag}_i c_i c^{\dag}_j c_j ) / Tr ( e^{- \beta H} )
!-------------------------------------------------------------------------
     do flvr=1,norbs-1
         do i=flvr+1,norbs
             raux1 = zero
             call sp_csr_mm_csr(        ncfgs, ncfgs, ncfgs, nzero, &
                             spm_s(2)%vv, spm_s(2)%jv, spm_s(2)%iv, &
              spm_m(flvr,i)%vv, spm_m(flvr,i)%jv, spm_m(flvr,i)%iv, &
                                      spm_t%vv, spm_t%jv, spm_t%iv )
             do j=1,ncfgs
                 raux1 = raux1 + sp_csr_cp_elm( j, j, ncfgs, nzero, &
                                      spm_t%vv, spm_t%jv, spm_t%iv )
             enddo ! over j={1,ncfgs} loop
             nnmat(flvr,i) = nnmat(flvr,i) + csign * raux1 / raux2

             raux1 = zero
             call sp_csr_mm_csr(        ncfgs, ncfgs, ncfgs, nzero, &
                             spm_s(2)%vv, spm_s(2)%jv, spm_s(2)%iv, &
              spm_m(i,flvr)%vv, spm_m(i,flvr)%jv, spm_m(i,flvr)%iv, &
                                      spm_t%vv, spm_t%jv, spm_t%iv )
             do j=1,ncfgs
                 raux1 = raux1 + sp_csr_cp_elm( j, j, ncfgs, nzero, &
                                      spm_t%vv, spm_t%jv, spm_t%iv )
             enddo ! over j={1,ncfgs} loop
             nnmat(i,flvr) = nnmat(i,flvr) + csign * raux1 / raux2
         enddo ! over i={flvr+1,norbs} loop
     enddo ! over flvr={1,norbs-1} loop
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate <N^2>
!-------------------------------------------------------------------------
     paux(6) = paux(6) + csign * ( sum(nvec) )**2
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate <N^1>
!-------------------------------------------------------------------------
     paux(5) = paux(5) + csign * sum(nvec)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate spin magnetization: < Sz >
!-------------------------------------------------------------------------
     do flvr=1,nband
         paux(4) = paux(4) + csign * ( nvec(flvr) - nvec(flvr+nband) )
     enddo ! over flvr={1,nband} loop
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
         paux(2) = paux(2) + csign * cprob(i) * ( eigs(i) + U )
     enddo ! over i={1,ncfgs} loop
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate total energy: etot
! equation : E_{tot} = < H_{loc} > - T < k > + \mu N
!-------------------------------------------------------------------------
     paux(1) = paux(2) + paux(3) + mune * sum(nvec)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! deallocate memory for the dummy sparse matrix
     call ctqmc_del_spmat(spm_t)

     return
  end subroutine ctqmc_record_nmat

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
  subroutine ctqmc_reduce_hist(hist_mpi)
     use constants, only : dp, zero
     use mmpi, only : mp_allreduce, mp_barrier

     use control, only : mkink
     use control, only : nprocs
     use context, only : hist

     implicit none

! external arguments
! histogram for perturbation expansion series
     real(dp), intent(out) :: hist_mpi(mkink)

! initialize hist_mpi
     hist_mpi = zero

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

!!>>> ctqmc_make_gtau: build imaginary green's function using normal
!!>>> representation
  subroutine ctqmc_make_gtau(gtau, gaux)
     use constants, only : dp, zero

     use control, only : norbs
     use control, only : ntime
     use control, only : beta

     implicit none

! external arguments
! impurity green's function
     real(dp), intent(in)  :: gtau(ntime,norbs,norbs)

! calculated impurity green's function
     real(dp), intent(out) :: gaux(ntime,norbs,norbs)

! local variables
! loop index
     integer  :: i
     integer  :: j

! dummy variables
     real(dp) :: raux

! initialize gaux
     gaux = zero

     raux = real(ntime) / (beta * beta)
     do i=1,norbs
         do j=1,ntime
             gaux(j,i,i) = gtau(j,i,i) * raux
         enddo ! over j={1,ntime} loop
     enddo ! over i={1,norbs} loop

     return
  end subroutine ctqmc_make_gtau

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
     use constants, only : dp, zero, one, czi, czero, eps6

     use control, only : norbs, ncfgs, nzero
     use control, only : mfreq
     use control, only : nfreq
     use control, only : mune
     use control, only : myid, master
     use context, only : rmesh
     use context, only : prob
     use context, only : eimp, eigs
     use context, only : op_d
     use context, only : grnf
     use context, only : hybf
     use context, only : sig2

     implicit none

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k
     integer  :: m

! dummy integer variables
     integer  :: start

! dummy real variables used to build F matrix
     real(dp) :: value

! dummy real variables, used to interpolate self-energy function
     real(dp) :: ob, oe
     real(dp) :: d0, d1

! dummy complex variables, used to interpolate self-energy function
     complex(dp) :: cb, ce
     complex(dp) :: sinf

! F matrix, < alpha | f_{n} | beta >
     integer  :: fcounter(norbs)
     integer  :: fa(nzero,norbs)
     integer  :: fb(nzero,norbs)

     real(dp) :: fv(nzero,norbs)

! dummy imurity green's function: G^{-1}
     complex(dp) :: gaux(norbs,norbs)

! atomic green's function and self-energy function in Hubbard-I approximation
     complex(dp) :: ghub(mfreq,norbs)
     complex(dp) :: shub(mfreq,norbs)

! build F matrix < alpha | f_{n} | beta >
! note 1: to save the memory and accelerate the computation, we only store
! the non-zero element of F matrix
! note 2: it is crucial to check whether the number of non-zero elements
! exceed limit (nzero)
     fcounter = 0
     alpha_loop: do i=1,ncfgs
         beta_loop: do j=1,ncfgs
             orbital_loop: do m=1,norbs
                 value = op_d(i,j,m)
                 if ( abs(value - zero) > eps6 ) then
                     fcounter(m) = fcounter(m) + 1
                     if ( fcounter(m) > nzero ) then
                         call s_print_error('ctqmc_make_hub1','non-zero elements exceed limit')
                     endif ! back if ( fcounter(m) > nzero ) block
                     fa(fcounter(m),m) = i
                     fb(fcounter(m),m) = j
                     fv(fcounter(m),m) = value
                 endif ! back if ( abs(value - zero) > eps6 ) block
             enddo orbital_loop ! over m={1,norbs} loop
         enddo beta_loop ! over j={1,ncfgs} loop
     enddo alpha_loop ! over i={1,ncfgs} loop

! calculate atomic green's function using Hubbard-I approximation
     do i=1,norbs
         do k=1,mfreq
             ghub(k,i) = czero
             do m=1,fcounter(i)
                 ob = fv(m,i) * fv(m,i) * ( prob(fa(m,i)) + prob(fb(m,i)) )
                 cb = czi * rmesh(k) + eigs(fa(m,i)) - eigs(fb(m,i))
                 ghub(k,i) = ghub(k,i) + ob / cb
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
         enddo ! over i={1,norbs} loop
         call s_inv_z(norbs, gaux)
         grnf(k,:,:) = gaux
     enddo ! over k={1,mfreq} loop

     return
  end subroutine ctqmc_make_hub1
