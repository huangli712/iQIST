!-------------------------------------------------------------------------
! project : lavender
! program : ctqmc_record_gtau
!           ctqmc_record_ftau
!           ctqmc_record_grnf
!           ctqmc_record_hist
!           ctqmc_record_nmat
!           ctqmc_record_schi
!           ctqmc_record_ochi
!           ctqmc_record_twop
!           ctqmc_record_vrtx
!           ctqmc_record_prob <<<---
!           ctqmc_reduce_gtau
!           ctqmc_reduce_ftau
!           ctqmc_reduce_grnf
!           ctqmc_reduce_hist
!           ctqmc_reduce_nmat
!           ctqmc_reduce_schi
!           ctqmc_reduce_ochi
!           ctqmc_reduce_twop
!           ctqmc_reduce_vrtx
!           ctqmc_reduce_prob <<<---
!           ctqmc_symm_nmat
!           ctqmc_symm_gtau
!           ctqmc_symm_grnf
!           ctqmc_smth_sigf   <<<---
!           ctqmc_make_gtau
!           ctqmc_make_ftau   <<<---
!           ctqmc_make_hub1
!           ctqmc_make_hub2   <<<---
! source  : ctqmc_record.f90
! type    : subroutine
! author  : li huang (email:huangli712@yahoo.com.cn)
! history : 09/16/2009 by li huang
!           09/18/2009 by li huang
!           09/20/2009 by li huang
!           09/25/2009 by li huang
!           09/27/2009 by li huang
!           10/29/2009 by li huang
!           11/01/2009 by li huang
!           11/03/2009 by li huang
!           11/10/2009 by li huang
!           11/19/2009 by li huang
!           11/30/2009 by li huang
!           12/06/2009 by li huang
!           12/09/2009 by li huang
!           12/18/2009 by li huang
!           12/22/2009 by li huang
!           12/26/2009 by li huang
!           12/29/2009 by li huang
!           01/14/2010 by li huang
!           02/01/2010 by li huang
!           02/24/2010 by li huang
!           02/27/2010 by li huang
!           09/29/2010 by li huang
! purpose : measure, record, and postprocess the key observables produced
!           by the hybridization expansion version continuous time quantum
!           Monte Carlo (CTQMC) quantum impurity solver
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> record the impurity green's function in imaginary time axis
  subroutine ctqmc_record_gtau()
     use constants
     use control
     use context

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

!>>> record impurity green's function using normal representation
  subroutine cat_record_gtau1()
     implicit none

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

!>>> record impurity green's function using legendre polynomial representation
  subroutine cat_record_gtau2()
     implicit none

! evaluate step at first
     step = real(legrd - 1) / two

     CTQMC_FLAVOR_LOOP: do flvr=1,norbs

! get imaginary time value for create and destroy operators
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

!>>> record impurity green's function using chebyshev polynomial representation
  subroutine cat_record_gtau3()
     implicit none

! evaluate step at first
     step = real(chgrd - 1) / two

     CTQMC_FLAVOR_LOOP: do flvr=1,norbs

! get imaginary time value for create and destroy operators
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

!>>> record the auxiliary correlation function in imaginary time axis
  subroutine ctqmc_record_ftau()
     use constants
     use control
     use context

     implicit none

     call ctqmc_print_error('ctqmc_record_ftau', 'this subroutine is not implemented')

     return
  end subroutine ctqmc_record_ftau

!>>> record the impurity green's function in matsubara frequency space
  subroutine ctqmc_record_grnf()
     use constants
     use control
     use context

     implicit none

! local variables
! loop index over matsubara frequencies
     integer :: ifrq

! loop index for flavor channel
     integer :: flvr

! note: only the first nfreq points of grnf are modified
! we enforce csign equal to 1, the influence of sign problem is unclear so far
!<     csign = 1
     do flvr=1,norbs
         do ifrq=1,nfreq
             grnf(ifrq, flvr, flvr) = grnf(ifrq, flvr, flvr) + csign * gmat(ifrq, flvr, flvr)
         enddo ! over ifrq={1,nfreq} loop
     enddo ! over flvr={1,norbs} loop

     return
  end subroutine ctqmc_record_grnf

!>>> record the histogram of perturbation expansion series
  subroutine ctqmc_record_hist()
     use context

     implicit none

! record current sign as a byproduct
     caves = caves + csign

! note: if ckink == 0, we record its count in hist(mkink)
     if ( ckink > 0 ) then
         hist(ckink) = hist(ckink) + 1
     else
         hist(mkink) = hist(mkink) + 1
     endif

     return
  end subroutine ctqmc_record_hist

!>>> record the occupation matrix, double occupation matrix, and auxiliary
! physical observables simulataneously
  subroutine ctqmc_record_nmat()
     use constants
     use control
     use context

     use sparse

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
     integer  :: sop_it(ncfgs+1)
     integer  :: sop_jt(nzero)
     real(dp) :: sop_t(nzero)

! evaluate cprob at first, it is current atomic propability
     do i=1,ncfgs
         cprob(i) = diag(i,2) / matrix_ptrace
     enddo ! over i={1,ncfgs} loop

! evaluate raux2, it is Tr ( e^{- \beta H} )
! i think it is equal to matrix_ptrace, to be checked
     raux2 = zero
     do i=1,ncfgs
         raux2 = raux2 + sparse_csr_cp_elm( i, i, ncfgs, nzero, sop_s(:,2), sop_js(:,2), sop_is(:,2) )
     enddo ! over i={1,ncfgs} loop

! check validity of raux2
!<     if ( abs(raux2) < epss ) then
!<         call ctqmc_print_exception('ctqmc_record_nmat()','Z trace is too small')
!<     endif

! evaluate occupation matrix: < n_i >
! equation : Tr ( e^{- \beta H} c^{\dag}_i c_i ) / Tr ( e^{- \beta H} )
!-------------------------------------------------------------------------
     do flvr=1,norbs
         call sparse_csr_mm_csr(             ncfgs, ncfgs, ncfgs, nzero, &
                             sop_s(:,2),    sop_js(:,2),    sop_is(:,2), &
                          sop_n(:,flvr), sop_jn(:,flvr), sop_in(:,flvr), &
                                  sop_t,         sop_jt,         sop_it )
         raux1 = zero
         do i=1,ncfgs
             raux1 = raux1 + sparse_csr_cp_elm( i, i, ncfgs, nzero, sop_t, sop_jt, sop_it )
         enddo ! over i={1,ncfgs} loop
         nvec(flvr) = raux1 / raux2
     enddo ! over flvr={1,norbs} loop

! update nmat
     nmat = nmat + nvec
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate double occupation matrix: < n_i n_j >
! equation : Tr ( e^{- \beta H} c^{\dag}_i c_i c^{\dag}_j c_j ) / Tr ( e^{- \beta H} )
!-------------------------------------------------------------------------
     do flvr=1,norbs-1
         do j=flvr+1,norbs
             call sparse_csr_mm_csr(         ncfgs, ncfgs, ncfgs, nzero, &
                         sop_s(:,2),      sop_js(:,2),      sop_is(:,2), &
                    sop_m(:,flvr,j), sop_jm(:,flvr,j), sop_im(:,flvr,j), &
                              sop_t,           sop_jt,           sop_it )

             raux1 = zero
             do i=1,ncfgs
                 raux1 = raux1 + sparse_csr_cp_elm( i, i, ncfgs, nzero, sop_t, sop_jt, sop_it )
             enddo ! over i={1,ncfgs} loop
             nnmat(flvr,j) = nnmat(flvr,j) + raux1 / raux2

             call sparse_csr_mm_csr(         ncfgs, ncfgs, ncfgs, nzero, &
                         sop_s(:,2),      sop_js(:,2),      sop_is(:,2), &
                    sop_m(:,j,flvr), sop_jm(:,j,flvr), sop_im(:,j,flvr), &
                              sop_t,           sop_jt,           sop_it )

             raux1 = zero
             do i=1,ncfgs
                 raux1 = raux1 + sparse_csr_cp_elm( i, i, ncfgs, nzero, sop_t, sop_jt, sop_it )
             enddo ! over i={1,ncfgs} loop
             nnmat(j,flvr) = nnmat(j,flvr) + raux1 / raux2
         enddo ! over j={flvr+1,norbs} loop
     enddo ! over flvr={1,norbs-1} loop
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate spin magnetization: < Sz >
!-------------------------------------------------------------------------
     do flvr=1,nband
         paux(4) = paux(4) + ( nvec(flvr) - nvec(flvr+nband) )
     enddo ! over flvr={1,nband} loop
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate kinetic energy: ekin
! equation : -T < k >
!-------------------------------------------------------------------------
     paux(3) = paux(3) - real(ckink * norbs) / beta
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate potential energy: epot
! it is < H_{loc} > in fact, not equal to the definition in azalea project
! equation : \sum_m P_m E_m
! note: here U denotes as energy zero point
!-------------------------------------------------------------------------
     do i=1,ncfgs
         paux(2) = paux(2) + cprob(i) * ( eigs(i) + U )
     enddo ! over i={1,ncfgs} loop
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate total energy: etot
! equation : E_{tot} = < H_{loc} > - T < k > + \mu N
!-------------------------------------------------------------------------
     paux(1) = paux(2) + paux(3) + mune * sum(nvec)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     return
  end subroutine ctqmc_record_nmat

!>>> record the spin-spin correlation function
  subroutine ctqmc_record_schi()
     use constants
     use control
     use context

     implicit none

     call ctqmc_print_error('ctqmc_record_schi', 'this subroutine is not implemented')

     return
  end subroutine ctqmc_record_schi

!>>> record the orbital-orbital correlation function
  subroutine ctqmc_record_ochi()
     use constants
     use control
     use context

     implicit none

     call ctqmc_print_error('ctqmc_record_ochi', 'this subroutine is not implemented')

     return
  end subroutine ctqmc_record_ochi

!>>> record the two-particle green's function
  subroutine ctqmc_record_twop()
     use constants
     use control
     use context

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

!>>> record the vertex function
  subroutine ctqmc_record_vrtx()
     use constants
     use control
     use context

     implicit none

     call ctqmc_print_error('ctqmc_record_vrtx', 'this subroutine is not implemented')

     return
  end subroutine ctqmc_record_vrtx

!>>> record the probability of atomic states
  subroutine ctqmc_record_prob()
     use constants
     use control
     use context

     implicit none

! local variables
! loop index
     integer :: i

     do i=1,ncfgs
         prob(i) = prob(i) + csign * diag(i,2) / matrix_ptrace
     enddo ! over i={1,ncfgs} loop

     return
  end subroutine ctqmc_record_prob

!>>> reduce the gtau from all children processes
  subroutine ctqmc_reduce_gtau(gtau_mpi)
     use constants
     use context

     use mmpi

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

!>>> reduce the ftau from all children processes
  subroutine ctqmc_reduce_ftau(ftau_mpi)
     use constants
     use context

     use mmpi

     implicit none

! external arguments
! auxiliary correlation function
     real(dp), intent(out) :: ftau_mpi(ntime,norbs,norbs)

! initialize ftau_mpi
     ftau_mpi = zero

! build ftau_mpi, collect data from all children processes
# if defined (MPI)

! collect data
     !call mp_allreduce(ftau, ftau_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     !ftau_mpi = ftau

# endif /* MPI */

! calculate the average
     ftau_mpi = ftau_mpi / real(nprocs)

     return
  end subroutine ctqmc_reduce_ftau

!>>> reduce the grnf from all children processes
  subroutine ctqmc_reduce_grnf(grnf_mpi)
     use constants
     use context

     use mmpi

     implicit none

! external arguments
! impurity green's function
     complex(dp), intent(out) :: grnf_mpi(mfreq,norbs,norbs)

! initialize grnf_mpi
     grnf_mpi = zero

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

!>>> reduce the hist from all children processes
! note: since hist_mpi and hist are integer (kind=4) type, it is important
! to avoid data overflow in them
  subroutine ctqmc_reduce_hist(hist_mpi)
     use constants
     use context

     use mmpi

     implicit none

! external arguments
! histogram for perturbation expansion series
     real(dp), intent(out) :: hist_mpi(mkink)

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

!>>> reduce the nmat and nnmat from all children processes
  subroutine ctqmc_reduce_nmat(nmat_mpi, nnmat_mpi)
     use constants
     use context

     use mmpi

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

!>>> reduce the schi and sschi from all children processes
  subroutine ctqmc_reduce_schi(schi_mpi, sschi_mpi)
     use constants
     use context

     use mmpi

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
     !call mp_allreduce(schi, schi_mpi)
     !call mp_allreduce(sschi, sschi_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     !schi_mpi = schi
     !sschi_mpi = sschi

# endif /* MPI */

! calculate the average
     schi_mpi = schi_mpi / real(nprocs)
     sschi_mpi = sschi_mpi / real(nprocs)

     return
  end subroutine ctqmc_reduce_schi

!>>> reduce the ochi and oochi from all children processes
  subroutine ctqmc_reduce_ochi(ochi_mpi, oochi_mpi)
     use constants
     use context

     use mmpi

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
     !call mp_allreduce(ochi, ochi_mpi)
     !call mp_allreduce(oochi, oochi_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     !ochi_mpi = ochi
     !oochi_mpi = oochi

# endif /* MPI */

! calculate the average
     ochi_mpi = ochi_mpi / real(nprocs)
     oochi_mpi = oochi_mpi / real(nprocs)

     return
  end subroutine ctqmc_reduce_ochi

!>>> reduce the g2_re_mpi and g2_im_mpi from all children processes
  subroutine ctqmc_reduce_twop(g2_re_mpi, g2_im_mpi)
     use constants
     use context

     use mmpi

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

!>>> reduce the h2_re_mpi and h2_im_mpi from all children processes
  subroutine ctqmc_reduce_vrtx(h2_re_mpi, h2_im_mpi)
     use constants
     use context

     use mmpi

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
     !call mp_allreduce(h2_re, h2_re_mpi)
     !call mp_allreduce(h2_im, h2_im_mpi)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     !h2_re_mpi = h2_re
     !h2_im_mpi = h2_im

# endif /* MPI */

! calculate the average
     h2_re_mpi = h2_re_mpi / real(nprocs)
     h2_im_mpi = h2_im_mpi / real(nprocs)

     return
  end subroutine ctqmc_reduce_vrtx

!>>> reduce the prob from all children processes
  subroutine ctqmc_reduce_prob(prob_mpi)
     use constants
     use context

     use mmpi

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

!>>> symmetrize the nmat according to symm vector
  subroutine ctqmc_symm_nmat(symm, nmat)
     use constants
     use control

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

!>>> symmetrize the gtau according to symm vector
! only the diagonal elements are taken into considerations
  subroutine ctqmc_symm_gtau(symm, gtau)
     use constants
     use control

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

!>>> symmetrize the grnf according to symm vector
! only the diagonal elements are taken into considerations
  subroutine ctqmc_symm_grnf(symm, grnf)
     use constants
     use control

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

!>>> smooth impurity self-energy function in low frequency region
  subroutine ctqmc_smth_sigf(sigf)
     use constants
     use control

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

!>>> build imaginary green's function using orthogonal polynomial representation
  subroutine ctqmc_make_gtau(tmesh, gtau, gaux)
     use constants
     use control
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

!>>> build the integral kernel function
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

!>>> build impurity green's function using normal representation
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

!>>> build impurity green's function using legendre polynomial representation
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

!>>> build impurity green's function using chebyshev polynomial representation
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

!>>> build auxiliary correlation function using orthogonal polynomial representation
  subroutine ctqmc_make_ftau(tmesh, ftau, faux)
     use constants
     use control
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

!>>> build atomic green's function and self-energy function using improved
! Hubbard-I approximation, and then make interpolation for self-energy
! function between low frequency QMC data and high frequency Hubbard-I
! approximation data, the full impurity green's function can be obtained by
! using dyson's equation finally
  subroutine ctqmc_make_hub1()
     use constants
     use control
     use context

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
     complex(dp) :: caux

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
                         call ctqmc_print_error('ctqmc_make_hub1','non-zero elements exceed limit')
                     endif
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
             caux = czero
             do m=1,fcounter(i)
                 ob = fv(m,i) * fv(m,i) * ( prob(fa(m,i)) + prob(fb(m,i)) )
                 cb = czi * rmesh(k) + eigs(fa(m,i)) - eigs(fb(m,i))
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
     endif

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

!>>> build atomic green's function and self-energy function using improved
! Hubbard-I approximation, and then make forward fourier transformation
! for impurity green's function and auxiliary correlation function. then
! the final self-energy function is obtained by analytical formula.
  subroutine ctqmc_make_hub2()
     use constants
     use control
     use context

     implicit none

     call ctqmc_print_error('ctqmc_make_hub2', 'this subroutine is not implemented')

     return
  end subroutine ctqmc_make_hub2
