!-------------------------------------------------------------------------
! project : begonia
! program : ctqmc_record_gtau
!           ctqmc_record_grnf
!           ctqmc_record_hist
!           ctqmc_record_nmat
!           ctqmc_record_prob <<<---
!           ctqmc_reduce_gtau
!           ctqmc_reduce_grnf
!           ctqmc_reduce_hist
!           ctqmc_reduce_nmat
!           ctqmc_reduce_prob <<<---
!           ctqmc_symm_nmat
!           ctqmc_symm_gtau
!           ctqmc_symm_grnf
!           ctqmc_smth_sigf   <<<---
!           ctqmc_make_hub1   <<<---
! source  : ctqmc_record.f90
! type    : subroutine
! author  : li huang (email:huangli712@yahoo.com.cn)
!           yilin wang (email: qhwyl2006@126.com)
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
!           07/19/2014 by yilin wang
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

! record gtau, we normalize gtau in ctqmc_dump_gtau() subroutine
                 gtau(curr, flvr, flvr) = gtau(curr, flvr, flvr) - maux

             enddo ! over ie={1,rank(flvr)} loop
         enddo ! over is={1,rank(flvr)} loop

     enddo CTQMC_FLAVOR_LOOP ! over flvr={1,norbs} loop

     return
  end subroutine ctqmc_record_gtau

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

     use m_sector

     implicit none

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k

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
     real(dp) :: tmp_mat(max_dim_sect, max_dim_sect)

! evaluate cprob at first, it is current atomic propability
     do i=1,ncfgs
         cprob(i) = ddmat(i,2) / matrix_ptrace
     enddo ! over i={1,ncfgs} loop

! evaluate raux2, it is Tr ( e^{- \beta H} )
! i think it is equal to matrix_ptrace, to be checked
     raux2 = zero
     do i=1, nsectors
         do j=1, sectors(i)%ndim
             raux2 = raux2 + sectors(i)%final_product(j, j, 2)
         enddo
     enddo

! check validity of raux2
!<     if ( abs(raux2) < epss ) then
!<         call ctqmc_print_exception('ctqmc_record_nmat()','Z trace is too small')
!<     endif

! evaluate occupation matrix: < n_i >
! equation : Tr ( e^{- \beta H} c^{\dag}_i c_i ) / Tr ( e^{- \beta H} )
!-------------------------------------------------------------------------
     tmp_mat = zero
     do flvr=1, norbs
         raux1 = zero
         do i=1, nsectors
             call ctqmc_dmat_gemm( sectors(i)%ndim, sectors(i)%ndim, sectors(i)%ndim, &
                          sectors(i)%final_product(:,:,2), sectors(i)%occu(:,:,flvr),& 
                          tmp_mat(1:sectors(i)%ndim, 1:sectors(i)%ndim) )
             do j=1, sectors(i)%ndim
                 raux1 = raux1 + tmp_mat(j,j)    
             enddo
         enddo 
         nvec(flvr) = raux1 / raux2
     enddo

! update nmat
     nmat = nmat + nvec
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! evaluate double occupation matrix: < n_i n_j >
! equation : Tr ( e^{- \beta H} c^{\dag}_i c_i c^{\dag}_j c_j ) / Tr ( e^{- \beta H} )
!-------------------------------------------------------------------------
     do flvr=1, norbs-1
         do i=flvr+1, norbs
             raux1 = zero
             do j=1, nsectors
                 call ctqmc_dmat_gemm( sectors(j)%ndim, sectors(j)%ndim, sectors(j)%ndim, &
                              sectors(j)%final_product(:,:,2), sectors(j)%double_occu(:,:,flvr,i),& 
                              tmp_mat(1:sectors(j)%ndim, 1:sectors(j)%ndim) )
                 do k=1, sectors(j)%ndim
                     raux1 = raux1 + tmp_mat(k,k)    
                 enddo
             enddo 
             nnmat(flvr,i) = nnmat(flvr,i) + raux1 / raux2

             raux1 = zero
             do j=1, nsectors
                 call ctqmc_dmat_gemm( sectors(j)%ndim, sectors(j)%ndim, sectors(j)%ndim, &
                              sectors(j)%final_product(:,:,2), sectors(j)%double_occu(:,:,i,flvr),& 
                              tmp_mat(1:sectors(j)%ndim, 1:sectors(j)%ndim) )
                 do k=1, sectors(j)%ndim
                     raux1 = raux1 + tmp_mat(k,k)    
                 enddo
             enddo 
             nnmat(i,flvr) = nnmat(i,flvr) + raux1 / raux2
         enddo
     enddo
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
         prob(i) = prob(i) + csign * ddmat(i,2) / matrix_ptrace
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

!>>> build atomic green's function and self-energy function using improved
! Hubbard-I approximation, and then make interpolation for self-energy
! function between low frequency QMC data and high frequency Hubbard-I
! approximation data, the full impurity green's function can be obtained by
! using dyson's equation finally
  subroutine ctqmc_make_hub1()
     use constants
     use control
     use context

     use m_sector

     implicit none

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k
     integer  :: l
     integer  :: m
     integer  :: kk

! start index of sectors
     integer  :: indx1
     integer  :: indx2

! dummy integer variables
     integer  :: start

! dummy real variables, used to interpolate self-energy function
     real(dp) :: ob, oe
     real(dp) :: d0, d1

! dummy complex variables, used to interpolate self-energy function
     complex(dp) :: cb, ce
     complex(dp) :: sinf
     complex(dp) :: caux

! dummy imurity green's function: G^{-1}
     complex(dp) :: gaux(norbs,norbs)

! atomic green's function and self-energy function in Hubbard-I approximation
     complex(dp) :: ghub(mfreq,norbs)
     complex(dp) :: shub(mfreq,norbs)

! calculate atomic green's function using Hubbard-I approximation
     do i=1,norbs
         do j=1, mfreq
             caux = czero
             do k=1, nsectors
                 kk = sectors(k)%next_sector(i,0)
                 if (kk == -1) cycle
                 indx1 = sectors(k)%istart
                 indx2 = sectors(kk)%istart
                 do l=1, sectors(k)%ndim
                     do m=1, sectors(kk)%ndim
                         ob = sectors(k)%myfmat(i,0)%item(m,l) ** 2 * (prob(indx2+m-1) + prob(indx1+l-1))    
                         cb = cmesh(j) + eigs(indx2+m-1) - eigs(indx1+l-1)
                         caux = caux + ob / cb
                     enddo 
                 enddo
             enddo ! over k={1,nsectors}
             ghub(j,i) = caux
         enddo ! over j={1,mfreq}
     enddo ! over i={1,norbs}

! calculate atomic self-energy function using dyson's equation
     do i=1,norbs
         do k=1,mfreq
             shub(k,i) = cmesh(k) + mune - eimp(i) - one / ghub(k,i)
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
         call ctqmc_zmat_inv(norbs, gaux)
         do i=1,norbs
             sig2(k,i,i) = cmesh(k) + mune - eimp(i) - gaux(i,i) - hybf(k,i,i)
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
             gaux(i,i) = cmesh(k) + mune - eimp(i) - sig2(k,i,i) - hybf(k,i,i)
         enddo ! over i={1,norbs} loop
         call ctqmc_zmat_inv(norbs, gaux)
         grnf(k,:,:) = gaux
     enddo ! over k={1,mfreq} loop

     return
  end subroutine ctqmc_make_hub1
