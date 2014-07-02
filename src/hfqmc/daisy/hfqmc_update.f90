!-------------------------------------------------------------------------
! project : daisy
! program : hfqmc_make_detrat
!           hfqmc_make_accept
!           hfqmc_make_vertex
!           cat_dirty_update
!           cat_clean_update
!           cat_delay_update
!           cat_clear_update
! source  : hfqmc_update.f90
! type    : subroutine
! author  : li huang (email:huangli712@yahoo.com.cn)
! history : 01/06/2006 by li huang
!           01/16/2008 by li huang
!           10/26/2008 by li huang
!           10/31/2008 by li huang
!           11/02/2008 by li huang
!           04/18/2009 by li huang
!           12/24/2009 by li huang
!           02/26/2010 by li huang
!           03/08/2010 by li huang
!           03/28/2010 by li huang
! purpose : to provide the core subroutines for Hirsch-Fye quantum Monte
!           Carlo quantum impurity solver.
!           to calculate the transition probability, to wrap the green's
!           function matrix, to update the green's function from weiss's
!           function by using traditional update algorithm, and recently
!           proposed delayed update algorithm.
! input   :
! output  :
! status  : unstable
! comment : need blas and lapack support
!-------------------------------------------------------------------------

!>>> to calculate the transition probability using heat-bath algorithm
  subroutine hfqmc_make_detrat(n, m, p)
     use constants
     use control
     use context

     implicit none

! external arguments
! current index of time slices
     integer, intent(in)   :: n

! current index of auxiliary ising-like fields
     integer, intent(in)   :: m

! transition probability for two successive configurations
     real(dp), intent(out) :: p

! local variables
! index of spin up orbitals
     integer  :: pup

! index of spin down orbitals
     integer  :: pdn

! transition ratio for spin up and spin down parts
     real(dp) :: rup
     real(dp) :: rdn

! real(dp) dummy variables
     real(dp) :: raux

! calculate intermediate variable
     raux = two * imat( n, m )

! evaluate index of spin up and spin down orbitals
     pup = pmat( m, 1 )
     pdn = pmat( m, 2 )

! calculate transition ratios for spin up and spin down part, respectively
     if ( mstep == 1 ) then ! traditional update algorithm
         rup = one + ( one - gmat(n,n,pup) ) * ( exp(-raux) - one )
         rdn = one + ( one - gmat(n,n,pdn) ) * ( exp( raux) - one )
     else                   ! delayed update algorithm
         rup = one + ( one - diag(n,  pup) ) * ( exp(-raux) - one )
         rdn = one + ( one - diag(n,  pdn) ) * ( exp( raux) - one )
     endif ! back if ( mstep == 1 ) block

! calculate total ratio between old and new determinants
     raux = rup * rdn

! calculate transition probability between two successive configurations
! using heat-bath dynamics
     p = raux / ( one + raux )

     return
  end subroutine hfqmc_make_detrat

!>>> core subroutine, to update the auxiliary ising-like fields and the
! green's function matrix
  subroutine hfqmc_make_accept(n, m, cstep)
     use constants
     use control
     use context

     implicit none

! external arguments
! current index of time slices
     integer, intent(in) :: n

! current index of auxiliary ising-like fields
     integer, intent(in) :: m

! current nsweep steps
     integer, intent(in) :: cstep

! local variables
! loop index
     integer :: j

! index of spin up orbitals
     integer :: pup

! index of spin down orbitals
     integer :: pdn

! evaluate index of spin up and spin down orbitals
     pup = pmat( m, 1 )
     pdn = pmat( m, 2 )

! flip current auxiliary ising-like fields
     imat(n,m) = -imat(n,m)

! applying clean update formula
! update green's function matrix from original weiss's function
     if ( mod(cstep, nclean) == 0 .and. mod(n, ntime/4) == 0 ) then
         call cat_clean_update( pup ) ! spin up
         call cat_clean_update( pdn ) ! spin down

! if delayed update algorithm is used, the diag must be update on time
         if ( mstep > 1 ) then
             do j=1,ntime
                 diag(j,pup) = gmat(j,j,pup)
                 diag(j,pdn) = gmat(j,j,pdn)
             enddo ! over j={1,ntime} loop
         endif ! back if ( mstep > 1 ) block

! if delayed update algorithm is used, the ktep must be update on time
         if ( mstep > 1 ) then
             ktep(pup) = 0
             ktep(pdn) = 0
         endif ! back if ( mstep > 1 ) block

! applying dirty update formula
! update green's function matrix from previous ising field configuration
     else
         if ( mstep > 1 ) then ! delayed update algorithm
             call cat_delay_update( n, m, pup, +1 ) ! spin up
             call cat_delay_update( n, m, pdn, -1 ) ! spin down
         else                  ! traditional update algorithm
             call cat_dirty_update( n, m, pup, +1 ) ! spin up
             call cat_dirty_update( n, m, pdn, -1 ) ! spin down
         endif ! back if ( mstep > 1 ) block

     endif ! back if ( mod(cstep, nclean) == 0 .and. mod(n, ntime/4) == 0 ) block

     return
  end subroutine hfqmc_make_accept

!>>> wrap green's function matrix to green's function
  subroutine hfqmc_make_vertex(pp, tmat, msum)
     use constants
     use control
     use context

     implicit none

! external arguments
! current index of orbitals
     integer, intent(in)   :: pp

! green's function matrix
     real(dp), intent(in)  :: tmat(ntime,ntime,norbs)

! wrapped green's function matrix
     real(dp), intent(out) :: msum(-ntime+1:ntime-1)

! local variables
! loop index
     integer  :: i
     integer  :: j

! real(dp) dummy wrapping matrix
     real(dp) :: gaux(-ntime+1:ntime-1)

! init msum matrix
     msum = zero

! init gaux matrix
     gaux = zero

! consider negative times (-) part
     do j=0,ntime-1
         do i=1,ntime-j
             gaux(-j) = gaux(-j) + tmat(i,i+j,pp)
         enddo ! over i={1,ntime-j} loop
     enddo ! over j={0,ntime-1} loop

! consider positive times (+) part
     do j=1,ntime-1
         do i=1,ntime-j
             gaux(+j) = gaux(+j) + tmat(i+j,i,pp)
         enddo ! over i={1,ntime-j} loop
     enddo ! over j={1,ntime-1} loop

! it is enforced the time antiperiodicity
     do i=1,ntime-1
         msum(i) = ( gaux(i) - gaux(i-ntime) ) / real(ntime)
         msum(i-ntime) = -msum(i)
     enddo ! over i={1,ntime-1} loop

! normalize properly
     msum(0) = gaux(0) / real(ntime)

     return
  end subroutine hfqmc_make_vertex

!>>> record changes of accepted move on the green's function matrix, by
! using traditional update algorithm
  subroutine cat_dirty_update(it, is, pp, ps)
     use constants
     use control
     use context

     implicit none

! external arguments
! current index of time slices
     integer, intent(in) :: it

! current index of auxiliary ising-like fields
     integer, intent(in) :: is

! index of spin up and spin down orbitals
     integer, intent(in) :: pp

! sign for auxiliary ising-like fields
     integer, intent(in) :: ps

! local variables
! loop index
     integer  :: i
     integer  :: j

! $\mathcal(A)$ for spin up and down
     real(dp) :: aa

! real(dp) dummy variables
     real(dp) :: dd

! real(dp) dummy matrix, used to accelerate the subroutine
     real(dp) :: tt(ntime)
     real(dp) :: uu(ntime)

! setup immediate variables
     dd = exp( ps * two * imat(it,is) ) - one

! setup $mathcal(A)$
     aa = dd / ( one + ( one - gmat(it,it,pp) ) * dd )

! prepare tt, uu
     tt = aa * gmat(it,:,pp)
     uu = gmat(:,it,pp) - unity(:,it)

! build new gmat using fast-update equation
     do j=1,ntime
         do i=1,ntime
             gmat(i,j,pp) = gmat(i,j,pp) + uu(i) * tt(j)
         enddo ! over i={1,ntime} loop
     enddo ! over j={1,ntime} loop

     return
  end subroutine cat_dirty_update

!>>> update green's function by orignial weiss's function
  subroutine cat_clean_update(pp)
     use constants
     use control
     use context

     implicit none

! external arguments
! index of spin up and spin down orbitals
     integer, intent(in)  :: pp

! local variables
! loop index
     integer  :: i

! status flag
     integer  :: ierr

! working space for lapack subroutines
     integer  :: ipiv(ntime)

! \exp(V) matrix
     real(dp) :: vmat(ntime)

! $\mathcal(A)$ matrix
     real(dp) :: amat(ntime,ntime)

! build vmat matrix
     do i=1,ntime
         vmat(i) = exp( dot_product( smat(pp,:), imat(i,:) ) ) - one
     enddo ! over i={1,ntime} loop

! build amat matrix
     do i=1,ntime
         amat(:,i) = -vmat(i) * wmat(:,i,pp)
     enddo ! over j={1,ntime} loop

     do i=1,ntime
         amat(i,i) = amat(i,i) + one + vmat(i)
     enddo ! over i={1,ntime} loop

! calculate A^{-1} G0 = G
! since A G = G0, we found that solving linear systems is much faster than
! inverse A matrix directly
! note: on input, gmat contain wmat (G0), on output, it contain G. amat is
! destroyed on output
     gmat(:,:,pp) = wmat(:,:,pp)
     call dgesv(ntime, ntime, amat, ntime, ipiv, gmat(:,:,pp), ntime, ierr)
     if ( ierr /= 0 ) then
         call hfqmc_print_error('cat_clean_update','error in lapack subroutine dgesv')
     endif

     return
  end subroutine cat_clean_update

!>>> record changes of accepted move on the green's function matrix, by
! using delayed update algorithm
  subroutine cat_delay_update(it, is, pp, ps)
     use constants
     use control
     use context

     implicit none

! external arguments
! current index of time slices
     integer, intent(in) :: it

! current index of auxiliary ising-like fields
     integer, intent(in) :: is

! index of spin up and spin down orbitals
     integer, intent(in) :: pp

! sign for auxiliary ising-like fields
     integer, intent(in) :: ps

! local variables
! loop index
     integer  :: i

! current cyclic steps
     integer  :: ks

! $\mathcal(A)$ for spin up and down
     real(dp) :: aa

! real(dp) dummy variables
     real(dp) :: dd

! increase cyclic steps for current orbitals
     ktep(pp) = ktep(pp) + 1
     ks = ktep(pp)

! setup immediate variables
     dd = exp( ps * two * imat(it,is) ) - one

! setup $mathcal(A)$
     aa = dd / ( one + ( one - diag(it,pp) ) * dd )

! evaluate the cyclic arrays needed by delayed update algorithm
! here atep is used to save the column elements of green's function matrix
! and  btep is used to save the row    elements of green's function matrix
     atep(:,ks,pp) = gmat(:,it,pp)
     btep(:,ks,pp) = gmat(it,:,pp)
     do i=1,ks-1
         atep(:,ks,pp) = atep(:,ks,pp) + btep(it,i,pp) * atep(:,i,pp)
         btep(:,ks,pp) = btep(:,ks,pp) + atep(it,i,pp) * btep(:,i,pp)
     enddo ! over i={1,ks-1} loop
     atep(:,ks,pp) = aa * ( atep(:,ks,pp) - unity(:,it) )

! evaluate the new diagonal elements of green's function matrix,
! it is useful for calculating the transition ratio
     do i=1,ntime
         diag(i,pp) = diag(i,pp) + atep(i,ks,pp) * btep(i,ks,pp)
     enddo ! over i={1,ntime} loop

! if cyclic steps arrive the maximun allowed value, we need to recalculate
! the green's function matrix using atep and btep, and reset everything
     if ( ktep(pp) == mstep ) then
         call dgemm('N','T', ntime, ntime, mstep, &
                        one, atep(:,:,pp), ntime, &
                             btep(:,:,pp), ntime, &
                        one, gmat(:,:,pp), ntime )

! update diag using new green's function matrix
         do i=1,ntime
             diag(i,pp) = gmat(i,i,pp)
         enddo ! over i={1,ntime} loop

! reset ktep for current orbital to zero
         ktep(pp) = 0
     endif ! back if ( ktep(pp) == mstep ) block

     return
  end subroutine cat_delay_update

!>>> record changes of accepted move on the green's function matrix, by
! using delayed update algorithm
  subroutine cat_clear_update()
     use constants
     use control
     use context

     implicit none

! local variables
! loop index
     integer :: i
     integer :: j

! sometimes, though cyclic steps do not arrive the maximun allowed value,
! we need to prepare gmat and related cyclic arrays for next QMC sweep.
! so it is necessary to recalculate the green's function matrix using atep
! and btep, and reset everything
     do i=1,norbs
         if ( ktep(i) /= 0 ) then
             call dgemm('N', 'T', ntime, ntime, ktep(i), &
                                one, atep(:,:,i), ntime, &
                                     btep(:,:,i), ntime, &
                                one, gmat(:,:,i), ntime )

! update diag using new green's function matrix
             do j=1,ntime
                 diag(j,i) = gmat(j,j,i)
             enddo ! over j={1,ntime} loop

! reset ktep for current orbital to zero
             ktep(i) = 0
         endif ! back if ( ktep(i) /= 0 ) block
     enddo ! over i={1,norbs} loop

     return
  end subroutine cat_clear_update
