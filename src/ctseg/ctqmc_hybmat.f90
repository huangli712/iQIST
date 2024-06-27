!!!-----------------------------------------------------------------------
!!! project : iqist @ narcissus
!!! program : cat_insert_matrix
!!!           cat_remove_matrix
!!!           cat_lshift_matrix
!!!           cat_rshift_matrix
!!!           cat_reflip_matrix
!!!           cat_reload_matrix
!!!           cat_insert_detrat
!!!           cat_remove_detrat
!!!           cat_lshift_detrat
!!!           cat_rshift_detrat
!!!           cat_reflip_detrat
!!!           cat_reload_detrat
!!! source  : ctqmc_hybmat.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli@caep.cn)
!!! history : 09/16/2009 by li huang (created)
!!!           07/03/2023 by li huang (last modified)
!!! purpose : offer basic infrastructure (elementary updating subroutines)
!!!           for hybridization expansion version continuous time quantum
!!!           Monte Carlo (CTQMC) quantum impurity solver. the following
!!!           subroutines mainly deal with the \mathscr{M} matrix: mmat,
!!!           and \mathscr{G} matrix: gmat.
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> service layer: update M and G matrices                           <<<
!!========================================================================

!!
!! @sub cat_insert_matrix
!!
!! update the mmat matrix and gmat matrix for inserting new creation and
!! annihilation operators
!!
  subroutine cat_insert_matrix(flvr, is, ie, tau_start, tau_end, deter_ratio)
     use constants, only : dp
     use constants, only : zero, one, czero

     use control, only : nfreq
     use control, only : beta

     use context, only : ckink
     use context, only : index_s, index_e
     use context, only : exp_s, exp_e
     use context, only : lspace, rspace
     use context, only : lsaves, rsaves
     use context, only : mmat, gmat

     implicit none

!! external arguments
     ! current flavor channel
     integer, intent(in)  :: flvr

     ! index address to insert new creation and annihilation operators
     integer, intent(in)  :: is
     integer, intent(in)  :: ie

     ! imaginary time \tau_s for new creation operator
     real(dp), intent(in) :: tau_start

     ! imaginary time \tau_e for new annihilation operator
     real(dp), intent(in) :: tau_end

     ! previous calculated determinant ratio
     real(dp), intent(in) :: deter_ratio

!! local variables
     ! loop index over operators
     integer  :: i
     integer  :: j

     ! loop index over frequencies
     integer  :: k

     ! dummy real(dp) variables
     real(dp) :: p

!! [body

     ! evaluate p at first
     p = one / deter_ratio

     ! shift lspace and rspace, and then supplement them with -1 at the end
     do i=ckink,ie,-1
         lspace(i+1, flvr) = lspace(i, flvr)
     enddo ! over i={ckink,ie,-1} loop
     lspace(ie, flvr) = -one

     do j=ckink,is,-1
         rspace(j+1, flvr) = rspace(j, flvr)
     enddo ! over j={ckink,is,-1} loop
     rspace(is, flvr) = -one

     ! scale lspace with p
     do i=1,ckink+1
         lspace(i, flvr) = lspace(i, flvr) * p
     enddo ! over i={1,ckink+1} loop

     ! shift mmat matrix
     do j=ckink,is,-1
         do i=ckink,ie,-1
             mmat(i+1, j+1, flvr) = mmat(i, j, flvr)
         enddo ! over i={ckink,ie,-1} loop
     enddo ! over j={ckink,is,-1} loop

     do j=ckink,is,-1
         do i=1,ie-1
             mmat(i, j+1, flvr) = mmat(i, j, flvr)
         enddo ! over i={1,ie-1} loop
     enddo ! over j={ckink,is,-1} loop

     do j=1,is-1
         do i=ckink,ie,-1
             mmat(i+1, j, flvr) = mmat(i, j, flvr)
         enddo ! over i={ckink,ie,-1} loop
     enddo ! over j={1,is-1} loop

     ! supplement mmat matrix with zero
     do i=1,ckink+1
         mmat(i, is, flvr) = zero
     enddo ! over i={1,ckink+1} loop

     do j=1,ckink+1
         mmat(ie, j, flvr) = zero
     enddo ! over j={1,ckink+1} loop

     ! finally evaluate mmat matrix
     do j=1,ckink+1
         do i=1,ckink+1
             mmat(i, j, flvr) = mmat(i, j, flvr) + lspace(i, flvr) * rspace(j, flvr)
         enddo ! over i={1,ckink+1} loop
     enddo ! over j={1,ckink+1} loop

     ! update the perturbation expansion series
     call cat_insert_colour(flvr, is, ie, tau_start, tau_end)

     ! update gmat matrix
     lsaves(:, flvr) = czero
     rsaves(:, flvr) = czero

     do i=1,ckink+1
         do k=1,nfreq
             lsaves(k, flvr) = lsaves(k, flvr) +         exp_e(k, index_e(i, flvr), flvr)   * lspace(i, flvr)
             rsaves(k, flvr) = rsaves(k, flvr) + dconjg( exp_s(k, index_s(i, flvr), flvr) ) * rspace(i, flvr)
         enddo ! over k={1,nfreq} loop
     enddo ! over i={1,ckink+1} loop

     p = one / beta ! we redefine p here
     do k=1,nfreq
         gmat(k, flvr, flvr) = gmat(k, flvr, flvr) - lsaves(k, flvr) * rsaves(k, flvr) * p
     enddo ! over k={1,nfreq} loop

!<   ! only for debug
!<   do i=1,ckink+1
!<       do j=1,ckink+1
!<           print *, 'M:', i, j, mmat(i, j, flvr)
!<       enddo ! over j={1,ckink+1} loop
!<   enddo ! over i={1,ckink+1} loop
!<
!<   print *, 'G1:', flvr, gmat(1, flvr, flvr)
!<   print *, 'G2:', flvr, gmat(2, flvr, flvr)
!<   print *, 'G3:', flvr, gmat(3, flvr, flvr)
!<   print *, 'Gn:', flvr, gmat(nfreq, flvr, flvr)

!! body]

     return
  end subroutine cat_insert_matrix

!!
!! @sub cat_remove_matrix
!!
!! update the mmat matrix and gmat matrix for removing old creation and
!! annihilation operators
!!
  subroutine cat_remove_matrix(flvr, is, ie)
     use constants, only : dp
     use constants, only : one, czero

     use control, only : nfreq
     use control, only : beta

     use context, only : ckink
     use context, only : index_s, index_e
     use context, only : exp_s, exp_e
     use context, only : lsaves, rsaves
     use context, only : mmat, gmat

     implicit none

!! external arguments
     ! current flavor channel
     integer, intent(in) :: flvr

     ! index address to remove old creation and annihilation operators
     integer, intent(in) :: is
     integer, intent(in) :: ie

!! local variables
     ! loop index over operators
     integer  :: i
     integer  :: j

     ! loop index over frequencies
     integer  :: k

     ! dummy real(dp) variables
     real(dp) :: p

!! [body

     ! update gmat matrix
     lsaves(:, flvr) = czero
     rsaves(:, flvr) = czero

     do i=1,ckink
         do k=1,nfreq
             lsaves(k, flvr) = lsaves(k, flvr) +         exp_e(k, index_e(i, flvr), flvr)   * mmat(i, is, flvr)
             rsaves(k, flvr) = rsaves(k, flvr) + dconjg( exp_s(k, index_s(i, flvr), flvr) ) * mmat(ie, i, flvr)
         enddo ! over k={1,nfreq} loop
     enddo ! over i={1,ckink} loop

     p = one / ( mmat(ie, is, flvr) * beta )
     do k=1,nfreq
         gmat(k, flvr, flvr) = gmat(k, flvr, flvr) + lsaves(k, flvr) * rsaves(k, flvr) * p
     enddo ! over k={1,nfreq} loop

     ! update mmat matrix
     p = one / mmat(ie, is, flvr) ! we redefine p here
     do j=1,ckink
         do i=1,ckink
             if ( i /= ie .and. j /= is ) then
                 mmat(i, j, flvr) = mmat(i, j, flvr) - mmat(i, is, flvr) * mmat(ie, j, flvr) * p
             endif ! back if ( i /= ie .and. j /= is ) block
         enddo ! over i={1,ckink} loop
     enddo ! over j={1,ckink} loop

     do j=is,ckink-1
         do i=ie,ckink-1
             mmat(i, j, flvr) = mmat(i+1, j+1, flvr)
         enddo ! over i={ie,ckink-1} loop
     enddo ! over j={is,ckink-1} loop

     do j=is,ckink-1
         do i=1,ie-1
             mmat(i, j, flvr) = mmat(i, j+1, flvr)
         enddo ! over i={1,ie-1} loop
     enddo ! over j={is,ckink-1} loop

     do j=1,is-1
         do i=ie,ckink-1
             mmat(i, j, flvr) = mmat(i+1, j, flvr)
         enddo ! over i={ie,ckink-1} loop
     enddo ! over j={1,is-1} loop

     ! update the perturbation expansion series
     call cat_remove_colour(flvr, is, ie)

!<   ! only for debug
!<   do i=1,ckink
!<       do j=1,ckink
!<           print *, 'M:', i, j, mmat(i, j, flvr)
!<       enddo ! over j={1,ckink} loop
!<   enddo ! over i={1,ckink} loop
!<
!<   print *, 'G1:', flvr, gmat(1, flvr, flvr)
!<   print *, 'G2:', flvr, gmat(2, flvr, flvr)
!<   print *, 'G3:', flvr, gmat(3, flvr, flvr)
!<   print *, 'Gn:', flvr, gmat(nfreq, flvr, flvr)

!! body]

     return
  end subroutine cat_remove_matrix

!!
!! @sub cat_lshift_matrix
!!
!! update the mmat matrix and gmat matrix for shifting creation operator
!!
  subroutine cat_lshift_matrix(flvr, iso, isn, tau_start1, tau_start2, deter_ratio)
     use constants, only : dp
     use constants, only : zero, czero

     use control, only : mkink
     use control, only : nfreq
     use control, only : beta

     use context, only : ckink
     use context, only : index_s, index_e
     use context, only : time_e
     use context, only : exp_s, exp_e
     use context, only : rmesh
     use context, only : lspace, rspace
     use context, only : lsaves, rsaves
     use context, only : mmat, gmat

     implicit none

!! external arguments
     ! current flavor channel
     integer, intent(in)  :: flvr

     ! index address to shift creation operator
     ! iso and isn are old and new indices, respectively
     integer, intent(in)  :: iso
     integer, intent(in)  :: isn

     ! imaginary time \tau_s for creation operator (the old one)
     real(dp), intent(in) :: tau_start1

     ! imaginary time \tau_s for creation operator (the new one)
     real(dp), intent(in) :: tau_start2

     ! previous calculated determinant ratio
     real(dp), intent(in) :: deter_ratio

!! external arguments
     ! used to interpolate the hybridization function
     procedure( real(dp) ) :: ctqmc_eval_htau

!! local variables
     ! loop index over operators
     integer  :: i
     integer  :: j

     ! loop index over frequencies
     integer  :: k

     ! used to store matrix element of mmat
     real(dp) :: md

     ! dummy real(dp) variables
     real(dp) :: xs
     real(dp) :: rs

     ! dummy real(dp) arrays, used to interpolate the hybridization function
     real(dp) :: lvec(mkink)
     real(dp) :: rvec(mkink)

     ! dummy complex(dp) arrays, used to calculate gmat matrix
     complex(dp) :: lexp(nfreq)
     complex(dp) :: gsum(nfreq)
     complex(dp) :: gdel(nfreq)

!! [body

     ! evaluate lexp
     lexp = czero
     do k=1,nfreq
         xs = tau_start2 * rmesh(k)
         lexp(k) = dconjg( exp_s(k, index_s(iso, flvr), flvr) ) - dcmplx( cos(xs), -sin(xs) )
     enddo ! over k={1,nfreq} loop
     lexp = lexp / beta

     ! evaluate gsum
     gsum = czero
     do i=1,ckink
         md = mmat(i, iso, flvr)
         do k=1,nfreq
             gsum(k) = gsum(k) + exp_e(k, index_e(i, flvr), flvr) * md
         enddo ! over k={1,nfreq} loop
     enddo ! over i={1,ckink} loop

     ! evaluate gdel, \delta G for gmat matrix
     gdel = czero
     do k=1,nfreq
         gdel(k) = gsum(k) * lexp(k)
     enddo ! over k={1,nfreq} loop

     ! calculate rvec by cubic spline interpolation
     do i=1,ckink
         if ( tau_start1 < time_e(index_e(i, flvr), flvr) ) then
             rvec(i) = -ctqmc_eval_htau(flvr, tau_start1 - time_e(index_e(i, flvr), flvr) + beta)
         else
             rvec(i) =  ctqmc_eval_htau(flvr, tau_start1 - time_e(index_e(i, flvr), flvr))
         endif ! back if ( tau_start1 < time_e(index_e(i, flvr), flvr) ) block
     enddo ! over i={1,ckink} loop

     ! calculate lvec by cubic spline interpolation
     do j=1,ckink
         if ( tau_start2 < time_e(index_e(j, flvr), flvr) ) then
             lvec(j) = -ctqmc_eval_htau(flvr, tau_start2 - time_e(index_e(j, flvr), flvr) + beta)
         else
             lvec(j) =  ctqmc_eval_htau(flvr, tau_start2 - time_e(index_e(j, flvr), flvr))
         endif ! back if ( tau_start2 < time_e(index_e(j, flvr), flvr) ) block
     enddo ! over j={1,ckink} loop

     ! adjust rvec
     do i=1,ckink
         rvec(i) = lvec(i) - rvec(i)
     enddo ! over i={1,ckink} loop

     ! prepare rspace
     do i=1,ckink
         rs = zero
         do j=1,ckink
             rs = rs + rvec(j) * mmat(j, i, flvr)
         enddo ! over j={1,ckink} loop
         rspace(i, flvr) = rs / deter_ratio
     enddo ! over i={1,ckink} loop

     ! prepare lspace
     do i=1,ckink
         lspace(i, flvr) = -mmat(i, iso, flvr)
     enddo ! over i={1,ckink} loop

     ! calculate mmat matrix
     do j=1,ckink
         do i=1,ckink
             mmat(i, j, flvr) = mmat(i, j, flvr) + lspace(i, flvr) * rspace(j, flvr)
         enddo ! over i={1,ckink} loop
     enddo ! over j={1,ckink} loop

     ! shuffle rows if time order changed because of move
     if ( isn /= iso ) then
         rs = rspace(iso, flvr)
         do i=1,ckink
             rvec(i) = mmat(i, iso, flvr)
         enddo ! over i={1,ckink} loop

         if ( isn < iso ) then
             do j=iso,isn+1,-1
                 do i=1,ckink
                     mmat(i, j, flvr) = mmat(i, j-1, flvr)
                 enddo ! over i={1,ckink} loop
                 rspace(j, flvr) = rspace(j-1, flvr)
             enddo ! over j={iso,isn+1,-1} loop
             do i=1,ckink
                 mmat(i, isn, flvr) = rvec(i)
             enddo ! over i={1,ckink} loop
             rspace(isn, flvr) = rs
         else
             do j=iso,isn-1
                 do i=1,ckink
                     mmat(i, j, flvr) = mmat(i, j+1, flvr)
                 enddo ! over i={1,ckink} loop
                 rspace(j, flvr) = rspace(j+1, flvr)
             enddo ! over j={iso,isn-1} loop
             do i=1,ckink
                 mmat(i, isn, flvr) = rvec(i)
             enddo ! over i={1,ckink} loop
             rspace(isn, flvr) = rs
         endif ! back if ( isn < iso ) block
     endif ! back if ( isn /= iso ) block

     ! update the perturbation expansion series
     call cat_lshift_colour(flvr, iso, isn, tau_start2)

     ! update gmat matrix
     lsaves(:, flvr) = czero
     rsaves(:, flvr) = czero

     do i=1,ckink
         do k=1,nfreq
             lsaves(k, flvr) = lsaves(k, flvr) +         exp_e(k, index_e(i, flvr), flvr)   * lspace(i, flvr)
             rsaves(k, flvr) = rsaves(k, flvr) + dconjg( exp_s(k, index_s(i, flvr), flvr) ) * rspace(i, flvr)
         enddo ! over k={1,nfreq} loop
     enddo ! over i={1,ckink} loop

     do k=1,nfreq
         gmat(k, flvr, flvr) = gmat(k, flvr, flvr) - lsaves(k, flvr) * rsaves(k, flvr) / beta + gdel(k)
     enddo ! over k={1,nfreq} loop

!<   ! only for debug
!<   do i=1,ckink
!<       do j=1,ckink
!<           print *,'M:',i, j, mmat(i, j, flvr)
!<       enddo ! over j={1,ckink} loop
!<   enddo ! over i={1,ckink} loop
!<
!<   print *, 'G1:', flvr, gmat(1, flvr, flvr)
!<   print *, 'G2:', flvr, gmat(2, flvr, flvr)
!<   print *, 'G3:', flvr, gmat(3, flvr, flvr)
!<   print *, 'Gn:', flvr, gmat(nfreq, flvr, flvr)

!! body]

     return
  end subroutine cat_lshift_matrix

!!
!! @sub cat_rshift_matrix
!!
!! update the mmat matrix and gmat matrix for shifting annihilation operator
!!
  subroutine cat_rshift_matrix(flvr, ieo, ien, tau_end1, tau_end2, deter_ratio)
     use constants, only : dp
     use constants, only : zero, czero

     use control, only : mkink
     use control, only : nfreq
     use control, only : beta

     use context, only : ckink
     use context, only : index_s, index_e
     use context, only : time_s
     use context, only : exp_s, exp_e
     use context, only : rmesh
     use context, only : lspace, rspace
     use context, only : lsaves, rsaves
     use context, only : mmat, gmat

     implicit none

!! external arguments
     ! current flavor channel
     integer, intent(in)  :: flvr

     ! index address to shift annihilation operator
     ! ieo and ien are old and new indices, respectively
     integer, intent(in)  :: ieo
     integer, intent(in)  :: ien

     ! imaginary time \tau_e for annihilation operator (the old one)
     real(dp), intent(in) :: tau_end1

     ! imaginary time \tau_e for annihilation operator (the new one)
     real(dp), intent(in) :: tau_end2

     ! previous calculated determinant ratio
     real(dp), intent(in) :: deter_ratio

!! external arguments
     ! used to interpolate the hybridization function
     procedure( real(dp) ) :: ctqmc_eval_htau

!! local variables
     ! loop index over operators
     integer  :: i
     integer  :: j

     ! loop index over frequencies
     integer  :: k

     ! used to store matrix element of mmat
     real(dp) :: md

     ! dummy real(dp) variables
     real(dp) :: xe
     real(dp) :: ls

     ! dummy real(dp) arrays, used to interpolate the hybridization function
     real(dp) :: lvec(mkink)
     real(dp) :: rvec(mkink)

     ! dummy complex(dp) arrays, used to calculate gmat matrix
     complex(dp) :: rexp(nfreq)
     complex(dp) :: gsum(nfreq)
     complex(dp) :: gdel(nfreq)

!! [body

     ! evaluate rexp
     rexp = czero
     do k=1,nfreq
         xe = tau_end2 * rmesh(k)
         rexp(k) = exp_e(k, index_e(ieo, flvr), flvr) - dcmplx( cos(xe), sin(xe) )
     enddo ! over k={1,nfreq} loop
     rexp = rexp / beta

     ! evaluate gsum
     gsum = czero
     do i=1,ckink
         md = mmat(ieo, i, flvr)
         do k=1,nfreq
             gsum(k) = gsum(k) + dconjg( exp_s(k, index_s(i, flvr), flvr) ) * md
         enddo ! over k={1,nfreq} loop
     enddo ! over i={1,ckink} loop

     ! evaluate gdel, \delta G for gmat matrix
     gdel = czero
     do k=1,nfreq
         gdel(k) = gsum(k) * rexp(k)
     enddo ! over k={1,nfreq} loop

     ! calculate lvec by cubic spline interpolation
     do i=1,ckink
         if ( time_s(index_s(i, flvr), flvr) < tau_end1 ) then
             lvec(i) = -ctqmc_eval_htau(flvr, time_s(index_s(i, flvr), flvr) - tau_end1 + beta)
         else
             lvec(i) =  ctqmc_eval_htau(flvr, time_s(index_s(i, flvr), flvr) - tau_end1)
         endif ! back if ( time_s(index_s(i, flvr), flvr) < tau_end1 ) block
     enddo ! over i={1,ckink} loop

     ! calculate rvec by cubic spline interpolation
     do j=1,ckink
         if ( time_s(index_s(j, flvr), flvr) < tau_end2 ) then
             rvec(j) = -ctqmc_eval_htau(flvr, time_s(index_s(j, flvr), flvr) - tau_end2 + beta)
         else
             rvec(j) =  ctqmc_eval_htau(flvr, time_s(index_s(j, flvr), flvr) - tau_end2)
         endif ! back if ( time_s(index_s(j, flvr), flvr) < tau_end2 ) block
     enddo ! over j={1,ckink} loop

     ! adjust lvec
     do i=1,ckink
         lvec(i) = rvec(i) - lvec(i)
     enddo ! over i={1,ckink} loop

     ! prepare lspace
     do i=1,ckink
         ls = zero
         do j=1,ckink
             ls = ls + mmat(i, j, flvr) * lvec(j)
         enddo ! over j={1,ckink} loop
         lspace(i, flvr) = ls / deter_ratio
     enddo ! over i={1,ckink} loop

     ! prepare rspace
     do i=1,ckink
         rspace(i, flvr) = -mmat(ieo, i, flvr)
     enddo ! over i={1,ckink} loop

     ! calculate mmat matrix
     do j=1,ckink
         do i=1,ckink
             mmat(i, j, flvr) = mmat(i, j, flvr) + lspace(i, flvr) * rspace(j, flvr)
         enddo ! over i={1,ckink} loop
     enddo ! over j={1,ckink} loop

     ! shuffle columns if time order changed because of move
     if ( ien /= ieo ) then
         ls = lspace(ieo, flvr)
         do i=1,ckink
             lvec(i) = mmat(ieo, i, flvr)
         enddo ! over i={1,ckink} loop

         if ( ien < ieo ) then
             do j=ieo,ien+1,-1
                 do i=1,ckink
                     mmat(j, i, flvr) = mmat(j-1, i, flvr)
                 enddo ! over i={1,ckink} loop
                 lspace(j, flvr) = lspace(j-1, flvr)
             enddo ! over j={ieo,ien+1,-1} loop
             do i=1,ckink
                 mmat(ien, i, flvr) = lvec(i)
             enddo ! over i={1,ckink} loop
             lspace(ien, flvr) = ls
         else
             do j=ieo,ien-1
                 do i=1,ckink
                     mmat(j, i, flvr) = mmat(j+1, i, flvr)
                 enddo ! over i={1,ckink} loop
                 lspace(j, flvr) = lspace(j+1, flvr)
             enddo ! over j={ieo,ien-1} loop
             do i=1,ckink
                 mmat(ien, i, flvr) = lvec(i)
             enddo ! over i={1,ckink} loop
             lspace(ien, flvr) = ls
         endif ! back if ( ien < ieo ) block
     endif ! back if ( ien /= ieo ) block

     ! update the perturbation expansion series
     call cat_rshift_colour(flvr, ieo, ien, tau_end2)

     ! update gmat matrix
     lsaves(:, flvr) = czero
     rsaves(:, flvr) = czero

     do i=1,ckink
         do k=1,nfreq
             lsaves(k, flvr) = lsaves(k, flvr) +         exp_e(k, index_e(i, flvr), flvr)   * lspace(i, flvr)
             rsaves(k, flvr) = rsaves(k, flvr) + dconjg( exp_s(k, index_s(i, flvr), flvr) ) * rspace(i, flvr)
         enddo ! over k={1,nfreq} loop
     enddo ! over i={1,ckink} loop

     do k=1,nfreq
         gmat(k, flvr, flvr) = gmat(k, flvr, flvr) - lsaves(k, flvr) * rsaves(k, flvr) / beta + gdel(k)
     enddo ! over k={1,nfreq} loop

!<   ! only for debug
!<   do i=1,ckink
!<       do j=1,ckink
!<           print *,'M:',i, j, mmat(i, j, flvr)
!<       enddo ! over j={1,ckink} loop
!<   enddo ! over i={1,ckink} loop
!<
!<   print *, 'G1:', flvr, gmat(1, flvr, flvr)
!<   print *, 'G2:', flvr, gmat(2, flvr, flvr)
!<   print *, 'G3:', flvr, gmat(3, flvr, flvr)
!<   print *, 'Gn:', flvr, gmat(nfreq, flvr, flvr)

!! body]

     return
  end subroutine cat_rshift_matrix

!!
!! @sub cat_reflip_matrix
!!
!! global flip the time_s, time_e, mmat matrix, gmat matrix, and other
!! related global variables between spin up and spin down states. it is
!! used to avoid trapped by unphysical phase
!!
  subroutine cat_reflip_matrix(fup, fdn, kmax)
     use stack, only : istack
     use stack, only : istack_create
     use stack, only : istack_copyer
     use stack, only : istack_destroy

     use control, only : mkink
     use control, only : nfreq

     use context, only : empty_s, empty_e
     use context, only : index_s, index_e
     use context, only : time_s, time_e
     use context, only : exp_s, exp_e
     use context, only : rank, stts
     use context, only : gmat

     implicit none

!! external arguments
     ! current flavor channel
     integer, intent(in) :: fup
     integer, intent(in) :: fdn

     ! maximum rank order for current flavor channel
     integer, intent(in) :: kmax

!! local variables
     ! maximum memory index accessed by index_s and index_e
     integer :: ismax
     integer :: iemax

     ! dummy copy for rank and stts
     integer :: Trank
     integer :: Tstts

     ! dummy copy for empty_s and empty_e
     type (istack) :: Tempty_s
     type (istack) :: Tempty_e

!! [body

     ! allocate memory for Tempty_s and Tempty_e
     call istack_create(Tempty_s, mkink)
     call istack_create(Tempty_e, mkink)

     ! swap empty_s and empty_e
     call istack_copyer(empty_s(fup), Tempty_s)
     call istack_copyer(empty_e(fup), Tempty_e)

     call istack_copyer(empty_s(fdn), empty_s(fup))
     call istack_copyer(empty_e(fdn), empty_e(fup))

     call istack_copyer(Tempty_s, empty_s(fdn))
     call istack_copyer(Tempty_e, empty_e(fdn))

     ! deallocate memory for Tempty_s and Tempty_e
     call istack_destroy(Tempty_s)
     call istack_destroy(Tempty_e)

     ! swap rank
     Trank = rank(fup)
     rank(fup) = rank(fdn)
     rank(fdn) = Trank

     ! swap stts
     Tstts = stts(fup)
     stts(fup) = stts(fdn)
     stts(fdn) = Tstts

     ! swap gmat matrix when needed
     call s_swap_z(nfreq, gmat(1:nfreq, fup, fup), gmat(1:nfreq, fdn, fdn))

     if ( kmax > 0 ) then

         ! determine ismax and iemax
         ismax = max( maxval( index_s(1:kmax, fup) ), maxval( index_s(1:kmax, fdn) ) )
         iemax = max( maxval( index_e(1:kmax, fup) ), maxval( index_e(1:kmax, fdn) ) )

         ! swap index_s and index_e
         call s_swap_i(kmax, index_s(1:kmax, fup), index_s(1:kmax, fdn))
         call s_swap_i(kmax, index_e(1:kmax, fup), index_e(1:kmax, fdn))

         ! swap time_s and time_e
         call s_swap_d(ismax, time_s(1:ismax, fup), time_s(1:ismax, fdn))
         call s_swap_d(iemax, time_e(1:iemax, fup), time_e(1:iemax, fdn))

         ! swap exp_s and exp_e
         call s_swap_z(nfreq*ismax, exp_s(1:nfreq, 1:ismax, fup), exp_s(1:nfreq, 1:ismax, fdn))
         call s_swap_z(nfreq*iemax, exp_e(1:nfreq, 1:iemax, fup), exp_e(1:nfreq, 1:iemax, fdn))

         ! update mmat and gmat matrix when needed
         if ( rank(fup) > 0 ) call cat_reload_matrix(fup)
         if ( rank(fdn) > 0 ) call cat_reload_matrix(fdn)

     endif ! back if ( kmax > 0 ) block

!! body]

     return
  end subroutine cat_reflip_matrix

!!
!! @sub cat_reload_matrix
!!
!! global update the mmat matrix and gmat matrix from scratch
!!
  subroutine cat_reload_matrix(flvr)
     use constants, only : dp
     use constants, only : zero, czero

     use control, only : nfreq
     use control, only : beta

     use context, only : index_s, index_e
     use context, only : time_s, time_e
     use context, only : exp_s, exp_e
     use context, only : rank
     use context, only : mmat, gmat

     implicit none

!! external arguments
     ! current flavor channel
     integer, intent(in) :: flvr

!! external functions
     ! used to interpolate the hybridization function
     procedure( real(dp) ) :: ctqmc_eval_htau

!! local variables
     ! loop index over operators
     integer  :: i
     integer  :: j

     ! loop index over frequencies
     integer  :: k

     ! dummy perturbation expansion order
     integer  :: kaux

     ! used to store matrix element of mmat
     real(dp) :: maux

     ! imaginary time for creation and annihilation operators
     real(dp) :: tau_start
     real(dp) :: tau_end

     ! dummy complex(dp) variables
     complex(dp) :: x_start
     complex(dp) :: x_end

!! [body

     ! evaluate kaux
     kaux = rank(flvr)

     ! reset mmat matrix
     mmat(1:kaux, 1:kaux, flvr) = zero

     ! recalculate mmat from scratch
     do j=1,kaux
         tau_end = time_e(index_e(j, flvr), flvr)
         do i=1,kaux
             tau_start = time_s(index_s(i, flvr), flvr)
             if ( tau_start < tau_end ) then
                 mmat(i, j, flvr) = -ctqmc_eval_htau(flvr, tau_start - tau_end + beta)
             else
                 mmat(i, j, flvr) =  ctqmc_eval_htau(flvr, tau_start - tau_end)
             endif ! back if ( tau_start < tau_end ) block
         enddo ! over i={1,kaux} loop
     enddo ! over j={1,kaux} loop

     ! now we obtain dmat matrix, while what we need is its inversion
     call s_inv_d(kaux, mmat(1:kaux, 1:kaux, flvr))

     ! reset gmat matrix
     gmat(:, flvr, flvr) = czero

     ! recalculate gmat from scratch
     do j=1,kaux
         do i=1,kaux
             maux = -mmat(i, j, flvr) / beta
             do k=1,nfreq
                 x_start = dconjg( exp_s(k, index_s(j, flvr), flvr) )
                 x_end   =         exp_e(k, index_e(i, flvr), flvr)
                 gmat(k, flvr, flvr) = gmat(k, flvr, flvr) + x_end * maux * x_start
             enddo ! over k={1,nfreq} loop
         enddo ! over i={1,kaux} loop
     enddo ! over j={1,kaux} loop

!! body]

     return
  end subroutine cat_reload_matrix

!!========================================================================
!!>>> service layer: evaluate the determinant ratio                    <<<
!!========================================================================

!!
!! @sub cat_insert_detrat
!!
!! calculate the determinant ratio for inserting new creation and
!! annihilation operators
!!
  subroutine cat_insert_detrat(flvr, tau_start, tau_end, deter_ratio)
     use constants, only : dp
     use constants, only : zero

     use control, only : mkink
     use control, only : beta

     use context, only : ckink
     use context, only : index_s, index_e
     use context, only : time_s, time_e
     use context, only : lspace, rspace
     use context, only : mmat

     implicit none

!! external arguments
     ! current flavor channel
     integer, intent(in)   :: flvr

     ! imaginary time \tau_s for new creation operator
     real(dp), intent(in)  :: tau_start

     ! imaginary time \tau_e for new annihilation operator
     real(dp), intent(in)  :: tau_end

     ! the desired determinant ratio
     real(dp), intent(out) :: deter_ratio

!! external arguments
     ! used to interpolate the hybridization function
     procedure( real(dp) ) :: ctqmc_eval_htau

!! local variables
     ! loop index over operators
     integer  :: i
     integer  :: j

     ! dummy real(dp) variables
     real(dp) :: sl
     real(dp) :: sr

     ! dummy real(dp) arrays, used to interpolate the hybridization function
     real(dp) :: lvec(mkink)
     real(dp) :: rvec(mkink)

!! [body

     ! calculate lvec by cubic spline interpolation
     do i=1,ckink
         if ( time_s(index_s(i, flvr), flvr) < tau_end   ) then
             lvec(i) = -ctqmc_eval_htau(flvr, time_s(index_s(i, flvr), flvr) - tau_end + beta)
         else
             lvec(i) =  ctqmc_eval_htau(flvr, time_s(index_s(i, flvr), flvr) - tau_end)
         endif ! back if ( time_s(index_s(i, flvr), flvr) < tau_end   ) block
     enddo ! over i={1,ckink} loop

     ! calculate rvec by cubic spline interpolation
     do j=1,ckink
         if ( tau_start < time_e(index_e(j, flvr), flvr) ) then
             rvec(j) = -ctqmc_eval_htau(flvr, tau_start - time_e(index_e(j, flvr), flvr) + beta)
         else
             rvec(j) =  ctqmc_eval_htau(flvr, tau_start - time_e(index_e(j, flvr), flvr))
         endif ! back if ( tau_start < time_e(index_e(j, flvr), flvr) ) block
     enddo ! over j={1,ckink} loop

     ! calculate deter_ratio by cubic spline interpolation
     if ( tau_start > tau_end ) then
         deter_ratio =  ctqmc_eval_htau(flvr, tau_start - tau_end)
     else
         deter_ratio = -ctqmc_eval_htau(flvr, tau_start - tau_end + beta)
     endif ! back if ( tau_start > tau_end ) block

     ! calculate lspace and rspace
     do i=1,ckink
         sl = zero
         sr = zero

         do j=1,ckink
             sl = sl + mmat(i, j, flvr) * lvec(j)
             sr = sr + rvec(j) * mmat(j, i, flvr)
         enddo ! over j={1,ckink} loop

         lspace(i, flvr) = sl
         rspace(i, flvr) = sr
     enddo ! over i={1,ckink} loop

     ! calculate final determinant ratio
     do i=1,ckink
         deter_ratio = deter_ratio - rvec(i) * lspace(i, flvr)
     enddo ! over i={1,ckink} loop

!! body]

     return
  end subroutine cat_insert_detrat

!!
!! @sub cat_remove_detrat
!!
!! calculate the determinant ratio for removing old creation and
!! annihilation operators
!!
  subroutine cat_remove_detrat(flvr, is, ie, deter_ratio)
     use constants, only : dp

     use context, only : mmat

     implicit none

!! external arguments
     ! current flavor channel
     integer, intent(in)   :: flvr

     ! index address to remove old creation and annihilation operators
     ! is and ie are for creation and annihilation operators, respectively
     integer, intent(in)   :: is
     integer, intent(in)   :: ie

     ! the desired determinant ratio
     real(dp), intent(out) :: deter_ratio

!! [body

     deter_ratio = mmat(ie, is, flvr)

!! body]

     return
  end subroutine cat_remove_detrat

!!
!! @sub cat_lshift_detrat
!!
!! calculate the determinant ratio for shifting creation operator
!!
  subroutine cat_lshift_detrat(flvr, addr, tau_start1, tau_start2, deter_ratio)
     use constants, only : dp
     use constants, only : one

     use control, only : mkink
     use control, only : beta

     use context, only : ckink
     use context, only : index_e
     use context, only : time_e
     use context, only : mmat

     implicit none

!! external arguments
     ! current flavor channel
     integer, intent(in)   :: flvr

     ! index address to shift creation operator (old index = iso)
     integer, intent(in)   :: addr

     ! imaginary time \tau_s for creation operator (the old one)
     real(dp), intent(in)  :: tau_start1

     ! imaginary time \tau_s for creation operator (the new one)
     real(dp), intent(in)  :: tau_start2

     ! the desired determinant ratio
     real(dp), intent(out) :: deter_ratio

!! external functions
     ! used to interpolate the hybridization function
     procedure( real(dp) ) :: ctqmc_eval_htau

!! local variables
     ! loop index over operators
     integer  :: i
     integer  :: j

     ! dummy real(dp) arrays, used to interpolate the hybridization function
     real(dp) :: lvec(mkink)
     real(dp) :: rvec(mkink)

!! [body

     ! calculate rvec by cubic spline interpolation
     do i=1,ckink
         if ( tau_start1 < time_e(index_e(i, flvr), flvr) ) then
             rvec(i) = -ctqmc_eval_htau(flvr, tau_start1 - time_e(index_e(i, flvr), flvr) + beta)
         else
             rvec(i) =  ctqmc_eval_htau(flvr, tau_start1 - time_e(index_e(i, flvr), flvr))
         endif ! back if ( tau_start1 < time_e(index_e(i, flvr), flvr) ) block
     enddo ! over i={1,ckink} loop

     ! calculate lvec by cubic spline interpolation
     do j=1,ckink
         if ( tau_start2 < time_e(index_e(j, flvr), flvr) ) then
             lvec(j) = -ctqmc_eval_htau(flvr, tau_start2 - time_e(index_e(j, flvr), flvr) + beta)
         else
             lvec(j) =  ctqmc_eval_htau(flvr, tau_start2 - time_e(index_e(j, flvr), flvr))
         endif ! back if ( tau_start2 < time_e(index_e(j, flvr), flvr) ) block
     enddo ! over j={1,ckink} loop

     ! adjust rvec
     do i=1,ckink
         rvec(i) = lvec(i) - rvec(i)
     enddo ! over i={1,ckink} loop

     ! calculate final determinant ratio
     deter_ratio = one
     do i=1,ckink
         deter_ratio = deter_ratio + rvec(i) * mmat(i, addr, flvr)
     enddo ! over i={1,ckink} loop

!! body]

     return
  end subroutine cat_lshift_detrat

!!
!! @sub cat_rshift_detrat
!!
!! calculate the determinant ratio for shifting annihilation operator
!!
  subroutine cat_rshift_detrat(flvr, addr, tau_end1, tau_end2, deter_ratio)
     use constants, only : dp
     use constants, only : one

     use control, only : mkink
     use control, only : beta

     use context, only : ckink
     use context, only : index_s
     use context, only : time_s
     use context, only : mmat

     implicit none

!! external arguments
     ! current flavor channel
     integer, intent(in)   :: flvr

     ! index address to shift annihilation operator (old index = ieo)
     integer, intent(in)   :: addr

     ! imaginary time \tau_e for annihilation operator (the old one)
     real(dp), intent(in)  :: tau_end1

     ! imaginary time \tau_e for annihilation operator (the new one)
     real(dp), intent(in)  :: tau_end2

     ! the desired determinant ratio
     real(dp), intent(out) :: deter_ratio

!! external functions
     ! used to interpolate the hybridization function
     procedure( real(dp) ) :: ctqmc_eval_htau

!! local variables
     ! loop index over operators
     integer  :: i
     integer  :: j

     ! dummy real(dp) arrays, used to interpolate the hybridization function
     real(dp) :: lvec(mkink)
     real(dp) :: rvec(mkink)

!! [body

     ! calculate lvec by cubic spline interpolation
     do i=1,ckink
         if ( time_s(index_s(i, flvr), flvr) < tau_end1 ) then
             lvec(i) = -ctqmc_eval_htau(flvr, time_s(index_s(i, flvr), flvr) - tau_end1 + beta)
         else
             lvec(i) =  ctqmc_eval_htau(flvr, time_s(index_s(i, flvr), flvr) - tau_end1)
         endif ! back if ( time_s(index_s(i, flvr), flvr) < tau_end1 ) block
     enddo ! over i={1,ckink} loop

     ! calculate rvec by cubic spline interpolation
     do j=1,ckink
         if ( time_s(index_s(j, flvr), flvr) < tau_end2 ) then
             rvec(j) = -ctqmc_eval_htau(flvr, time_s(index_s(j, flvr), flvr) - tau_end2 + beta)
         else
             rvec(j) =  ctqmc_eval_htau(flvr, time_s(index_s(j, flvr), flvr) - tau_end2)
         endif ! back if ( time_s(index_s(j, flvr), flvr) < tau_end2 ) block
     enddo ! over j={1,ckink} loop

     ! adjust lvec
     do i=1,ckink
         lvec(i) = rvec(i) - lvec(i)
     enddo ! over i={1,ckink} loop

     ! calculate final determinant ratio
     deter_ratio = one
     do i=1,ckink
         deter_ratio = deter_ratio + mmat(addr, i, flvr) * lvec(i)
     enddo ! over i={1,ckink} loop

!! body]

     return
  end subroutine cat_rshift_detrat

!!
!! @sub cat_reflip_detrat
!!
!! calculate the determinant ratio for global spin flip
!!
  subroutine cat_reflip_detrat(up, dn, ratio)
     use constants, only : dp
     use constants, only : zero, one

     use control, only : beta

     use context, only : index_s, index_e
     use context, only : time_s, time_e
     use context, only : rank
     use context, only : mmat

     implicit none

!! external arguments
     ! band index for spin up state
     integer, intent(in)   :: up

     ! band index for spin dn state
     integer, intent(in)   :: dn

     ! the desired determinant ratio
     real(dp), intent(out) :: ratio

!! external functions
     ! used to interpolate the hybridization function
     procedure( real(dp) ) :: ctqmc_eval_htau

!! local variables
     ! loop index over operators
     integer  :: i
     integer  :: j

     ! dummy perturbation expansion order
     integer  :: kaux

     ! status flag
     integer  :: istat

     ! imaginary time for creation and annihilation operators
     real(dp) :: tau_start
     real(dp) :: tau_end

     ! dummy mmat matrix
     real(dp), allocatable :: Dmm(:,:)
     real(dp), allocatable :: Tmm(:,:)

!! [body

     ! evaluate kaux
     kaux = rank(up)

     ! check the status of kaux, if there does not exist any operators
     ! in up state ( kaux == 0 ), we need to return immediately and the
     ! ratio is one
     if ( kaux == 0 ) then
         ratio = one; RETURN
     endif ! back if ( kaux == 0 ) block

     ! allocate memory
     allocate(Dmm(kaux,kaux), stat=istat)
     allocate(Tmm(kaux,kaux), stat=istat)

     if ( istat /= 0 ) then
         call s_print_error('cat_reflip_detrat','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

     ! init Dmm and Tmm matrix
     Dmm = zero
     Tmm = zero

     ! recalculate Dmm from scratch
     do j=1,kaux
         tau_end = time_e(index_e(j, up), up)
         do i=1,kaux
             tau_start = time_s(index_s(i, up), up)
             if ( tau_start < tau_end ) then
                 Dmm(i, j) = -ctqmc_eval_htau(dn, tau_start - tau_end + beta)
             else
                 Dmm(i, j) =  ctqmc_eval_htau(dn, tau_start - tau_end)
             endif ! back if ( tau_start < tau_end ) block
         enddo ! over i={1,kaux} loop
     enddo ! over j={1,kaux} loop

     ! calculate Tmm: Tmm = Dmm . Mmat
     call dgemm('N', 'N', kaux, kaux, kaux, one, Dmm, kaux, mmat(1:kaux, 1:kaux, up), kaux, zero, Tmm, kaux)

     ! calculate the determinant of Tmm, it is the desired ratio
     call s_det_d(kaux, Tmm, ratio)

     ! deallocate memory
     deallocate(Dmm)
     deallocate(Tmm)

!! body]

     return
  end subroutine cat_reflip_detrat

!!
!! @sub cat_reload_detrat
!!
!! to do nothing, it is an empty subroutine
!!
  subroutine cat_reload_detrat()
     implicit none

!! [body

     CONTINUE

!! body]

     return
  end subroutine cat_reload_detrat
