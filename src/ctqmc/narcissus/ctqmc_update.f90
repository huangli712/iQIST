!!!-----------------------------------------------------------------------
!!! project : narcissus
!!! program : ctqmc_warming
!!!           ctqmc_walking
!!!           ctqmc_warning
!!!           ctqmc_insert_kink
!!!           ctqmc_remove_kink
!!!           ctqmc_lshift_kink
!!!           ctqmc_rshift_kink
!!!           ctqmc_reflip_kink
!!!           ctqmc_reload_kink
!!! source  : ctqmc_update.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           04/25/2017 by li huang (last modified)
!!! purpose : basic update actions for the hybridization expansion version
!!!           continuous time quantum Monte Carlo (CTQMC) quantum impurity
!!!           solver. they are called by ctqmc_impurity_solver()
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> service layer: updating perturbation expansion series 1          <<<
!!========================================================================

!!
!! @sub ctqmc_warming
!!
!! perform thermalization or warmup on the perturbation expansion series
!! to achieve thermodynamics stable equilibrium state
!!
  subroutine ctqmc_warming()
     use constants, only : zero

     use control, only : ntherm
     use context, only : ins_t, ins_a, ins_r
     use context, only : rmv_t, rmv_a, rmv_r
     use context, only : lsh_t, lsh_a, lsh_r
     use context, only : rsh_t, rsh_a, rsh_r
     use context, only : rfl_t, rfl_a, rfl_r

     implicit none

! local variables
! loop index
     integer :: i

! warm up the diagram series
     do i=1,ntherm
         call ctqmc_walking(i)
     enddo ! over i={1,ntherm} loop

! reinit statistics variables
     ins_t = zero; ins_a = zero; ins_r = zero
     rmv_t = zero; rmv_a = zero; rmv_r = zero
     lsh_t = zero; lsh_a = zero; lsh_r = zero
     rsh_t = zero; rsh_a = zero; rsh_r = zero
     rfl_t = zero; rfl_a = zero; rfl_r = zero

     return
  end subroutine ctqmc_warming

!!
!! @sub ctqmc_walking
!!
!! visit the perturbation expansion diagrams randomly
!!
  subroutine ctqmc_walking(cstep)
     use constants, only : dp
     use spring, only : spring_sfmt_stream

     use control, only : nflip, nclean

     implicit none

! external arguments
! current QMC sweep steps
     integer, intent(in) :: cstep

! change the order of perturbation expansion series
     if ( spring_sfmt_stream() < 0.9_dp ) then
         if ( spring_sfmt_stream() > 0.5_dp ) then
             call ctqmc_insert_kink()  ! insert one new kink
         else
             call ctqmc_remove_kink()  ! remove one old kink
         endif ! back if ( spring_sfmt_stream() > 0.5_dp ) block
! do not change the order of perturbation expansion series
     else
         if ( spring_sfmt_stream() > 0.5_dp ) then
             call ctqmc_lshift_kink()  ! shift the left  endpoints
         else
             call ctqmc_rshift_kink()  ! shift the right endpoints
         endif ! back if ( spring_sfmt_stream() > 0.5_dp ) block
     endif ! back if ( spring_sfmt_stream() < 0.9_dp ) block

! numerical trick: perform global spin flip periodically
     if ( nflip > 0  .and. mod(cstep, +nflip) == 0 ) then
         if ( spring_sfmt_stream() < 0.9_dp ) then
             call ctqmc_reflip_kink(1) ! flip intra-orbital spins one by one
         else
             call ctqmc_reflip_kink(2) ! flip intra-orbital spins globally
         endif ! back if ( spring_sfmt_stream() < 0.9_dp ) block
     endif ! back if ( nflip > 0  .and. mod(cstep, +nflip) == 0 ) block

     if ( nflip < 0  .and. mod(cstep, -nflip) == 0 ) then
         if ( spring_sfmt_stream() > 0.9_dp ) then
             call ctqmc_reflip_kink(1) ! flip intra-orbital spins one by one
         else
             call ctqmc_reflip_kink(2) ! flip intra-orbital spins globally
         endif ! back if ( spring_sfmt_stream() > 0.9_dp ) block
     endif ! back if ( nflip < 0  .and. mod(cstep, -nflip) == 0 ) block

! numerical trick: perform global update periodically
     if ( nclean > 0 .and. mod(cstep, nclean) == 0 ) then
         call ctqmc_reload_kink()
     endif ! back if ( nclean > 0 .and. mod(cstep, nclean) == 0 ) block

     return
  end subroutine ctqmc_walking

!!
!! @sub ctqmc_warning
!!
!! checking whether the quantum impurity solver is consistent internally
!!
  subroutine ctqmc_warning(cflag)
     use constants, only : mystd

     use control, only : norbs
     use control, only : myid, master
     use context, only : index_s, index_e
     use context, only : time_s, time_e
     use context, only : rank, stts

     implicit none

! external arguments
! control flag
     integer, intent(inout) :: cflag

! local variables
! loop index
     integer :: i
     integer :: j

     if ( cflag == 1 ) then ! if cflag /= 1, nothing to do

! check perturbation expansion order
         do i=1,norbs
             if ( stts(i) == 0 ) then
                 if ( rank(i) /= 0 ) cflag = 99
             endif ! back if ( stts(i) == 0 ) block

             if ( stts(i) == 1 ) then
                 if ( rank(i) == 0 ) cflag = 99
             endif ! back if ( stts(i) == 1 ) block

             if ( stts(i) == 2 ) then
                 if ( rank(i) == 0 ) cflag = 99
             endif ! back if ( stts(i) == 2 ) block

             if ( stts(i) == 3 ) then
                 if ( rank(i) /= 0 ) cflag = 99
             endif ! back if ( stts(i) == 3 ) block
         enddo ! over i={1,norbs} loop

! check time order of operators
         do i=1,norbs
             do j=1,rank(i)-1
                 if ( time_s( index_s(j, i), i ) > time_s( index_s(j+1, i), i ) ) then
                     cflag = 99
                 endif ! back if ( time_s( index_s(j, i), i ) > time_s( index_s(j+1, i), i ) ) block
                 if ( time_e( index_e(j, i), i ) > time_e( index_e(j+1, i), i ) ) then
                     cflag = 99
                 endif ! back if ( time_e( index_e(j, i), i ) > time_e( index_e(j+1, i), i ) ) block
             enddo ! over j={1,rank(i)-1} loop
         enddo ! over i={1,norbs} loop

! check segment and anti-segment
         do i=1,norbs
             if ( stts(i) == 1 ) then
                 if ( time_s( index_s(1, i), i ) > time_e( index_e(1, i), i ) ) then
                     cflag = 99
                 endif ! back if ( time_s( index_s(1, i), i ) > time_e( index_e(1, i), i ) ) block
             endif ! back if ( stts(i) == 1 ) block

             if ( stts(i) == 2 ) then
                 if ( time_s( index_s(1, i), i ) < time_e( index_e(1, i), i ) ) then
                     cflag = 99
                 endif ! back if ( time_s( index_s(1, i), i ) < time_e( index_e(1, i), i ) ) block
             endif ! back if ( stts(i) == 2 ) block
         enddo ! over i={1,norbs} loop

! write the results, only master node can do it
         if ( myid == master ) then
             if ( cflag == 99 ) then
                 write(mystd,'(4X,a)') '>>> quantum impurity solver status: fatel error'
                 write(mystd,'(4X,a)') '>>> please check the status file: solver.status.dat'
                 call ctqmc_save_status()
                 call s_print_error('ctqmc_warning','unknown fatal error occur')
             else
                 write(mystd,'(4X,a)') '>>> quantum impurity solver status: normal'
             endif ! back if ( cflag == 99 ) block
         endif ! back if ( myid == master ) block

     endif ! back if ( cflag == 1 ) block

     return
  end subroutine ctqmc_warning

!!========================================================================
!!>>> service layer: updating perturbation expansion series 2          <<<
!!========================================================================

!!
!! @sub ctqmc_insert_kink
!!
!! insert new segment or anti-segment in the perturbation expansion series
!!
  subroutine ctqmc_insert_kink()
     use constants, only : dp, zero, one
     use spring, only : spring_sfmt_stream

     use control, only : norbs
     use control, only : mkink
     use control, only : beta
     use context, only : ckink, cstat
     use context, only : ins_t, ins_a, ins_r
     use context, only : rank, stts

     implicit none

! local variables
! whether the new segment or anti-segment can be inserted diagrammatically
     logical  :: ladd

! whether it is an anti-segment
! anti = .true., anti-segment
! anti = .false., segment
     logical  :: anti

! whether the update operation is accepted
     logical  :: pass

! current flavor channel for both band and spin
     integer  :: flvr

! index address to insert new segment or anti-segment
! is and ie are for start and end points, respectively
     integer  :: is, ie

! transition probability for insert new segment or anti-segment
     real(dp) :: p

! \tau_s, start point of the new segment
     real(dp) :: tau_start

! \tau_e, end   point of the new segment
     real(dp) :: tau_end

! the possible maximum length of the new segment
     real(dp) :: tau_max

! ratio between old and new configurations, the local trace part
     real(dp) :: trace_ratio

! ratio between old and new configurations, the determinant part
     real(dp) :: deter_ratio

! initialize logical variables
     ladd = .false.
     anti = .false.
     pass = .false.

! select the flavor channel randomly among 1 ~ norbs
     flvr = ceiling( spring_sfmt_stream() * norbs )

! get the perturbation expansion order ( number of existing segments or
! anti-segments ) for current flavor channel
     ckink = rank(flvr)
     if ( ckink == mkink ) then
!<         call s_print_exception('ctqmc_insert_kink','can not insert any segments')
         ins_t = ins_t + one
         ins_r = ins_r + one
         RETURN
     endif ! back if ( ckink == mkink ) block

! randomly choose anti and tau_start, and then check whether tau_start is
! valid, if tau_start is valid, then determine tau_end, tau_max, is, and
! ie consistently, and set ladd to .true., if tau_start is not valid, then
! set ladd to .false.
     call cat_insert_flavor(flvr, is, ie, anti, ladd, tau_start, tau_end, tau_max)

! calculate the transition ratio between old and new configurations,
! for the local trace part
     if ( ladd .eqv. .true. ) then
         call cat_insert_ztrace(flvr, anti, tau_start, tau_end, trace_ratio)
     else
         trace_ratio = zero
     endif ! back if ( ladd .eqv. .true. ) block

! calculate the transition ratio between old and new configurations,
! for the determinant part
     if ( ladd .eqv. .true. ) then
         call cat_insert_detrat(flvr, tau_start, tau_end, deter_ratio)
     else
         deter_ratio = zero
     endif ! back if ( ladd .eqv. .true. ) block

! calculate the transition probability for insert new segment or anti-segment
     p = deter_ratio * trace_ratio * ( beta * tau_max / real( ckink + 1 ) )

! determine pass, using important sampling algorithm (metropolis algorithm)
     pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if the update action is accepted
     if ( pass .eqv. .true. ) then

! update the mmat matrix and gmat matrix, respectively,
! cat_insert_action() subroutine is invoked internally to update the
! perturbation expansion series
         call cat_insert_matrix(flvr, is, ie, tau_start, tau_end, deter_ratio)

! update ckink for current flavor channel
         ckink = ckink + 1

! update stts for current flavor channel
         stts(flvr) = cstat

! update rank for current flavor channel
         rank(flvr) = rank(flvr) + 1

     endif ! back if ( pass .eqv. .true. ) block

! update the insert statistics
     ins_t = ins_t + one
     if ( pass .eqv. .true. ) then
         ins_a = ins_a + one
     else
         ins_r = ins_r + one
     endif ! back if ( pass .eqv. .true. ) block

     return
  end subroutine ctqmc_insert_kink

!!
!! @sub ctqmc_remove_kink
!!
!! remove old segment or anti-segment in the perturbation expansion series
!!
  subroutine ctqmc_remove_kink()
     use constants, only : dp, one
     use spring, only : spring_sfmt_stream

     use control, only : norbs
     use control, only : beta
     use context, only : ckink, cstat
     use context, only : rmv_t, rmv_a, rmv_r
     use context, only : rank, stts

     implicit none

! local variables
! whether it is an anti-segment
! anti = .true., anti-segment
! anti = .false., segment
     logical  :: anti

! whether the update operation is accepted
     logical  :: pass

! current flavor channel for both band and spin
     integer  :: flvr

! index address to remove old segment or anti-segment
! is and ie are for start and end points, respectively
     integer  :: is, ie

! transition probability for remove old segment or anti-segment
     real(dp) :: p

! \tau_s, start point of the old segment
     real(dp) :: tau_start

! \tau_e, end   point of the old segment
     real(dp) :: tau_end

! the possible maximum length of the old segment
     real(dp) :: tau_max

! ratio between old and new configurations, the local trace part
     real(dp) :: trace_ratio

! ratio between old and new configurations, the determinant part
     real(dp) :: deter_ratio

! initialize logical variables
     anti = .false.
     pass = .false.

! select the flavor channel randomly among 1 ~ norbs
     flvr = ceiling( spring_sfmt_stream() * norbs )

! get the perturbation expansion order ( number of existing segments or
! anti-segments ) for current flavor channel
     ckink = rank(flvr)
     if ( ckink == 0 ) then
!<         call s_print_exception('ctqmc_remove_kink','can not remove any segments')
         rmv_t = rmv_t + one
         rmv_r = rmv_r + one
         RETURN
     endif ! back if ( ckink == 0 ) block

! at first determine anti and is randomly, then tau_start is obtained by
! is. and then ie, tau_end, and tau_max are evaluated carefully according
! to is and ie
     call cat_remove_flavor(flvr, is, ie, anti, tau_start, tau_end, tau_max)

! calculate the transition ratio between old and new configurations,
! for the local trace part
     call cat_remove_ztrace(flvr, anti, tau_start, tau_end, trace_ratio)

! calculate the transition ratio between old and new configurations,
! for the determinant part
     call cat_remove_detrat(flvr, is, ie, deter_ratio)

! calculate the transition probability for remove old segment or anti-segment
     p = deter_ratio * trace_ratio * real( ckink ) / ( beta * tau_max )

! determine pass, using important sampling algorithm (metropolis algorithm)
     pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if update action is accepted
     if ( pass .eqv. .true. ) then

! update the mmat matrix and gmat matrix, respectively,
! cat_remove_action() subroutine is invoked internally to update the
! perturbation expansion series
         call cat_remove_matrix(flvr, is, ie)

! decrease ckink for current flavor channel
         ckink = ckink - 1

! update stts for current flavor channel
         stts(flvr) = cstat

! update rank for current flavor channel
         rank(flvr) = rank(flvr) - 1

     endif ! back if ( pass .eqv. .true. ) block

! update the remove statistics
     rmv_t = rmv_t + one
     if ( pass .eqv. .true. ) then
         rmv_a = rmv_a + one
     else
         rmv_r = rmv_r + one
     endif ! back if ( pass .eqv. .true. ) block

     return
  end subroutine ctqmc_remove_kink

!!
!! @sub ctqmc_lshift_kink
!!
!! left shift old segment or anti-segment in the perturbation expansion series
!!
  subroutine ctqmc_lshift_kink()
     use constants, only : dp, one
     use spring, only : spring_sfmt_stream

     use control, only : norbs
     use context, only : ckink, cstat
     use context, only : lsh_t, lsh_a, lsh_r
     use context, only : rank, stts

     implicit none

! local variables
! whether the update operation winds around the circle
     logical  :: ring

! whether the update operation is accepted
     logical  :: pass

! current flavor channel for both band and spin
     integer  :: flvr

! index address to left shift old segment or anti-segment
! iso and isn are for old and new indices, respectively
     integer  :: iso, isn

! transition probability for left shift old segment or anti-segment
     real(dp) :: p

! \tau_s, start point of the old segment (old point)
     real(dp) :: tau_start1

! \tau_s, start point of the old segment (new point)
     real(dp) :: tau_start2

! ratio between old and new configurations, the local trace part
     real(dp) :: trace_ratio

! ratio between old and new configurations, the determinant part
     real(dp) :: deter_ratio

! initialize logical variables
     ring = .false.
     pass = .false.

! select the flavor channel randomly among 1 ~ norbs
     flvr = ceiling( spring_sfmt_stream() * norbs )

! get the perturbation expansion order ( number of existing segments or
! anti-segments ) for current flavor channel
     ckink = rank(flvr)
     if ( ckink == 0 ) then
!<         call s_print_exception('ctqmc_lshift_kink','can not lshift any segments')
         lsh_t = lsh_t + one
         lsh_r = lsh_r + one
         RETURN
     endif ! back if ( ckink == 0 ) block

! at first, we select iso randomly, and then obtain tau_start1. according
! to the existing segments, we determine tau_start2 and related index isn,
! finally ring is evaluated.
     call cat_lshift_flavor(flvr, iso, isn, ring, tau_start1, tau_start2)

! calculate the transition ratio between old and new configurations,
! for the local trace part
     call cat_lshift_ztrace(flvr, ring, tau_start1, tau_start2, trace_ratio)

! calculate the transition ratio between old and new configurations,
! for the determinant part
     call cat_lshift_detrat(flvr, iso, tau_start1, tau_start2, deter_ratio)

! calculate the transition probability for left shift old segment or anti-segment
     p = deter_ratio * trace_ratio

! determine pass, using important sampling algorithm (metropolis algorithm)
     pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if update action is accepted
     if ( pass .eqv. .true. ) then

! update the mmat matrix and gmat matrix, respectively,
! cat_lshift_action() subroutine is invoked internally to update the
! perturbation expansion series
         call cat_lshift_matrix(flvr, iso, isn, tau_start1, tau_start2, deter_ratio)

! update stts for current flavor channel
         stts(flvr) = cstat

     endif ! back if ( pass .eqv. .true. ) block

! update the lshift statistics
     lsh_t = lsh_t + one
     if ( pass .eqv. .true. ) then
         lsh_a = lsh_a + one
     else
         lsh_r = lsh_r + one
     endif ! back if ( pass .eqv. .true. ) block

     return
  end subroutine ctqmc_lshift_kink

!!
!! @sub ctqmc_rshift_kink
!!
!! right shift old segment or anti-segment in the perturbation expansion series
!!
  subroutine ctqmc_rshift_kink()
     use constants, only : dp, one
     use spring, only : spring_sfmt_stream

     use control, only : norbs
     use context, only : ckink, cstat
     use context, only : rsh_t, rsh_a, rsh_r
     use context, only : rank, stts

     implicit none

! local variables
! whether the update operation winds around the circle
     logical  :: ring

! whether the update operation is accepted
     logical  :: pass

! current flavor channel for both band and spin
     integer  :: flvr

! index address to right shift old segment or anti-segment
! ieo and ien are for old and new indices, respectively
     integer  :: ieo, ien

! transition probability for right shift old segment or anti-segment
     real(dp) :: p

! \tau_e, end   point of the old segment (old point)
     real(dp) :: tau_end1

! \tau_e, end   point of the old segment (new point)
     real(dp) :: tau_end2

! ratio between old and new configurations, the local trace part
     real(dp) :: trace_ratio

! ratio between old and new configurations, the determinant part
     real(dp) :: deter_ratio

! initialize logical variables
     ring = .false.
     pass = .false.

! select the flavor channel randomly among 1 ~ norbs
     flvr = ceiling( spring_sfmt_stream() * norbs )

! get the perturbation expansion order ( number of existing segments or
! anti-segments ) for current flavor channel
     ckink = rank(flvr)
     if ( ckink == 0 ) then
!<         call s_print_exception('ctqmc_rshift_kink','can not rshift any segments')
         rsh_t = rsh_t + one
         rsh_r = rsh_r + one
         RETURN
     endif ! back if ( ckink == 0 ) block

! at first, we select ieo randomly, and then obtain tau_end1. according
! to the existing segments, we determine tau_end2 and related index ien,
! finally ring is evaluated.
     call cat_rshift_flavor(flvr, ieo, ien, ring, tau_end1, tau_end2)

! calculate the transition ratio between old and new configurations,
! for the local trace part
     call cat_rshift_ztrace(flvr, ring, tau_end1, tau_end2, trace_ratio)

! calculate the transition ratio between old and new configurations,
! for the determinant part
     call cat_rshift_detrat(flvr, ieo, tau_end1, tau_end2, deter_ratio)

! calculate the transition probability for right shift old segment or anti-segment
     p = deter_ratio * trace_ratio

! determine pass, using important sampling algorithm (metropolis algorithm)
     pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if update action is accepted
     if ( pass .eqv. .true. ) then

! update the mmat matrix and gmat matrix, respectively,
! cat_rshift_action() subroutine is invoked internally to update the
! perturbation expansion series
         call cat_rshift_matrix(flvr, ieo, ien, tau_end1, tau_end2, deter_ratio)

! update stts for current flavor channel
         stts(flvr) = cstat

     endif ! back if ( pass .eqv. .true. ) block

! update the rshift statistics
     rsh_t = rsh_t + one
     if ( pass .eqv. .true. ) then
         rsh_a = rsh_a + one
     else
         rsh_r = rsh_r + one
     endif ! back if ( pass .eqv. .true. ) block

     return
  end subroutine ctqmc_rshift_kink

!!
!! @sub ctqmc_reflip_kink
!!
!! perform a global update, exchange the states between spin up and spin
!! down, it maybe useful for magnetic systems
!!
  subroutine ctqmc_reflip_kink(cflip)
     use constants, only : dp, one
     use spring, only : spring_sfmt_stream

     use control, only : nband, norbs
     use context, only : rfl_t, rfl_a, rfl_r
     use context, only : rank, symm

     implicit none

! external arguments
! control flag
! if cflip = 1, flip inter-orbital spins randomly;
! if cflip = 2, flip intra-orbital spins one by one;
! if cflip = 3, flip intra-orbital spins globally
     integer, intent(in) :: cflip

! local variables
! whether the update operation is accepted
     logical  :: pass

! selected flavor pairs
     integer  :: fup, fdn

! loop index for flavor channel
     integer  :: flvr

! maximum rank order
     integer  :: kmax

! transition probability for global spin flip
     real(dp) :: p

! global flip determinant ratio
     real(dp) :: ratup
     real(dp) :: ratdn

! initialize logical variables
     pass = .false.

! initialize transition probability
     p = one

     if ( cflip == 1 ) then
! determine fup and fdn, and fup /= fdn
         fup = ceiling( spring_sfmt_stream() * norbs )
         do while ( .true. )
             fdn = ceiling( spring_sfmt_stream() * norbs )
             if ( fup /= fdn .and. symm(fup) == symm(fdn) ) EXIT
         enddo ! over do while loop

! calculate the transition ratio between old and new configurations,
! for the determinant part, spin up case
         call cat_reflip_detrat(fup, fdn, ratup)

! calculate the transition ratio between old and new configurations,
! for the determinant part, spin dn case
         call cat_reflip_detrat(fdn, fup, ratdn)

! calculate the transition probability for global spin flip
         p = ratup * ratdn

! determine pass, using important sampling algorithm (metropolis algorithm)
         pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if update action is accepted
         if ( pass .eqv. .true. ) then

! get maximum rank order in spin up and spin down states
             kmax = max( rank(fup), rank(fdn) )

! swap global variables between spin up and spin down states
             call cat_reflip_matrix(fup, fdn, kmax)

         endif ! back if ( pass .eqv. .true. ) block

! update the reflip statistics
         reflip_tcount = reflip_tcount + one
         if ( pass .eqv. .true. ) then
             reflip_accept = reflip_accept + one
         else
             reflip_reject = reflip_reject + one
         endif ! back if ( pass .eqv. .true. ) block

     else if ( cflip == 2 ) then ! cflip = 2, local flip
         do flvr=1,nband

! get fup and fdn
             fup = flvr; fdn = flvr + nband

! calculate the transition ratio between old and new configurations,
! for the determinant part, spin up case
             call cat_reflip_detrat(fup, fdn, ratup)

! calculate the transition ratio between old and new configurations,
! for the determinant part, spin dn case
             call cat_reflip_detrat(fdn, fup, ratdn)

! calculate the transition probability for global spin flip
             p = ratup * ratdn

! determine pass, using important sampling algorithm (metropolis algorithm)
             pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if update action is accepted
             if ( pass .eqv. .true. ) then

! get maximum rank order in spin up and spin down states
                 kmax = max( rank(fup), rank(fdn) )

! swap global variables between spin up and spin down states
                 call cat_reflip_matrix(fup, fdn, kmax)

             endif ! back if ( pass .eqv. .true. ) block

! update the reflip statistics
             reflip_tcount = reflip_tcount + one
             if ( pass .eqv. .true. ) then
                 reflip_accept = reflip_accept + one
             else
                 reflip_reject = reflip_reject + one
             endif ! back if ( pass .eqv. .true. ) block

         enddo ! over flvr={1,nband} loop

     else if ( cflip == 3 ) then ! cflip = 3, global flip
         do flvr=1,nband

! get fup and fdn
             fup = flvr; fdn = flvr + nband

! calculate the transition ratio between old and new configurations,
! for the determinant part, spin up case
             call cat_reflip_detrat(fup, fdn, ratup)

! calculate the transition ratio between old and new configurations,
! for the determinant part, spin dn case
             call cat_reflip_detrat(fdn, fup, ratdn)

! calculate the transition probability for global spin flip
             p = p * ( ratup * ratdn )

         enddo ! over flvr={1,nband} loop

! determine pass, using important sampling algorithm (metropolis algorithm)
         pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if update action is accepted
         if ( pass .eqv. .true. ) then

             do flvr=1,nband

! get fup and fdn
                 fup = flvr; fdn = flvr + nband

! get maximum rank order in spin up and spin down states
                 kmax = max( rank(fup), rank(fdn) )

! swap global variables between spin up and spin down states
                 call cat_reflip_matrix(fup, fdn, kmax)

             enddo ! over flvr={1,nband} loop

         endif ! back if ( pass .eqv. .true. ) block

! update the reflip statistics
         reflip_tcount = reflip_tcount + one
         if ( pass .eqv. .true. ) then
             reflip_accept = reflip_accept + one
         else
             reflip_reject = reflip_reject + one
         endif ! back if ( pass .eqv. .true. ) block

     endif ! back if ( cflip == 1 ) block

     return
  end subroutine ctqmc_reflip_kink

!!
!! @sub ctqmc_reload_kink
!!
!! global update all segments or anti-segments in the perturbation
!! expansion series
!!
  subroutine ctqmc_reload_kink()
     use control, only : norbs
     use context, only : rank

     implicit none

! local variables
! loop index for flavor channel
     integer :: flvr

     do flvr=1,norbs

! check the perturbation expansion order ( number of existing segments or
! anti-segments ) for current flavor channel
         if ( rank(flvr) == 0 ) CYCLE

! regenerate the mmat matrix and gmat matrix from scratch
         call cat_reload_matrix(flvr)

     enddo ! over flvr={1,norbs} loop

     return
  end subroutine ctqmc_reload_kink
