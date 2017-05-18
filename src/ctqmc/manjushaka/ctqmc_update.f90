!!!-----------------------------------------------------------------------
!!! project : manjushaka
!!! program : ctqmc_try_warming
!!!           ctqmc_try_walking
!!!           ctqmc_try_warning <<<---
!!!           ctqmc_insert_kink
!!!           ctqmc_remove_kink
!!!           ctqmc_lshift_kink
!!!           ctqmc_rshift_kink
!!!           ctqmc_reflip_kink
!!!           ctqmc_reload_kink <<<---
!!! source  : ctqmc_update.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           05/18/2017 by li huang (last modified)
!!! purpose : basic update actions for the hybridization expansion version
!!!           continuous time quantum Monte Carlo (CTQMC) quantum impurity
!!!           solver. they are called by ctqmc_impurity_solver().
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!========================================================================
!!>>> service layer: updating perturbation expansion series 1          <<<
!!========================================================================

!!
!! @sub ctqmc_try_warming
!!
!! perform thermalization or warmup on the perturbation expansion series
!! to achieve thermodynamics stable equilibrium state
!!
  subroutine ctqmc_try_warming()
     use constants, only : zero

     use control, only : ntherm

     use context, only : cnegs, caves
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
         call ctqmc_try_walking(i)
     enddo ! over i={1,ntherm} loop

! reset cnegs
     cnegs = 0

! reset caves
     caves = 0

! reinit statistics variables
     ins_t = zero
     ins_a = zero
     ins_r = zero

     rmv_t = zero
     rmv_a = zero
     rmv_r = zero

     lsh_t = zero
     lsh_a = zero
     lsh_r = zero

     rsh_t = zero
     rsh_a = zero
     rsh_r = zero

     rfl_t = zero
     reflip_accept = zero
     reflip_reject = zero

     return
  end subroutine ctqmc_try_warming

  subroutine ctqmc_try_walking(cstep)
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
             call ctqmc_lshift_kink()  ! shift the create  operators
         else
             call ctqmc_rshift_kink()  ! shift the destroy operators
         endif ! back if ( spring_sfmt_stream() > 0.5_dp ) block
     endif ! back if ( spring_sfmt_stream() < 0.9_dp ) block

! numerical trick: perform global spin flip periodically
     if ( nflip > 0  .and. mod(cstep, +nflip) == 0 ) then
         if ( spring_sfmt_stream() < 0.8_dp ) then
             call ctqmc_reflip_kink(2) ! flip intra-orbital spins one by one
         else
             call ctqmc_reflip_kink(3) ! flip intra-orbital spins globally
         endif ! back if ( spring_sfmt_stream() < 0.8_dp ) block
     endif ! back if ( nflip > 0  .and. mod(cstep, +nflip) == 0 ) block

     if ( nflip < 0  .and. mod(cstep, -nflip) == 0 ) then
         if ( spring_sfmt_stream() < 0.8_dp ) then
             call ctqmc_reflip_kink(1) ! flip inter-orbital spins randomly
         else
             call ctqmc_reflip_kink(3) ! flip intra-orbital spins globally
         endif ! back if ( spring_sfmt_stream() < 0.8_dp ) block
     endif ! back if ( nflip < 0  .and. mod(cstep, -nflip) == 0 ) block

! numerical trick: perform global update periodically
     if ( nclean > 0 .and. mod(cstep, nclean) == 0 ) then
         call ctqmc_reload_kink()
     endif ! back if ( nclean > 0 .and. mod(cstep, nclean) == 0 ) block

     return
  end subroutine ctqmc_try_walking

  subroutine ctqmc_try_warning(cflag)
     use constants, only : mystd

     use control, only : norbs
     use control, only : myid, master
     use context, only : index_s, index_e, time_s, time_e
     use context, only : index_v, time_v
     use context, only : rank

     implicit none

! external arguments
! control flag
     integer, intent(inout) :: cflag

! local variables
! loop index
     integer :: i
     integer :: j

     if ( cflag == 1 ) then

! check time order of operators in colour part
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

! check time order of operators in flavor part
         do j=1,2*sum(rank)-1
             if ( index_v(j) <= 0 .or. index_v(j+1) <= 0 ) then
                 cflag = 99
             endif ! back if ( index_v(j) <= 0 .or. index_v(j+1) <= 0 ) block
         enddo ! over j={1,2*sum(rank)-1} loop

         do j=1,2*sum(rank)-1
             if ( time_v( index_v(j) ) > time_v( index_v(j+1) ) ) then
                 cflag = 99
             endif ! back if ( time_v( index_v(j) ) > time_v( index_v(j+1) ) ) block
         enddo ! over j={1,2*sum(rank)-1} loop

! write the results, only master node can do it
         if ( myid == master ) then
             if ( cflag == 99 ) then
                 write(mystd,'(4X,a)') '>>> quantum impurity solver status: error?'
                 write(mystd,'(4X,a)') '>>> please check the status file: solver.status.dat'
                 call ctqmc_save_status()
                 call s_print_error('ctqmc_try_warning','unknown fatal error occur')
             else
                 write(mystd,'(4X,a)') '>>> quantum impurity solver status: normal'
             endif ! back if ( cflag == 99 ) block
         endif ! back if ( myid == master ) block

     endif ! back if ( cflag == 1 ) block

     return
  end subroutine ctqmc_try_warning

!!========================================================================
!!>>> driver layer: updating perturbation expansion series             <<<
!!========================================================================

!!>>> ctqmc_insert_kink: insert new create and destroy operators in the
!!>>> perturbation expansion series
  subroutine ctqmc_insert_kink()
     use constants, only : dp, zero, one
     use spring, only : spring_sfmt_stream

     use control, only : norbs
     use control, only : mkink
     use control, only : beta
     use context, only : ckink, csign, cnegs
     use context, only : ins_t, ins_a, ins_r
     use context, only : rank

     implicit none

! local variables
! whether the new create and destroy operators can be inserted diagrammatically
     logical  :: ladd

! whether the update operation is accepted
     logical  :: pass

! current flavor channel for both band and spin
     integer  :: flvr

! index address to insert new create and destroy operators
! is and ie are for create and destroy operators, respectively
! cis and cie are for the colour part, and fis and fie are for the flavor part
     integer  :: cis, cie
     integer  :: fis, fie

! transition probability for insert new create and destroy operators
     real(dp) :: p

! random number
     real(dp) :: r

! \tau_s, imaginary time point of the create operator
     real(dp) :: tau_start

! \tau_e, imaginary time point of the destroy operator
     real(dp) :: tau_end

! ratio between old and new configurations, the local trace part
     real(dp) :: trace_ratio

! ratio between old and new configurations, the determinant part
     real(dp) :: deter_ratio

! initialize logical variables
     ladd = .false.
     pass = .false.

! select the flavor channel randomly among 1 ~ norbs
     flvr = ceiling( spring_sfmt_stream() * norbs )

! get the perturbation expansion order ( number of existing create or
! destroy operators ) for current flavor channel
     ckink = rank(flvr)
     if ( ckink == mkink ) then
!<         call s_print_exception('ctqmc_insert_kink','can not insert any operators')
         ins_t = ins_t + one
         ins_r = ins_r + one
         if ( csign < 0 )  cnegs = cnegs + 1
         RETURN
     endif ! back if ( ckink == mkink ) block

! randomly generate tau_start and tau_end at selected flvr channel, and
! then determine index address cis and cie for them
     call try_insert_colour(flvr, cis, cie, tau_start, tau_end)

! fast look up the flavor part of perturbation expansion series, determine
! corresponding fis and fie, and determine whether the operators trace is
! not equal to zero
     call try_insert_flavor(flvr, fis, fie, tau_start, tau_end, ladd)

! calculate the transition ratio between old and new configurations,
! for the local trace part
     if ( ladd .eqv. .true. ) then
         call cat_insert_ztrace(flvr, fis, fie, tau_start, tau_end, trace_ratio)
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

! we will determine the pass by lazy trace evalution
! if ladd is false, we set the pass as false immediately
     r = spring_sfmt_stream()
     trace_ratio = deter_ratio * ( beta / real( ckink + 1 ) ) ** 2
     if ( ladd .eqv. .true. ) then
         call ctqmc_lazy_ztrace( 1, 2*sum(rank) + 2, trace_ratio, tau_start, tau_end, r, p, pass )
     else
         pass = .false.
     endif ! back if ( ladd .eqv. .true. ) block

! if the update action is accepted
     if ( pass .eqv. .true. ) then

! update the mmat matrix and gmat matrix, respectively,
! cat_insert_colour() subroutine is invoked internally to update the colour
! part of perturbation expansion series
         call cat_insert_matrix(flvr, cis, cie, tau_start, tau_end, deter_ratio)

! update the flavor part of perturbation expansion series
         call cat_insert_flavor(flvr, fis, fie, tau_start, tau_end)

! update the operators trace
         call ctqmc_make_evolve()

! update ckink for current flavor channel
         ckink = ckink + 1

! update rank for current flavor channel
         rank(flvr) = rank(flvr) + 1

! determine the sign, TO BE CHECKED
         csign = csign * int ( sign(one, p) )

     endif ! back if ( pass .eqv. .true. ) block

! record negative sign
     if ( csign < 0 ) then
         cnegs = cnegs + 1
!<         call s_print_exception('ctqmc_insert_kink','csign is negative')
     endif ! back if ( csign < 0 ) block

! update the insert statistics
     ins_t = ins_t + one
     if ( pass .eqv. .true. ) then
         ins_a = ins_a + one
     else
         ins_r = ins_r + one
     endif ! back if ( pass .eqv. .true. ) block

     return
  end subroutine ctqmc_insert_kink

!!>>> ctqmc_remove_kink: remove old create and destroy operators in the
!!>>> perturbation expansion series
  subroutine ctqmc_remove_kink()
     use constants, only : dp, zero, one
     use spring, only : spring_sfmt_stream

     use control, only : norbs
     use control, only : beta
     use context, only : ckink, csign, cnegs
     use context, only : rmv_t, rmv_a, rmv_r
     use context, only : rank

     implicit none

! local variables
! whether the old create and destroy operators can be removed diagrammatically
     logical  :: lrmv

! whether the update operation is accepted
     logical  :: pass

! current flavor channel for both band and spin
     integer  :: flvr

! index address to remove old create and destroy operators
! is and ie are for create and destroy operators, respectively
! cis and cie are for the colour part, and fis and fie are for the flavor part
     integer  :: cis, cie
     integer  :: fis, fie

! transition probability for remove old create and destroy operators
     real(dp) :: p

! random number
     real(dp) :: r

! \tau_s, imaginary time point of the create operator
     real(dp) :: tau_start

! \tau_e, imaginary time point of the destroy operator
     real(dp) :: tau_end

! ratio between old and new configurations, the local trace part
     real(dp) :: trace_ratio

! ratio between old and new configurations, the determinant part
     real(dp) :: deter_ratio

! initialize logical variables
     lrmv = .false.
     pass = .false.

! select the flavor channel randomly among 1 ~ norbs
     flvr = ceiling( spring_sfmt_stream() * norbs )

! get the perturbation expansion order ( number of existing create or
! destroy operators ) for current flavor channel
     ckink = rank(flvr)
     if ( ckink == 0 ) then
!<         call s_print_exception('ctqmc_remove_kink','can not remove any operators')
         rmv_t = rmv_t + one
         rmv_r = rmv_r + one
         if ( csign < 0 )  cnegs = cnegs + 1
         RETURN
     endif ! back if ( ckink == 0 ) block

! randomly generate cis and cie at selected flvr channel, and then determine
! tau_start and tau_end for them
     call try_remove_colour(flvr, cis, cie, tau_start, tau_end)

! fast look up the flavor part of perturbation expansion series, determine
! corresponding fis and fie, and determine whether the operators trace is
! not equal to zero
     call try_remove_flavor(fis, fie, tau_start, tau_end, lrmv)

! calculate the transition ratio between old and new configurations,
! for the local trace part
     if ( lrmv .eqv. .true. ) then
         call cat_remove_ztrace(fis, fie, tau_start, tau_end, trace_ratio)
     else
         trace_ratio = zero
     endif ! back if ( lrmv .eqv. .true. ) block

! calculate the transition ratio between old and new configurations,
! for the determinant part
     if ( lrmv .eqv. .true. ) then
         call cat_remove_detrat(flvr, cis, cie, deter_ratio)
     else
         deter_ratio = zero
     endif ! back if ( lrmv .eqv. .true. ) block

! we will determine the pass by lazy trace evalution
! if lrmv is false, we set the pass as false immediately
     r = spring_sfmt_stream()
     trace_ratio = deter_ratio * ( real( ckink ) / beta ) ** 2
     if ( lrmv .eqv. .true. ) then
         call ctqmc_lazy_ztrace( 1, 2*sum(rank) - 2, trace_ratio, tau_start, tau_end, r, p, pass )
     else
         pass = .false.
     endif ! back if ( lrmv .eqv. .true. ) block

! if update action is accepted
     if ( pass .eqv. .true. ) then

! update the flavor part of perturbation expansion series
         call cat_remove_flavor(fis, fie, tau_start, tau_end)

! update the mmat matrix and gmat matrix, respectively,
! cat_remove_colour() subroutine is invoked internally to update the colour
! part of perturbation expansion series
         call cat_remove_matrix(flvr, cis, cie)

! update the operators trace
         call ctqmc_make_evolve()

! decrease ckink for current flavor channel
         ckink = ckink - 1

! update rank for current flavor channel
         rank(flvr) = rank(flvr) - 1

! determine the sign, TO BE CHECKED
         csign = csign * int ( sign(one, p) )

     endif ! back if ( pass .eqv. .true. ) block

! record negative sign
     if ( csign < 0 ) then
         cnegs = cnegs + 1
!<         call s_print_exception('ctqmc_remove_kink','csign is negative')
     endif ! back if ( csign < 0 ) block

! update the remove statistics
     rmv_t = rmv_t + one
     if ( pass .eqv. .true. ) then
         rmv_a = rmv_a + one
     else
         rmv_r = rmv_r + one
     endif ! back if ( pass .eqv. .true. ) block

     return
  end subroutine ctqmc_remove_kink

!!>>> ctqmc_lshift_kink: shift old create operators in the perturbation
!!>>> expansion series
  subroutine ctqmc_lshift_kink()
     use constants, only : dp, zero, one
     use spring, only : spring_sfmt_stream

     use control, only : norbs
     use context, only : ckink, csign, cnegs
     use context, only : lsh_t, lsh_a, lsh_r
     use context, only : rank

     implicit none

! local variables
! whether the old create operators can be shifted diagrammatically
     logical  :: lshf

! whether the update operation is accepted
     logical  :: pass

! current flavor channel for both band and spin
     integer  :: flvr

! index address to shift old create operators
! iso and isn are for old and new indices, respectively
! ciso and cisn are for the colour part, and fiso and fisn are for the flavor part
     integer  :: ciso, cisn
     integer  :: fiso, fisn

! transition probability for shift old create operators
     real(dp) :: p

! random number
     real(dp) :: r

! \tau_s, imaginary time point of the old create operator
     real(dp) :: tau_start1

! \tau_s, imaginary time point of the new create operator
     real(dp) :: tau_start2

! ratio between old and new configurations, the local trace part
     real(dp) :: trace_ratio

! ratio between old and new configurations, the determinant part
     real(dp) :: deter_ratio

! initialize logical variables
     lshf = .false.
     pass = .false.

! select the flavor channel randomly among 1 ~ norbs
     flvr = ceiling( spring_sfmt_stream() * norbs )

! get the perturbation expansion order ( number of existing create or
! destroy operators ) for current flavor channel
     ckink = rank(flvr)
     if ( ckink == 0 ) then
!<         call s_print_exception('ctqmc_lshift_kink','can not lshift any operators')
         lsh_t = lsh_t + one
         lsh_r = lsh_r + one
         if ( csign < 0 )  cnegs = cnegs + 1
         RETURN
     endif ! back if ( ckink == 0 ) block

! at first, we select ciso randomly, and then obtain tau_start1. according
! to the existing operators, we determine tau_start2 and related index cisn
     call try_lshift_colour(flvr, ciso, cisn, tau_start1, tau_start2)

! fast look up the flavor part of perturbation expansion series, determine
! corresponding fiso and fisn, and determine whether the operators trace is
! not equal to zero
     call try_lshift_flavor(flvr, fiso, fisn, tau_start1, tau_start2, lshf)

! calculate the transition ratio between old and new configurations,
! for the local trace part
     if ( lshf .eqv. .true. ) then
         call cat_lshift_ztrace(flvr, fiso, fisn, tau_start1, tau_start2, trace_ratio)
     else
         trace_ratio = zero
     endif ! back if ( lshf .eqv. .true. ) block

! calculate the transition ratio between old and new configurations,
! for the determinant part
     if ( lshf .eqv. .true. ) then
         call cat_lshift_detrat(flvr, ciso, tau_start1, tau_start2, deter_ratio)
     else
         deter_ratio = zero
     endif ! back if ( lshf .eqv. .true. ) block

! we will determine the pass by lazy trace evalution
! if lshf is false, we set the pass as false immediately
     r = spring_sfmt_stream()
     trace_ratio = deter_ratio * one
     if ( lshf .eqv. .true. ) then
         call ctqmc_lazy_ztrace( 1, 2*sum(rank), trace_ratio, tau_start1, tau_start2, r, p, pass )
     else
         pass = .false.
     endif ! back if ( lshf .eqv. .true. ) block

! if update action is accepted
     if ( pass .eqv. .true. ) then

! update the flavor part of perturbation expansion series
         call cat_lshift_flavor(flvr, fiso, fisn, tau_start2)

! update the mmat matrix and gmat matrix, respectively,
! cat_lshift_colour() subroutine is invoked internally to update the colour
! part of perturbation expansion series
         call cat_lshift_matrix(flvr, ciso, cisn, tau_start1, tau_start2, deter_ratio)

! update the operators trace
         call ctqmc_make_evolve()

! determine the sign, TO BE CHECKED
         csign = csign * int ( sign(one, p) )

     endif ! back if ( pass .eqv. .true. ) block

! record negative sign
     if ( csign < 0 ) then
         cnegs = cnegs + 1
!<         call s_print_exception('ctqmc_lshift_kink','csign is negative')
     endif ! back if ( csign < 0 ) block

! update the lshift statistics
     lsh_t = lsh_t + one
     if ( pass .eqv. .true. ) then
         lsh_a = lsh_a + one
     else
         lsh_r = lsh_r + one
     endif ! back if ( pass .eqv. .true. ) block

     return
  end subroutine ctqmc_lshift_kink

!!>>> ctqmc_rshift_kink: shift old destroy operators in the perturbation
!!>>> expansion series
  subroutine ctqmc_rshift_kink()
     use constants, only : dp, zero, one
     use spring, only : spring_sfmt_stream

     use control, only : norbs
     use context, only : ckink, csign, cnegs
     use context, only : rsh_t, rsh_a, rsh_r
     use context, only : rank

     implicit none

! local variables
! whether the old destroy operators can be shifted diagrammatically
     logical  :: rshf

! whether the update operation is accepted
     logical  :: pass

! current flavor channel for both band and spin
     integer  :: flvr

! index address to shift old destroy operators
! ieo and ien are for old and new indices, respectively
! cieo and cien are for the colour part, and fieo and fien are for the flavor part
     integer  :: cieo, cien
     integer  :: fieo, fien

! transition probability for shift old destroy operators
     real(dp) :: p

! random number
     real(dp) :: r

! \tau_e, imaginary time point of the old destroy operator
     real(dp) :: tau_end1

! \tau_e, imaginary time point of the new destroy operator
     real(dp) :: tau_end2

! ratio between old and new configurations, the local trace part
     real(dp) :: trace_ratio

! ratio between old and new configurations, the determinant part
     real(dp) :: deter_ratio

! initialize logical variables
     rshf = .false.
     pass = .false.

! select the flavor channel randomly among 1 ~ norbs
     flvr = ceiling( spring_sfmt_stream() * norbs )

! get the perturbation expansion order ( number of existing create or
! destroy operators ) for current flavor channel
     ckink = rank(flvr)
     if ( ckink == 0 ) then
!<         call s_print_exception('ctqmc_rshift_kink','can not rshift any operators')
         rsh_t = rsh_t + one
         rsh_r = rsh_r + one
         if ( csign < 0 )  cnegs = cnegs + 1
         RETURN
     endif ! back if ( ckink == 0 ) block

! at first, we select cieo randomly, and then obtain tau_end1. according
! to the existing operators, we determine tau_end2 and related index cien
     call try_rshift_colour(flvr, cieo, cien, tau_end1, tau_end2)

! fast look up the flavor part of perturbation expansion series, determine
! corresponding fieo and fien, and determine whether the operators trace is
! not equal to zero
     call try_rshift_flavor(flvr, fieo, fien, tau_end1, tau_end2, rshf)

! calculate the transition ratio between old and new configurations,
! for the local trace part
     if ( rshf .eqv. .true. ) then
         call cat_rshift_ztrace(flvr, fieo, fien, tau_end1, tau_end2, trace_ratio)
     else
         trace_ratio = zero
     endif ! back if ( rshf .eqv. .true. ) block

! calculate the transition ratio between old and new configurations,
! for the determinant part
     if ( rshf .eqv. .true. ) then
         call cat_rshift_detrat(flvr, cieo, tau_end1, tau_end2, deter_ratio)
     else
         deter_ratio = zero
     endif ! back if ( rshf .eqv. .true. ) block

! we will determine the pass by lazy trace evalution
! if rshf is false, we set the pass as false immediately
     r = spring_sfmt_stream()
     trace_ratio = deter_ratio * one
     if ( rshf .eqv. .true. ) then
         call ctqmc_lazy_ztrace( 1, 2*sum(rank), trace_ratio, tau_end1, tau_end2, r, p, pass )
     else
         rshf = .false.
     endif ! back if ( rshf .eqv. .true. ) block

! if update action is accepted
     if ( pass .eqv. .true. ) then

! update the flavor part of perturbation expansion series
         call cat_rshift_flavor(flvr, fieo, fien, tau_end2)

! update the mmat matrix and gmat matrix, respectively,
! cat_rshift_colour() subroutine is invoked internally to update the colour
! part of perturbation expansion series
         call cat_rshift_matrix(flvr, cieo, cien, tau_end1, tau_end2, deter_ratio)

! update the operators trace
         call ctqmc_make_evolve()

! determine the sign, TO BE CHECKED
         csign = csign * int ( sign(one, p) )

     endif ! back if ( pass .eqv. .true. ) block

! record negative sign
     if ( csign < 0 ) then
         cnegs = cnegs + 1
!<         call s_print_exception('ctqmc_rshift_kink','csign is negative')
     endif ! back if ( csign < 0 ) block

! update the rshift statistics
     rsh_t = rsh_t + one
     if ( pass .eqv. .true. ) then
         rsh_a = rsh_a + one
     else
         rsh_r = rsh_r + one
     endif ! back if ( pass .eqv. .true. ) block

     return
  end subroutine ctqmc_rshift_kink

!!>>> ctqmc_reflip_kink: perform a global update, exchange the states
!!>>> between spin up and spin down, it maybe useful for magnetic systems
  subroutine ctqmc_reflip_kink(cflip)
     use constants, only : dp, zero, one
     use spring, only : spring_sfmt_stream
     use stack, only : istack_getrest

     use control, only : nband, norbs
     use context, only : rfl_t, reflip_accept, reflip_reject
     use context, only : empty_v, index_t, index_v, flvr_v
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

! loop index
     integer  :: i

! selected flavor pairs
     integer  :: fup, fdn

! loop index for flavor channel
     integer  :: flvr

! maximum rank order
     integer  :: kmax

! total number of operators
     integer  :: nsize

! transition probability for global spin flip
     real(dp) :: p

! random number
     real(dp) :: r

! global flip determinant ratio
     real(dp) :: ratup
     real(dp) :: ratdn

! initialize logical variables
     pass = .false.

! initialize transition probability
     p = one

! get total number of operators
     nsize = istack_getrest( empty_v )

! not need to perform global flip if there are no operators at all
     if ( nsize == 0 ) then
!<         call s_print_exception('ctqmc_reflip_kink','can not reflip any operators')
         rfl_t = rfl_t + one
         reflip_reject = reflip_reject + one
         RETURN
     endif ! back if ( nsize == 0 ) block

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

! make a trial swap for flvr_v
         do i=1,nsize
             if ( flvr_v ( index_v(i) ) == fup ) then
                 flvr_v ( index_v(i) ) = fdn; CYCLE
             endif ! back if ( flvr_v ( index_v(i) ) == fup ) block
             if ( flvr_v ( index_v(i) ) == fdn ) then
                 flvr_v ( index_v(i) ) = fup; CYCLE
             endif ! back if ( flvr_v ( index_v(i) ) == fdn ) block
         enddo ! over i={1,nsize} loop

! make a copy of index_v, index_t is need by ctqmc_make_ztrace()
         do i=1,nsize
             index_t(i) = index_v(i)
         enddo ! over i={1,nsize} loop

! calculate the transition ratio between old and new configurations,
! for the local trace part, by lazy trace evaluation
         r = spring_sfmt_stream()
         ratup = p
         call ctqmc_lazy_ztrace( 3, nsize, ratup, zero, zero, r, p, pass )

! if update action is accepted
         if ( pass .eqv. .true. ) then

! get maximum rank order in spin up and spin down states
             kmax = max( rank(fup), rank(fdn) )

! swap global variables between spin up and spin down states
             call cat_reflip_matrix(fup, fdn, kmax)

! update the operators trace
             call ctqmc_make_evolve()

! if this update action can not be accepted, reset it
         else

! recover the original status of flvr_v
             do i=1,nsize
                 if ( flvr_v ( index_v(i) ) == fup ) then
                     flvr_v ( index_v(i) ) = fdn; CYCLE
                 endif ! back if ( flvr_v ( index_v(i) ) == fup ) block
                 if ( flvr_v ( index_v(i) ) == fdn ) then
                     flvr_v ( index_v(i) ) = fup; CYCLE
                 endif ! back if ( flvr_v ( index_v(i) ) == fdn ) block
             enddo ! over i={1,nsize} loop

! print exception information
!<             call s_print_exception('ctqmc_reflip_kink','quantum impurity solver refuse to reflip')

         endif ! back if ( pass .eqv. .true. ) block

! update the reflip statistics
         rfl_t = rfl_t + one
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

! make a trial swap for flvr_v
             do i=1,nsize
                 if ( flvr_v ( index_v(i) ) == fup ) then
                     flvr_v ( index_v(i) ) = fdn; CYCLE
                 endif ! back if ( flvr_v ( index_v(i) ) == fup ) block
                 if ( flvr_v ( index_v(i) ) == fdn ) then
                     flvr_v ( index_v(i) ) = fup; CYCLE
                 endif ! back if ( flvr_v ( index_v(i) ) == fdn ) block
             enddo ! over i={1,nsize} loop

! make a copy of index_v, index_t is need by ctqmc_make_ztrace()
             do i=1,nsize
                 index_t(i) = index_v(i)
             enddo ! over i={1,nsize} loop

! calculate the transition ratio between old and new configurations,
! for the local trace part, by lazy trace evaluation
             r = spring_sfmt_stream()
             ratup = p
             call ctqmc_lazy_ztrace( 3, nsize, ratup, zero, zero, r, p, pass )

! if update action is accepted
             if ( pass .eqv. .true. ) then

! get maximum rank order in spin up and spin down states
                 kmax = max( rank(fup), rank(fdn) )

! swap global variables between spin up and spin down states
                 call cat_reflip_matrix(fup, fdn, kmax)

! update the operators trace
                 call ctqmc_make_evolve()

! if this update action can not be accepted, reset it
             else

! recover the original status of flvr_v
                 do i=1,nsize
                     if ( flvr_v ( index_v(i) ) == fup ) then
                         flvr_v ( index_v(i) ) = fdn; CYCLE
                     endif ! back if ( flvr_v ( index_v(i) ) == fup ) block
                     if ( flvr_v ( index_v(i) ) == fdn ) then
                         flvr_v ( index_v(i) ) = fup; CYCLE
                     endif ! back if ( flvr_v ( index_v(i) ) == fdn ) block
                 enddo ! over i={1,nsize} loop

! print exception information
!<                 call s_print_exception('ctqmc_reflip_kink','quantum impurity solver refuse to reflip')

             endif ! back if ( pass .eqv. .true. ) block

! update the reflip statistics
             rfl_t = rfl_t + one
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

! make a trial swap for flvr_v
         do i=1,nsize
             flvr = flvr_v( index_v(i) )
             if ( flvr <= nband ) then
                 flvr_v ( index_v(i) ) = flvr + nband
             else
                 flvr_v ( index_v(i) ) = flvr - nband
             endif ! back if ( flvr <= nband ) block
         enddo ! over i={1,nsize} loop

! make a copy of index_v, index_t is need by ctqmc_make_ztrace()
         do i=1,nsize
             index_t(i) = index_v(i)
         enddo ! over i={1,nsize} loop

! calculate the transition ratio between old and new configurations,
! for the local trace part, by lazy trace evaluation
         r = spring_sfmt_stream()
         ratup = p
         call ctqmc_lazy_ztrace( 3, nsize, ratup, zero, zero, r, p, pass )

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

! update the operators trace
             call ctqmc_make_evolve()

! if this update action can not be accepted, reset it
         else

! recover the original status of flvr_v
             do i=1,nsize
                 flvr = flvr_v( index_v(i) )
                 if ( flvr <= nband ) then
                     flvr_v ( index_v(i) ) = flvr + nband
                 else
                     flvr_v ( index_v(i) ) = flvr - nband
                 endif ! back if ( flvr <= nband ) block
             enddo ! over i={1,nsize} loop

! print exception information
!<             call s_print_exception('ctqmc_reflip_kink','quantum impurity solver refuse to reflip')

         endif ! back if ( pass .eqv. .true. ) block

! update the reflip statistics
         rfl_t = rfl_t + one
         if ( pass .eqv. .true. ) then
             reflip_accept = reflip_accept + one
         else
             reflip_reject = reflip_reject + one
         endif ! back if ( pass .eqv. .true. ) block

     endif ! back if ( cflip == 1 ) block

     return
  end subroutine ctqmc_reflip_kink

!!>>> ctqmc_reload_kink: global update all operators in the perturbation
!!>>> expansion series
  subroutine ctqmc_reload_kink()
     use control, only : norbs
     use context, only : rank

     implicit none

! local variables
! loop index for flavor channel
     integer :: flvr

     do flvr=1,norbs

! check the perturbation expansion order ( number of existing create or
! destroy operators ) for current flavor channel
         if ( rank(flvr) == 0 ) CYCLE

! regenerate the mmat matrix and gmat matrix from scratch
         call cat_reload_matrix(flvr)

     enddo ! over flvr={1,norbs} loop

     return
  end subroutine ctqmc_reload_kink
