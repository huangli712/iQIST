!!!-----------------------------------------------------------------------
!!! project : narcissus
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
!!!           08/14/2017 by li huang (last modified)
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
!! to achieve thermodynamics stable equilibrium state. we also calculate
!! the autocorrelation function (for the total occupation number), and
!! then adjust the nmonte parameter here
!!
  subroutine ctqmc_try_warming()
     use constants, only : dp
     use constants, only : zero
     use constants, only : mystd

     use control, only : ntime
     use control, only : ntherm
     use control, only : nmonte
     use control, only : myid, master

     use context, only : ins_t, ins_a, ins_r
     use context, only : rmv_t, rmv_a, rmv_r
     use context, only : lsh_t, lsh_a, lsh_r
     use context, only : rsh_t, rsh_a, rsh_r
     use context, only : rfl_t, rfl_a, rfl_r
     use context, only : ac_f

     implicit none

! local parameters
! possible (predefined) values for the nmonte parameter
     integer, parameter :: P_NMONTE(10) = (/10, 20, 40, 50, 100, 200, 400, 500, 1000, 2000/)

! local variables
! loop index
     integer  :: i

! estimated autocorrelation time
     real(dp) :: ac_t

! autocorrelation function
     real(dp) :: ac_f_mpi(ntime + 2)
     real(dp) :: ac_f_err(ntime + 2)

! warm up the diagram series at first
     do i=1,ntherm
         call ctqmc_try_walking(i)
     enddo ! over i={1,ntherm} loop

! and then measure the autocorrelation function
     do i=1,ntherm
         call ctqmc_try_walking(i)
         call ctqmc_record_ac_f()
     enddo ! over i={1,ntherm} loop

! reduce the autocorrelation function, ac_f -> ac_f_mpi
     call ctqmc_reduce_ac_f(ac_f_mpi, ac_f_err); ac_f = ac_f_mpi

! normalize the autocorrelation function
     ac_f(1:ntime) = ac_f(1:ntime) / float( ntherm - ntime )
     ac_f(ntime + 1) = ac_f(ntime + 1) / float( ntherm )
     ac_f(ntime + 2) = ac_f(ntime + 2) / float( ntherm )

! calculate the autocorrelation function (numerator part)
     ac_f(1:ntime) = ac_f(1:ntime) - ac_f(ntime + 1)**2

! normalize the autocorrelation function again
     ac_f(1:ntime) = ac_f(1:ntime) / ( ac_f(ntime + 2) - ac_f(ntime + 1)**2 )

! evaluate the integrated autocorrelation time
     ac_t = 0.5_dp
     do i=2,ntime
         if ( ac_f(i) > zero ) then
             ac_t = ac_t + ac_f(i)
         else
             EXIT
         endif ! back if ( ac_f(i) > zero ) block
     enddo ! over i={2,ntime} loop

! update nmonte parameter to reduce autocorrelation (old algorithm)
!<     do while ( nmonte < ac_t )
!<         nmonte = nmonte * 10
!<     enddo ! over do while loop
!<     if ( nmonte > ac_t * 5 ) nmonte = nmonte / 5 ! adjust nmonte further

! update nmonte parameter to reduce autocorrelation (new algorithm)
     nmonte = P_NMONTE( count ( real(P_NMONTE, dp) < ac_t ) + 1 )

! write the autocorrelation function, only master node can do it
     if ( myid == master ) then
         write(mystd,'(4X,a,f11.4)',advance='no') 'ac_t:', ac_t
         write(mystd,'(1X,a,i4,a)') '(nmonte ->', nmonte, ')'
         call ctqmc_dump_ac_f(ac_f)
     endif ! back if ( myid == master ) block

! reset statistics variables
     ins_t = zero; ins_a = zero; ins_r = zero
     rmv_t = zero; rmv_a = zero; rmv_r = zero
     lsh_t = zero; lsh_a = zero; lsh_r = zero
     rsh_t = zero; rsh_a = zero; rsh_r = zero
     rfl_t = zero; rfl_a = zero; rfl_r = zero

     return
  end subroutine ctqmc_try_warming

!!
!! @sub ctqmc_try_walking
!!
!! visit the perturbation expansion diagrams randomly
!!
  subroutine ctqmc_try_walking(cstep)
     use constants, only : dp

     use spring, only : spring_sfmt_stream

     use control, only : nflip
     use control, only : nclean

     implicit none

! external arguments
! current QMC sweep steps
     integer, intent(in) :: cstep

! random walking in C_Z space
!-------------------------------------------------------------------------
     C_Z_SPACE: BLOCK

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
                 call ctqmc_lshift_kink()  ! shift the creation operators
             else
                 call ctqmc_rshift_kink()  ! shift the annihilation operators
             endif ! back if ( spring_sfmt_stream() > 0.5_dp ) block
         endif ! back if ( spring_sfmt_stream() < 0.9_dp ) block

     END BLOCK C_Z_SPACE

! numerical trick: perform global spin flip periodically
!-------------------------------------------------------------------------
     GLOBAL_REFLIP: BLOCK

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

     END BLOCK GLOBAL_REFLIP

! numerical trick: perform global update periodically
!-------------------------------------------------------------------------
     GLOBAL_RELOAD: BLOCK

         if ( nclean > 0 .and. mod(cstep, nclean) == 0 ) then
             call ctqmc_reload_kink()
         endif ! back if ( nclean > 0 .and. mod(cstep, nclean) == 0 ) block

     END BLOCK GLOBAL_RELOAD

     return
  end subroutine ctqmc_try_walking

!!
!! @sub ctqmc_try_warning
!!
!! checking whether the quantum impurity solver is consistent internally
!!
  subroutine ctqmc_try_warning(cflag)
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

! check time order of operators
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
                 write(mystd,'(4X,a)') '>>> quantum impurity solver status: fatal error'
                 write(mystd,'(4X,a)') '>>> please check the status file: solver.status.dat'
                 call ctqmc_save_status()
                 call s_print_error('ctqmc_try_warning','unknown fatal error occur')
             else
                 write(mystd,'(4X,a)') '>>> quantum impurity solver status: healthy'
             endif ! back if ( cflag == 99 ) block
         endif ! back if ( myid == master ) block

     endif ! back if ( cflag == 1 ) block

     return
  end subroutine ctqmc_try_warning

!!========================================================================
!!>>> service layer: updating perturbation expansion series 2          <<<
!!========================================================================

!!
!! @sub ctqmc_insert_kink
!!
!! insert new creation and annihilation operators in the perturbation
!! expansion series
!!
  subroutine ctqmc_insert_kink()
     use constants, only : dp
     use constants, only : zero, one

     use spring, only : spring_sfmt_stream

     use control, only : norbs
     use control, only : mkink
     use control, only : beta

     use context, only : ckink, cstat
     use context, only : ins_t, ins_a, ins_r
     use context, only : rank, stts

     implicit none

! local variables
! whether the new creation and annihilation operators can be inserted
     logical  :: ladd

! whether it is an anti-segment
! if anti = .true., anti-segment
! if anti = .false., segment
     logical  :: anti

! whether the update operation is accepted
     logical  :: pass

! current flavor channel
     integer  :: flvr

! index address to insert new creation and annihilation operators
! is and ie are for creation and annihilation operators, respectively
     integer  :: is, ie

! transition probability
     real(dp) :: p

! \tau_s, imaginary time point of the creation operator
     real(dp) :: tau_start

! \tau_e, imaginary time point of the annihilation operator
     real(dp) :: tau_end

! possible maximum length of the new segment
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

! get the perturbation expansion order for current flavor channel
     ckink = rank(flvr)
     if ( ckink == mkink ) then
         ins_t = ins_t + one
         ins_r = ins_r + one
         RETURN
     endif ! back if ( ckink == mkink ) block

! try to generate new configuration
! (1) randomly choose anti and tau_start
! (2) check whether tau_start is valid
! (3) if tau_start is valid, determine tau_end, tau_max, is, and ie
!     consistently, and set ladd to .true.
! (4) if tau_start is not valid, set ladd to .false.
     call try_insert_colour(flvr, is, ie, anti, ladd, tau_start, tau_end, tau_max)

! calculate the transition ratio for the local trace part
     if ( ladd .eqv. .true. ) then
         call cat_insert_ztrace(flvr, anti, tau_start, tau_end, trace_ratio)
     else
         trace_ratio = zero
     endif ! back if ( ladd .eqv. .true. ) block

! calculate the transition ratio for the determinant part
     if ( ladd .eqv. .true. ) then
         call cat_insert_detrat(flvr, tau_start, tau_end, deter_ratio)
     else
         deter_ratio = zero
     endif ! back if ( ladd .eqv. .true. ) block

! calculate the transition probability
     p = deter_ratio * trace_ratio * ( beta * tau_max / real( ckink + 1 ) )

! determine pass using metropolis algorithm
     pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if the update action is accepted
     if ( pass .eqv. .true. ) then

! update the mmat matrix and gmat matrix, respectively
! the perturbation expansion series are updated as well
         call cat_insert_matrix(flvr, is, ie, tau_start, tau_end, deter_ratio)

! update ckink for current flavor channel
         ckink = ckink + 1

! update stts for current flavor channel
         stts(flvr) = cstat

! update rank for current flavor channel
         rank(flvr) = rank(flvr) + 1

     endif ! back if ( pass .eqv. .true. ) block

! update monte carlo statistics
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
!! remove old creation and annihilation operators in the perturbation
!! expansion series
!!
  subroutine ctqmc_remove_kink()
     use constants, only : dp
     use constants, only : one

     use spring, only : spring_sfmt_stream

     use control, only : norbs
     use control, only : beta

     use context, only : ckink, cstat
     use context, only : rmv_t, rmv_a, rmv_r
     use context, only : rank, stts

     implicit none

! local variables
! whether it is an anti-segment
! if anti = .true., anti-segment
! if anti = .false., segment
     logical  :: anti

! whether the update operation is accepted
     logical  :: pass

! current flavor channel
     integer  :: flvr

! index address to remove old creation and annihilation operators
! is and ie are for creation and annihilation operators, respectively
     integer  :: is, ie

! transition probability
     real(dp) :: p

! \tau_s, imaginary time point of the creation operator
     real(dp) :: tau_start

! \tau_e, imaginary time point of the annihilation operator
     real(dp) :: tau_end

! possible maximum length of the old segment
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

! get the perturbation expansion order for current flavor channel
     ckink = rank(flvr)
     if ( ckink == 0 ) then
         rmv_t = rmv_t + one
         rmv_r = rmv_r + one
         RETURN
     endif ! back if ( ckink == 0 ) block

! try to generate new configuration
! (1) determine anti and is randomly
! (2) tau_start is obtained by is
! (3) determine ie
! (4) tau_end and tau_max are evaluated carefully according to is and ie
     call try_remove_colour(flvr, is, ie, anti, tau_start, tau_end, tau_max)

! calculate the transition ratio for the local trace part
     call cat_remove_ztrace(flvr, anti, tau_start, tau_end, trace_ratio)

! calculate the transition ratio for the determinant part
     call cat_remove_detrat(flvr, is, ie, deter_ratio)

! calculate the transition probability
     p = deter_ratio * trace_ratio * real( ckink ) / ( beta * tau_max )

! determine pass using metropolis algorithm
     pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if update action is accepted
     if ( pass .eqv. .true. ) then

! update the mmat matrix and gmat matrix, respectively
! the perturbation expansion series are updated as well
         call cat_remove_matrix(flvr, is, ie)

! update ckink for current flavor channel
         ckink = ckink - 1

! update stts for current flavor channel
         stts(flvr) = cstat

! update rank for current flavor channel
         rank(flvr) = rank(flvr) - 1

     endif ! back if ( pass .eqv. .true. ) block

! update monte carlo statistics
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
!! shift old creation operator in the perturbation expansion series
!!
  subroutine ctqmc_lshift_kink()
     use constants, only : dp
     use constants, only : one

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

! current flavor channel
     integer  :: flvr

! index address to shift old creation operator
! iso and isn are for old and new indices, respectively
     integer  :: iso, isn

! transition probability
     real(dp) :: p

! \tau_s, imaginary time point of the old creation operator
     real(dp) :: tau_start1

! \tau_s, imaginary time point of the new creation operator
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

! get the perturbation expansion order for current flavor channel
     ckink = rank(flvr)
     if ( ckink == 0 ) then
         lsh_t = lsh_t + one
         lsh_r = lsh_r + one
         RETURN
     endif ! back if ( ckink == 0 ) block

! try to generate new configuration
! (1) at first, we select iso randomly
! (2) obtain tau_start1 according to iso
! (3) based on the existing configuration, we determine tau_start2 and
!     related index isn
! (4) finally ring is evaluated
     call try_lshift_colour(flvr, iso, isn, ring, tau_start1, tau_start2)

! calculate the transition ratio for the local trace part
     call cat_lshift_ztrace(flvr, ring, tau_start1, tau_start2, trace_ratio)

! calculate the transition ratio for the determinant part
     call cat_lshift_detrat(flvr, iso, tau_start1, tau_start2, deter_ratio)

! calculate the transition probability
     p = deter_ratio * trace_ratio

! determine pass using metropolis algorithm
     pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if update action is accepted
     if ( pass .eqv. .true. ) then

! update the mmat matrix and gmat matrix, respectively
! the perturbation expansion series are updated as well
         call cat_lshift_matrix(flvr, iso, isn, tau_start1, tau_start2, deter_ratio)

! update stts for current flavor channel
         stts(flvr) = cstat

     endif ! back if ( pass .eqv. .true. ) block

! update monte carlo statistics
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
!! shift old annihilation operator in the perturbation expansion series
!!
  subroutine ctqmc_rshift_kink()
     use constants, only : dp
     use constants, only : one

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

! current flavor channel
     integer  :: flvr

! index address to shift old annihilation operator
! ieo and ien are for old and new indices, respectively
     integer  :: ieo, ien

! transition probability
     real(dp) :: p

! \tau_e, imaginary time point of the old annihilation operator
     real(dp) :: tau_end1

! \tau_e, imaginary time point of the new annihilation operator
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

! get the perturbation expansion order for current flavor channel
     ckink = rank(flvr)
     if ( ckink == 0 ) then
         rsh_t = rsh_t + one
         rsh_r = rsh_r + one
         RETURN
     endif ! back if ( ckink == 0 ) block

! try to generate new configuration
! (1) at first, we select ieo randomly
! (2) obtain tau_end1 according to ieo
! (3) based on the existing configuration, we determine tau_end2 and
!     related index ien
! (4) finally ring is evaluated
     call try_rshift_colour(flvr, ieo, ien, ring, tau_end1, tau_end2)

! calculate the transition ratio for the local trace part
     call cat_rshift_ztrace(flvr, ring, tau_end1, tau_end2, trace_ratio)

! calculate the transition ratio for the determinant part
     call cat_rshift_detrat(flvr, ieo, tau_end1, tau_end2, deter_ratio)

! calculate the transition probability
     p = deter_ratio * trace_ratio

! determine pass using metropolis algorithm
     pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if update action is accepted
     if ( pass .eqv. .true. ) then

! update the mmat matrix and gmat matrix, respectively
! the perturbation expansion series are updated as well
         call cat_rshift_matrix(flvr, ieo, ien, tau_end1, tau_end2, deter_ratio)

! update stts for current flavor channel
         stts(flvr) = cstat

     endif ! back if ( pass .eqv. .true. ) block

! update monte carlo statistics
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
     use constants, only : dp
     use constants, only : one

     use spring, only : spring_sfmt_stream

     use control, only : nband

     use context, only : rfl_t, rfl_a, rfl_r
     use context, only : rank

     implicit none

! external arguments
! control flag
! if cflip = 1, flip intra-orbital spins one by one
! if cflip = 2, flip intra-orbital spins globally
     integer, intent(in) :: cflip

! local variables
! whether the update operation is accepted
     logical  :: pass

! selected flavor pairs
     integer  :: fup
     integer  :: fdn

! loop index for flavor channel
     integer  :: flvr

! maximum rank order
     integer  :: kmax

! transition probability
     real(dp) :: p

! global flip determinant ratio
     real(dp) :: ratup
     real(dp) :: ratdn

! initialize logical variables
     pass = .false.

! initialize transition probability
     p = one

! case 1: cflip = 1, local flip
!-------------------------------------------------------------------------
     if ( cflip == 1 ) then
         do flvr=1,nband

! get fup and fdn
             fup = flvr; fdn = flvr + nband

! calculate the transition ratio for the determinant part, spin up case
             call cat_reflip_detrat(fup, fdn, ratup)

! calculate the transition ratio for the determinant part, spin dn case
             call cat_reflip_detrat(fdn, fup, ratdn)

! calculate the transition probability
             p = ratup * ratdn

! determine pass using metropolis algorithm
             pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if update action is accepted
             if ( pass .eqv. .true. ) then

! get maximum rank order in spin up and spin down states
                 kmax = max( rank(fup), rank(fdn) )

! exchange global variables between spin up and spin down states
                 call cat_reflip_matrix(fup, fdn, kmax)

             endif ! back if ( pass .eqv. .true. ) block

! update monte carlo statistics
             rfl_t = rfl_t + one
             if ( pass .eqv. .true. ) then
                 rfl_a = rfl_a + one
             else
                 rfl_r = rfl_r + one
             endif ! back if ( pass .eqv. .true. ) block

         enddo ! over flvr={1,nband} loop
     endif ! back if ( cflip == 1 ) block

! case 2: cflip = 2, global flip
!-------------------------------------------------------------------------
     if ( cflip == 2 ) then
         do flvr=1,nband

! get fup and fdn
             fup = flvr; fdn = flvr + nband

! calculate the transition ratio for the determinant part, spin up case
             call cat_reflip_detrat(fup, fdn, ratup)

! calculate the transition ratio for the determinant part, spin dn case
             call cat_reflip_detrat(fdn, fup, ratdn)

! calculate the transition probability
             p = p * ( ratup * ratdn )

         enddo ! over flvr={1,nband} loop

! determine pass using metropolis algorithm
         pass = ( min( one, abs(p) ) > spring_sfmt_stream() )

! if update action is accepted
         if ( pass .eqv. .true. ) then

             do flvr=1,nband

! get fup and fdn
                 fup = flvr; fdn = flvr + nband

! get maximum rank order in spin up and spin down states
                 kmax = max( rank(fup), rank(fdn) )

! exchange global variables between spin up and spin down states
                 call cat_reflip_matrix(fup, fdn, kmax)

             enddo ! over flvr={1,nband} loop

         endif ! back if ( pass .eqv. .true. ) block

! update monte carlo statistics
         rfl_t = rfl_t + one
         if ( pass .eqv. .true. ) then
             rfl_a = rfl_a + one
         else
             rfl_r = rfl_r + one
         endif ! back if ( pass .eqv. .true. ) block
     endif ! back if ( cflip == 2 ) block

     return
  end subroutine ctqmc_reflip_kink

!!
!! @sub ctqmc_reload_kink
!!
!! reload all creation and annihilation operators in the perturbation
!! expansion series, then rebuild all related global matrices from scratch
!!
  subroutine ctqmc_reload_kink()
     use control, only : norbs

     use context, only : rank

     implicit none

! local variables
! loop index for flavor channel
     integer :: flvr

     do flvr=1,norbs

! check the perturbation expansion order for current flavor channel
         if ( rank(flvr) == 0 ) CYCLE

! generate the mmat matrix and gmat matrix from scratch
         call cat_reload_matrix(flvr)

     enddo ! over flvr={1,norbs} loop

     return
  end subroutine ctqmc_reload_kink
