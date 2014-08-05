!=========+=========+=========+=========+=========+=========+=========+>>>
! build spectral function from imaginary-time green's function using the !
! well-known maximum entropy method. in principle, it solves the laplace !
! transformation                                                         !
!     G(\tau) = \int kernel A(\omega) d\omega                            !
! where                                                                  !
!     kernel = \frac{ \exp{-\tau\omega} }{1.0+\exp{-\beta\omega}}        !
! for details of maximum entropy method, please refer to:                !
!     Physics Reports 269 (1996) 133-195                                 !
! author  : li huang                                                     !
! version : v2011.08.18T                                                 !
! status  : WARNING: IN TESTING STAGE, USE IT IN YOUR RISK               !
! comment : the code is written by Anders W. Sandvik (Akademi University,!
!           Finland, email:asandvik@ra.abo.fi) originally, and modified  !
!           by li huang using fortran 90 language                        !
!           any question, please contact with huangli712@yahoo.com.cn    !
!=========+=========+=========+=========+=========+=========+=========+>>>

  program entropy_main
     use context

     use mmpi

     implicit none

! local variables
! loop index for bands
     integer :: i

! initialize mpi envirnoment
# if defined (MPI)

! initialize the mpi execution environment
     call mp_init()

! determines the rank of the calling process in the communicator
     call mp_comm_rank(myid)

! determines the size of the group associated with a communicator
     call mp_comm_size(nprocs)

# endif  /* MPI */

! print the running header for classic maximum entropy method code
     if ( myid == master ) then ! only master node can do it
         call entropy_print_header()
     endif

! setup the important parameters for classic maximum entropy method code
     call entropy_config()

! print out runtime parameters in summary, only for check
     if ( myid == master ) then
         call entropy_print_summary()
     endif

! allocate memory and initialize
     call entropy_allocate_memory()

! read in imaginary time green's function and related imaginary time mesh
     call entropy_make_init1(tmesh, G_qmc, G_dev)

! prepare necessary data for classic maximum entropy method code
     call entropy_make_init2(tmesh, wmesh, model, fnorm, fkern)

     ENTROPY_BAND_LOOP: do i=1,norbs

! write out helpful information
         if ( myid == master ) then ! only master node can do it
             write(mystd,'(2X,2(a,i2))') 'HIBISCUS >>> data index: ', i, ' in ', norbs
         endif

! perform the classic maximum entropy algorithm
         call entropy_make_image(G_qmc(:,i), G_dev(:,i), model, fnorm(:,1), fkern, image(:,i))

! calculate the sum-rules values for image function
         call entropy_make_srule(fnorm, image(:,i), srule(:,i))

     enddo ENTROPY_BAND_LOOP ! over i={1,norbs} loop

! write out the final spectrum data
     if ( myid == master ) then ! only master node can do it
         call entropy_dump_image(wmesh, image)
     endif

! write out check data for sum-rules
     if ( myid == master ) then ! only master node can do it
         call entropy_dump_srule(srule)
     endif

! deallocate memory and finalize
     call entropy_deallocate_memory()

! print the footer for classic maximum entropy method code
     if ( myid == master ) then ! only master node can do it
         call entropy_print_footer()
     endif

! finalize mpi envirnoment
# if defined (MPI)

! blocks until all processes have reached this routine
     call mp_barrier()

! terminates mpi execution environment
     call mp_finalize()

# endif  /* MPI */

  end program entropy_main
