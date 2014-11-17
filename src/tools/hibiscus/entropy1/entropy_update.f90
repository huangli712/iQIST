!!!-----------------------------------------------------------------------
!!! project : hibiscus
!!! program : entropy_make_image
!!!           entropy_make_sampling
!!!           entropy_make_updating
!!! source  : entropy_update.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 01/09/2011 by li huang
!!!           01/26/2011 by li huang
!!!           11/18/2014 by li huang
!!! purpose : to provide basic subroutines (in other words, elementary
!!!           updating subroutines) for classic maximum entropy method
!!!           code
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> entropy_make_image: the driver subroutine of classic maximum entropy
!!>>> program
  subroutine entropy_make_image(G_qmc, G_dev, model, fnorm, fkern, image)
     use constants, only : dp, zero, one, two, half, mystd
     use mmpi, only : mp_allreduce, mp_barrier
     use spring, only : spring_sfmt_stream

     use control, only : ntime, nwmax, niter, ntune
     use control, only : ainit
     use control, only : nprocs, myid, master

     implicit none

! external arguments
! imaginary time green's function
     real(dp), intent(in)  :: G_qmc(ntime)

! error bar data
     real(dp), intent(in)  :: G_dev(ntime)

! reference model function
     real(dp), intent(in)  :: model(-nwmax:nwmax)

! normalization function
     real(dp), intent(in)  :: fnorm(-nwmax:nwmax)

! fermion kernel function
     real(dp), intent(in)  :: fkern(-nwmax:nwmax,ntime)

! spectrum function
     real(dp), intent(out) :: image(-nwmax:nwmax)

! local variables
! iteration number
     integer  :: iter

! variables to control the annealing procedure
     real(dp) :: rfac
     real(dp) :: tfac

! alpha factor
     real(dp) :: alpha

! entropy term
     real(dp) :: sterm

! trace of classic maximum entropy method linear algebra equation
     real(dp) :: trace

! dummy image function, used in mpi case
     real(dp) :: image_t(-nwmax:nwmax)

! see entropy_make_ckern() subroutine, important immediate arrays
     real(dp) :: ckern(2*nwmax+1,2*nwmax+1)

! calculate 2-order differential of Likelihood function
     call entropy_make_ckern(G_dev, fkern, ckern)

! evaluate initial image function using random number generator
     call random_number(image); image_t = zero

! and then normalize it
     call entropy_make_normal(one, fnorm, image)

! initialize some key parameters
     tfac  = 10.0_dp
     rfac  = 1.00_dp

     sterm = 2.0_dp
     trace = 1.0_dp

! setup initial alpha parameter
     alpha = ainit

! initialize the iteration number
     iter  = 0

! main loop
     ENTROPY_MAIN_LOOP: do while ( abs( sterm / trace - one ) > 0.0005_dp .and. iter <= niter )

! evaluate new image function
         call entropy_make_sampling(rfac, tfac, alpha, G_qmc, G_dev, model, fkern, image)

! calculate new entropy term
         call entropy_make_sterm(sterm, image, model)

! calculate trace
         call entropy_make_trace(trace, alpha, image, ckern)

! note: we implement the classical maximum entropy algorithm
! i.e. the following equation should be satisfied:
!    -2*\alpha*S = trace \Lambda * [\alpha * I + \Lambda ]^{-1}
! please refer to equation (4.28) in the reference
         sterm = - two * sterm * alpha

! adjust alpha parameter dynamically
         if ( trace / sterm < 0.05_dp ) then
             alpha = alpha * 0.05_dp
         else
             alpha = alpha * ( trace / sterm ) * ( one + 0.001_dp * ( spring_sfmt_stream() - half ) )
         endif ! back if ( trace / sterm < 0.05_dp ) block

! increase iteration number
         iter = iter + 1

! scaling the parameters
         tfac = 0.001_dp
         rfac = 0.050_dp

! write out the necessary information, only master node can do it
         if ( myid == master ) then
             write(mystd,'(4X,a,i3)')    '>>> iteration ', iter
             write(mystd,'(4X,a,f12.6)') '-2*alpha*S = ' , sterm
             write(mystd,'(4X,a,f12.6)') '     trace = ' , trace
             write(mystd,'(4X,a,f12.6)') '     alpha = ' , alpha
             write(mystd,'(4X,a,f12.6)') '     delta = ' , abs(trace - sterm)
             write(mystd,*)
         endif ! back if ( myid == master ) block

     enddo ENTROPY_MAIN_LOOP ! over do while loop

! write out the necessary information, only master node can do it
     if ( myid == master ) then
         write(mystd,'(4X,a)') '>>> waiting ...'
         write(mystd,'(4X,a)',advance='no') '>>> smoothing '
     endif ! back if ( myid == master ) block

! smooth the image function
     ENTROPY_TUNE_LOOP: do iter=1,ntune

! smoothing it
         call entropy_make_smooth(3, image)

! normalizing it
         call entropy_make_normal(one, fnorm, image)

! change the image function slightly
         tfac = 0.005_dp
         rfac = 0.005_dp
         call entropy_make_sampling(rfac, tfac, alpha, G_qmc, G_dev, model, fkern, image)

! write out the necessary information, only master node can do it
         if ( myid == master ) then
             write(mystd,'(a)',advance='no') '.'
         endif ! back if ( myid == master ) block

     enddo ENTROPY_TUNE_LOOP ! over iter={1,ntune} loop

! smoothing it again
     call entropy_make_smooth(3, image)

! write out the necessary information, only master node can do it
     if ( myid == master ) then
         write(mystd,'(a)') ' done '
         write(mystd,*)
     endif ! back if ( myid == master ) block

! record image to image_t, collect data from all children processes
# if defined (MPI)

! collect data
     call mp_allreduce(image, image_t)

! block until all processes have reached here
     call mp_barrier()

# else  /* MPI */

     image_t = image

# endif /* MPI */

! calculate the average image function
     image = image_t / real(nprocs)

     return
  end subroutine entropy_make_image

!!>>> entropy_make_sampling: implement the kernel algorithm of classic
!!>>> maximum entropy method
  subroutine entropy_make_sampling(rfac, tfac, alpha, G_qmc, G_dev, model, fkern, image)
     use constants, only : dp, zero

     use control, only : ntime, nwmax, nstep

     implicit none

! external arguments
! parameters to control the monte carlo annealing steps
     real(dp), intent(inout) :: rfac
     real(dp), intent(inout) :: tfac

! spectrum function
     real(dp), intent(inout) :: image(-nwmax:nwmax)

! alpha parameter
     real(dp), intent(in) :: alpha

! imaginary time green's function
     real(dp), intent(in) :: G_qmc(ntime)

! error bar data
     real(dp), intent(in) :: G_dev(ntime)

! reference model function
     real(dp), intent(in) :: model(-nwmax:nwmax)

! fermion kernel function
     real(dp), intent(in) :: fkern(-nwmax:nwmax,ntime)

! local variables
! loop index
     integer  :: i
     integer  :: j

! accept statistics
     real(dp) :: acc

! try statistics
     real(dp) :: try

! maximum value in image function
     real(dp) :: amax

! \delta image function
     real(dp) :: dimg(-nwmax:nwmax)

! local parameters
! how often to adjust the annealing parameters
     integer, parameter :: nfast = 100

     entropy_annealing_loop: do i=1,nstep,nfast

! reset try and acc to zero
         try = zero
         acc = zero

! calculate amax and dimg
         amax = maxval(image)
         do j=-nwmax,nwmax
             dimg(j) = rfac * ( image(j) / 10.0_dp + 0.01_dp * amax)
         enddo ! over j={-nwmax,nwmax} loop

         do j=1,nfast
             call entropy_make_updating(acc, try, tfac, amax, alpha, G_qmc, G_dev, dimg, model, fkern, image)
         enddo ! over j={1,nfast} loop

! adjust rfac factor every nfast steps
         if ( acc / try > 0.1_dp ) then
             if ( rfac < 0.01_dp ) then
                 rfac = rfac * 1.5_dp
             endif ! back if ( rfac < 0.01_dp ) block
         else
             if ( rfac > 0.001_dp ) then
                 rfac = rfac / 1.5_dp
             endif ! back if ( rfac > 0.001_dp ) block
         endif ! back if ( arat > 0.1_dp ) block

! scaling tfac
         tfac = tfac / 1.5_dp

     enddo entropy_annealing_loop ! over i={1,nstep} loop

     return
  end subroutine entropy_make_sampling

!!>>> entropy_make_updating: update the image function by monte
!!>>> carlo procedure
  subroutine entropy_make_updating(acc, try, tfac, amax, alpha, G_qmc, G_dev, delta, model, fkern, image)
     use constants, only : dp, zero, one, two, half, epss
     use spring, only : spring_sfmt_stream

     use control, only : ntime, nwmax
     use control, only : wstep

     implicit none

! external arguments
! accept statistics
     real(dp), intent(inout) :: acc

! try statistics
     real(dp), intent(inout) :: try

! image function
     real(dp), intent(inout) :: image(-nwmax:nwmax)

! parameters to control the monte carlo annealing steps
     real(dp), intent(in) :: tfac

! maximum value in image function
     real(dp), intent(in) :: amax

! alpha parameter
     real(dp), intent(in) :: alpha

! imaginary time green's function
     real(dp), intent(in) :: G_qmc(ntime)

! error bar data
     real(dp), intent(in) :: G_dev(ntime)

! \delta image function
     real(dp), intent(in) :: delta(-nwmax:nwmax)

! reference model function
     real(dp), intent(in) :: model(-nwmax:nwmax)

! fermion kernel function
     real(dp), intent(in) :: fkern(-nwmax:nwmax,ntime)

! local variables
! loop index
     integer  :: i
     integer  :: j

! selected index of frequency mesh
     integer  :: j1
     integer  :: j2

! transition probability
     real(dp) :: p

! \delta entropy term
     real(dp) :: del_ent

! \delta image on j1 and j2
     real(dp) :: dj1, del_img1
     real(dp) :: dj2, del_img2

! old and new chihc ( = \chi^{2} )
     real(dp) :: old_chihc
     real(dp) :: new_chihc

! old and new akern ( = K . A )
     real(dp) :: old_akern(ntime)
     real(dp) :: new_akern(ntime)

! calculate akern ( = K * A )
     call entropy_make_akern(old_akern, image, fkern)

! calculate chihc ( = \chi^{2} )
     call entropy_make_chihc(old_chihc, old_akern, G_qmc, G_dev)

     entropy_mesh_loop: do j=1,2*nwmax+1

! generate j1, j2, dj1, and dj2
         entropy_j1j2_loop: do while ( .true. )

! choose j1 and j2 between [-nwmax,nwmax], and j1 /= j2
             do while ( .true. )
                 j1 = min( int( spring_sfmt_stream() * ( 2 * nwmax + 1 ) ), 2 * nwmax ) - nwmax
                 j2 = min( int( spring_sfmt_stream() * ( 2 * nwmax + 1 ) ), 2 * nwmax ) - nwmax
                 if ( j1 /= j2 ) EXIT
             enddo ! over do while loop

! calculate dj1
             do while ( .true. )
                 dj1 = delta(j1) * ( spring_sfmt_stream() - half )
                 if ( image(j1) + dj1 > zero ) EXIT
             enddo ! over do while loop

! calculate dj2
             dj2 = -dj1 + 0.05_dp * delta(j2) * ( spring_sfmt_stream() - half )
             if ( image(j2) + dj2 > zero ) EXIT entropy_j1j2_loop

         enddo entropy_j1j2_loop ! over do while loop

! calculate new akern ( = K * delta A )
         do i=1,ntime
             new_akern(i) = old_akern(i) + dj1 * fkern(j1,i) + dj2 * fkern(j2,i)
         enddo ! over i={1,ntime} loop

! calculate new \chi^{2} using new akern
         call entropy_make_chihc(new_chihc, new_akern, G_qmc, G_dev)

! calculate \delta image on j1 and j2
         del_img1 = image(j1) + dj1
         del_img2 = image(j2) + dj2

! calculate \delta entropy term
         del_ent = zero

         if ( del_img1 > epss ) then
             del_ent = del_ent - wstep * del_img1 * log( del_img1 / model(j1) )
         endif ! back if ( del_img1 > epss ) block

         if ( del_img2 > epss ) then
             del_ent = del_ent - wstep * del_img2 * log( del_img2 / model(j2) )
         endif ! back if ( del_img2 > epss ) block

         if ( image(j1) > epss ) then
             del_ent = del_ent + wstep * image(j1) * log( image(j1) / model(j1) )
         endif ! back if ( image(j1) > epss ) block

         if ( image(j2) > epss ) then
             del_ent = del_ent + wstep * image(j2) * log( image(j2) / model(j2) )
         endif ! back if ( image(j2) > epss ) block

! calculate weight
         p = ( ( old_chihc - new_chihc ) / two + alpha * del_ent ) / tfac

! weight > 0.0
         if ( p > zero ) then
             p = one
! weight > -100.0 and weight < 0.0
         else if ( p > -100.0_dp ) then
             p = exp(p)
! weight < -100.0
         else
             p = zero
         endif ! back if ( p > zero ) block

! update the image and \chi^{2} if it is accepted
         if ( spring_sfmt_stream() < p ) then

! record accept statustics
             if ( image(j1) > 0.1_dp * amax ) then
                 acc = acc + one
                 try = try + one
             endif ! back if ( image(j1) > 0.1_dp * amax ) block

! update related variables
             image(j1) = del_img1
             image(j2) = del_img2
             old_akern = new_akern
             old_chihc = new_chihc

         else

! record try statistics
             if ( image(j1) > 0.1_dp * amax ) then
                 try = try + one
             endif ! back if ( image(j1) > 0.1_dp * amax ) block

         endif ! back if ( spring_sfmt_stream() < p ) block

     enddo entropy_mesh_loop ! over j={1,2*nwmax+1} loop

     return
  end subroutine entropy_make_updating
