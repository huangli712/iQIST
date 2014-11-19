!-------------------------------------------------------------------------
! project : hibiscus
! program : context    module
! source  : sai_context.f90
! type    : module
! author  : li huang (email:huangli712@yahoo.com.cn)
! history : 01/08/2011 by li huang
!           01/09/2011 by li huang
!           01/10/2011 by li huang
! purpose : define the key data structure and global arrays/variables for
!           stochastic analytic continuation code
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

  module context
     use constants
     use control

     implicit none

! monte carlo configuration: a_{\gamma}
     integer, public, save, allocatable  :: igamm(:,:)

! monte carlo configuration: r_{\gamma}
     real(dp), public, save, allocatable :: rgamm(:,:)

! fermion kernel function
     real(dp), public, save, allocatable :: fkern(:,:)

! legendre polynomial defined on [-1,1]
     real(dp), public, save, allocatable :: ppleg(:,:)

! precalculated delta function, to facilite the calculation of image function
     real(dp), public, save, allocatable :: delta(:,:)

! alpha-resolved image function, it is what we need
     real(dp), public, save, allocatable :: image(:,:)

! alpha-resolved image function, it is what we need, used in mpi case
     real(dp), public, save, allocatable :: immpi(:,:)

! \phi(\omega) function
     real(dp), public, save, allocatable :: F_phi(:)

! default model D(\omega)
     real(dp), public, save, allocatable :: model(:)

! real frequency mesh
     real(dp), public, save, allocatable :: wmesh(:)

! imaginary time mesh
     real(dp), public, save, allocatable :: tmesh(:)

! interval [-1,1] on which legendre polynomial is defined
     real(dp), public, save, allocatable :: pmesh(:)

! very dense grid on [0,1]
     real(dp), public, save, allocatable :: xgrid(:)

! very dense grid on [ -wstep*nwmax, +wstep*nwmax ]
     real(dp), public, save, allocatable :: wgrid(:)

! imaginary time green's function from QMC quantum impurity solver: G(t)
     real(dp), public, save, allocatable :: G_qmc(:)

! imaginary time green's function rescaled by the variance: G(t)/sigma(t)
     real(dp), public, save, allocatable :: G_tau(:)

! standard deviation of G(t) from QMC data: \sigma(t)
     real(dp), public, save, allocatable :: G_dev(:)

! list of alpha parameters used in parallel tempering
     real(dp), public, save, allocatable :: alpha(:)

! hamiltonian for different alpha parameters
     real(dp), public, save, allocatable :: hamil(:)

! move action statistics of monte carlo procedure
     real(dp), public, save, allocatable :: move_accept(:)
     real(dp), public, save, allocatable :: move_reject(:)
     real(dp), public, save, allocatable :: move_tcount(:)

! swap action statistics of monte carlo procedure
     real(dp), public, save, allocatable :: swap_accept(:)
     real(dp), public, save, allocatable :: swap_reject(:)
     real(dp), public, save, allocatable :: swap_tcount(:)

  contains

!>>> allocate module memory
     subroutine sai_allocate_memory()
         implicit none

! status flag
         integer :: istat

! allocate memory
         allocate(igamm(nalph,ngamm),        stat=istat)
         allocate(rgamm(nalph,ngamm),        stat=istat)
         allocate(fkern(ntime,ngrid),        stat=istat)
         allocate(ppleg(legrd,lemax),        stat=istat)

         allocate(delta(-nwmax:nwmax,ngrid), stat=istat)
         allocate(image(nalph,-nwmax:nwmax), stat=istat)
         allocate(immpi(nalph,-nwmax:nwmax), stat=istat)

         allocate(F_phi(-nwmax:nwmax),       stat=istat)
         allocate(model(-nwmax:nwmax),       stat=istat)

         allocate(wmesh(-nwmax:nwmax),       stat=istat)
         allocate(tmesh(ntime),              stat=istat)
         allocate(pmesh(legrd),              stat=istat)
         allocate(xgrid(ngrid),              stat=istat)
         allocate(wgrid(ngrid),              stat=istat)

         allocate(G_qmc(ntime),              stat=istat)
         allocate(G_tau(ntime),              stat=istat)
         allocate(G_dev(ntime),              stat=istat)

         allocate(alpha(nalph),              stat=istat)
         allocate(hamil(nalph),              stat=istat)

         allocate(move_accept(nalph),        stat=istat)
         allocate(move_reject(nalph),        stat=istat)
         allocate(move_tcount(nalph),        stat=istat)

         allocate(swap_accept(nalph),        stat=istat)
         allocate(swap_reject(nalph),        stat=istat)
         allocate(swap_tcount(nalph),        stat=istat)


! check the status
         if ( istat /= 0 ) then
             call sai_print_error('sai_allocate_memory','can not allocate enough memory')
         endif

! initialize them
         igamm = 0
         rgamm = zero
         fkern = zero
         ppleg = zero

         delta = zero
         image = zero
         immpi = zero

         F_phi = zero
         model = zero

         wmesh = zero
         tmesh = zero
         pmesh = zero
         xgrid = zero
         wgrid = zero

         G_qmc = zero
         G_tau = zero
         G_dev = zero

         alpha = zero
         hamil = zero

         move_accept = zero
         move_reject = zero
         move_tcount = zero

         swap_accept = zero
         swap_reject = zero
         swap_tcount = zero

         return
     end subroutine sai_allocate_memory

!>>> deallocate module memory
     subroutine sai_deallocate_memory()
         implicit none

         if ( allocated(igamm) )  deallocate(igamm)
         if ( allocated(rgamm) )  deallocate(rgamm)
         if ( allocated(fkern) )  deallocate(fkern)
         if ( allocated(ppleg) )  deallocate(ppleg)

         if ( allocated(delta) )  deallocate(delta)
         if ( allocated(image) )  deallocate(image)
         if ( allocated(immpi) )  deallocate(immpi)

         if ( allocated(F_phi) )  deallocate(F_phi)
         if ( allocated(model) )  deallocate(model)

         if ( allocated(wmesh) )  deallocate(wmesh)
         if ( allocated(tmesh) )  deallocate(tmesh)
         if ( allocated(pmesh) )  deallocate(pmesh)
         if ( allocated(xgrid) )  deallocate(xgrid)
         if ( allocated(wgrid) )  deallocate(wgrid)

         if ( allocated(G_qmc) )  deallocate(G_qmc)
         if ( allocated(G_tau) )  deallocate(G_tau)
         if ( allocated(G_dev) )  deallocate(G_dev)

         if ( allocated(alpha) )  deallocate(alpha)
         if ( allocated(hamil) )  deallocate(hamil)

         if ( allocated(move_accept) ) deallocate(move_accept)
         if ( allocated(move_reject) ) deallocate(move_reject)
         if ( allocated(move_tcount) ) deallocate(move_tcount)

         if ( allocated(swap_accept) ) deallocate(swap_accept)
         if ( allocated(swap_reject) ) deallocate(swap_reject)
         if ( allocated(swap_tcount) ) deallocate(swap_tcount)

         return
     end subroutine sai_deallocate_memory

  end module context
