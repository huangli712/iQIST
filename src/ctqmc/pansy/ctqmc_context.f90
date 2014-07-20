!-------------------------------------------------------------------------
! project : pansy
! program : ctqmc_core module
!           ctqmc_clur module
!           ctqmc_flvr module
!           ctqmc_umat module
!           ctqmc_fmat module
!           ctqmc_mmat module
!           ctqmc_gmat module
!           ctqmc_wmat module
!           ctqmc_smat module
!           ctqmc_sect module
!           context    module
! source  : ctqmc_context.f90
! type    : module
! author  : li huang (email:huangli712@yahoo.com.cn)
!         : yilin wang (qhwyl2006@126.com)
! history : 09/16/2009 by li huang
!           09/17/2009 by li huang
!           09/19/2009 by li huang
!           09/20/2009 by li huang
!           09/21/2009 by li huang
!           09/22/2009 by li huang
!           09/27/2009 by li huang
!           11/01/2009 by li huang
!           11/10/2009 by li huang
!           11/18/2009 by li huang
!           12/01/2009 by li huang
!           12/05/2009 by li huang
!           02/21/2010 by li huang
!           02/23/2010 by li huang
!           06/08/2010 by li huang
!           07/16/2014 by yilin wang
!           07/19/2014 by yilin wang
! purpose : define the key data structure and global arrays/variables for
!           hybridization expansion version continuous time quantum Monte
!           Carlo (CTQMC) quantum impurity solver and dynamical mean field
!           theory (DMFT) self-consistent engine
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!=========================================================================
!>>> module ctqmc_core                                                 <<<
!=========================================================================
!>>> containing core (internal) variables used by continuous time quantum
! Monte Carlo quantum impurity solver
  module ctqmc_core
     use constants, only : dp, zero

     implicit none

! current perturbation expansion order
     integer, public, save  :: ckink = 0

! sign change related with current diagram update operation
     integer, public, save  :: csign = 1

! counter for negative sign, used to measure the sign problem
     integer, public, save  :: cnegs = 0

! averaged sign values, used to measure the sign problem
     integer, public, save  :: caves = 0
!-------------------------------------------------------------------------
!::: core variables: real, matrix trace                                :::
!-------------------------------------------------------------------------

! matrix trace of flavor part, current value
     real(dp), public, save :: matrix_ptrace = zero

! matrix trace of flavor part, proposed value
     real(dp), public, save :: matrix_ntrace = zero

!-------------------------------------------------------------------------
!::: core variables: real, insert action counter                       :::
!-------------------------------------------------------------------------

! insert kink (operators pair) statistics: total insert count
     real(dp), public, save :: insert_tcount = zero

! insert kink (operators pair) statistics: total accepted insert count
     real(dp), public, save :: insert_accept = zero

! insert kink (operators pair) statistics: total rejected insert count
     real(dp), public, save :: insert_reject = zero

!-------------------------------------------------------------------------
!::: core variables: real, remove action counter                       :::
!-------------------------------------------------------------------------

! remove kink (operators pair) statistics: total remove count
     real(dp), public, save :: remove_tcount = zero

! remove kink (operators pair) statistics: total accepted remove count
     real(dp), public, save :: remove_accept = zero

! remove kink (operators pair) statistics: total rejected remove count
     real(dp), public, save :: remove_reject = zero

!-------------------------------------------------------------------------
!::: core variables: real, lshift action counter                       :::
!-------------------------------------------------------------------------

! lshift kink (operators pair) statistics: total lshift count
     real(dp), public, save :: lshift_tcount = zero

! lshift kink (operators pair) statistics: total accepted lshift count
     real(dp), public, save :: lshift_accept = zero

! lshift kink (operators pair) statistics: total rejected lshift count
     real(dp), public, save :: lshift_reject = zero

!-------------------------------------------------------------------------
!::: core variables: real, rshift action counter                       :::
!-------------------------------------------------------------------------

! rshift kink (operators pair) statistics: total rshift count
     real(dp), public, save :: rshift_tcount = zero

! rshift kink (operators pair) statistics: total accepted rshift count
     real(dp), public, save :: rshift_accept = zero

! rshift kink (operators pair) statistics: total rejected rshift count
     real(dp), public, save :: rshift_reject = zero

!-------------------------------------------------------------------------
!::: core variables: real, reflip action counter                       :::
!-------------------------------------------------------------------------

! reflip kink (operators pair) statistics: total reflip count
     real(dp), public, save :: reflip_tcount = zero

! reflip kink (operators pair) statistics: total accepted reflip count
     real(dp), public, save :: reflip_accept = zero

! reflip kink (operators pair) statistics: total rejected reflip count
     real(dp), public, save :: reflip_reject = zero

  end module ctqmc_core

!=========================================================================
!>>> module ctqmc_clur                                                 <<<
!=========================================================================
!>>> containing perturbation expansion series related arrays (colour part)
! used by continuous time quantum Monte Carlo quantum impurity solver
  module ctqmc_clur
     use constants, only : dp
     use stack

     implicit none

! memory address index for the imaginary time \tau_s
     integer, public, save, allocatable :: index_s(:,:)

! memory address index for the imaginary time \tau_e
     integer, public, save, allocatable :: index_e(:,:)

! imaginary time \tau_s of create  operators
     real(dp), public, save, allocatable :: time_s(:,:)

! imaginary time \tau_e of destroy operators
     real(dp), public, save, allocatable :: time_e(:,:)

! exp(i\omega t), s means create  operators
     complex(dp), public, save, allocatable :: exp_s(:,:,:)

! exp(i\omega t), e means destroy operators
     complex(dp), public, save, allocatable :: exp_e(:,:,:)

! container for the empty (unoccupied) memory address index
     type (istack), public, save, allocatable :: empty_s(:)

! container for the empty (unoccupied) memory address index
     type (istack), public, save, allocatable :: empty_e(:)

  end module ctqmc_clur

!=========================================================================
!>>> module ctqmc_flvr                                                 <<<
!=========================================================================
!>>> containing perturbation expansion series related arrays (flavor part)
! used by continuous time quantum Monte Carlo quantum impurity solver
  module ctqmc_flvr
     use constants, only : dp
     use stack

     implicit none

! container for the empty (unoccupied) memory address index of operators
     type (istack), public, save :: empty_v

! memory address index for the imaginary time \tau (auxiliary)
     integer, public, save, allocatable  :: index_t(:)

! memory address index for the imaginary time \tau
     integer, public, save, allocatable  :: index_v(:)

! to record type of operators, 1 means create operators, 0 means destroy operators
     integer, public, save, allocatable  :: type_v(:)

! to record flavor of operators, from 1 to norbs
     integer, public, save, allocatable  :: flvr_v(:)

! imaginary time \tau for create and destroy operators
     real(dp), public, save, allocatable :: time_v(:)

! exp(-H\tau), exponent matrix for local hamiltonian multiply \tau (the last point)
     real(dp), public, save, allocatable :: expt_t(:,:)

! exp(-H\tau), exponent matrix for local hamiltonian multiply \tau
     real(dp), public, save, allocatable :: expt_v(:,:)

  end module ctqmc_flvr

!=========================================================================
!>>> module ctqmc_umat                                                 <<<
!=========================================================================
!>>> containing util-matrix related arrays used by continuous time quantum
! Monte Carlo quantum impurity solver
  module ctqmc_umat
     use constants, only : dp

     implicit none

! histogram for perturbation expansion series
     integer,  public, save, allocatable :: hist(:)

! current perturbation expansion order for different flavor channel
     integer,  public, save, allocatable :: rank(:)

!-------------------------------------------------------------------------
!::: input data variables                                              :::
!-------------------------------------------------------------------------
! symmetry properties for correlated orbitals
     integer,  public, save, allocatable :: symm(:)

! impurity level for correlated orbitals
     real(dp), public, save, allocatable :: eimp(:)

! eigenvalues for local hamiltonian matrix
     real(dp), public, save, allocatable :: eigs(:)

! occupation number for the eigenstates of local hamiltonian matrix
     real(dp), public, save, allocatable :: naux(:)

! total spin for the eigenstates of local hamiltonian matrix
     real(dp), public, save, allocatable :: saux(:)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-------------------------------------------------------------------------
!::: physical observables                                              :::
!-------------------------------------------------------------------------
! probability of eigenstates of local hamiltonian matrix
     real(dp), public, save, allocatable :: prob(:)

! auxiliary physical observables
! paux(1) : total energy, Etot
! paux(2) : potential engrgy, Epot
! paux(3) : kinetic energy, Ekin
! paux(4) : magnetic moment, < Sz >
     real(dp), public, save, allocatable :: paux(:)

! impurity occupation number, < n_i >
     real(dp), public, save, allocatable :: nmat(:)

! impurity double occupation number matrix, < n_i n_j >
     real(dp), public, save, allocatable :: nnmat(:,:)

! diagonal elements of current matrix product of flavor part
! it is used to calculate the probability of eigenstates
     real(dp), public, save, allocatable :: ddmat(:,:)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!-------------------------------------------------------------------------
!::: mesh data variables                                               :::
!-------------------------------------------------------------------------
! imaginary time mesh
     real(dp), public, save, allocatable :: tmesh(:)

! real matsubara frequency mesh
     real(dp), public, save, allocatable :: rmesh(:)

! complex matsubara frequency mesh
     complex(dp), public, save, allocatable :: cmesh(:)
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

! identity matrix
     complex(dp), public, save, allocatable :: unity(:,:)

  end module ctqmc_umat

!=========================================================================
!>>> module ctqmc_mmat                                                 <<<
!=========================================================================
!>>> containing M-matrix and G-matrix related arrays used by continuous
! time quantum Monte Carlo quantum impurity solver
  module ctqmc_mmat
     use constants, only : dp

     implicit none

! helper matrix for evaluating M & G matrices
     real(dp), public, save, allocatable    :: lspace(:,:)

! helper matrix for evaluating M & G matrices
     real(dp), public, save, allocatable    :: rspace(:,:)

! M matrix, $ \mathscr{M} $
     real(dp), public, save, allocatable    :: mmat(:,:,:)

! helper matrix for evaluating G matrix
     complex(dp), public, save, allocatable :: lsaves(:,:)

! helper matrix for evaluating G matrix
     complex(dp), public, save, allocatable :: rsaves(:,:)

! G matrix, $ \mathscr{G} $
     complex(dp), public, save, allocatable :: gmat(:,:,:)

  end module ctqmc_mmat

!=========================================================================
!>>> module ctqmc_gmat                                                 <<<
!=========================================================================
!>>> containing green's function matrix related arrays used by continuous
! time quantum Monte Carlo quantum impurity solver
  module ctqmc_gmat
     use constants, only : dp

     implicit none

! impurity green's function, in imaginary time axis, matrix form
     real(dp), public, save, allocatable    :: gtau(:,:,:)

! impurity green's function, in matsubara frequency axis, matrix form
     complex(dp), public, save, allocatable :: grnf(:,:,:)

  end module ctqmc_gmat

!=========================================================================
!>>> module ctqmc_wmat                                                 <<<
!=========================================================================
!>>> containing weiss's function and hybridization function matrix related
! arrays used by continuous time quantum Monte Carlo quantum impurity solver
  module ctqmc_wmat
     use constants, only : dp

     implicit none

! bath weiss's function, in imaginary time axis, matrix form
     real(dp), public, save, allocatable    :: wtau(:,:,:)

! bath weiss's function, in matsubara frequency axis, matrix form
     complex(dp), public, save, allocatable :: wssf(:,:,:)

! hybridization function, in imaginary time axis, matrix form
     real(dp), public, save, allocatable    :: htau(:,:,:)

! hybridization function, in matsubara frequency axis, matrix form
     complex(dp), public, save, allocatable :: hybf(:,:,:)

! second order derivates for hybridization function, used to interpolate htau
     real(dp), public, save, allocatable    :: hsed(:,:,:)

  end module ctqmc_wmat

!=========================================================================
!>>> module ctqmc_smat                                                 <<<
!=========================================================================
!>>> containing self-energy function matrix related arrays used by
! continuous time quantum Monte Carlo quantum impurity solver
  module ctqmc_smat
     use constants, only : dp

     implicit none

! self-energy function, in matsubara frequency axis, matrix form
     complex(dp), public, save, allocatable :: sig1(:,:,:)

! self-energy function, in matsubara frequency axis, matrix form
     complex(dp), public, save, allocatable :: sig2(:,:,:)

  end module ctqmc_smat

!=========================================================================
!>>> module ctqmc_sect                                                 <<<
!=========================================================================
!>>> containing the sector information for good quantum number algorithm
  module ctqmc_sect
     use constants, only: dp
     use m_sector
 
     implicit none

! the total number of sectors
     integer, public, save :: nsectors

! the max dimension of the sectors
     integer, public, save :: max_dim_sect

! the average dimension of the sectors
     real(dp), public, save :: ave_dim_sect

! the array contains all the sectors
     type(t_sector), public, save, allocatable :: sectors(:)

! whether this sector forming a string
     logical, public, save, allocatable :: is_string(:,:)

! how to treat each part when calculate trace
     integer, public, save, allocatable :: is_save(:,:)

! part index
     integer, public, save, allocatable :: part_indx(:,:)

! nop, ops, ope
     integer, public, save, allocatable :: nop(:)
     integer, public, save, allocatable :: ops(:)
     integer, public, save, allocatable :: ope(:)

! saved parts of matrices product, previous configuration 
     real(dp), public, save, allocatable :: saved_a(:,:,:,:)

! start and end index of sectors for saved parts of matrices product, 
! previous configuration
     integer,  public, save, allocatable :: saved_a_nm(:,:,:)

! saved parts of matrices product, current configuration
     real(dp), public, save, allocatable :: saved_b(:,:,:,:)

! start and end index of sectors for saved parts of matrices product, 
! current configuration
     integer,  public, save, allocatable :: saved_b_nm(:,:,:)

  end module ctqmc_sect

!=========================================================================
!>>> module context                                                    <<<
!=========================================================================
!>>> containing memory management subroutines and define global variables
  module context
     use constants
     use control

     use ctqmc_core
     use ctqmc_clur
     use ctqmc_flvr

     use ctqmc_umat
     use ctqmc_mmat

     use ctqmc_gmat
     use ctqmc_wmat
     use ctqmc_smat

     use ctqmc_sect
     implicit none

! status flag
     integer, private :: istat

! declaration of module procedures: allocate memory
     public :: ctqmc_allocate_memory_clur
     public :: ctqmc_allocate_memory_flvr
     public :: ctqmc_allocate_memory_umat
     public :: ctqmc_allocate_memory_mmat
     public :: ctqmc_allocate_memory_gmat
     public :: ctqmc_allocate_memory_wmat
     public :: ctqmc_allocate_memory_smat
     public :: ctqmc_allocate_memory_sect

! declaration of module procedures: deallocate memory
     public :: ctqmc_deallocate_memory_clur
     public :: ctqmc_deallocate_memory_flvr
     public :: ctqmc_deallocate_memory_umat
     public :: ctqmc_deallocate_memory_mmat
     public :: ctqmc_deallocate_memory_gmat
     public :: ctqmc_deallocate_memory_wmat
     public :: ctqmc_deallocate_memory_smat
     public :: ctqmc_deallocate_memory_sect

     contains

!=========================================================================
!>>> allocate memory subroutines                                       <<<
!=========================================================================

!>>> allocate memory for clur-related variables
     subroutine ctqmc_allocate_memory_clur()
         implicit none

! loop index
         integer :: i

! allocate memory
         allocate(index_s(mkink,norbs),     stat=istat)
         allocate(index_e(mkink,norbs),     stat=istat)

         allocate(time_s(mkink,norbs),      stat=istat)
         allocate(time_e(mkink,norbs),      stat=istat)

         allocate(exp_s(nfreq,mkink,norbs), stat=istat)
         allocate(exp_e(nfreq,mkink,norbs), stat=istat)

         allocate(empty_s(norbs),           stat=istat)
         allocate(empty_e(norbs),           stat=istat)

! check the status
         if ( istat /= 0 ) then
             call ctqmc_print_error('ctqmc_allocate_memory_clur','can not allocate enough memory')
         endif

! initialize them
         index_s = 0
         index_e = 0

         time_s  = zero
         time_e  = zero

         exp_s   = czero
         exp_e   = czero

         do i=1,norbs
             empty_s(i) = istack_create(mkink)
             empty_e(i) = istack_create(mkink)
         enddo ! over i={1,norbs} loop

         return
     end subroutine ctqmc_allocate_memory_clur

!>>> allocate memory for flvr-related variables
     subroutine ctqmc_allocate_memory_flvr()
         implicit none

! allocate memory
         allocate(index_t(mkink),      stat=istat)
         allocate(index_v(mkink),      stat=istat)

         allocate(type_v(mkink),       stat=istat)
         allocate(flvr_v(mkink),       stat=istat)

         allocate(time_v(mkink),       stat=istat)

         allocate(expt_t(ncfgs,  4  ), stat=istat)
         allocate(expt_v(ncfgs,mkink), stat=istat)

! check the status
         if ( istat /= 0 ) then
             call ctqmc_print_error('ctqmc_allocate_memory_flvr','can not allocate enough memory')
         endif

! initialize them
         index_t = 0
         index_v = 0

         type_v  = 1
         flvr_v  = 1

         time_v  = zero

         expt_t  = zero
         expt_v  = zero

         empty_v = istack_create(mkink)

         return
     end subroutine ctqmc_allocate_memory_flvr

!>>> allocate memory for umat-related variables
     subroutine ctqmc_allocate_memory_umat()
         implicit none

! allocate memory
         allocate(hist(mkink),        stat=istat)
         allocate(rank(norbs),        stat=istat)

         allocate(symm(norbs),        stat=istat)

         allocate(eimp(norbs),        stat=istat)
         allocate(eigs(ncfgs),        stat=istat)
         allocate(naux(ncfgs),        stat=istat)
         allocate(saux(ncfgs),        stat=istat)

         allocate(prob(ncfgs),        stat=istat)
         allocate(paux(  4  ),        stat=istat)
         allocate(nmat(norbs),        stat=istat)

         allocate(nnmat(norbs,norbs), stat=istat)
         allocate(ddmat(ncfgs,  2  ), stat=istat)

         allocate(tmesh(ntime),       stat=istat)
         allocate(rmesh(mfreq),       stat=istat)

         allocate(cmesh(mfreq),       stat=istat)

         allocate(unity(norbs,norbs), stat=istat)

! check the status
         if ( istat /= 0 ) then
             call ctqmc_print_error('ctqmc_allocate_memory_umat','can not allocate enough memory')
         endif

! initialize them
         hist  = 0
         rank  = 0

         symm  = 0

         eimp  = zero
         eigs  = zero
         naux  = zero
         saux  = zero

         prob  = zero
         paux  = zero
         nmat  = zero

         nnmat = zero
         ddmat = zero

         tmesh = zero
         rmesh = zero

         cmesh = czero

         unity = czero

         return
     end subroutine ctqmc_allocate_memory_umat

!>>> allocate memory for mmat-related variables
     subroutine ctqmc_allocate_memory_mmat()
         implicit none

! allocate memory
         allocate(lspace(mkink,norbs),     stat=istat)
         allocate(rspace(mkink,norbs),     stat=istat)

         allocate(mmat(mkink,mkink,norbs), stat=istat)

         allocate(lsaves(nfreq,norbs),     stat=istat)
         allocate(rsaves(nfreq,norbs),     stat=istat)

         allocate(gmat(nfreq,norbs,norbs), stat=istat)

! check the status
         if ( istat /= 0 ) then
             call ctqmc_print_error('ctqmc_allocate_memory_mmat','can not allocate enough memory')
         endif

! initialize them
         lspace = zero
         rspace = zero

         mmat   = zero

         lsaves = czero
         rsaves = czero

         gmat   = czero

         return
     end subroutine ctqmc_allocate_memory_mmat

!>>> allocate memory for gmat-related variables
     subroutine ctqmc_allocate_memory_gmat()
         implicit none

! allocate memory
         allocate(gtau(ntime,norbs,norbs), stat=istat)

         allocate(grnf(mfreq,norbs,norbs), stat=istat)

! check the status
         if ( istat /= 0 ) then
             call ctqmc_print_error('ctqmc_allocate_memory_gmat','can not allocate enough memory')
         endif

! initialize them
         gtau = zero

         grnf = czero

         return
     end subroutine ctqmc_allocate_memory_gmat

!>>> allocate memory for wmat-related variables
     subroutine ctqmc_allocate_memory_wmat()
         implicit none

! allocate memory
         allocate(wtau(ntime,norbs,norbs), stat=istat)
         allocate(htau(ntime,norbs,norbs), stat=istat)
         allocate(hsed(ntime,norbs,norbs), stat=istat)

         allocate(wssf(mfreq,norbs,norbs), stat=istat)
         allocate(hybf(mfreq,norbs,norbs), stat=istat)

! check the status
         if ( istat /= 0 ) then
             call ctqmc_print_error('ctqmc_allocate_memory_wmat','can not allocate enough memory')
         endif

! initialize them
         wtau = zero
         htau = zero
         hsed = zero

         wssf = czero
         hybf = czero

         return
     end subroutine ctqmc_allocate_memory_wmat

!>>> allocate memory for smat-related variables
     subroutine ctqmc_allocate_memory_smat()
         implicit none

! allocate memory
         allocate(sig1(mfreq,norbs,norbs), stat=istat)
         allocate(sig2(mfreq,norbs,norbs), stat=istat)

! check the status
         if ( istat /= 0 ) then
             call ctqmc_print_error('ctqmc_allocate_memory_smat','can not allocate enough memory')
         endif

! initialize them
         sig1 = czero
         sig2 = czero

         return
     end subroutine ctqmc_allocate_memory_smat

!>>> allocate memory for sect-related variables
     subroutine ctqmc_allocate_memory_sect()
         implicit none

! local variables
         integer :: i

! allocate memory
         allocate(sectors(nsectors),              stat=istat)
         allocate(is_string(nsectors,2),          stat=istat)
         allocate(is_save(npart, nsectors),       stat=istat)
         allocate(part_indx(npart, nsectors),     stat=istat)
         allocate(nop(npart),                     stat=istat)
         allocate(ops(npart),                     stat=istat)
         allocate(ope(npart),                     stat=istat)
         allocate(saved_a_nm(2, npart, nsectors), stat=istat)
         allocate(saved_b_nm(2, npart, nsectors), stat=istat)
         allocate(saved_a(max_dim_sect, max_dim_sect, npart, nsectors), stat=istat)
         allocate(saved_b(max_dim_sect, max_dim_sect, npart, nsectors), stat=istat)

! check the status
         if ( istat /= 0 ) then
             call ctqmc_print_error('ctqmc_allocate_memory_sect','can not allocate enough memory')
         endif

! initialize them
         do i=1, nsectors
             sectors(i)%ndim = 0
             sectors(i)%nelectron = 0
             sectors(i)%nops = norbs
             sectors(i)%istart = 0
             call nullify_one_sector(sectors(i))
         enddo 

         is_string = .false.
         is_save = 1
         part_indx = -1
         nop = 0
         ops = 0
         ope = 0
         saved_a = zero
         saved_b = zero
         saved_a_nm = 0
         saved_b_nm = 0

         return
     end subroutine ctqmc_allocate_memory_sect

!=========================================================================
!>>> deallocate memory subroutines                                     <<<
!=========================================================================

!>>> deallocate memory for clur-related variables
     subroutine ctqmc_deallocate_memory_clur()
         implicit none

! loop index
         integer :: i

         do i=1,norbs
             call istack_destroy(empty_s(i))
             call istack_destroy(empty_e(i))
         enddo ! over i={1,norbs} loop

         if ( allocated(index_s) ) deallocate(index_s)
         if ( allocated(index_e) ) deallocate(index_e)

         if ( allocated(time_s)  ) deallocate(time_s )
         if ( allocated(time_e)  ) deallocate(time_e )

         if ( allocated(exp_s)   ) deallocate(exp_s  )
         if ( allocated(exp_e)   ) deallocate(exp_e  )

         if ( allocated(empty_s) ) deallocate(empty_s)
         if ( allocated(empty_e) ) deallocate(empty_e)

         return
     end subroutine ctqmc_deallocate_memory_clur

!>>> deallocate memory for flvr-related variables
     subroutine ctqmc_deallocate_memory_flvr()
         implicit none

         call istack_destroy(empty_v)

         if ( allocated(index_t) ) deallocate(index_t)
         if ( allocated(index_v) ) deallocate(index_v)

         if ( allocated(type_v)  ) deallocate(type_v )
         if ( allocated(flvr_v)  ) deallocate(flvr_v )

         if ( allocated(time_v)  ) deallocate(time_v )

         if ( allocated(expt_t)  ) deallocate(expt_t )
         if ( allocated(expt_v)  ) deallocate(expt_v )

         return
     end subroutine ctqmc_deallocate_memory_flvr

!>>> deallocate memory for umat-related variables
     subroutine ctqmc_deallocate_memory_umat()
         implicit none

         if ( allocated(hist)  )   deallocate(hist )
         if ( allocated(rank)  )   deallocate(rank )

         if ( allocated(symm)  )   deallocate(symm )

         if ( allocated(eimp)  )   deallocate(eimp )
         if ( allocated(eigs)  )   deallocate(eigs )
         if ( allocated(naux)  )   deallocate(naux )
         if ( allocated(saux)  )   deallocate(saux )

         if ( allocated(prob)  )   deallocate(prob )
         if ( allocated(paux)  )   deallocate(paux )
         if ( allocated(nmat)  )   deallocate(nmat )

         if ( allocated(nnmat) )   deallocate(nnmat)
         if ( allocated(ddmat) )   deallocate(ddmat)

         if ( allocated(tmesh) )   deallocate(tmesh)
         if ( allocated(rmesh) )   deallocate(rmesh)

         if ( allocated(cmesh) )   deallocate(cmesh)

         if ( allocated(unity) )   deallocate(unity)

         return
     end subroutine ctqmc_deallocate_memory_umat


!>>> deallocate memory for mmat-related variables
     subroutine ctqmc_deallocate_memory_mmat()
         implicit none

         if ( allocated(lspace) )  deallocate(lspace)
         if ( allocated(rspace) )  deallocate(rspace)

         if ( allocated(mmat)   )  deallocate(mmat  )

         if ( allocated(lsaves) )  deallocate(lsaves)
         if ( allocated(rsaves) )  deallocate(rsaves)

         if ( allocated(gmat)   )  deallocate(gmat  )

         return
     end subroutine ctqmc_deallocate_memory_mmat

!>>> deallocate memory for gmat-related variables
     subroutine ctqmc_deallocate_memory_gmat()
         implicit none

         if ( allocated(gtau) )    deallocate(gtau)

         if ( allocated(grnf) )    deallocate(grnf)

         return
     end subroutine ctqmc_deallocate_memory_gmat

!>>> deallocate memory for wmat-related variables
     subroutine ctqmc_deallocate_memory_wmat()
         implicit none

         if ( allocated(wtau) )    deallocate(wtau)
         if ( allocated(htau) )    deallocate(htau)
         if ( allocated(hsed) )    deallocate(hsed)

         if ( allocated(wssf) )    deallocate(wssf)
         if ( allocated(hybf) )    deallocate(hybf)

         return
     end subroutine ctqmc_deallocate_memory_wmat

!>>> deallocate memory for smat-related variables
     subroutine ctqmc_deallocate_memory_smat()
         implicit none

         if ( allocated(sig1) )    deallocate(sig1)
         if ( allocated(sig2) )    deallocate(sig2)

         return
     end subroutine ctqmc_deallocate_memory_smat

!>>> deallocate memory for sect-related variables
     subroutine ctqmc_deallocate_memory_sect()
         implicit none

! local variables
         integer :: i

         if ( allocated(sectors) ) then
! first, loop over all the sectors and deallocate their memory
             do i=1, nsectors
                 call dealloc_one_sector(sectors(i))
             enddo
! then, deallocate memory of sect
             deallocate(sectors)
         endif

         if ( allocated(is_string) )    deallocate(is_string)
         if ( allocated(is_save) )      deallocate(is_save)
         if ( allocated(nop) )          deallocate(nop)
         if ( allocated(ops) )          deallocate(ops)
         if ( allocated(ope) )          deallocate(ope)
         if ( allocated(part_indx) )    deallocate(part_indx)
         if ( allocated(saved_a) )      deallocate(saved_a)
         if ( allocated(saved_a_nm) )   deallocate(saved_a_nm)
         if ( allocated(saved_b) )      deallocate(saved_b)
         if ( allocated(saved_b_nm) )   deallocate(saved_b_nm)
        
         return
     end subroutine ctqmc_deallocate_memory_sect

  end module context
