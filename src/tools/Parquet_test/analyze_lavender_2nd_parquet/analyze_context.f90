!-------------------------------------------------------------------------
! project : analyze_lavender-soc
! program : analyze_gmat module
!           context    module
! source  : analyze_context.f90
! type    : module
! author  : li huang (email:huangli712@yahoo.com.cn)
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
!           10/27/2013 by Kuang-Shing Chen 
! purpose : define the key data structure and global arrays/variables for
!           hybridization expansion version continuous time quantum Monte
!           Carlo (analyze) quantum impurity solver and dynamical mean field
!           theory (DMFT) self-consistent engine
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------


!=========================================================================
!>>> module analyze_gmat                                                 <<<
!=========================================================================
!>>> containing green's function matrix related arrays used by continuous
! time quantum Monte Carlo quantum impurity solver
  module analyze_gmat
     use constants, only : dp

     implicit none

! real matsubara frequency mesh
     real(dp), public, save, allocatable :: rmesh(:)

! complex matsubara frequency mesh
     complex(dp), public, save, allocatable :: cmesh(:)

! impurity green's function, in matsubara frequency axis, matrix form
     complex(dp), public, save, allocatable :: grnf(:,:)
     complex(dp), public, save, allocatable :: sigf(:,:)

! lattice green's function, in matsubara frequency axis, matrix form
     complex(dp), public, save, allocatable :: grnf_lattice(:,:,:)

  end module analyze_gmat

!=========================================================================
!>>> module analyze_two_ptle_mat                                       <<<
 !=========================================================================
!>>> containing 2-tple green's function matrix related arrays 
  module analyze_two_ptle_mat
     use constants, only : dp

     implicit none

!   two-particle green's function, s:singlet; t:triplet
     complex(dp), public, save, allocatable :: g2_pp_s(:,:,:)
!!     complex(sp), public, save, allocatable :: g2_pp_t(:,:,:)
!!     complex(dp), public, save, allocatable :: g2_ph_c(:,:,:)
!!     complex(dp), public, save, allocatable :: g2_ph_m(:,:,:)

!   two-particle vertex function, s:singlet; t:triplet
!!     complex(dp), public, save, allocatable :: Fc_pp_s(:,:,:)
!!     complex(sp), public, save, allocatable :: Fc_pp_t(:,:,:)
     complex(dp), public, save, allocatable :: Fc_ph_c(:,:,:)
     complex(dp), public, save, allocatable :: Fc_ph_m(:,:,:)
     complex(dp), public, save, allocatable :: gamma_ph_c(:,:,:)
     complex(dp), public, save, allocatable :: gamma_ph_m(:,:,:)
     complex(dp), public, save, allocatable :: gamma_pp_s(:,:,:)
     complex(dp), public, save, allocatable :: gamma_pp_t(:,:,:)


! lattice bare bubble, ph, in matsubara frequency axis, matrix form
!     complex(dp), public, save, allocatable :: chi0tilde_dm(:,:,:,:)
     complex(dp), public, save, allocatable :: chi0_dm(:,:,:,:)
     complex(dp), public, save, allocatable :: chi0_st(:,:,:,:)
     
! vertex ladder, ph, in momentum and matsubara frequency axis, matrix form
     complex(dp), public, save, allocatable :: phi_local_c(:,:,:)
     complex(dp), public, save, allocatable :: phi_local_m(:,:,:)
     complex(dp), public, save, allocatable :: phi_local_s(:,:,:)
     complex(dp), public, save, allocatable :: phi_local_t(:,:,:)
     complex(dp), public, save, allocatable :: phi_ph_c(:,:,:,:,:)
     complex(dp), public, save, allocatable :: phi_ph_m(:,:,:,:,:)
     complex(dp), public, save, allocatable :: phi_pp_s(:,:,:,:,:)
     complex(dp), public, save, allocatable :: phi_pp_t(:,:,:,:,:)

! Fully irreducible vertex lambda
     complex(dp), public, save, allocatable :: lambda_s(:,:)
     complex(dp), public, save, allocatable :: lambda_t(:,:)
     complex(dp), public, save, allocatable :: lambda_c(:,:)
     complex(dp), public, save, allocatable :: lambda_m(:,:)

! irreducible singlet vertex Gamma
     complex(dp), public, save, allocatable :: gamat_s(:,:)

! irreducible density/magnetic vertex Gamma_m, 
     complex(dp), public, save, allocatable :: gamat_m(:,:)
     complex(dp), public, save, allocatable :: gamat_c(:,:)

  end module analyze_two_ptle_mat


!=========================================================================
!>>> module context                                                    <<<
!=========================================================================
!>>> containing memory management subroutines and define global variables
  module context
     use constants
     use control

!     use analyze_core

!     use analyze_umat
!     use analyze_fmat
!     use analyze_mmat

     use analyze_gmat
     use analyze_two_ptle_mat
!     use analyze_wmat
!     use analyze_smat

     implicit none

! status flag
     integer, private :: istat

! declaration of module procedures: allocate memory
!     public :: analyze_allocate_memory_umat
!     public :: analyze_allocate_memory_fmat
!     public :: analyze_allocate_memory_mmat
     public :: analyze_allocate_memory_gmat
     public :: analyze_allocate_memory_two_mat_local
     public :: analyze_allocate_memory_two_mat_non_local
!     public :: analyze_allocate_memory_wmat
!     public :: analyze_allocate_memory_smat

! declaration of module procedures: deallocate memory
!     public :: analyze_deallocate_memory_umat
!     public :: analyze_deallocate_memory_fmat
!     public :: analyze_deallocate_memory_mmat
!!     public :: analyze_deallocate_useless_array0
!!     public :: analyze_deallocate_useless_array1
     public :: analyze_deallocate_memory_gmat
     public :: analyze_deallocate_memory_two_mat_local
     public :: analyze_deallocate_memory_two_mat_non_local
!     public :: analyze_deallocate_memory_wmat
!     public :: analyze_deallocate_memory_smat

     contains




!>>> allocate memory for gmat-related variables
     subroutine analyze_allocate_memory_gmat()
         implicit none
! Single particle frequency allocation example:
! nffrq =8, nbfrq=8:
! n:          -3 -2 -1  0 1 2 3 4 
! w_n/(pi*T): -7 -5 -3 -1 1 3 5 7
! m:           -7  -6  -5 -4 -3 -2 -1 0 1 2 3 4  5  6  7
! nu_m/(pi*T):-14 -12 -10 -8 -6 -4 -2 0 2 4 6 8 10 12 14
! Range of n+m is from (-nffrq/2 +1)-(nbfrq-1) to nffrq/2+(nbfrq-1), -10 to 11.
! allocate memory
         allocate(rmesh(-nffrq/2+1-(nbfrq-1):nffrq/2+(nbfrq-1)), stat=istat)
         allocate(cmesh(1:nffrq/2+(nbfrq-1)), stat=istat)
         allocate(grnf(-nffrq/2+1-(nbfrq-1):nffrq/2+(nbfrq-1),2), stat=istat)
         allocate(sigf(-nffrq/2+1-(nbfrq-1):nffrq/2+(nbfrq-1),2), stat=istat)
         allocate(grnf_lattice(-ktime+1:ktime,-ktime+1:ktime, &
                  -nffrq/2+1-(nbfrq-1):nffrq/2+(nbfrq-1)), stat=istat)
! check the status
         if ( istat /= 0 ) then
             call analyze_print_error('analyze_allocate_memory_gmat','can not allocate enough memory')
         endif

! initialize them
         rmesh = zero
         cmesh = czero
         grnf = czero
         sigf = czero
         grnf_lattice = czero

         return
     end subroutine analyze_allocate_memory_gmat

!>>> allocate memory for two-ptle-related variables
     subroutine analyze_allocate_memory_two_mat_local()
     implicit none

! allocate memory
         allocate(g2_pp_s(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2,-nbfrq+1:nbfrq-1), stat=istat)
!!         allocate(g2_pp_t(nffrq,nffrq,nbfrq), stat=istat)
!!         allocate(g2_ph_c(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2,-nbfrq+1:nbfrq-1), stat=istat)
!!         allocate(g2_ph_m(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2,-nbfrq+1:nbfrq-1), stat=istat)
!!        allocate(Fc_pp_s(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2,-nbfrq+1:nbfrq-1), stat=istat)
!!         allocate(Fc_pp_t(nffrq,nffrq,nbfrq), stat=istat)
         allocate(Fc_ph_c(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2,-nbfrq+1:nbfrq-1), stat=istat)
         allocate(Fc_ph_m(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2,-nbfrq+1:nbfrq-1), stat=istat)
         allocate(gamma_ph_c(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2,-nbfrq+1:nbfrq-1), stat=istat)
         allocate(gamma_ph_m(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2,-nbfrq+1:nbfrq-1), stat=istat)
         allocate(gamma_pp_s(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2,-nbfrq+1:nbfrq-1), stat=istat)
         allocate(gamma_pp_t(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2,-nbfrq+1:nbfrq-1), stat=istat)

	allocate(phi_local_c(-nbfrq+1:nbfrq-1,-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2), stat=istat)
	allocate(phi_local_m(-nbfrq+1:nbfrq-1,-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2), stat=istat)
	allocate(phi_local_s(-nbfrq+1:nbfrq-1,-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2), stat=istat)
	allocate(phi_local_t(-nbfrq+1:nbfrq-1,-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2), stat=istat)
	
	allocate(lambda_s(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2), stat=istat)
	allocate(lambda_t(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2), stat=istat)
	allocate(lambda_c(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2), stat=istat)
	allocate(lambda_m(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2), stat=istat)
	
	allocate(chi0_dm(-ktime+1:ktime,-ktime+1:ktime, &
           -nffrq/2+1:nffrq/2,-nbfrq+1:nbfrq-1), stat=istat)
        allocate(chi0_st(-ktime+1:ktime,-ktime+1:ktime, &
           -nffrq/2+1:nffrq/2,-nbfrq+1:nbfrq-1), stat=istat) 
           
         allocate(phi_ph_c(-ktime+1:ktime,-ktime+1:ktime, &
          -nbfrq+1:nbfrq-1,-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2), stat=istat)
         allocate(phi_ph_m(-ktime+1:ktime,-ktime+1:ktime, &
          -nbfrq+1:nbfrq-1,-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2), stat=istat)
	 allocate(phi_pp_s(-ktime+1:ktime,-ktime+1:ktime, &
          -nbfrq+1:nbfrq-1,-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2), stat=istat)
         allocate(phi_pp_t(-ktime+1:ktime,-ktime+1:ktime, &
          -nbfrq+1:nbfrq-1,-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2), stat=istat)

! check the status
         if ( istat /= 0 ) then
             write(*,*) "istat=",istat
             call analyze_print_error('analyze_allocate_memory_two_mat','can not allocate enough memory')
         endif

! initialize them
         g2_pp_s = czero
!!         g2_pp_t = czero
!!         g2_ph_c = czero
!!         g2_ph_m = czero
!!         Fc_pp_s = czero
!!         Fc_pp_t = czero
         Fc_ph_c = czero
         Fc_ph_m = czero
         gamma_ph_c = czero
         gamma_ph_m = czero
         gamma_pp_s = czero
         gamma_pp_t = czero
         
         phi_local_c = czero
         phi_local_m = czero
         phi_local_s = czero
         phi_local_t = czero
   
         lambda_s = czero
         lambda_t = czero
         lambda_c = czero
         lambda_m = czero

         chi0_dm = czero
         chi0_st = czero
         phi_ph_c = czero
         phi_ph_m = czero
         phi_pp_s = czero
         phi_pp_t = czero

     return
     end subroutine analyze_allocate_memory_two_mat_local
     
     
!>>> allocate memory for non-local two-ptle-related variables
     subroutine analyze_allocate_memory_two_mat_non_local()
     implicit none
  
	 integer :: nt_2p_q
! allocate memory
!!	 nt_2p_q = 4*ktime*ktime*(2*nbfrq-1)
!!	we only focus on the Q=(q,v=0)
!!	 nt_2p_q = 4*ktime*ktime
	 allocate(gamat_m(4*ktime*ktime*nffrq,4*ktime*ktime*nffrq), stat=istat)
	 allocate(gamat_c(4*ktime*ktime*nffrq,4*ktime*ktime*nffrq), stat=istat)
         allocate(gamat_s(4*ktime*ktime*nffrq,4*ktime*ktime*nffrq), stat=istat)

! check the status
         if ( istat /= 0 ) then
             write(*,*) "istat=",istat
             call analyze_print_error('analyze_allocate_memory_two_mat_non_local','can not allocate enough memory')
         endif

! initialize them
         gamat_m = czero
         gamat_c = czero
         gamat_s = czero

     return
     end subroutine analyze_allocate_memory_two_mat_non_local
!=========================================================================
!>>> deallocate memory subroutines                                     <<<
!=========================================================================

!>>> deallocate memory for gmat-related variables
     subroutine analyze_deallocate_memory_gmat()
         implicit none

         if ( allocated(rmesh) )   deallocate(rmesh)  
         if ( allocated(cmesh) )   deallocate(cmesh)
         if ( allocated(grnf) )    deallocate(grnf)
         if ( allocated(sigf) )    deallocate(sigf)
	 if ( allocated(grnf_lattice)) deallocate(grnf_lattice)

         return
     end subroutine analyze_deallocate_memory_gmat

!>>> deallocate memory for gmat-related variables
     subroutine analyze_deallocate_memory_two_mat_local()
         implicit none
 
         if ( allocated(g2_pp_s) )    deallocate(g2_pp_s)
!!         if ( allocated(g2_pp_t) )    deallocate(g2_pp_t)
!!         if ( allocated(g2_ph_c) )    deallocate(g2_ph_c)
!!         if ( allocated(g2_ph_m) )    deallocate(g2_ph_m)
!!         if ( allocated(Fc_pp_s) )    deallocate(Fc_pp_s)
!!         if ( allocated(Fc_pp_t) )    deallocate(Fc_pp_t)
         if ( allocated(Fc_ph_c) )    deallocate(Fc_ph_c)
         if ( allocated(Fc_ph_m) )    deallocate(Fc_ph_m)
         if ( allocated(gamma_ph_c) )    deallocate(gamma_ph_c)
         if ( allocated(gamma_ph_m) )    deallocate(gamma_ph_m)
         if ( allocated(gamma_pp_s) )    deallocate(gamma_pp_s)
         if ( allocated(gamma_pp_t) )    deallocate(gamma_pp_t)
         if ( allocated(phi_local_c))    deallocate(phi_local_c)
	 if ( allocated(phi_local_m))    deallocate(phi_local_m)
	 if ( allocated(phi_local_s))    deallocate(phi_local_s)
	 if ( allocated(phi_local_t))    deallocate(phi_local_t)
	  if ( allocated(chi0_st))      deallocate(chi0_st)
 
         return
     end subroutine analyze_deallocate_memory_two_mat_local
     
     !>>> deallocate memory for gmat-related variables
     subroutine analyze_deallocate_memory_two_mat_non_local()
         implicit none
 
	 if ( allocated(phi_ph_c))    deallocate(phi_ph_c)
	 if ( allocated(phi_ph_m))    deallocate(phi_ph_m)
	 if ( allocated(phi_pp_s))    deallocate(phi_pp_s)
	 if ( allocated(phi_pp_t))    deallocate(phi_pp_t)
	 if ( allocated(lambda_s))    deallocate(lambda_s)
	 if ( allocated(lambda_t))    deallocate(lambda_t)
	 if ( allocated(lambda_c))    deallocate(lambda_c)
	 if ( allocated(lambda_m))    deallocate(lambda_m)
	 if ( allocated(chi0_dm))      deallocate(chi0_dm)
	 if ( allocated(gamat_m))    deallocate(gamat_m)
	 if ( allocated(gamat_c))    deallocate(gamat_c)
	 if ( allocated(gamat_s))    deallocate(gamat_s)
 
         return
     end subroutine analyze_deallocate_memory_two_mat_non_local

  end module context
