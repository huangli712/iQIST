

!!========================================================================
!!>>> module m_sect                                                    <<<
!!========================================================================

!!>>> define the data structure for good quantum numbers (GQNs) algorithm
  module m_sect
     use constants, only : dp, zero, eps6, mystd, mytmp
     use mmpi, only : mp_bcast, mp_barrier

     use control, only : norbs, ncfgs
     use control, only : mkink
     use control, only : myid, master
     use context, only : type_v, flvr_v

     implicit none

!!========================================================================
!!>>> declare global structures                                        <<<
!!========================================================================

! data structure for one F-matrix
!-------------------------------------------------------------------------
     public :: t_fmat
     type t_fmat

! the dimension, n x m
         integer :: n
         integer :: m

! the memory space for the matrix
         real(dp), allocatable :: val(:,:)

     end type t_fmat

! data structure for one sector
!-------------------------------------------------------------------------
     public :: t_sector
     type t_sector

! number of states in this sector
         integer :: ndim

! number of fermion operators, it should be equal to norbs
         integer :: nops

! start index of this sector
         integer :: istart

! total number of electrons
         integer :: nele

! z component of spin: Sz
         integer :: sz

! z component of spin-orbit momentum: Jz
         integer :: jz

! PS good quantum number
         integer :: ps

! the next sector when a fermion operator acts on the sector
! next(nops,0) for annihilation and next(nops,1) for creation operators
! -1 means it is outside the Hilbert space,
! otherwise, it is the index of next sector
         integer, allocatable  :: next(:,:)

! the eigenvalues
         real(dp), allocatable :: eval(:)

! final products of matrices
         real(dp), allocatable :: prod(:)

! the F-matrix between this sector and all other sectors
! fmat(nops,0) for annihilation and fmat(nops,1) for creation operators
! if this sector doesn't point to some other sectors, it is not allocated
         type (t_fmat), allocatable :: fmat(:,:)

     end type t_sector

!!========================================================================
!!>>> declare global variables                                         <<<
!!========================================================================

! total number of sectors
     integer, public, save  :: nsect

! maximal dimension of the sectors
     integer, public, save  :: max_dim_sect

! average dimension of the sectors
     real(dp), public, save :: ave_dim_sect

! which sectors should be truncated?
     logical, public, save, allocatable :: sectoff(:)

! array of t_sector contains all the sectors
     type (t_sector), public, save, allocatable :: sectors(:)

!!========================================================================
!!>>> declare private variables                                        <<<
!!========================================================================

! status flag
     integer, private :: istat

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     public :: cat_alloc_one_fmat
     public :: cat_alloc_one_sect
     public :: cat_alloc_sect

     public :: cat_free_one_fmat
     public :: cat_free_one_sect
     public :: cat_free_sect

     public :: cat_make_string
     public :: cat_trun_sector

  contains ! encapsulated functionality

!!========================================================================
!!>>> allocate memory subroutines                                      <<<
!!========================================================================

!!>>> cat_alloc_one_fmat: allocate memory for one F-matrix
  subroutine cat_alloc_one_fmat(mat)
     implicit none

! external variables
! F-matrix structure
     type (t_fmat), intent(inout) :: mat

! allocate memory
     allocate(mat%val(mat%n,mat%m), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_one_fmat','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize it
     mat%val = zero

     return
  end subroutine cat_alloc_one_fmat

!!>>> cat_alloc_one_sect: allocate memory for one sector
  subroutine cat_alloc_one_sect(sect)
     implicit none

! external variables
! sector structure
     type (t_sector), intent(inout) :: sect

! local variables
! loop index
     integer :: i
     integer :: j

! allocate memory
     allocate(sect%next(sect%nops,0:1), stat=istat)

     allocate(sect%eval(sect%ndim),     stat=istat)
     allocate(sect%prod(sect%ndim),     stat=istat)

     allocate(sect%fmat(sect%nops,0:1), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_one_sect','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     sect%next = 0

     sect%eval = zero
     sect%prod = zero

! initialize fmat one by one
     do i=1,sect%nops
         do j=0,1
             sect%fmat(i,j)%n = 0
             sect%fmat(i,j)%m = 0
         enddo ! over j={0,1} loop
     enddo ! over i={1,sect%nops} loop

     return
  end subroutine cat_alloc_one_sect

!!>>> cat_alloc_sect: allocate memory for sector related variables
  subroutine cat_alloc_sect()
     implicit none

! local variables
! loop index
     integer :: i

! allocate memory
     allocate(sectors(nsect), stat=istat)
     allocate(sectoff(nsect), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_sect','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     do i=1,nsect
         sectors(i)%ndim   = 0
         sectors(i)%nops   = norbs
         sectors(i)%istart = 0
         sectors(i)%nele   = 0
         sectors(i)%sz     = 0
         sectors(i)%jz     = 0
         sectors(i)%ps     = 0
     enddo ! over i={1,nsect} loop
     sectoff = .false.

     return
  end subroutine cat_alloc_sect

!!========================================================================
!!>>> deallocate memory subroutines                                    <<<
!!========================================================================

!!>>> cat_free_one_fmat: deallocate memory for one F-matrix
  subroutine cat_free_one_fmat(mat)
     implicit none

! external variables
! F-matrix structure
     type (t_fmat), intent(inout) :: mat

     if ( allocated(mat%val) ) deallocate(mat%val)

     return
  end subroutine cat_free_one_fmat

!!>>> cat_free_one_sect: deallocate memory for one sector
  subroutine cat_free_one_sect(sect)
     implicit none

! external variables
! sector structure
     type (t_sector), intent(inout) :: sect

! local variables
! loop index
     integer :: i
     integer :: j

     if ( allocated(sect%next) ) deallocate(sect%next)

     if ( allocated(sect%eval) ) deallocate(sect%eval)
     if ( allocated(sect%prod) ) deallocate(sect%prod)

! deallocate fmat one by one
     if ( allocated(sect%fmat) ) then
         do i=1,sect%nops
             do j=0,1
                 call cat_free_one_fmat(sect%fmat(i,j))
             enddo ! over j={0,1} loop
         enddo ! over i={1,sect%nops} loop
         deallocate(sect%fmat)
     endif ! back if ( allocated(sect%fmat) ) block

     return
  end subroutine cat_free_one_sect

!!>>> cat_free_sect: deallocate memory for sector
!!>>> related variables
  subroutine cat_free_sect()
     implicit none

! local variables
! loop index
     integer :: i

! first, loop over all the sectors and deallocate their component's memory
! then, deallocate memory of the sectors itself
     if ( allocated(sectors) ) then
         do i=1,nsect
             call cat_free_one_sect(sectors(i))
         enddo ! over i={1,nsect} loop
         deallocate(sectors)
     endif ! back if ( allocated(sectors) ) block

     if ( allocated(sectoff) )  deallocate(sectoff)

     return
  end subroutine cat_free_sect

!!========================================================================
!!>>> core service subroutines                                         <<<
!!========================================================================

!!>>> cat_make_string: it is used to build a time evolution string
  subroutine cat_make_string(csize, vindex, string)
     implicit none

! external variables
! number of fermion operators for the current diagram
     integer, intent(in)  :: csize

! memory address index of fermion operators
     integer, intent(in)  :: vindex(mkink)

! time evolution string, i.e., sequence of sector index
! if it is not a valid string, then all of its values should be -1
     integer, intent(out) :: string(csize+1,nsect)

! local variables
! loop index
     integer :: i
     integer :: j

! flavor and type of fermion operators
     integer :: vf
     integer :: vt

! current sector index and next sector index
     integer :: curr_sect
     integer :: next_sect

! init return array, we assume all of strings are invalid
     string = -1

! we try to build a string from left to right, that is, 0 -> \beta
! we assume the sectors are S1, S2, S3, ..., SM, and the fermion
! operators are F1, F2, F3, F4, .... FN. here, F1 is in \tau_1, F2
! is in \tau_2, F3 is in \tau_3, and so on, and 
!     0 < \tau_1 < \tau_2 < \tau_3 < ... < \beta
! is always guaranteed. then a typical (and also valid) string must
! look like this:
!     F1       F2       F3       F4       F5        FN
! S1 ----> S2 ----> S3 ----> S4 ----> S5 ----> ... ----> S1
! then the sequence of sector indices is the so-called string. if some
! Si are -1 (null sector), this string is invalid. we will enforce all
! elements in it to be -1. it is easy to speculate that if the number
! of fermion operators is csize, the length of string must be csize + 1
     SECTOR_SCAN_LOOP: do i=1,nsect
! setup starting sector
         curr_sect = i
         string(1,i) = curr_sect
         OPERATOR_SCAN_LOOP: do j=1,csize
! determine the type and flavor of current operator
             vt = type_v( vindex(j) )
             vf = flvr_v( vindex(j) )
! get the next sector
             next_sect = sectors(curr_sect)%next(vf,vt)
! meet null sector, it is an invalid string. we will try another
! new string
             if ( next_sect == -1 ) then
                 string(:,i) = -1; EXIT OPERATOR_SCAN_LOOP
! the string is still alive, we record the sector, and set it to
! the current sector
             else
                 string(j+1,i) = next_sect
                 curr_sect = next_sect
             endif ! back if ( next_sect == -1 ) block
         enddo OPERATOR_SCAN_LOOP ! over j={1,csize} loop
! we have to ensure that the first sector is the same with the last
! sector in this string, or else it is invalid
         if ( string(1,i) /= string(csize+1,i) ) then
             string(:,i) = -1
         endif ! back if ( string(1,i) /= string(csize+1,i) ) block
     enddo SECTOR_SCAN_LOOP ! over i={1,nsect} loop

     return
  end subroutine cat_make_string

!!>>> cat_trun_sector: it is used to truncate the Hilbert space
!!>>> of H_{loc} according to the probatility of atomic states
  subroutine cat_trun_sector()
     implicit none

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k
     integer  :: m

! used to calculate the averaged dimension of sectors
     integer  :: sum_dim

! file status
     logical  :: exists

! number of sectors after truncation
     integer  :: nsect_t

! maximal dimension of sectors after truncation
     integer  :: max_dim_sect_t

! averaged dimension of sectors after truncation
     real(dp) :: ave_dim_sect_t

! real(dp) dummy variable
     real(dp) :: rt

! probability for sector, used to do truncation
     real(dp) :: sprob(nsect)

! read file solver.prob.dat, only master node can do it
     if ( myid == master ) then
         sectoff = .false.
         sprob = zero
         inquire (file = 'solver.prob.dat', exist = exists)
! solver.prob.dat is available, we read the sector probability data
         if ( exists .eqv. .true.) then
             open(mytmp, file='solver.prob.dat', form='formatted', status='unknown')
             do i=1,ncfgs+2
                 read(mytmp,*) ! skip header
             enddo ! over i={1,ncfgs+2} loop
             do i=1,nsect
                 read(mytmp,*) m, rt, sprob(i)
             enddo ! over i={1,nsect} loop
             close(mytmp)
! determine which sector should be truncated
             do i=1,nsect
                 if ( sprob(i) < eps6 ) sectoff(i) = .true.
             enddo ! over i={1,nsect} loop
! solver.prob.dat is not available, we can not do truncation
         else
             call s_print_exception('cat_trun_sector','sector probability data are unavailable')
         endif ! back if ( exists .eqv. .true. ) block
     endif ! back if ( myid == master ) block

# if defined (MPI)

! broadcast data
     call mp_bcast(sectoff, master)

! block until all processes have reached here
     call mp_barrier()

# endif  /* MPI */

! make truncation for the sectors whose sector probabilities are too low
     max_dim_sect_t = -1
     nsect_t = 0
     sum_dim = 0
     do i=1,nsect
         if ( sectoff(i) .eqv. .true. ) then
             sectors(i)%next = -1
         else
             if ( max_dim_sect_t < sectors(i)%ndim ) then
                 max_dim_sect_t = sectors(i)%ndim
             endif ! back if ( max_dim_sect_t < sectors(i)%ndim ) block
             sum_dim = sum_dim + sectors(i)%ndim
             nsect_t = nsect_t + 1
             do j=1,sectors(i)%nops
                 do k=0,1
                     m = sectors(i)%next(j,k)
                     if ( m == -1 ) CYCLE
                     if ( sectoff(m) .eqv. .true. ) then
                         sectors(i)%next(j,k) = -1
                     endif ! back if ( sectoff(m) .eqv. .true. ) block
                 enddo ! over k={0,1} loop
             enddo ! over j={1,sectors(i)%nops} loop
         endif ! back if ( sectoff(i) .eqv. .true. ) block
     enddo ! over i={1,nsect} loop

! calculate ave_dim_sect_t
     ave_dim_sect_t = real(sum_dim) / real(nsect_t)

! print summary of sectors after truncation
     if ( myid == master ) then
         write(mystd,'(4X,a)') 'WARNING: TRUNCATION APPROXIMATION IS USED...'
         write(mystd,'(4X,a)') 'BEFORE TRUNCATION:'
         write(mystd,'(4X,a,i8)')    'tot_num_sect: ', nsect
         write(mystd,'(4X,a,i8)')    'max_dim_sect: ', max_dim_sect
         write(mystd,'(4X,a,f8.1)')  'ave_dim_sect: ', ave_dim_sect
         write(mystd,'(4X,a)') 'AFTER  TRUNCATION:'
         write(mystd,'(4X,a,i8)')    'tot_num_sect: ', nsect_t
         write(mystd,'(4X,a,i8)')    'max_dim_sect: ', max_dim_sect_t
         write(mystd,'(4X,a,f8.1)')  'ave_dim_sect: ', ave_dim_sect_t
         write(mystd,'(4X,a)') 'TRUNCATED SECTORS:'
         do i=1,nsect
             write(mystd,'(4X,a,i4,2X,a,L2)') 'index: ', i, 'status:', .not. sectoff(i)
         enddo ! over i={1,nsect} loop
         write(mystd,*)
     endif ! back if ( myid == master ) block

     return
  end subroutine cat_trun_sector

  end module m_sect

!!========================================================================
!!>>> module m_part                                                    <<<
!!========================================================================

!!>>> contains some key global variables and subroutines for divide and
!!>>> conquer algorithm to speed up the trace evaluation
  module m_part
     use constants, only : dp, zero, one

     use control, only : ncfgs
     use control, only : mkink
     use control, only : npart
     use control, only : beta
     use context, only : type_v, flvr_v, time_v, expt_v

     use m_sect, only : nsect, max_dim_sect
     use m_sect, only : sectors

     implicit none

!!========================================================================
!!>>> declare global variables                                         <<<
!!========================================================================

! number of operators for each part
     integer, public, save, allocatable  :: nop(:)

! start index of operators for each part
     integer, public, save, allocatable  :: ops(:)

! end index of operators for each part
     integer, public, save, allocatable  :: ope(:)

! how to treat each part when calculating trace
! 0: matrices product for this part has been calculated previously
! 1: this part should be recalculated, and the result must be
!    stored in saved_p, if this Monte Caro move has been accepted
     integer, public, save, allocatable  :: renew(:)

! determine which parts of saved_p are unsafe or invalid (we just call
! it asynchronization), and have to be updated (or synchronized) for
! future trace calculations
! 0: synchronous, this part of saved_p is OK
! 1: asynchronous, this part of saved_p is invalid
! Q: why is renew not enough? why do we need async and is_cp?
! A: because string is not always valid. string broken is possible. at
! that time, even renew(j) is 1, some sectors in this part will be not
! updated successfully. of course, saved_p for them will be not updated
! as well. so we have to mark the corresponding saved_p as wrong value.
! this is the role of async. due to the same reason, we cann't use renew
! to control which parts of saved_p should be updated with saved_n only.
! so we need is_cp as well.
     integer, public, save, allocatable  :: async(:,:)

! determine which parts of saved_p should be updated by the corresponding
! parts of saved_n
! 0: do nothing, saved_n and saved_p have the same values, or saved_n is
!    unavailable, we can not use it to update saved_p
! 1: saved_n will be copied to saved_p in ctqmc_make_evolve() subroutine
     integer, public, save, allocatable  :: is_cp(:,:)

! number of columns to be copied, in order to save copy time
     integer, public, save, allocatable  :: nc_cp(:,:)

! saved parts of matrices product, for previous accepted configuration
     real(dp), public, save, allocatable :: saved_p(:,:,:,:)

! saved parts of matrices product, for new proposed configuration
     real(dp), public, save, allocatable :: saved_n(:,:,:,:)

!!========================================================================
!!>>> declare private variables                                        <<<
!!========================================================================

! status flag
     integer, private :: istat

!!========================================================================
!!>>> declare accessibility for module routines                        <<<
!!========================================================================

     public :: cat_alloc_part
     public :: cat_free_part

     public :: cat_make_npart
     public :: cat_make_trace

  contains ! encapsulated functionality

!!========================================================================
!!>>> allocate memory subroutines                                      <<<
!!========================================================================

!!>>> cat_alloc_part: allocate memory for part related variables
  subroutine cat_alloc_part()
     implicit none

! allocate memory
     allocate(nop(npart),         stat=istat)
     allocate(ops(npart),         stat=istat)
     allocate(ope(npart),         stat=istat)

     allocate(renew(npart),       stat=istat)
     allocate(async(npart,nsect), stat=istat)
     allocate(is_cp(npart,nsect), stat=istat)
     allocate(nc_cp(npart,nsect), stat=istat)

     allocate(saved_p(max_dim_sect,max_dim_sect,npart,nsect), stat=istat)
     allocate(saved_n(max_dim_sect,max_dim_sect,npart,nsect), stat=istat)

! check the status
     if ( istat /= 0 ) then
         call s_print_error('cat_alloc_part','can not allocate enough memory')
     endif ! back if ( istat /= 0 ) block

! initialize them
     nop   = 0
     ops   = 0
     ope   = 0

     renew = 0
     async = 0
     is_cp = 0
     nc_cp = 0

     saved_p = zero
     saved_n = zero

     return
  end subroutine cat_alloc_part

!!========================================================================
!!>>> deallocate memory subroutines                                    <<<
!!========================================================================

!!>>> cat_free_part: deallocate memory for part related variables
  subroutine cat_free_part()
     implicit none

     if ( allocated(nop)     ) deallocate(nop    )
     if ( allocated(ops)     ) deallocate(ops    )
     if ( allocated(ope)     ) deallocate(ope    )

     if ( allocated(renew)   ) deallocate(renew  )
     if ( allocated(async)   ) deallocate(async  )
     if ( allocated(is_cp)   ) deallocate(is_cp  )
     if ( allocated(nc_cp)   ) deallocate(nc_cp  )

     if ( allocated(saved_p) ) deallocate(saved_p)
     if ( allocated(saved_n) ) deallocate(saved_n)

     return
  end subroutine cat_free_part

!!========================================================================
!!>>> core service subroutines                                         <<<
!!========================================================================

!!>>> cat_make_npart: it is used to determine renew, which parts should
!!>>> be recalculated, is_cp is also reseted in this subroutine
  subroutine cat_make_npart(cmode, csize, index_loc, tau_s, tau_e)
     implicit none

! external arguments
! mode for different Monte Carlo moves
     integer, intent(in)  :: cmode

! total number of operators for current diagram
     integer, intent(in)  :: csize

! local version of index_t
     integer, intent(in)  :: index_loc(mkink)

! imaginary time value of operator A, only valid in cmode = 1 or 2
     real(dp), intent(in) :: tau_s

! imaginary time value of operator B, only valid in cmode = 1 or 2
     real(dp), intent(in) :: tau_e

! local variables
! loop index
     integer  :: i
     integer  :: j

! position of the operator A and operator B, index of part
     integer  :: tis
     integer  :: tie
     integer  :: tip

! length in imaginary time axis for each part
     real(dp) :: interval

! evaluate interval at first
     interval = beta / real(npart)

! init key arrays
     nop = 0
     ops = 0
     ope = 0

! init global arrays (renew and is_cp)
     renew = 0
     is_cp = 0

! calculate number of operators for each part
     do i=1,csize
         j = ceiling( time_v( index_loc(i) ) / interval )
         nop(j) = nop(j) + 1
     enddo ! over i={1,csize} loop

! calculate the start and end index of operators for each part
     do i=1,npart
         if ( nop(i) > 0 ) then
             ops(i) = 1
             do j=1,i-1
                 ops(i) = ops(i) + nop(j)
             enddo ! over j={1,i-1} loop
             ope(i) = ops(i) + nop(i) - 1
         endif ! back if ( nop(i) > 0 ) block
     enddo ! over i={1,npart} loop

! next we have to figure out which parts should be updated
! case 1: only some parts need to be updated
     if ( cmode == 1 .or. cmode == 2 ) then

! get the position of operator A and operator B
         tis = ceiling( tau_s / interval )
         tie = ceiling( tau_e / interval )

! determine the influence of operator A, which part should be recalculated
         renew(tis) = 1
! special attention: if operator A is on the left or right boundary, then
! the neighbour part should be recalculated as well
         if ( nop(tis) > 0 ) then
             if ( tau_s >= time_v( index_loc( ope(tis) ) ) ) then
                 tip = tis + 1
                 do while ( tip <= npart )
                     if ( nop(tip) > 0 ) then
                         renew(tip) = 1; EXIT
                     endif ! back if ( nop(tip) > 0 ) block
                     tip = tip + 1
                 enddo ! over do while loop
             endif ! back if ( tau_s >= time_v( index_t( ope(tis) ) ) ) block
         else
             tip = tis + 1
             do while ( tip <= npart )
                 if ( nop(tip) > 0 ) then
                     renew(tip) = 1; EXIT
                 endif ! back if ( nop(tip) > 0 ) block
                 tip = tip + 1
             enddo ! over do while loop
         endif ! back if ( nop(tis) > 0 ) block

! determine the influence of operator B, which part should be recalculated
         renew(tie) = 1
! special attention: if operator B is on the left or right boundary, then
! the neighbour part should be recalculated as well
         if ( nop(tie) > 0 ) then
             if ( tau_e >= time_v( index_loc( ope(tie) ) ) ) then
                 tip = tie + 1
                 do while ( tip <= npart )
                     if ( nop(tip) > 0 ) then
                         renew(tip) = 1; EXIT
                     endif ! back if ( nop(tip) > 0 ) block
                     tip = tip + 1
                 enddo ! over do while loop
             endif ! back if ( tau_e >= time_v( index_t( ope(tie) ) ) ) block
         else
             tip = tie + 1
             do while ( tip <= npart )
                 if ( nop(tip) > 0 ) then
                     renew(tip) = 1; EXIT
                 endif ! back if ( nop(tip) > 0 ) block
                 tip = tip + 1
             enddo ! over do while loop
         endif ! back if ( nop(tie) > 0 ) block

! case 2: all parts should be updated
     else
         renew = 1
     endif

     return
  end subroutine cat_make_npart

!!>>> cat_make_trace: calculate the contribution to final trace for
!!>>> a given string
  subroutine cat_make_trace(csize, string, index_loc, expt_loc, trace)
     implicit none

! external variables
! number of total fermion operators
     integer, intent(in)   :: csize

! evolution string for this sector
     integer, intent(in)   :: string(csize+1)

! memory address index of fermion operators
     integer, intent(in)   :: index_loc(mkink)

! diagonal elements of last time-evolution matrices
     real(dp), intent(in)  :: expt_loc(ncfgs)

! the calculated trace of this sector
     real(dp), intent(out) :: trace

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k
     integer  :: l

! type for current operator
     integer  :: vt

! flavor channel for current operator
     integer  :: vf

! start index of this sector
     integer  :: indx

! dimension for the sectors
     integer  :: dim1
     integer  :: dim2
     integer  :: dim3
     integer  :: dim4

! index for sectors
     integer  :: isect
     integer  :: sect1
     integer  :: sect2

! the first part with non-zero fermion operators
     integer  :: fpart

! counter for fermion operators
     integer  :: counter

! real(dp) dummy matrices
     real(dp) :: mat_r(max_dim_sect,max_dim_sect)
     real(dp) :: mat_t(max_dim_sect,max_dim_sect)

! initialize dummy arrays
     mat_r = zero
     mat_t = zero

! select the first sector in the string
     isect = string(1)
     dim1  = sectors( string(1) )%ndim

! determine fpart
     fpart = 0
     do i=1,npart
         if ( nop(i) > 0 ) then
             fpart = i; EXIT
         endif ! back if ( nop(i) > 0 ) block
     enddo ! over i={1,npart} loop

! next we perform time evolution from left to right: 0 -> \beta
! loop over all the parts
     do i=1,npart

! empty part, we just skip it
         if ( nop(i) == 0 ) CYCLE

! this part should be recalcuated
         if ( renew(i) == 1 .or. async(i,isect) == 1 ) then
             sect1 = string(ope(i)+1)
             sect2 = string(ops(i))
             dim4 = sectors(sect2)%ndim
             saved_n(:,:,i,isect) = zero

! set its copy status
             is_cp(i,isect) = 1
             nc_cp(i,isect) = dim4

! loop over all the fermion operators in this part
             counter = 0
             do j=ops(i),ope(i)
                 counter = counter + 1
                 indx = sectors( string(j)   )%istart
                 dim2 = sectors( string(j+1) )%ndim
                 dim3 = sectors( string(j)   )%ndim

! multiply the diagonal matrix of time evolution operator
                 if ( counter > 1 ) then
                     do l=1,dim4
                         do k=1,dim3
                             mat_t(k,l) = saved_n(k,l,i,isect) * expt_v(indx+k-1,index_loc(j))
                         enddo ! over k={1,dim3} loop
                     enddo ! over l={1,dim4} loop
                 else
                     mat_t = zero
                     do k=1,dim3
                         mat_t(k,k) = expt_v(indx+k-1,index_loc(j))
                     enddo ! over k={1,dim3} loop
                 endif ! back if ( counter > 1 ) block

! multiply the matrix of fermion operator
                 vt = type_v( index_loc(j) )
                 vf = flvr_v( index_loc(j) )
                 call dgemm( 'N', 'N', dim2, dim4, dim3, &
                                                    one, &
                   sectors( string(j) )%fmat(vf,vt)%val, &
                                            dim2, mat_t, &
                                           max_dim_sect, &
                             zero, saved_n(:,:,i,isect), &
                                           max_dim_sect )
             enddo ! over j={ops(i),ope(i)} loop

! multiply this part with the rest parts
             if ( i > fpart ) then
                 call dgemm( 'N', 'N', dim2, dim1, dim4, &
                              one, saved_n(:,:,i,isect), &
                                           max_dim_sect, &
                                                  mat_r, &
                                           max_dim_sect, &
                                            zero, mat_t, &
                                           max_dim_sect )
                 mat_r(:,1:dim1) = mat_t(:,1:dim1)
             else
                 mat_r(:,1:dim1) = saved_n(:,1:dim1,i,isect)
             endif ! back if ( i > fpart ) block

! this part has been calculated previously, just use its results
         else
             sect1 = string(ope(i)+1)
             sect2 = string(ops(i))
             dim2 = sectors(sect1)%ndim
             dim3 = sectors(sect2)%ndim
             if ( i > fpart ) then
                 call dgemm( 'N', 'N', dim2, dim1, dim3, &
                              one, saved_p(:,:,i,isect), &
                                           max_dim_sect, &
                                                  mat_r, &
                                           max_dim_sect, &
                              zero, mat_t, max_dim_sect )
                 mat_r(:,1:dim1) = mat_t(:,1:dim1)
             else
                 mat_r(:,1:dim1) = saved_p(:,1:dim1,i,isect)
             endif ! back if ( i > fpart ) block

         endif ! back if ( renew(i) == 1 .or. async(i,isect) == 1 )  block

! setup the start sector for next part
         isect = sect1
     enddo ! over i={1,npart} loop

! special treatment of the last time evolution operator
     indx = sectors( string(1) )%istart

! no fermion operators
     if ( csize == 0 ) then
         do k=1,dim1
             mat_r(k,k) = expt_loc(indx+k-1)
         enddo ! over k={1,dim1} loop
! multiply the last time evolution operator
     else
         do l=1,dim1
             do k=1,dim1
                 mat_r(k,l) = mat_r(k,l) * expt_loc(indx+k-1)
             enddo ! over k={1,dim1} loop
         enddo ! over l={1,dim1} loop
     endif ! back if ( csize == 0 ) block

! calculate the trace and store the final product
     trace = zero
     do j=1,sectors( string(1) )%ndim
         trace = trace + mat_r(j,j)
         sectors( string(1) )%prod(j) = mat_r(j,j)
     enddo ! over j={1,sectors( string(1) )%ndim} loop

     return
  end subroutine cat_make_trace

  end module m_part

!!========================================================================
!!>>> service layer: utility subroutines to calculate trace            <<<
!!========================================================================

!!>>> ctqmc_lazy_ztrace: core subroutine of manjushaka
!!>>> (1) use good quantum numbers (GQNs) algorithm, split the total
!!>>>     Hibert space to small subspace, the dimension of F-matrix will
!!>>>     be smaller.
!!>>> (2) use divide and conqure algorithm, split the imaginary time axis
!!>>>     into several parts, save the matrices products of each part,
!!>>>     which may be used by next Monte Carlo move.
!!>>> (3) use lazy trace algorithm to reject some proposed moves immediately.
!!>>> (4) truncate the Hilbert space according to the total occupancy and
!!>>>     the probability of atomic eigenstates.
!!>>> note: you should carefully choose npart in order to obtain the
!!>>> best speedup.
  subroutine ctqmc_lazy_ztrace(cmode, csize, ratio, tau_s, tau_e, r, p, pass)
     use constants, only : dp, zero, one

     use control, only : ncfgs
     use control, only : mkink
     use context, only : c_mtr, n_mtr
     use context, only : index_t, index_v, expt_t, expt_v
     use context, only : diag

     use m_sect, only : nsect
     use m_sect, only : sectors
     use m_sect, only : cat_make_string
     use m_part, only : cat_make_npart
     use m_part, only : cat_make_trace

     implicit none

! external arguments
! different type of Monte Carlo moves
! if cmode = 1, partly-trial calculation, useful for ctqmc_insert_ztrace() etc
! if cmode = 2, partly-normal calculation, not used by now
! if cmode = 3, fully-trial calculation, useful for ctqmc_reflip_kink()
! if cmode = 4, fully-normal calculation, useful for ctqmc_retrieve_status()
     integer,  intent(in)  :: cmode

! total number of operators for current diagram
     integer,  intent(in)  :: csize

! the calculated determinant ratio and prefactor
     real(dp), intent(in)  :: ratio

! imaginary time value of operator A, only needed in cmode = 1 or 2
     real(dp), intent(in)  :: tau_s

! imaginary time value of operator B, only needed in cmode = 1 or 2
     real(dp), intent(in)  :: tau_e

! random number
     real(dp), intent(in)  :: r

! the final transition probability
     real(dp), intent(out) :: p

! whether accept this move
     logical, intent(out)  :: pass

! local variables
! loop index
     integer  :: i
     integer  :: j

! start index of a sector
     integer  :: indx

! number of alive sectors
     integer  :: nlive

! maximum and minimum bounds of acceptance ratio
     real(dp) :: pmax
     real(dp) :: pmin

! sum of btrace (trace boundary)
     real(dp) :: sbound

! sum of absolute value of trace
     real(dp) :: cumsum

! index of the living sector
     integer  :: living(nsect)

! minimum dimension of the sectors in valid strings
     integer  :: mindim(nsect)

! sector index of a string
     integer  :: string(csize+1,nsect)

! local version of index_t
     integer  :: index_loc(mkink)

! local version of expt_t
     real(dp) :: expt_loc(ncfgs)

! trace boundary for sectors
     real(dp) :: btrace(nsect)

! trace for each sector
     real(dp) :: strace(nsect)

! copy data from index_t or index_v to index_loc
! copy data from expt_t to expt_loc
     select case (cmode)

         case (1)
             index_loc = index_t
             expt_loc = expt_t(:,1)

         case (2)
             index_loc = index_v
             expt_loc = expt_t(:,2)

         case (3)
             index_loc = index_t
             expt_loc = expt_t(:,2)

         case (4)
             index_loc = index_v
             expt_loc = expt_t(:,2)

     end select

! build all possible strings for all the sectors. if one string may be
! invalid, then all of its elements must be -1
     call cat_make_string(csize, index_loc, string)

! determine which part should be recalculated (global variables renew
! and is_cp will be updated in this subroutine)
     call cat_make_npart(cmode, csize, index_loc, tau_s, tau_e)

! we can verify string here to see whether this diagram can survive?
! if not, return immediately.
     pass = .true.
     if ( all( string == -1 ) ) then
         pass = .false.; p = zero; RETURN
     endif ! back if ( all( string == -1 ) ) block

! calculate the trace bounds for each sector and determine the
! number of sectors which actually contribute to the total trace
     nlive  = 0
     living = -1
     mindim = 0
     btrace = zero
     do i=1,nsect
! if the string is invalid, we just skip it
         if ( string(1,i) == -1 ) then
             CYCLE
! find valid string which may contribute to the total trace
         else
! increase the counter
             nlive = nlive + 1

! record its index
             living(nlive) = i

! calculate its trace bound and determine the minimal dimension for
! each alive string
             mindim(i) = sectors(i)%ndim
             btrace(nlive) = one
             do j=1,csize
                 if ( mindim(i) > sectors( string(j,i) )%ndim ) then
                     mindim(i) = sectors( string(j,i) )%ndim
                 endif ! back if ( mindim(i) > sectors( string(j,i) )%ndim ) block
                 indx = sectors(string(j,i))%istart
                 btrace(nlive) = btrace(nlive) * expt_v(indx, index_loc(j))
             enddo ! over j={1,csize} loop

! special treatment for the last time evolution operator
             indx = sectors(string(1,i))%istart
             btrace(nlive) = btrace(nlive) * expt_loc(indx) * mindim(i)
         endif ! back if ( string(1,i) == -1 ) block
     enddo ! over i={1,nsect} loop

! calculate the summmation of trace bounds and the maximum bound of the
! acceptance ratio, and then we check whether pmax < r. if it is true,
! reject this move immediately
     sbound = sum( btrace(1:nlive) )
     pmax = abs(ratio) * abs(sbound / c_mtr)
     if ( pmax < r ) then
         pass = .false.; p = zero; RETURN
     endif ! back if ( pmax < r ) block

! sort the btrace to speed up the refining process. here, we use simple
! bubble sort algorithm, because nalive_sect is usually small
     call s_sorter2( nlive, btrace(1:nlive), living(1:nlive) )

! begin to refine the trace bounds
     pass = .false.
     cumsum = zero
     strace = zero
     do i=1,nlive
! calculate the trace for one sector, this call will consume a lot of
! time if the dimension of fmat and expansion order is large, so we
! should carefully optimize it.
         call cat_make_trace(csize, string(:,living(i)), index_loc, expt_loc, strace(i))
! if this move is not accepted, refine the trace bound to see whether
! we can reject it before calculating the trace of all of the sectors
         if ( .not. pass ) then
             cumsum = cumsum + abs( strace(i) )
             sbound = sbound - btrace(i)
! calculate pmax and pmin
             pmax = abs(ratio) * abs( (cumsum + sbound) / c_mtr )
             pmin = abs(ratio) * abs( (cumsum - sbound) / c_mtr )
! check whether pmax < r
             if ( pmax < r ) then
                 pass = .false.; p = zero; RETURN
             endif ! back if ( pmax < r ) block
! this move is accepted, stop refining process, calculate the trace of
! remaining sectors to get the final result of trace.
             if ( pmin > r ) then
                 pass = .true.
             endif ! back if ( pmin > r ) block
         endif ! back if ( .not. pass ) block
     enddo ! over i={1,nlive} loop

! if we arrive here, two cases
! case 1: pass == .false., we haven't determined the pass
! case 2: pass == .true. we have determined the pass
! anyway, we have to calculate the final transition probability (p), and
! update n_mtr and pass.
     n_mtr = sum(strace(1:nlive))
     p = ratio * (n_mtr / c_mtr)
     pass = ( min(one, abs(p)) > r )

! store the diagonal elements of final product in diag(:,1)
     diag(:,1) = zero
     do i=1,nlive
         indx = sectors( living(i) )%istart
         do j=1,sectors( living(i) )%ndim
             diag(indx+j-1,1) = sectors( living(i) )%prod(j)
         enddo ! over j={1,sectors( living(i) )%ndim} loop
     enddo ! over i={1,nlive} loop

     return
  end subroutine ctqmc_lazy_ztrace

!!>>> ctqmc_retrieve_ztrace: calculate the trace for retrieve status
  subroutine ctqmc_retrieve_ztrace(csize, trace)
     use constants, only : dp, zero

     use control, only : ncfgs
     use control, only : mkink
     use context, only : index_v, expt_t
     use context, only : diag

     use m_sect, only : nsect
     use m_sect, only : sectors
     use m_sect, only : cat_make_string
     use m_part, only : cat_make_npart
     use m_part, only : cat_make_trace

     implicit none

! external arguments
! total number of operators for current diagram
     integer,  intent(in)  :: csize

! the calculated trace
     real(dp), intent(out) :: trace

! local variables
! loop index
     integer  :: i
     integer  :: j

! start index of a sector
     integer  :: indx

! sector index of a string
     integer  :: string(csize+1,nsect)

! local version of index_v
     integer  :: index_loc(mkink)

! local version of expt_t
     real(dp) :: expt_loc(ncfgs)

! trace for each sector
     real(dp) :: strace(nsect)

! copy data from index_v to index_loc
! copy data from expt_t to expt_loc
     index_loc = index_v
     expt_loc = expt_t(:,2)

! build all possible strings for all the sectors. if one string may be
! invalid, then all of its elements must be -1
     call cat_make_string(csize, index_loc, string)

! determine which part should be recalculated (global variables renew
! and is_cp will be updated in this subroutine)
     call cat_make_npart(4, csize, index_loc, zero, zero)

! calculate the trace of each sector one by one
     strace = zero
     do i=1,nsect
! invalid string, its contribution is neglected
! note: here we only check the first element of this string. it is enough
         if ( string(1,i) == -1 ) then
             strace(i) = zero
             sectors(i)%prod = zero
! valid string, we have to calculate its contribution to trace
         else
             call cat_make_trace(csize, string(:,i), index_loc, expt_loc, strace(i))
         endif ! back if ( string(1,i) == -1 ) block
     enddo ! over i={1,nsect} loop
     trace = sum(strace)

! store the diagonal elements of final product in diag(:,1), which can be
! used to calculate the atomic probability
     do i=1,nsect
         indx = sectors(i)%istart
         do j=1,sectors(i)%ndim
             diag(indx+j-1,1) = sectors(i)%prod(j)
         enddo ! over j={1,sectors(i)%ndim} loop
     enddo ! over i={1,nsect} loop

     return
  end subroutine ctqmc_retrieve_ztrace

!!>>> ctqmc_make_evolve: used to update the operator traces of the
!!>>> modified part
  subroutine ctqmc_make_evolve()
     use control, only : npart
     use context, only : c_mtr, n_mtr
     use context, only : diag

     use m_sect, only : nsect
     use m_part, only : renew, async, is_cp, nc_cp, saved_p, saved_n

     implicit none

! local variables
! loop index
     integer :: i
     integer :: j

! update the operator traces
     c_mtr = n_mtr

! update diag for the calculation of atomic state probability
     diag(:,2) = diag(:,1)

! even if renew(j) is 1, not all of the sectors in this part (the j-th
! part) will be renewed. there are many reasons. one of them is the
! broken string. anyway, at this time, we have to remind the solver that
! the matrix products for these sectors in j-th part is unsafe. so it is
! necessary to update async here
     do i=1,nsect
         do j=1,npart
             if ( renew(j) == 1 .and. is_cp(j,i) == 0 ) then
                 async(j,i) = 1
             endif ! back if ( renew(j) == 1 .and. is_cp(j,i) == 0 ) block
         enddo ! over j={1,npart} loop
     enddo ! over i={1,nsect} loop

! if we used the divide-and-conquer algorithm, then we had to save the
! change matrices products when proposed moves were accepted. and sine
! the matrices products are updated, we also update the corresponding
! async variable to tell the impurity solver that these saved_p is OK
     do i=1,nsect
         do j=1,npart
             if ( is_cp(j,i) == 1 ) then
                 saved_p(:,1:nc_cp(j,i),j,i) = saved_n(:,1:nc_cp(j,i),j,i)
                 async(j,i) = 0
             endif ! back if ( is_cp(j,i) == 1 ) block
         enddo ! over j={1,npart} loop
     enddo ! over i={1,nsect} loop

     return
  end subroutine ctqmc_make_evolve
