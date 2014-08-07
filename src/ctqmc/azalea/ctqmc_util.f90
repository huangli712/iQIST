!-------------------------------------------------------------------------
! project : azalea
! program : ctqmc_make_uumat
!           ctqmc_make_state
! source  : ctqmc_util.f90
! type    : functions & subroutines
! author  : li huang (email:huangli712@gmail.com)
! history : 10/01/2008 by li huang
!           02/08/2009 by li huang
!           09/23/2009 by li huang
!           09/26/2009 by li huang
!           11/17/2009 by li huang
!           11/21/2009 by li huang
!           12/18/2009 by li huang
!           12/22/2009 by li huang
!           12/29/2009 by li huang
!           01/12/2010 by li huang
!           02/27/2010 by li huang
!           06/08/2010 by li huang
!           06/22/2010 by li huang
! purpose : to provide utility functions and subroutines for hybridization
!           expansion version continuous time quantum Monte Carlo (CTQMC)
!           quantum impurity solver
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> to build general U interaction matrix: uumat, using my own style
! note: do not support spin-flip and pair-hopping term so far
! note: only Uc and Jz are need, the other Coulomb interaction parameters
! are used as backup
  subroutine ctqmc_make_uumat(uumat)
     use constants
     use control

     implicit none

! external arguments
! Coulomb interaction matrix
     real(dp), intent(out) :: uumat(norbs, norbs)

! local variables
! loop index
     integer  :: i
     integer  :: j

     integer  :: k
     integer  :: m

! dummy u vector
     real(dp) :: ut(nband*(norbs-1))

! initialize it
     uumat = zero

! calculate it
     k = 0
     do i=1,norbs-1
         do j=i+1,norbs
             k = k + 1
             if ( i <= nband .and. j > nband ) then
                 m = j - nband
                 if ( m == i ) then
                     ut(k) = Uc
                 else
                     ut(k) = Uc - 2.0_dp * Jz
                 endif
             else
                 ut(k) = Uc - 3.0_dp * Jz
             endif

             uumat(i,j) = ut(k)
             uumat(j,i) = ut(k)
         enddo ! over j={i+1,norbs} loop
     enddo ! over i={1,norbs-1} loop

     return
  end subroutine ctqmc_make_uumat

!>>> convert current atomic state array into a decimal number (state index)
  subroutine ctqmc_make_state(norbs, pstat, state)
     implicit none

! external arguments
! index of atomic state
     integer, intent(out) :: pstat

! number of orbitals
     integer, intent(in)  :: norbs

! atomic state array
     integer, intent(in)  :: state(norbs)

! local variables
! loop index
     integer :: i

! init pstat
     pstat = 1

! evaluate pstat, for example, 0101 = 0*2^0 + 1*2^1 + 0*2^2 + 1*2^3 = 10
     do i=1,norbs
         if ( state(i) > 0 ) pstat = pstat + ishft(1, i-1)
     enddo ! over i={1,norbs} loop

     return
  end subroutine ctqmc_make_state
