!-------------------------------------------------------------------------!
! project : clematis
! program : atomic_config
! history : Apr 27, 2011
! author  : duliang (duleung@gmail.com)
! purpose : setup atomic Hamiltonian parameters
! comment : 
!-------------------------------------------------------------------------!
!>>> atomic hamiltonian parameters
  subroutine atomic_config()
     use constants
     use control

     implicit none

! local variables
! check whether the input file (dft.atom.in) exist
     logical :: exists

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

     nband = 3            ! number of bands
     nspin = 2            ! number of spins
     norbs = 6            ! number of orbits
     ntots = 3            ! number of total electrons
     ncfgs = 20           ! number of many-body configurations

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

     Uc     = 4.00_dp      ! intraorbital Coulomb interaction
     Uv     = 2.00_dp      ! interorbital Coulomb interaction
     Jz     = 1.00_dp      ! Hund's exchange interaction
     Js     = 1.00_dp      ! spin-flip interaction
     Jp     = 1.00_dp      ! pair-hopping interaction
     lamb   = 0.00_dp      ! spin-orbit coupling parameter
     nmin   = 0            ! the minimal occupancy number
     nmax   = 6            ! the maximal occupancy number
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

! read from input file if it exists
     exists = .false.

! inquire the input file status
     inquire( file="dft.atom.in", exist=exists )

! read parameters from dft.atom.in
     if ( exists .eqv. .true. ) then
         open( mytmp, file="dft.atom.in", status="old" )

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
         read(mytmp, *)
         read(mytmp, *)
         read(mytmp, *)

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
         read(mytmp, *) nband
         read(mytmp, *) nspin
         read(mytmp, *) norbs
         read(mytmp, *) ntots
         read(mytmp, *) ncfgs

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
         read(mytmp, *)

         read(mytmp, *) Uc
         read(mytmp, *) Uv
         read(mytmp, *) Jz
         read(mytmp, *) Js
         read(mytmp, *) Jp
         read(mytmp, *) lamb

         read(mytmp, *) nmin
         read(mytmp, *) nmax
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

         close(mytmp)
     endif ! back if ( exists .eqv. .true. ) block

     return
  end subroutine atomic_config


