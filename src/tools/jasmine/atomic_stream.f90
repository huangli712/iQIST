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

!>>> allocate memory and initialize them
  subroutine atomic_setup_array()
     use constants
     use control
     use mod_global

     implicit none

! allocate memory for basis sheet
     call atomic_allocate_memory_basis()

! allocate memory for Hamiltonian matrix
     call atomic_allocate_memory_hmtrx()

     return
  end subroutine atomic_setup_array

!>>> deallocate memory and finalize them
  subroutine atomic_final_array()
     use constants
     use control
     use mod_global

     implicit none

! deallocate memory for basis sheet
     call atomic_deallocate_memory_basis()

! deallocate memory for hamiltonian matrix
     call atomic_deallocate_memory_hmtrx()

     return
  end subroutine atomic_final_array

!>>> calculate combination algebra 
  function state_pick(ntiny, nlarg) result(value)
     implicit none

! external variables
     integer, intent(in) :: ntiny
     integer, intent(in) :: nlarg

! local variables
     integer :: i

! auxiliary integer variable
     integer :: nlow

! numberator of the combination algebra
     real(8) :: numer

! denominator of the combination algebra
     real(8) :: denom

! result value of the combination algebra
     integer :: value

! transform the combination algebra
     nlow = min(ntiny, nlarg-ntiny)

! numerator in combination algebra
     numer = 1.0D0
     do i=nlarg-nlow+1,nlarg
        numer = numer * dble(i)
     enddo ! over i={nlarg-nlow+1,nlarg} loop

! denominator in combination algebra
     denom = 1.0D0
     do i=1,nlow
        denom = denom * dble(i)
     enddo ! over i={1,nlow} loop

! result value
     value = nint(numer / denom)

     return
  end function state_pick

!-------------------------------------------------------------------------!
!>>> build key matrix in atomic problem
!-------------------------------------------------------------------------!
  subroutine atomic_jovian(nstat, cemat, somat, cumat)
     use constants
     use control

     implicit none

! external arguments
! external functions
     integer, external :: state_pick

! number of configuration in each subspace
     integer, intent(out) :: nstat(0:norbs)

! impurity energy level
     complex(dp), intent(out) :: cemat(norbs, norbs)

! spin-orbit coupling matrix in orginal single particle basis
     complex(dp), intent(out) :: somat(norbs, norbs)

! general interaction coefficents in orginal single particle basis
     complex(dp), intent(out) :: cumat(norbs, norbs, norbs, norbs)

! local variables
! loop index over good quantum number
     integer :: ibit

! loop index over orbits
     integer :: iorb
     integer :: jorb

! initialize dimension of each subspace
     write(mystd, "(2X, A)") "-----------------------------------------------"
     write(mystd, "(2X, A)") ">>> number of configuation in each subspace <<<"
     write(mystd, "(2X, A)") "-----------------------------------------------"
     do ibit=0,norbs
         nstat(ibit) = state_pick(ibit, norbs)
         write(mystd, "(2I8)") ibit, nstat(ibit)
     enddo ! over ibit={0,norbs} loop 

     somat = czero
! build impurity energy level matrix
     call atomic_make_cemat(norbs, cemat)

! build spin-orbit coupling matrix in orginal single particle basis
!     call atomic_make_somat(norbs, lamb, somat)

! build general interaction coefficents matrix in orginal basis
     call atomic_make_cumat(norbs, Uc, Uv, Jz, Js, Jp, cumat)

     return
  end subroutine atomic_jovian

!>>> impurity energy level
  subroutine atomic_make_cemat(norbs, cemat)
     use constants

     implicit none

! external arguments
! number of orbits
     integer, intent(in) :: norbs

! impurity energy level
     complex(dp), intent(out) :: cemat(norbs, norbs)

! local variable
! check whether the input file (dft.cemat.in) exist
     logical :: exists

! status while reading data
     integer :: ierr

! orbital index
     integer :: iorb
     integer :: jorb

! auxiliary real(dp) variables
     real(dp) :: dtmpa
     real(dp) :: dtmpb

! initialize eimp to be czero
     cemat = dcmplx(0.0D0, 0.0D0)

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
! read from  input file if it exists
     exists = .false.

! inquire the input file status
     inquire( file="dft.cemat.in", exist=exists )

! read from rtgw.eimp.in
     if ( exists .eqv. .true. ) then
         open( mytmp, file="dft.cemat.in", status="old" )

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
         do while ( .true. )
             read(mytmp, *, iostat=ierr) iorb, jorb, dtmpa, dtmpb
             if (ierr /= 0) exit; cemat(iorb, jorb) = dcmplx(dtmpa, dtmpb)
         enddo ! over while ( .true. ) loop

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

         close(mytmp)
     endif ! back if ( exists .eqv. .true. ) block

# if defined (duliang)
     open(mytmp, file='solver.cemat.out')
     call zmat_dump2(mytmp, norbs, norbs, cemat)
     close(mytmp)
# endif /* duliang */
     return
  end subroutine atomic_make_cemat

!>>>  $U_{\alpha\betta\delta\gamma}$ coefficent matrix for Uij  <<<!
!--------------------three band for example------------------------!
!> $f_{\alpha}^{\dagger}f_{\betta}^{\dagger}f_{\delta}f_{\gamma}$ <!
! norbs   bandindex(?band)    spinindex(?spin)    representation   !
!   1          1               1->\uparrow             1\up        !
!   2          1               0->\doarrow             1\do        !
!   3          2               1->\uparrow             2\up        !
!   4          2               0->\doarrow             2\do        !
!   5          3               1->\uparrow             3\up        !
!   6          3               0->\doarrow             3\do        !
!------------------------------------------------------------------!
  subroutine atomic_make_cumat( norbs, Uc, Uv, Jz, Js, Jp, cumat )
     use constants

     implicit none

! external arguments
! number of orbits
     integer, intent(in) :: norbs

! general coulomb interaction parameters
     real(dp), intent(in) :: Uc
     real(dp), intent(in) :: Uv
     real(dp), intent(in) :: Jz
     real(dp), intent(in) :: Js
     real(dp), intent(in) :: Jp

! general coulomb interaction matrix
     complex(dp), intent(out) :: cumat(norbs, norbs, norbs, norbs)

! local varibales
! loop index over orbits
     integer :: alpha, betta
     integer :: delta, gamma

! band index and spin index
! band index of alpha and betta
     integer :: aband, bband 

! band index of delta and gamma
     integer :: dband, gband 

! spin index of alpha and betta
     integer :: aspin, bspin 

! spin index of delta and gamma
     integer :: dspin, gspin 

! auxiliary variables
     real(dp) :: dtmp

! initialize cumat to zero
     cumat = czero

! loop for creation operator
     alphaloop: do alpha=1,norbs-1
     bettaloop: do betta=alpha+1,norbs

! loop for destroy operator
        gammaloop: do gamma=1,norbs-1
        deltaloop: do delta=gamma+1,norbs
            aband = (alpha+1)/2; aspin = mod(alpha,2)
            bband = (betta+1)/2; bspin = mod(betta,2)
            gband = (gamma+1)/2; gspin = mod(gamma,2)
            dband = (delta+1)/2; dspin = mod(delta,2)

! here we use "res" due to overlap between "Uv and Jz"
            dtmp = 0.0D0

! intraorbital Coulomb interaction
            if ((alpha.eq.gamma) .and. (betta.eq.delta)) then
                if ((aband.eq.bband) .and. (aspin.ne.bspin)) then
                    dtmp = dtmp + Uc
!#                  write(mystd, "(A3,4I4)") " Uc", alpha, betta, delta, gamma
                endif
            endif

! interorbital Coulomb interaction
            if ((alpha.eq.gamma) .and. (betta.eq.delta)) then
                if (aband .ne. bband) then
                    dtmp = dtmp + Uv
!#                  write(mystd, "(A3,4I4)") " Uv", alpha, betta, delta, gamma
                endif
            endif

! Hund's exchange interaction 
            if ((alpha.eq.gamma) .and. (betta.eq.delta)) then
                if ((aband.ne.bband) .and. (aspin.eq.bspin)) then
                    dtmp = dtmp - Jz
!#                  write(mystd, "(A3,4I4)") " Jz", alpha, betta, delta, gamma
                endif
            endif
           
! spin flip term
            if ((aband.eq.gband) .and. (bband.eq.dband)) then
                if ((aspin.ne.gspin) .and. (bspin.ne.dspin) .and. (aspin.ne.bspin)) then
                    dtmp = dtmp - Js
!#                  write(mystd, "(A3,4I4)") " Js", alpha, betta, delta, gamma
                endif
            endif
          
! pair hopping term
            if ((aband.eq.bband) .and. (dband.eq.gband) .and. (aband.ne.dband)) then
                if ((aspin.ne.bspin) .and. (dspin.ne.gspin) .and. (aspin.eq.gspin)) then
                    dtmp = dtmp + Jp
!#                  write(mystd, "(A3,4I4)") " Jp", alpha, betta, delta, gamma
                endif
            endif
                 
            cumat(alpha, betta, delta, gamma) = dtmp

        enddo deltaloop ! over delta={gamma+1,norbs} loop
        enddo gammaloop ! over gamma={1,norbs-1} loop
     enddo bettaloop ! over betta={alpha+1,norbs} loop
     enddo alphaloop ! over alpha={1,norbs-1} loop

# if defined (duliang)
     open(mytmp, file='solver.cumat.out')
     call zmat_dump4(mytmp, norbs, norbs, norbs, norbs, cumat)
     close(mytmp)
# endif /* duliang */
     return
  end subroutine atomic_make_cumat

!-------------------------------------------------------------------------!
!>>> setup spin-orbit coupling matrix in orginal single particle basis <<<!
!>>> norbs=02 => l=0; norbs=06 => l=1; norbs=10 => l=2; norbs=14 => l=3 <<!
!-------------------------------------------------------------------------!
  subroutine atomic_make_somat(norbs, lamb, somat)
     use constants

     implicit none

! external arguments
! number of orbits
     integer, intent(in) :: norbs

! spin-orbital coupling parameter
     real(8), intent(in) :: lamb

! spin-orbit coupling matrix in {lz, sz} single particle basis
     complex(8), intent(out) :: somat(norbs, norbs)

! local variables
! loop index over orbits
     integer :: iorb
     integer :: jorb

! status while allocating memory
     integer :: istat

! status while reading from file
     integer :: ierr

! auxiliary logical variables
     logical :: exists

! auxiliary real(8) variables
     real(8) :: dtmpa, dtmpb

! initialize somat to be zero
     somat = dcmplx(0.0D0, 0.0D0)

! spin orbit coupling in real orbit {px,up; px,dn; py,up; py,dn; pz,up; pz,dn}
! for p-orbitals. for d-orbitals is {yz,up; yz,dn; zx,up; zx,dn; xy,up; xy,dn}
     somat(1, 3) = + dcmplx(0.0D0, 1.0D0); somat(1, 6) = - dcmplx(1.0D0, 0.0D0)
     somat(3, 1) = - dcmplx(0.0D0, 1.0D0); somat(3, 6) = + dcmplx(0.0D0, 1.0D0)
     somat(5, 2) = + dcmplx(1.0D0, 0.0D0); somat(5, 4) = - dcmplx(0.0D0, 1.0D0)
     somat(2, 5) = + dcmplx(1.0D0, 0.0D0); somat(2, 4) = - dcmplx(0.0D0, 1.0D0)
     somat(4, 5) = + dcmplx(0.0D0, 1.0D0); somat(4, 2) = + dcmplx(0.0D0, 1.0D0)
     somat(6, 1) = - dcmplx(1.0D0, 0.0D0); somat(6, 3) = - dcmplx(0.0D0, 1.0D0)

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
! read from  input file if it exists
     exists = .false.

! inquire the input file status
     inquire( file="dft.somat.in", exist=exists )
     if (exists .eqv. .true.) somat = dcmplx(0.0D0, 0.0D0)

! read from rtgw.eimp.in
     if ( exists .eqv. .true. ) then
         open( mytmp, file="dft.somat.in", status="old" )

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

         do while ( .true. )
             read(mytmp, *, iostat=ierr) iorb, jorb, dtmpa, dtmpb
             if (ierr /= 0) exit; somat(iorb, jorb) = dcmplx(dtmpa, dtmpb)
         enddo ! over while ( .true. ) loop

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!

         close(mytmp)
     endif ! back if ( exists .eqv. .true. ) block

     somat = somat * lamb / 2.0D0

# if defined (duliang)
     open(mytmp, file='solver.somat.out')
     call zmat_dump2(mytmp, norbs, norbs, somat)
     close(mytmp)
# endif /* duliang */
     return
  end subroutine atomic_make_somat

  subroutine atomic_tran_cumat(norbs, amtrx, cumat, cumat_t)
     use constants

     implicit none

! number of orbits
     integer, intent(in) :: norbs

! transformation matrix from orginal basis to natural basis
     complex(dp), intent(in) :: amtrx(norbs, norbs)

! coefficents matrix for generalized interaction U in orginal basis
     complex(dp), intent(in) :: cumat(norbs, norbs, norbs, norbs)
 
! coefficents matrix for generalized interaction U in natural basis
     complex(dp), intent(out) :: cumat_t(norbs, norbs, norbs, norbs)

! local varoables
! loop index over orbits in orginal single particle basis
     integer :: alpha1, alpha2
     integer :: alpha3, alpha4

! loop index over orbits in natural single particle basis
     integer :: sigma1, sigma2
     integer :: sigma3, sigma4

! auxiliary complex(dp) variables
     complex(dp) :: ctmp

! initialize cumat_t to be zero
     cumat_t = dcmplx(0.0D0, 0.0D0)

     sigma1loop: do sigma1=1,norbs
     sigma2loop: do sigma2=1,norbs
     sigma3loop: do sigma3=1,norbs
     sigma4loop: do sigma4=1,norbs
         ctmp = dcmplx(0.0D0, 0.0D0)

         alpha1loop: do alpha1=1,norbs
         alpha2loop: do alpha2=1,norbs
         alpha3loop: do alpha3=1,norbs
         alpha4loop: do alpha4=1,norbs
             if (abs(cumat(alpha1, alpha2, alpha3, alpha4)) .lt. eps6) cycle
             ctmp = ctmp + cumat(alpha1, alpha2, alpha3, alpha4)          &
                  * conjg(amtrx(alpha1, sigma1)) * amtrx(alpha3, sigma3)  &
                  * conjg(amtrx(alpha2, sigma2)) * amtrx(alpha4, sigma4)
         enddo alpha4loop ! over alpha4={1,norbs} loop
         enddo alpha3loop ! over alpha3={1,norbs} loop
         enddo alpha2loop ! over alpha2={1,norbs} loop
         enddo alpha1loop ! over alpha1={1,norbs} loop

         cumat_t(sigma1, sigma2, sigma3, sigma4) = ctmp
     enddo sigma4loop ! over sigma4={1,norbs} loop
     enddo sigma3loop ! over sigma3={1,norbs} loop
     enddo sigma2loop ! over sigma2={1,norbs} loop
     enddo sigma1loop ! over sigma1={1,norbs} loop

# if defined (duliang)
     open(mytmp, file='solver.tumat.out')
     call zmat_dump4(mytmp, norbs, norbs, norbs, norbs, cumat_t)
     close(mytmp)
# endif /* duliang */
     return
  end subroutine atomic_tran_cumat

  subroutine atomic_make_amtrx(norbs, amtrx)
     use constants

     implicit none

! external arguments
! number of orbits
     integer, intent(in) :: norbs

! transformation matrix from orginal to natural single particle basis
     complex(dp), intent(out) :: amtrx(norbs, norbs)


! local variable
! check whether the input file (dft.cemat.in) exist
     logical :: exists

! status while reading data
     integer :: ierr

! orbital index
     integer :: iorb
     integer :: jorb

! auxiliary real(dp) variables
     real(dp) :: dtmpa
     real(dp) :: dtmpb


! initialize some variables
     amtrx = dcmplx(0.0D0, 0.0D0)

     open( mytmp, file="dft.amat.in", status="old" )

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
         do while ( .true. )
             read(mytmp, *, iostat=ierr) iorb, jorb, dtmpa, dtmpb
             if (ierr /= 0) exit; amtrx(iorb, jorb) = dcmplx(dtmpa, dtmpb)
         enddo ! over while ( .true. ) loop

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
     close(mytmp)

     return
  end subroutine atomic_make_amtrx
