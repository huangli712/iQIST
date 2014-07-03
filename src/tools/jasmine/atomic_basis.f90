!========================================================================!
! project : clematis
! program : atomic_build_basis
! history : Apr 29, 2011
! authors : duliang (emial:duleung@gmail.com)
! purpose : build basis in "ntots" subspace or global space
! comment : 
!=========================================================================!
  subroutine atomic_make_basis(norbs, ncfgs, ntots, nstat, basis, invcd, invsn)
     use constants

     implicit none

! external variables
! number of orbits
     integer, intent(in) :: norbs

! number of Fock basis
     integer, intent(in) :: ncfgs

! number of total electrons
     integer, intent(in) :: ntots

! number of configuration
     integer, intent(in) :: nstat(0:norbs)

! decimal representation of Fock basis
     integer, intent(out) :: basis(ncfgs)

! serial number of decimal represented Fock basis
     integer, intent(out) :: invsn(0:2**norbs-1)

! binary representation of Fock basis
     integer, intent(out) :: invcd(norbs, ncfgs)

! local variables
! loop index over decimal number
     integer :: inum

! loop index over fock basis
     integer :: ibas

! loop index over orbits
     integer :: ipos

! basis counter in global space
     integer :: icnt
     integer :: jcnt

! number of total electrons
     integer :: itot

! number of total electrons denoted by a decical number
     integer :: nbit

! auxiliary integer variables
     integer :: casea, caseb
     integer :: ntiny, nlarg

! initialize some variables
     icnt = 0; basis = 0
     jcnt = 0; invsn = 0; invcd = 0

! casea => construct basis sheet in only itot subspace (default)
! caseb => construct basis sheet in global state space (special)
     casea = nstat(ntots); caseb = sum(nstat(0:norbs))

     ! case I : many body basis of a subspace
     if (ncfgs .eq. casea) then
         ntiny = ntots; nlarg = ntots
     endif ! back if (ncfgs .eq. casea) block

     ! case II: many body basis of global space
     if (ncfgs .eq. caseb) then
         ntiny = 0    ; nlarg = norbs
     endif ! back if (ncfgs .eq. caseb) block

     ! case III: severe error case
     if ((ncfgs.ne.casea) .and. (ncfgs.ne.caseb)) then
         stop "number of configurations ncfgs is wrong, modify dft.atom.in first"
     endif ! back if ((ncfgs.ne.casea) .and. (ncfgs.eq.caseb)) block

! construct basis(decimal) and "serial number(sn)" in Fock space
     itotloop: do itot=ntiny,nlarg

         do inum=0,2**norbs-1
             call verif( inum, norbs, nbit )
             if ( nbit .eq. itot ) then
                 icnt = icnt + 1; jcnt = jcnt + 1
                 basis(icnt) = inum; invsn(inum) = icnt
             endif ! back if ( nbit .eq. itot ) block
         enddo ! over inum={0,2**norbs-1} loop

!# very important, check number of configurations in itot subspace
         if (jcnt .ne. nstat(itot)) then
             stop "error happened in subroutine atomic_build_basis"
         endif ! back if ( icnt .ne. nstat(itot) ) block
         jcnt = 0 ! reinitialize the counter for next subspace
        
     enddo itotloop ! over itot={ntiny,nlarg} loop

! construct inverse binary code from a decimal number 
     do ibas=1,ncfgs
         do ipos=1,norbs
             if( btest(basis(ibas), ipos-1) ) invcd(ipos, ibas) = 1
         enddo ! over ipos={1,norbs} loop
     enddo ! over ibas={1,ncfgs} loop

! dump atomic configurations to file "atom.basis.in"
     call atomic_dump_basis(norbs, ncfgs, basis, invsn, invcd)

     return
  end subroutine atomic_make_basis

!>>> electron number denoted by a decimal number
  subroutine verif( inum, norbs, nbits )
     implicit none

! external arguments
! decimal number
     integer, intent(in)  :: inum

! number of orbits
     integer, intent(in)  :: norbs

! number of total electrons
     integer, intent(out) :: nbits

! local variables
! local index over orbits
     integer :: ipos

! number of total electrons
     nbits = 0
     do ipos=1,norbs
         if( btest(inum, ipos-1) ) nbits = nbits + 1
     enddo ! over ipos={1,norbs} loop

     return
  end subroutine verif
