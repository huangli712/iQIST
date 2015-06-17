!!!=========+=========+=========+=========+=========+=========+=========+!
!!! HIBISCUS/toolbox/makekra @ iQIST                                     !
!!!                                                                      !
!!! This tool is used to perform kramers-kronig transformation for the   !
!!! imaginary part of matsubara green's function                         !
!!! author  : Li Huang (at IOP/CAS & SPCLab/CAEP & UNIFR)                !
!!! version : v2015.01.06T                                               !
!!! status  : WARNING: IN TESTING STAGE, USE IT IN YOUR RISK             !
!!! comment : any question, please contact with lihuang.dmft@gmail.com   !
!!!=========+=========+=========+=========+=========+=========+=========+!

!!
!!
!! Introduction
!! ============
!!
!! The makekra code reads the density of states data, and then calculate
!! the imaginary part of matsubara green's function. And then using the
!! kramers-kronig transformation, we can calculate the real part of the
!! matsubara green's function easily. So the complete green's function
!! is obtained, which can be used to calculate the self-energy function
!! on real axis by the invert hilbert transformation.
!!
!! Now this code is interfaced with hibiscus/entropy code merely. it
!! can read the mem.dos.dat file as input data. While to interface it
!! with the hibiscus/stoch code is very simple. What you need to do is
!! to rename sac.imsum.dat to mem.dos.dat file, and then supplement
!! the data for different orbitals.
!!
!! Usage
!! =====
!!
!! # ./mkra or bin/mkra.x
!!
!! Input
!! =====
!!
!! mem.dos.dat (necessary)
!!
!! Output
!! ======
!!
!! kra.grn.dat
!!
!! Documents
!! =========
!!
!! For more details, please go to iqist/doc/manual directory.
!!
!!

  program makekra
     use constants, only : dp, zero, one, half, pi, mystd, mytmp

     implicit none

! local control parameters
! number of frequency grid points on half plane (total number = 2*nw + 1)
     integer :: nw = 400

! number of orbitals, include the spin degree of freedom
     integer :: nq = 2

! local variables
! loop index for frequency grid
     integer :: iw

! loop index for orbitals
     integer :: iq

! status flag
     integer :: istat

! logical file exist flag
     logical :: exists

! frequency grid
     real(dp), allocatable :: w(:)
     real(dp), allocatable :: dh(:)
     real(dp), allocatable :: logf(:)
     real(dp), allocatable :: delta(:)

! density of states data
     real(dp), allocatable :: dos(:,:)

! imaginary part of green's function
     real(dp), allocatable :: img(:,:)

! real part of green''s function
     real(dp), allocatable :: reg(:,:)

! print program header
     write(mystd,'(2X,a)') 'HIBISCUS/toolbox/makekra'
     write(mystd,'(2X,a)') '>>> Making kramer-kronig transformation for matsubara green''s function'
     write(mystd,*) ! print blank line

     write(mystd,'(2X,a)') 'Version: 2015.01.06T '//'(built at '//__TIME__//" "//__DATE__//')'
     write(mystd,'(2X,a)') 'Develop: by li huang (at IOP/CAS & SPCLab/CAEP & UNIFR)'
     write(mystd,'(2X,a)') 'Support: lihuang.dmft@gmail.com'
     write(mystd,'(2X,a)') 'License: GNU General Public License version 3'
     write(mystd,*) ! print blank line

! setup necessary parameters
     write(mystd,'(2X,a)')   'Number of orbitals (default = 2):'
     write(mystd,'(2X,a,$)') '>>> '
     read (*,*) nq
     write(mystd,*)

     write(mystd,'(2X,a)')   'Number of frequency points (default = 400):'
     write(mystd,'(2X,a,$)') '>>> '
     read (*,*) nw
     write(mystd,*)

! check the parameters
     call s_assert2( nq > 0 .and. nq < 15, 'wrong number of orbitals' )
     call s_assert2( nw > 0, 'wrong number of frequency points in half plane' )

! allocate memory
     allocate(w(-nw:nw),      stat=istat)
     allocate(dh(-nw:nw),     stat=istat)
     allocate(logf(-nw:nw),   stat=istat)
     allocate(delta(-nw:nw),  stat=istat)

     allocate(dos(-nw:nw,nq), stat=istat)
     allocate(img(-nw:nw,nq), stat=istat)
     allocate(reg(-nw:nw,nq), stat=istat)
     if ( istat /= 0 ) then
         call s_print_error('makekra','can not allocate enough memory')
     endif ! back if ( istat / = 0 ) block

! initialize arrays
     dos = zero
     img = zero
     reg = zero

! inquire data file
     inquire(file = 'mem.dos.dat', exist = exists)
     if ( exists .eqv. .false. ) then
         call s_print_error('makekra','file mem.dos.dat does not exist')
     endif ! back if ( exists .eqv. .false. ) block

! open density of states file
     write(mystd,'(2X,a)') 'Reading mem.dos.dat ...'
     open(mytmp, file='mem.dos.dat', form='formatted', status='unknown')

! read in density of states data
     do iq=1,nq
         do iw=-nw,nw
             read(mytmp,*) w(iw), dos(iw,iq)
         enddo ! over iw={-nw,nw} loop
         read(mytmp,*) ! skip two blank lines
         read(mytmp,*)
     enddo ! over iq={1,nq} loop

! close density of states file
     close(mytmp)
     write(mystd,'(2X,a)') '>>> status: OK'
     write(mystd,*)

! calculate imaginary part of green's function using the equation:
!    A(\omega) = -\frac{1}{\pi} IM G (\omega)
     img = -dos * pi

! build dh, logf, and delta meshes
     dh(-nw) = half * ( w(-nw+1) - w(-nw) )
     logf(-nw) = zero
     delta(-nw) = one / ( w(-nw+1) - w(-nw) )

     do iw=-nw+1,nw-1
         dh(iw) = half * ( w(iw+1) - w(iw-1) )
         logf(iw) = log( ( w(nw) - w(iw) ) / ( w(iw) - w(-nw) ) )
         delta(iw) = one / ( w(iw+1) - w(iw) )
     enddo ! over iw={-nw+1,nw-1} loop

     dh(nw) = half * ( w(nw) - w(nw-1) )
     logf(nw) = zero
     delta(nw) = zero

! perform the kramers-kronig transformations
     do iq=1,nq
         write(mystd,'(2X,a,i2)') 'Doing kramers-kronig transformation #', iq
         call kramers(nw, img(:,iq), reg(:,iq), w, dh, logf, delta)
         write(mystd,'(2X,a)')    '>>> status: OK'
         write(mystd,*)
     enddo ! over iq={1,nq} loop

! write out the green's function
     write(mystd,'(2X,a)') 'Writing kra.grn.dat ...'
     open(mytmp, file='kra.grn.dat', form='formatted', status='unknown')

     do iq=1,nq
         do iw=-nw,nw
             write(mytmp,'(3f16.8)') w(iw), reg(iw,iq), img(iw,iq)
         enddo ! over iw={-nw,nw} loop
         write(mytmp,*) ! write two blank lines
         write(mytmp,*)
     enddo ! over iq={1,nq} loop

     close(mytmp)
     write(mystd,'(2X,a)') '>>> status: OK'
     write(mystd,*)

! deallocate memory
     deallocate(w)
     deallocate(dh)
     deallocate(logf)
     deallocate(delta)

     deallocate(dos)
     deallocate(img)
     deallocate(reg)

  end program makekra

!!>>> kramers: implement the kramers-kronig transformation
  subroutine kramers(nw, img, reg, w, dh, logf, delta)
     use constants, only : dp, zero, half, pi

     implicit none

! external arguments
! number of frequency point
     integer, intent(in)   :: nw

! frequency grid
     real(dp), intent(in)  :: w(-nw:nw)
     real(dp), intent(in)  :: dh(-nw:nw)
     real(dp), intent(in)  :: logf(-nw:nw)
     real(dp), intent(in)  :: delta(-nw:nw)

! imaginary part of green''s function
     real(dp), intent(in)  :: img(-nw:nw)

! real part of green''s function
     real(dp), intent(out) :: reg(-nw:nw)

! local variables
! loop index
     integer  :: i
     integer  :: j

! integer dummy variables
     integer  :: ip1
     integer  :: im1

! real(dp) dummy variables
     real(dp) :: om1
     real(dp) :: om2
     real(dp) :: summ

     MESH1: do i=-nw,nw
! determine ip1 and im1
         if ( i < +nw ) then
             ip1 = i + 1
         else
             ip1 = i
         endif ! back if ( i < +nw ) block

         if ( i > -nw ) then
             im1 = i - 1
         else
             im1 = i
         endif ! back if ( i > -nw ) block

! determine om1 and om2
         if ( i > -nw ) then
             om1 = delta(i-1)
         else
             om1 = zero
         endif ! back if ( i > -nw ) block
         om2 = half * ( delta(i) * ( img(ip1) - img(i) ) + om1 * ( img(i) - img(im1) ) )

         summ = zero
         MESH2: do j=-nw,nw
             if ( i /= j ) then
                 summ = summ + ( img(j) - img(i) ) * dh(j) / ( w(j) - w(i) )
             else
                 summ = summ + om2 * dh(j)
             endif ! back if ( i /= j ) block
         enddo MESH2 ! over j={-nw,nw} loop
         reg(i) = ( summ + img(i) * logf(i) ) / pi
     enddo MESH1 ! over i={-nw,nw} loop

     return
  end subroutine kramers
