!-------------------------------------------------------------------------
! project : begonia
! program : ctqmc_dump_gtau
!           ctqmc_dump_wtau
!           ctqmc_dump_htau
!           ctqmc_dump_gbin
!           ctqmc_dump_grnf
!           ctqmc_dump_wssf
!           ctqmc_dump_hybf
!           ctqmc_dump_sigf
!           ctqmc_dump_hub1
!           ctqmc_dump_hist
!           ctqmc_dump_nmat
!           ctqmc_dump_prob
! source  : ctqmc_dump.f90
! type    : subroutine
! author  : li huang (email:huangli712@yahoo.com.cn)
! history : 09/16/2009 by li huang
!           09/17/2009 by li huang
!           09/18/2009 by li huang
!           09/20/2009 by li huang
!           09/22/2009 by li huang
!           10/25/2009 by li huang
!           11/01/2009 by li huang
!           11/30/2009 by li huang
!           12/01/2009 by li huang
!           12/04/2009 by li huang
!           12/09/2009 by li huang
!           12/26/2009 by li huang
!           12/30/2009 by li huang
!           02/28/2010 by li huang
!           03/04/2010 by li huang
!           08/23/2010 by li huang
! purpose : dump key observables produced by the hybridization expansion
!           version continuous time quantum Monte Carlo (CTQMC) quantum
!           impurity solver and dynamical mean field theory (DMFT) self
!           -consistent engine to disk files
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> write out impurity green's function in imaginary time space
  subroutine ctqmc_dump_gtau(tmesh, gtau)
     use constants
     use control

     implicit none

! external arguments
! imaginary time mesh
     real(dp), intent(in) :: tmesh(ntime)

! impurity green's function
     real(dp), intent(in) :: gtau(ntime,norbs,norbs)

! local variables
! loop index
     integer  :: i
     integer  :: j

! dummy variables
     real(dp) :: raux

! scaled impurity green's function
     real(dp) :: gaux(ntime,norbs,norbs)

! evaluate gaux first
     raux = real(ntime) / (beta * beta)
     do i=1,norbs
         do j=1,ntime
             gaux(j,i,i) = gtau(j,i,i) * raux
         enddo ! over j={1,ntime} loop
     enddo ! over i={1,norbs} loop

! open data file: solver.green.dat
     open(mytmp, file='solver.green.dat', form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,ntime
             write(mytmp,'(2i5,3f12.6)') i, j, tmesh(j), gaux(j,i,i), gaux(j,i+nband,i+nband)
         enddo ! over j={1,ntime} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_gtau

!>>> write out bath weiss's function in imaginary time space
  subroutine ctqmc_dump_wtau(tmesh, wtau)
     use constants
     use control

     implicit none

! external arguments
! imaginary time mesh
     real(dp), intent(in) :: tmesh(ntime)

! bath weiss's function
     real(dp), intent(in) :: wtau(ntime,norbs,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: solver.weiss.dat
     open(mytmp, file='solver.weiss.dat', form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,ntime
             write(mytmp,'(2i5,3f12.6)') i, j, tmesh(j), wtau(j,i,i), wtau(j,i+nband,i+nband)
         enddo ! over j={1,ntime} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_wtau

!>>> write out hybridization function in imaginary time space
  subroutine ctqmc_dump_htau(tmesh, htau)
     use constants
     use control

     implicit none

! external arguments
! imaginary time mesh
     real(dp), intent(in) :: tmesh(ntime)

! hybridization function
     real(dp), intent(in) :: htau(ntime,norbs,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: solver.hybri.dat
     open(mytmp, file='solver.hybri.dat', form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,ntime
             write(mytmp,'(2i5,3f12.6)') i, j, tmesh(j), htau(j,i,i), htau(j,i+nband,i+nband)
         enddo ! over j={1,ntime} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_htau

!>>> write out impurity green's function in imaginary time space (binning mode)
  subroutine ctqmc_dump_gbin(ibin, tmesh, gtau)
     use constants
     use control

     implicit none

! external arguments
! current bin index, integer representation
     integer, intent(in)  :: ibin

! imaginary time mesh
     real(dp), intent(in) :: tmesh(ntime)

! impurity green's function
     real(dp), intent(in) :: gtau(ntime,norbs,norbs)

! local variables
! loop index
     integer  :: i
     integer  :: j

! dummy variables
     real(dp) :: raux

! scaled impurity green's function
     real(dp) :: gaux(ntime,norbs,norbs)

! current bin index, string representation
     character(len=10) :: sbin

! evaluate gaux first
     raux = real(ntime) / (beta * beta)
     do i=1,norbs
         do j=1,ntime
             gaux(j,i,i) = gtau(j,i,i) * raux
         enddo ! over j={1,ntime} loop
     enddo ! over i={1,norbs} loop

! open data file: solver.green.bin.x
     write(sbin,'(i10)') ibin ! convert ibin to sbin
     open(mytmp, file='solver.green.bin.'//trim(adjustl(sbin)), form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,ntime
             write(mytmp,'(2i5,3f12.6)') i, j, tmesh(j), gaux(j,i,i), gaux(j,i+nband,i+nband)
         enddo ! over j={1,ntime} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_gbin

!>>> write out impurity green's function in matsubara frequency space
  subroutine ctqmc_dump_grnf(rmesh, grnf)
     use constants
     use control

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in) :: rmesh(mfreq)

! impurity green's function
     complex(dp), intent(in) :: grnf(mfreq,norbs,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: solver.grn.dat
     open(mytmp, file='solver.grn.dat', form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,mfreq
             write(mytmp,'(i5,5f16.8)') i, rmesh(j), &
                                  real(grnf(j,i,i)), &
                                 aimag(grnf(j,i,i)), &
                      real(grnf(j,i+nband,i+nband)), &
                     aimag(grnf(j,i+nband,i+nband))
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_grnf

!>>> write out bath weiss's function in matsubara frequency space
  subroutine ctqmc_dump_wssf(rmesh, wssf)
     use constants
     use control

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in) :: rmesh(mfreq)

! bath weiss's function
     complex(dp), intent(in) :: wssf(mfreq,norbs,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: solver.wss.dat
     open(mytmp, file='solver.wss.dat', form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,mfreq
             write(mytmp,'(i5,5f16.8)') i, rmesh(j), &
                                  real(wssf(j,i,i)), &
                                 aimag(wssf(j,i,i)), &
                      real(wssf(j,i+nband,i+nband)), &
                     aimag(wssf(j,i+nband,i+nband))
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_wssf

!>>> write out hybridization function in matsubara frequency space
  subroutine ctqmc_dump_hybf(rmesh, hybf)
     use constants
     use control

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in) :: rmesh(mfreq)

! hybridization function
     complex(dp), intent(in) :: hybf(mfreq,norbs,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: solver.hyb.dat
     open(mytmp, file='solver.hyb.dat', form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,mfreq
             write(mytmp,'(i5,5f16.8)') i, rmesh(j), &
                                  real(hybf(j,i,i)), &
                                 aimag(hybf(j,i,i)), &
                      real(hybf(j,i+nband,i+nband)), &
                     aimag(hybf(j,i+nband,i+nband))
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_hybf

!>>> write out self-energy function in matsubara frequency space
  subroutine ctqmc_dump_sigf(rmesh, sigf)
     use constants
     use control

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in) :: rmesh(mfreq)

! self-energy function
     complex(dp), intent(in) :: sigf(mfreq,norbs,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: solver.sgm.dat
     open(mytmp, file='solver.sgm.dat', form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,mfreq
             write(mytmp,'(i5,5f16.8)') i, rmesh(j), &
                                  real(sigf(j,i,i)), &
                                 aimag(sigf(j,i,i)), &
                      real(sigf(j,i+nband,i+nband)), &
                     aimag(sigf(j,i+nband,i+nband))
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_sigf

!>>> write out impurity green's function and self-energy function obtained
! by hubbard-I approximation in matsubara frequency space
  subroutine ctqmc_dump_hub1(rmesh, ghub, shub)
     use constants
     use control

     implicit none

! external arguments
! matsubara frequency mesh
     real(dp), intent(in)    :: rmesh(mfreq)

! impurity green's function by hubbard-I approximation
     complex(dp), intent(in) :: ghub(mfreq,norbs)

! self-energy function by hubbard-I approximation
     complex(dp), intent(in) :: shub(mfreq,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: solver.hub.dat
     open(mytmp, file='solver.hub.dat', form='formatted', status='unknown')

! write it
     do i=1,nband
         do j=1,mfreq
             write(mytmp,'(i5,5f16.8)') i, rmesh(j), &
                                    real(ghub(j,i)), &
                                   aimag(ghub(j,i)), &
                                    real(shub(j,i)), &
                                   aimag(shub(j,i))
         enddo ! over j={1,mfreq} loop
         write(mytmp,*) ! write empty lines
         write(mytmp,*)
     enddo ! over i={1,nband} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_hub1

!>>> write out the Monte Carlo sampling histogram for perturbation expansion series
  subroutine ctqmc_dump_hist(hist)
     use constants
     use control

     implicit none

! external arguments
! histogram data
     integer, intent(in) :: hist(mkink)

! local variables
! loop index
     integer  :: i

! dummy variables
     real(dp) :: raux

! scaled histogram data
     real(dp) :: haux(mkink)

! evaluate haux at first
     raux = real( sum(hist) )
     do i=1,mkink
         haux(i) = real( hist(i) ) / raux
     enddo ! over i={1,mkink} loop

! open data file: solver.hist.dat
     open(mytmp, file='solver.hist.dat', form='formatted', status='unknown')

! write it
     write(mytmp,'(a)') '# histogram: order | count | percent'
     do i=1,mkink
         write(mytmp,'(i5,i12,f12.6)') i, hist(i), haux(i)
     enddo ! over i={1,mkink} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_hist

!>>> write out the occupation matrix and double occupation matrix
  subroutine ctqmc_dump_nmat(nmat, nnmat)
     use constants
     use control

     implicit none

! external arguments
! occupation matrix data
     real(dp), intent(in) :: nmat(norbs)

! double occupation matrix data
     real(dp), intent(in) :: nnmat(norbs,norbs)

! local variables
! loop index
     integer :: i
     integer :: j

! open data file: solver.nmat.dat
     open(mytmp, file='solver.nmat.dat', form='formatted', status='unknown')

! write it
     write(mytmp,'(a)') '  < n_i >   data:'
     do i=1,norbs
         write(mytmp,'(i5,f12.6)') i, nmat(i)
     enddo ! over i={1,norbs} loop
     write(mytmp,'(a5,f12.6)') 'sup', sum( nmat(1:nband) )
     write(mytmp,'(a5,f12.6)') 'sdn', sum( nmat(nband+1:norbs) )
     write(mytmp,'(a5,f12.6)') 'sum', sum( nmat(1:norbs) )

     write(mytmp,'(a)') '< n_i n_j > data:'
     do i=1,norbs
         do j=1,norbs
             write(mytmp,'(2i5,f12.6)') i, j, nnmat(i,j)
         enddo ! over j={1,norbs} loop
     enddo ! over i={1,norbs} loop

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_nmat

!>>> write out the probability of eigenstates of local hamiltonian matrix
  subroutine ctqmc_dump_prob(prob, naux, saux)
     use constants
     use control

     implicit none

! external arguments
! probability data of eigenstates
     real(dp), intent(in) :: prob(ncfgs)

! occupation number of eigenstates
     real(dp), intent(in) :: naux(ncfgs)

! net spin of eigenstates
     real(dp), intent(in) :: saux(ncfgs)

! local variables
! loop index
     integer  :: i
     integer  :: j

! number of individual spin values
     integer  :: ns

! probability of occupation number distribution
     real(dp) :: oprob(0:norbs)

! dummy arrays, used to store spin of eigenstates
     real(dp) :: stmp1(ncfgs)
     real(dp) :: stmp2(ncfgs)

! probability of net spin distribution
     real(dp) :: sprob(ncfgs)

! evaluate oprob
     oprob = zero
     do i=1,ncfgs
         j = int( naux(i) )
         oprob(j) = oprob(j) + prob(i)
     enddo ! over i={1,ncfgs} loop

! sort all the spin values
     stmp1 = saux
     call ctqmc_time_qsorter(ncfgs, stmp1)

! find out individual spin values, and store them into stmp2
     ns = 1
     stmp2 = zero
     stmp2(1) = stmp1(1)
     do i=2,ncfgs
         if ( stmp2(ns) < stmp1(i) ) then
             ns = ns + 1
             stmp2(ns) = stmp1(i)
         endif
     enddo ! over i={2,ncfgs} loop

! evaluate sprob
     sprob = zero
     do i=1,ncfgs
         do j=1,ns
             if ( abs( stmp2(j) - saux(i) ) < eps6 ) then
                 sprob(j) = sprob(j) + prob(i)
                 EXIT
             endif
         enddo ! over j={1,ns} loop
     enddo ! over i={1,ncfgs} loop

! open data file: solver.prob.dat
     open(mytmp, file='solver.prob.dat', form='formatted', status='unknown')

! write it
     write(mytmp,'(a)') '# state probability: index | prob | occupy | spin'
     do i=1,ncfgs
         write(mytmp,'(i5,3f12.6)') i, prob(i), naux(i), saux(i)
     enddo ! over i={1,ncfgs} loop

     write(mytmp,'(a)') '# orbital probability: index | occupy | prob'
     do i=0,norbs
         write(mytmp,'(i5,2f12.6)') i+1, real(i), oprob(i)
     enddo ! over i={0,norbs} loop
     write(mytmp,'(a5,12X,f12.6)') 'sum', sum(oprob)

     write(mytmp,'(a)') '# spin probability: index | spin | prob'
     do i=1,ns
         write(mytmp,'(i5,2f12.6)') i, stmp2(i), sprob(i)
     enddo ! over i={1,ns} loop
     write(mytmp,'(a5,12X,f12.6)') 'sum', sum(sprob)

! close data file
     close(mytmp)

     return
  end subroutine ctqmc_dump_prob

!>>> dump the probability of sectors
  subroutine ctqmc_dump_psect()
     use constants
     use control
     use context

     implicit none

! local variables
! probability of sectors
     real(dp) :: psect(nsectors)

! loop index 
     integer :: i,j

! start index of sectors
     integer :: indx

     psect = zero
  
     do i=1, nsectors
         indx = sectors(i)%istart
         do j=1, sectors(i)%ndim
             psect(i) = psect(i) + prob(indx+j-1)
         enddo
     enddo      

! open file solver.psect.dat to write 
     open(mytmp, file='solver.psect.dat', form='unformatted', status='unknown')
     write(mytmp, '(a)') '#sector | probability | nelectron |'
     do i=1, nsectors
         write(mytmp, '(I10, F20.10, I10)') i, psect(i), sectors(i)%nelectron
     enddo
     close(mytmp)

     return
  end subroutine ctqmc_dump_psect
