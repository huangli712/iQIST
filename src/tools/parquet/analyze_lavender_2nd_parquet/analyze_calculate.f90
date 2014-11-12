!     	analyze_symmetrize_pairing_chi()
!       analyze_calculate_Gamma_local()
!       analyze_calculate_ladder_local()
!       analyze_calculate_lambda()
!       analyze_calculate_Fc()
!       analyze_symmetrize_pairing_Fc()

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine analyze_symmetrize_pairing_Gamma()
        use constants
        use control
        use context
        implicit none
! loop index for frequencies
        integer :: i,j,inu,ischeck=1

        complex(dp),allocatable :: g2_pp_tmp(:,:), gamma_pp(:,:,:)

        allocate(g2_pp_tmp(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2))
        allocate(gamma_pp(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2,-nbfrq+1:nbfrq-1))

! 1.1 gamma_pp(w,w') = 1/chi0(w,w+nu) - 1/chi(w,w+nu)
	    gamma_pp(:,:,:) = czero
	    
	do inu=-nbfrq+1,nbfrq-1	!-nbfrq+1,nbfrq-1

            g2_pp_tmp(:,:) = g2_pp_s(:,:,inu)
	    call analyze_zmat_inv(nffrq,g2_pp_tmp)
	    do i=-nffrq/2+1, nffrq/2 ! iw
	    do j=-nffrq/2+1, nffrq/2 ! iw'
            gamma_pp(i,j,inu) = - g2_pp_tmp(i,j)
              	if (i == j) then
                gamma_pp(i,j,inu) = gamma_pp_s(i,j,inu) + cone/((grnf(j+inu,2)*grnf(-j+1,2))*beta)
                end if ! i==j
	    end do ! over i={-nffrq/2+1, nffrq/2} loop
	    end do ! over j={-nffrq/2+1, nffrq/2} loop	
        end do ! over inu={-nbfrq+1,nbfrq-1} loop

! 1.2 gamma_pp_s(w,w',v) = gamma_pp(w,w',v) + gamma_pp(w,-w'-v,v)
!     gamma_pp_t(w,w',v) = gamma_pp(w,w',v) - gamma_pp(w,-w'-v,v)
	do inu=-nbfrq+1,nbfrq-1
            do i=-nffrq/2+1, nffrq/2 ! iw
            do j=-nffrq/2+1, nffrq/2 ! iw'
!!	prepare the singlet gamma_local
	    if (-j+1-inu.ge.-nffrq/2+1.and.-j+1-inu.le.nffrq/2) then
            gamma_pp_s(i,j,inu) = (gamma_pp(i,j,inu)+gamma_pp(i,-j+1-inu,inu))*beta*beta
!!	prepare the triplet gamma_local
	    gamma_pp_t(i,j,inu) = (gamma_pp(i,j,inu)-gamma_pp(i,-j+1-inu,inu))*beta*beta
	    end if
            end do ! over i={-nffrq/2+1, nffrq/2} loop
            end do ! over j={-nffrq/2+1, nffrq/2} loop  
         end do ! over inu={-nbfrq+1,nbfrq-1} loop

        deallocate(g2_pp_tmp)
        deallocate(gamma_pp)

!    print out gamma_pp_s for checking:
     if(ischeck.eq.1) then
     open(mytmp, file='analyze.gamma.pp.s.t.dat', form='formatted', status='unknown')
	  do inu=-nbfrq+1,nbfrq-1
	      write(mytmp,'(i5)') inu
	      do i=-nffrq/2+1, nffrq/2 ! iw
	      do j=-nffrq/2+1, nffrq/2 ! iw'
		write(mytmp,'(2i5,4f16.8)') 2*i-1, 2*j-1, gamma_pp_s(i,j,inu), gamma_pp_t(i,j,inu) 
	      end do
	      end do
		write(mytmp,*)
	  end do
     close(mytmp)
     end if! ischeck.eq.1

      return
      end subroutine analyze_symmetrize_pairing_Gamma
	
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine analyze_calculate_Gamma_local()
		use constants
		use control
		use context
		implicit none

! local variables
! loop index for frequencies
	integer :: i
	integer :: j
	integer :: k
	integer :: inu
	
! status flag
	integer :: istat
		
! dummy two-particle, intermediate, chi_magnet, chi_charge, chi_singlet
	complex(dp), allocatable :: chi_magnet(:,:)
	complex(dp), allocatable :: chi_charge(:,:)
!	complex(dp), allocatable :: chi_singlet(:,:)
	
! allocate memory
	allocate(chi_magnet(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2), stat=istat)
	allocate(chi_charge(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2), stat=istat)
!	allocate(chi_singlet(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2), stat=istat)
	if ( istat /= 0 ) then
	    call analyze_print_error('analyze_calculate_Gamma_local','can not allocate enough memory')
	endif
		
	chi_magnet(:,:) = czero
	chi_charge(:,:) = czero
!	chi_singlet(:,:) = czero
	
! open data files
	open(mytmp,file='analyze.twop.ph.charge.Gamma.local.dat',form='formatted',status='unknown')	
	open(mytmp1, file='analyze.twop.ph.magnet.Gamma.local.dat',form='formatted',status='unknown')
!	open(1002, file='analyze.twop.pp.singlet.Gamma.local.dat',form='formatted',status='unknown')

	    gamma_ph_c(:,:,:) = czero
	    gamma_ph_m(:,:,:) = czero
	
! here we consider positive and negative bosonic matsubara frequencies
	do inu=-nbfrq+1,nbfrq-1

	  if(inu.eq.0) write(mytmp,'(i5)') inu
	  if(inu.eq.0) write(mytmp1,'(i5)') inu
!!	   write(mytmp,'(i5)') inu
!!	   write(mytmp1,'(i5)') inu
!	  write(1002,'(i5)') inu
!  Step 1: Prepare chi_{w,w'}
	    chi_magnet(:,:) = czero
	    chi_charge(:,:) = czero
	    do i=-nffrq/2+1, nffrq/2 ! iw
	    do j=-nffrq/2+1, nffrq/2 ! iw'
	    
		chi_charge(i,j) = -Fc_ph_c(i,j,inu)*grnf(j+inu,2)*grnf(j,2)/beta
	   	chi_magnet(i,j) = -Fc_ph_m(i,j,inu)*grnf(j+inu,2)*grnf(j,2)/beta
!		chi_singlet(i,j) = -Fc_pp_s(i,j,inu)*grnf(j+inu,2)*grnf(-j+1,2)/(2.d0*beta)
		
		if (i == j) then
		    chi_charge(i,j) = cone + chi_charge(i,j)
		    chi_magnet(i,j) = cone + chi_magnet(i,j)
!		    chi_singlet(i,j) = cone + chi_singlet(i,j)
		end if 

	    end do ! over i={-nffrq/2+1, nffrq/2} loop
	    end do ! over j={-nffrq/2+1, nffrq/2} loop	
	    
	    call analyze_zmat_inv(nffrq,chi_charge)
	    call analyze_zmat_inv(nffrq,chi_magnet)
!	    call analyze_zmat_inv(nffrq,chi_singlet(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2))
	    
!  Step 2: sum_w'' {Step 1}_{w,w''} * F(inu)_{w'',w'}
	    
	    do i=-nffrq/2+1, nffrq/2 ! iw
	    do j=-nffrq/2+1, nffrq/2 ! iw'
	    do k=-nffrq/2+1, nffrq/2 ! iw'', sum over
		gamma_ph_c(i,j,inu) = gamma_ph_c(i,j,inu) + chi_charge(i,k)*Fc_ph_c(k,j,inu)
		gamma_ph_m(i,j,inu) = gamma_ph_m(i,j,inu) + chi_magnet(i,k)*Fc_ph_m(k,j,inu)
!		gamma_pp_s(i,j,inu) = gamma_pp_s(i,j,inu) + chi_singlet(i,k)*Fc_pp_s(k,j,inu)
	    end do ! over k={-nffrq/2+1, nffrq/2} loop
	    end do ! over i={-nffrq/2+1, nffrq/2} loop
	    end do ! over j={-nffrq/2+1, nffrq/2} loop
	    
!  Step 3: output local Gamma
	    do i=-nffrq/2+1, nffrq/2 ! iw
	    do j=-nffrq/2+1, nffrq/2 ! iw'
		if(inu.eq.0) write(mytmp,'(2i5,2f16.8)') 2*i-1, 2*j-1, gamma_ph_c(i,j,inu)
		if(inu.eq.0) write(mytmp1,'(2i5,2f16.8)') 2*i-1, 2*j-1, gamma_ph_m(i,j,inu)
!!		write(mytmp,'(2i5,2f16.8)') 2*i-1, 2*j-1, gamma_ph_c(i,j,inu)
!!		write(mytmp1,'(2i5,2f16.8)') 2*i-1, 2*j-1, gamma_ph_m(i,j,inu)
!		write(1002,'(2i5,2f16.8)') 2*i-1, 2*j-1, gamma_pp_s(i,j,inu)
	    end do
	    end do
	write(mytmp,*) ! write empty lines
	write(mytmp1,*) ! write empty lines
!	write(1002,*) ! write empty lines
	end do ! over inu={-nbfrq+1, nbfrq-1} loop		
	close(mytmp)	
	close(mytmp1)
!	close(1002)

        deallocate(chi_charge)
        deallocate(chi_magnet)
	  
	return
	end subroutine analyze_calculate_Gamma_local
	
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       subroutine analyze_calculate_ladder_local()
		use constants
		use control
		use context
		implicit none
! loop index for frequencies
	integer :: i,j
	integer :: iw1,iw2
	integer :: inu
	
! status flag
	integer :: istat
	complex(dp), allocatable :: ph_charge(:,:)
	complex(dp), allocatable :: ph_magnet(:,:)
	complex(dp), allocatable :: ph_singlet(:,:)
	complex(dp), allocatable :: ph_triplet(:,:)
! allocate memory
	allocate(ph_charge(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2), stat=istat)
	allocate(ph_magnet(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2), stat=istat)
	allocate(ph_singlet(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2), stat=istat)
	allocate(ph_triplet(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2), stat=istat)

	if ( istat /= 0 ) then
	    call analyze_print_error('analyze_calculate_ladder_local','can not allocate enough memory')
	endif

	ph_magnet = czero
	ph_charge = czero
	ph_singlet = czero
	ph_triplet = czero
!       For each fixed inu, calculate the ladder phi_local_X(inu;iw,iw'), where X = c and m.
        do inu=-nbfrq+1,nbfrq-1  ! inu, both positive and negative bosonic frequencies

!  Step 1: chi_X = { -1/GG - gamma_ph_X }^(-1)_{w,w'}
	    do i=-nffrq/2+1, nffrq/2 ! iw
	    do j=-nffrq/2+1, nffrq/2 ! iw'
		
		ph_charge(i,j) = -gamma_ph_c(i,j,inu)
		ph_magnet(i,j) = -gamma_ph_m(i,j,inu)
		ph_singlet(i,j) = -gamma_pp_s(i,j,inu)
		ph_triplet(i,j) = -gamma_pp_t(i,j,inu)
		
		if (i == j) then
		    ph_charge(i,j) = -cone*beta/(grnf(j+inu,2)*grnf(j,2)) + ph_charge(i,j)
		    ph_magnet(i,j) = -cone*beta/(grnf(j+inu,2)*grnf(j,2)) + ph_magnet(i,j)
		    ph_singlet(i,j) = 2.d0*beta/(grnf(j+inu,2)*grnf(-j+1,2)) + ph_singlet(i,j)
		    ph_triplet(i,j) = 2.d0*beta/(grnf(j+inu,2)*grnf(-j+1,2)) + ph_triplet(i,j)
		end if
		
	    end do ! over i={-nffrq/2+1, nffrq/2} loop
	    end do ! over j={-nffrq/2+1, nffrq/2} loop

	    call analyze_zmat_inv(nffrq,ph_magnet)
	    call analyze_zmat_inv(nffrq,ph_charge)
	    call analyze_zmat_inv(nffrq,ph_singlet)
	    call analyze_zmat_inv(nffrq,ph_triplet)
!   Now we have local chi_c, chi_m, chi_s and chi_t.

            do i=-nffrq/2+1, nffrq/2 ! iw
	    do j=-nffrq/2+1, nffrq/2 ! iw'
	    do iw1=-nffrq/2+1, nffrq/2 ! iw1
	    do iw2=-nffrq/2+1, nffrq/2 ! iw2
		phi_local_c(inu,i,j) = phi_local_c(inu,i,j) + gamma_ph_c(i,iw1,inu)*ph_charge(iw1,iw2)*gamma_ph_c(iw2,j,inu)
		phi_local_m(inu,i,j) = phi_local_m(inu,i,j) + gamma_ph_m(i,iw1,inu)*ph_magnet(iw1,iw2)*gamma_ph_m(iw2,j,inu)
		phi_local_s(inu,i,j) = phi_local_s(inu,i,j) + gamma_pp_s(i,iw1,inu)*ph_singlet(iw1,iw2)*gamma_pp_s(iw2,j,inu)
		phi_local_t(inu,i,j) = phi_local_t(inu,i,j) + gamma_pp_t(i,iw1,inu)*ph_triplet(iw1,iw2)*gamma_pp_t(iw2,j,inu)
	    end do ! over iw2={-nffrq/2+1, nffrq/2} loop
	    end do ! over iw1={-nffrq/2+1, nffrq/2} loop
	    end do ! over j={-nffrq/2+1, nffrq/2} loop
	    end do ! over i={-nffrq/2+1, nffrq/2} loop
        end do ! over inu={-nbfrq+1,nbfrq-1} loop
        
        phi_local_c(:,:,:) = phi_local_c(:,:,:)/(nffrq*nffrq)
        phi_local_m(:,:,:) = phi_local_m(:,:,:)/(nffrq*nffrq)
        phi_local_s(:,:,:) = phi_local_s(:,:,:)/(nffrq*nffrq)
        phi_local_t(:,:,:) = phi_local_t(:,:,:)/(nffrq*nffrq)

        deallocate(ph_charge)
        deallocate(ph_magnet)
        deallocate(ph_singlet)
        deallocate(ph_triplet)
       return
       end subroutine analyze_calculate_ladder_local


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       subroutine analyze_calculate_lambda
		use constants
		use control
		use context
		implicit none
! loop index for frequencies
	integer :: i,j,l,m
	integer :: inu,ischeck=1
	 
       do inu=-nbfrq+1,nbfrq-1 ! nu = (-7,-6,...,-1,0,1,...,6,7)

          do l=-nffrq/2+1,nffrq/2   ! w = (-3,-2,...,3,4)
          do m=-nffrq/2+1,nffrq/2   ! w'= (-3,-2,...,3,4)

!!	calculate fully irreducible vertex lambda(w,w')
          if(m-l.eq.inu) then         ! Phi(d/m)1(w'-w,k'-k)_-w',w+v ! no restriction in momentum
          lambda_s(l,m) = gamma_pp_s(l,m,0) + 0.5d0*phi_local_c(inu,-m+1,l) - 1.5d0*phi_local_m(inu,-m+1,l)
          lambda_t(l,m) = gamma_pp_t(l,m,0) + 0.5d0*phi_local_c(inu,-m+1,l+inu) + 0.5d0*phi_local_m(inu,-m+1,l+inu)
          lambda_m(l,m) = gamma_ph_m(l,m,0) - 0.5d0*phi_local_c(inu,l,l) + 0.5d0*phi_local_m(inu,l,l)
          lambda_c(l,m) = gamma_ph_c(l,m,0) - 0.5d0*phi_local_c(inu,l,l+inu) - 1.5d0*phi_local_m(inu,l,l+inu)
          else if(m+l-1.eq.inu) then  ! Phi(d/m)2(w'+w+v,k'+k+q)_-w',-w ! Note that the index of w+w'+v is (l+m-1+inu)
          lambda_s(l,m) = gamma_pp_s(l,m,0) + 0.5d0*phi_local_c(inu,-m+1,-l+1) - 1.5d0*phi_local_m(inu,-m+1,-l+1)
          lambda_t(l,m) = gamma_pp_t(l,m,0) - 0.5d0*phi_local_c(inu,-m+1,-l+1) - 0.5d0*phi_local_m(inu,-m+1,-l+1)
          lambda_m(l,m) = gamma_ph_m(l,m,0) - 0.5d0*phi_local_s(inu,-l+1,-l+1) + 0.5d0*phi_local_t(inu,-l+1,-l+1)
          lambda_c(l,m) = gamma_ph_c(l,m,0) + 0.5d0*phi_local_s(inu,-l+1,-l+1) + 1.5d0*phi_local_t(inu,-l+1,-l+1)
          end if 
          
!!	calculate fully irreducible vertex lambda(w,w'): triplet
!!          if(m-l.eq.inu) then         ! Phi(d/m)1(w'-w,k'-k)_-w',w ! no restriction in momentum
!!          lambdaF(l,m) = gamma_pp_s(l,m,0) + 0.5d0*phi_local_c(inu,-m+1,l) + 0.5d0*phi_local_m(inu,-m+1,l) 
!!          else if(m+l-1.eq.inu) then  ! Phi(d/m)2(w'+w,k'+k)_-w',-w ! Note that the index of w+w' is (l+m-1)
!!          lambdaF(l,m) = gamma_pp_s(l,m,0) - 0.5d0*phi_local_c(inu,-m+1,-l+1) - 0.5d0*phi_local_m(inu,-m+1,-l+1)  
!!          end if 

          end do ! over m={-nffrq/2+1,nffrq/2} loop
          end do ! over l={-nffrq/2+1,nffrq/2} loop

       end do ! over inu={-nbfrq+1,nbfrq-1} loop

!    print out lambda for checking:
     if(ischeck.eq.1) then
     open(mytmp, file='analyze.lambda.pp.s.t.dat', form='formatted', status='unknown')
	    do i=-nffrq/2+1, nffrq/2 ! iw
	    do j=-nffrq/2+1, nffrq/2 ! iw'
		write(mytmp,'(2i5,4f16.8)') 2*i-1, 2*j-1, lambda_s(i,j), lambda_t(i,j)
	    end do
	    end do  
     close(mytmp)
     end if! ischeck.eq.1

!    print out lambda.ph for checking:
     if(ischeck.eq.1) then
     open(mytmp, file='analyze.lambda.ph.c.m.dat', form='formatted', status='unknown')
	    do i=-nffrq/2+1, nffrq/2 ! iw
	    do j=-nffrq/2+1, nffrq/2 ! iw'
		write(mytmp,'(2i5,4f16.8)') 2*i-1, 2*j-1, lambda_c(i,j), lambda_m(i,j)
	    end do ! over i={-nffrq/2+1,nffrq/2} loop
	    end do ! over j={-nffrq/2+1,nffrq/2} loop
!!                write(mytmp,*) 
!!     end do  ! over inu={-nbfrq+1,nbfrq-1}
     close(mytmp)
     end if! ischeck.eq.1
     
       return
       end subroutine analyze_calculate_lambda


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!     subroutine analyze_calculate_Fc()
!!	use constants
!!	use control
!!	use context
!!	implicit none

! local variables
! loop index for frequencies
!!	integer :: i
!!	integer :: j
!!	integer :: k

! dummy complex(dp) variable, two-particle green's function, pp, singlet
!!	complex(dp) :: chit_pp_singlet_fc

! dummy complex(dp) variable, bare two-particle green's function
!!	complex(dp) :: chi0_pp_singlet

! open data files
!!	open(mytmp,file='analyze.twop.pp.pairing.Fc.dat',form='formatted', status='unknown')

! here we consider positive and negative bosonic matsubara frequencies
! Calculate Fc_pp_singlet from chi_pp_singlet
!!	do k=-nbfrq+1,nbfrq-1
!!	  write(mytmp,'(i5)') k
!!	do i=-nffrq/2+1,nffrq/2
!!	do j=-nffrq/2+1,nffrq/2
!!		if (k == 0) chi0 = chi0 - beta * grnf(j+k,2) * grnf(-j+1,2) 
!!		chi0_pp_singlet = czero
!!		if (i == j) chi0_pp_singlet = chi0_pp_singlet + beta * grnf(j+k,2) * grnf(-j+1,2)  ! Check chi0
!!		chit_pp_singlet_fc = g2_pp_s(i,j,k) - chi0_pp_singlet

!!		Fc_pp_s(i,j,k) = chit_pp_singlet_fc/(grnf(j+k,2)*grnf(-j+1,2)*grnf(-i+1,2)*grnf(i+k,2))!Check


!!		write(mytmp,'(2i5,2f16.8)') 2*i-1, 2*j-1, Fc_pp_s(i,j,k)
!!	end do ! over i={1,nffrq} loop
!!	end do ! over j={1,nffrq} loop
!!	write(mytmp,*) ! write empty lines
!!	end do ! over k={1,nbfrq} loop
	
!!	close(mytmp)

! Output Fc for charge and magnet
!!	open(mytmp,file='analyze.twop.ph.charge.Fc.dat', form='formatted', status='unknown')
!!	open(mytmp1,file='analyze.twop.ph.magnet.Fc.dat', form='formatted', status='unknown')
!!	do k=-nbfrq+1,nbfrq-1
!!	  write(mytmp,'(i5)') k
!!	  write(mytmp1,'(i5)') k
!!	do i=-nffrq/2+1,nffrq/2
!!	do j=-nffrq/2+1,nffrq/2
!!		write(mytmp,'(2i5,2f16.8)') 2*i-1, 2*j-1, Fc_ph_c(i,j,k)
!!		write(mytmp1,'(2i5,2f16.8)') 2*i-1, 2*j-1, Fc_ph_m(i,j,k)
!!	enddo ! over i={-nffrq/2+1,nffrq/2} loop
!!	enddo ! over j={-nffrq/2+1,nffrq/2} loop
!!	write(mytmp,*) ! write empty lines
!!	write(mytmp1,*) ! write empty lines
!!	enddo ! over k={-nbfrq+1,nbfrq-1} loop

!!	close(mytmp)
!!	close(mytmp1)
	
!!	return
!!	end subroutine analyze_calculate_Fc


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!      subroutine analyze_symmetrize_pairing_Fc()
!!        use constants
!!        use control
!!        use context
!!        implicit none
! loop index for frequencies
!!        integer :: i,j

!!        complex(dp),allocatable :: Fc_pp_tmp(:,:)

!!        allocate(Fc_pp_tmp(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2))

!!        Fc_pp_tmp(:,:) = Fc_pp_s(:,:,0)

!! Symmetrize Fc_pp_s(w,w',inu=0), note that we only symmetrize inu=0.
!!            do i=-nffrq/2+1, nffrq/2 ! iw
!!            do j=-nffrq/2+1, nffrq/2 ! iw'
!! prepare singlet full vertex
!!            Fc_pp_s(i,j,0) = (Fc_pp_tmp(i,j)+Fc_pp_tmp(i,-j+1))
!! prepera triplet full vertex
!!	    Fc_pp_s(i,j,0) = (Fc_pp_tmp(i,j)-Fc_pp_tmp(i,-j+1))
!!            end do ! over i={-nf_inh+1, nf_inh} loop
!!            end do ! over j={-nf_inh+1, nf_inh} loop    
        
!!        deallocate(Fc_pp_tmp)

!!      return
!!      end subroutine analyze_symmetrize_pairing_Fc


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine analyze_calculate_Gl_chi0tilde()
		use constants
		use control
		use context
		implicit none

! local variables
! loop index over orbitals, ktimes is the lattice momentum resoulution
	integer :: i, j, iw, inu, l, m, p, q
! function add_k_1D
        integer :: add_k_1D
	real(dp) :: rDk, kx, ky

! hopping parameters
	real(dp) :: t2, t3, t4, t5
	real(dp) :: Exy, Eyz, Ezx, E1D
	real(dp) :: E11, E13, E16, E22, E24, E25, E31, E33, E36
	real(dp) :: E42, E44, E45, E52, E54, E55, E61, E63, E66

	complex(dp) :: c11, c13, c16, c22, c24, c25, c31, c33, c36
	complex(dp) :: c42, c44, c45, c52, c54, c55, c61, c63, c66

! build matsubara frequency mesh: cmesh
	do iw=1,nffrq/2+nbfrq-1
	 cmesh(iw) = czi * (two * float(iw) - one) * (pi/beta)
	end do ! over iw={1,mfreq} loop
	
	grnf_lattice(:,:,:) = czero
	rDk = 1.0d0/real(ktime,dp)

! setup parameters
	Exy=0.0; Eyz=0.0; Ezx=0.0; E1D=0.0;
!!	this is the parameter set of the one-band square lattice
!!	t1=0.25d0; t2=-0.3*t1
!!	this is the parameter set of the typeII VHS paper
!!	t1=0.25d0; t2=-0.4d0*t1; t3=0.1d0*t1
!!	this is the parameter set of the t2g three bangs model
	t1=0.36; t2=0.18; t3=0.09; t4=0.37; t5=0.06;
	E11=0.0; E13=0.0; E16=0.0; E22=0.0; E24=0.0; E25=0.0; E31=0.0; E33=0.0; E36=0.0;
	E42=0.0; E44=0.0; E45=0.0; E52=0.0; E54=0.0; E55=0.0; E61=0.0; E63=0.0; E66=0.0;
	c11=czero; c13=czero; c16=czero; c22=czero; c24=czero; c25=czero; c31=czero; c33=czero; c36=czero;
	c42=czero; c44=czero; c45=czero; c52=czero; c54=czero; c55=czero; c61=czero; c63=czero; c66=czero;
	
	open(mytmp, file='analyze.grn_lattice.dat', form='formatted', status='unknown')
!!	open(mytmp1, file='analyze.grn_lattice_free.dat', form='formatted', status='unknown')
	do i=-ktime+1,ktime
	do j=-ktime+1,ktime
		kx = rDk*pi*(real(i,dp)-half)
		ky = rDk*pi*(real(j,dp)-half)
!! Here we use the 2D square lattice dispersion
!!		Exy = -2.d0*t1*(cos(kx)+cos(ky))
!!		Exy = -2.d0*t1*(cos(kx)+cos(ky))-2.d0*t2*(cos(kx+ky)+cos(kx-ky))

!! Here we use the typeII VHS paper dispersion
!!		Exy = -2.d0*t1*(cos(kx)+cos(ky))-2.d0*t2*(cos(kx+ky)+cos(kx-ky))-2.d0*t3*(cos(2.d0*kx)+cos(2.d0*ky))

!! Here we also have the 3 bands t2g dispersion
		Exy = -2.0d0*t1*(cos(kx)+cos(ky))-4.0d0*t2*cos(kx)*cos(ky)-2.0d0*t3*(cos(2*kx)+cos(2*ky)) 
		Ezx = -2.0d0*t4*cos(kx)-2.0d0*t5*cos(ky)
		Eyz = -2.0d0*t4*cos(ky)-2.0d0*t5*cos(kx)
		E55 = (Eyz+Ezx+Exy)/3.0d0 + lambda
		E66 = E55
		do iw=1,nffrq/2+nbfrq-1
		grnf_lattice(i,j,iw) = cone/(cmesh(iw)+mune-E55-sigf(iw,2))
		if(iw.eq.1) then
		write(mytmp,'(4f16.8)') kx, ky, grnf_lattice(i,j,1)
!!		write(mytmp1,'(4f16.8)') kx, ky, cone/(cmesh(k)-E55+mune)
		end if
		end do ! over iw={1,nffrq} loop
	end do ! over j={-ktime+1,ktime} loop
!!	write(mytmp,*)
!!	write(mytmp1,*)
	end do ! over i={-ktime+1,ktime} loop
	close(mytmp)
!!	close(mytmp1)

!! give the negative frequencies to grnf_lattice
	do i=-ktime+1,ktime
	do j=-ktime+1,ktime
	do iw = -nffrq/2+1-(nbfrq-1), 0
	    grnf_lattice(i,j,iw) = dconjg(grnf_lattice(i,j,-iw+1))
	end do
	end do
	end do
	
!! prepare the chi0_dm
	open(mytmp, file='analyze.chi0.ph.pp.dat', form='formatted', status='unknown')
	chi0_dm(:,:,:,:) = czero
	chi0_st(:,:,:,:) = czero
	do l=-ktime+1,ktime ! kx
	do m=-ktime+1,ktime ! ky
	   do iw=-nffrq/2+1,nffrq/2
	   do inu=-nbfrq+1,nbfrq-1 
	      do i=-ktime+1,ktime ! qx
	      do j=-ktime+1,ktime ! qy
		p = add_k_1D(l+i)
		q = add_k_1D(m+j)
		    chi0_dm(l,m,iw,inu) = chi0_dm(l,m,iw,inu) + grnf_lattice(p,q,iw+inu)*grnf_lattice(i,j,iw)
		    chi0_st(l,m,iw,inu) = chi0_st(l,m,iw,inu) + grnf_lattice(p,q,iw+inu)*grnf_lattice(-i+1,-j+1,-iw+1)
	      end do ! over j={-ktime+1,ktime} loop
	      end do ! over i={-ktime+1,ktime} loop
	    end do ! over inu={-nbfrq+1,nbfrq-1} loop
	    end do ! over iw ={1,nffrq} loop
	 end do ! over m={-ktime+1,ktime} loop
	 end do ! over l={-ktime+1,ktime} loop
	 
	 chi0_dm(:,:,:,:) = -chi0_dm(:,:,:,:)/(beta*real((2*ktime)**2,dp)) ! coarse-grain over k=(i,j)
	 chi0_st(:,:,:,:) = chi0_st(:,:,:,:)/(2*beta*real((2*ktime)**2,dp))  ! coarse-grain over k=(i,j)
		
	do l=-ktime+1,ktime
	do m=-ktime+1,ktime
		kx = rDk*pi*(real(l,dp)-half)
		ky = rDk*pi*(real(m,dp)-half)
		write(mytmp,'(6f16.8)') kx, ky, chi0_dm(l,m,1,0), chi0_st(l,m,1,0)
	end do ! over m={-ktime,ktime} loop
	end do ! over l={-ktime,ktime} loop
	close(mytmp)
	

	return
	end subroutine analyze_calculate_Gl_chi0tilde

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine analyze_calculate_vertexladder()
		use constants
		use control
		use context
		implicit none
!
! phi_ph_X(qx,qy,nu;w,w') = sum_{w1,w2} { gamma_ph_X(w,w1,nu)*chi_ph_X(qx,qy,nu,w1,w2)*gamma_ph_X(w2,w',nu) }
! where chi_ph_X(qx,qy,nu,w,w') = (1/chi0_dm(qx,qy,nu,w) - gamma_ph_X(w,w',nu))^-1
! phi_pp_X(qx,qy,nu;w,w') = sum_{w1,w2} { gamma_pp_X(w,w1,nu)*chi_pp_X(qx,qy,nu,w1,w2)*gamma_pp_X(w2,w',nu) }
! where chi_pp_X(qx,qy,nu,w,w') = (1/chi0_st(qx,qy,nu,w) - gamma_pp_X(w,w',nu))^-1
!
! local variables
! loop index over momentum, frequency, orbital spaces
		integer :: i, j, l, m, inu, iw1, iw2
! status flag
		integer :: istat
		real(dp) :: kx, ky, rDk
		
! dummy two-particle, intermediate, ph_magnet, ph_charge
		complex(dp), allocatable :: ph_magnet(:,:)
		complex(dp), allocatable :: ph_charge(:,:)
		complex(dp), allocatable :: ph_singlet(:,:)
		complex(dp), allocatable :: ph_triplet(:,:)
		
! allocate memory
	allocate(ph_magnet(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2), stat=istat)
	allocate(ph_charge(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2), stat=istat)
	allocate(ph_singlet(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2), stat=istat)
	allocate(ph_triplet(-nffrq/2+1:nffrq/2,-nffrq/2+1:nffrq/2), stat=istat)
	
	if ( istat /= 0 ) then
	    call analyze_print_error('analyze_calculate_vertexladder','can not allocate enough memory')
	endif
		
	phi_ph_c(:,:,:,:,:) = czero
	phi_ph_m(:,:,:,:,:) = czero
	phi_pp_s(:,:,:,:,:) = czero
	phi_pp_t(:,:,:,:,:) = czero

	    ph_magnet(:,:) = czero
	    ph_charge(:,:) = czero
	    ph_singlet(:,:) = czero
	    ph_triplet(:,:) = czero
	    
!! prepare the momentum-dependence of the chi(q)
	rDk = 1.0d0/real(ktime,dp)
	open(mytmp, file='analyze.twop.chi.ph.c.m.dat', form='formatted', status='unknown')
	open(mytmp1,file='analyze.twop.chi.pp.s.t.dat', form='formatted', status='unknown')
!       For each fixed (qx,qy) and inu, calculate the ladder phi_ph_X(qx,qy,inu;iw,iw'), where X = c and m.
!       For each fixed (qx,qy) and inu, calculate the ladder phi_pp_X(qx,qy,inu;iw,iw'), where X = s and t.
	do l=-ktime+1,ktime ! qx
	do m=-ktime+1,ktime ! qy
	  kx = rDk*pi*(real(l,dp)-half)
	  ky = rDk*pi*(real(m,dp)-half)
        do inu=-nbfrq+1,nbfrq-1  ! inu, both positive and negative bosonic frequencies
	
!  Step 1: Prepare 1/chi_ph_X(w,w') = 1/chi0_dm(qx,qy,w,nu) - gamma_ph_X(w,w',nu)
!		   1/chi_pp_X(w,w') = 1/chi0_st(qx,qy,w,nu) - gamma_pp_X(w,w',nu)
	    do i=-nffrq/2+1, nffrq/2 ! iw
	    do j=-nffrq/2+1, nffrq/2 ! iw'

		ph_charge(i,j) = -gamma_ph_c(i,j,inu)
		ph_magnet(i,j) = -gamma_ph_m(i,j,inu)
		ph_singlet(i,j) = -gamma_pp_s(i,j,inu)
		ph_triplet(i,j) = -gamma_pp_t(i,j,inu)
		
		if (i == j) then
		    ph_charge(i,j) = cone/chi0_dm(l,m,i,inu) + ph_charge(i,j)
		    ph_magnet(i,j) = cone/chi0_dm(l,m,i,inu) + ph_magnet(i,j)
		    ph_singlet(i,j) = cone/chi0_st(l,m,i,inu) + ph_singlet(i,j)
		    ph_triplet(i,j) = cone/chi0_st(l,m,i,inu) + ph_triplet(i,j)
		end if
		
	    end do ! over i={-nffrq/2+1, nffrq/2} loop
	    end do ! over j={-nffrq/2+1, nffrq/2} loop

	    call analyze_zmat_inv(nffrq,ph_magnet)
	    call analyze_zmat_inv(nffrq,ph_charge)
	    call analyze_zmat_inv(nffrq,ph_singlet)
	    call analyze_zmat_inv(nffrq,ph_triplet)
	    	    
	    if (inu == 0) then
	      write(mytmp,'(6f16.8)') kx, ky, ph_charge(1,1), ph_magnet(1,1)
	      write(mytmp1,'(6f16.8)')kx, ky, ph_singlet(1,1), ph_triplet(1,1)
	    end if

!  Step 2: sum_{w1,w2} gamma_ph_X{w,w1} * chi_ph_X(inu)_{qx,qy,w1,w2} * gamma_ph_X(w2,w')
	    do i=-nffrq/2+1, nffrq/2 ! iw
	    do j=-nffrq/2+1, nffrq/2 ! iw'
	    do iw1=-nffrq/2+1, nffrq/2 ! iw1
	    do iw2=-nffrq/2+1, nffrq/2 ! iw2
		phi_ph_c(l,m,inu,i,j) = phi_ph_c(l,m,inu,i,j) - gamma_ph_c(i,iw1,inu)*ph_charge(iw1,iw2)*gamma_ph_c(iw2,j,inu)
		phi_ph_m(l,m,inu,i,j) = phi_ph_m(l,m,inu,i,j) - gamma_ph_m(i,iw1,inu)*ph_magnet(iw1,iw2)*gamma_ph_m(iw2,j,inu)
		phi_pp_s(l,m,inu,i,j) = phi_pp_s(l,m,inu,i,j) - gamma_pp_s(i,iw1,inu)*ph_singlet(iw1,iw2)*gamma_pp_s(iw2,j,inu)
		phi_pp_t(l,m,inu,i,j) = phi_pp_t(l,m,inu,i,j) - gamma_pp_t(i,iw1,inu)*ph_triplet(iw1,iw2)*gamma_pp_t(iw2,j,inu)
	    end do ! over iw2={-nffrq/2+1, nffrq/2} loop
	    end do ! over iw1={-nffrq/2+1, nffrq/2} loop
	    end do ! over j={-nffrq/2+1, nffrq/2} loop
	    end do ! over i={-nffrq/2+1, nffrq/2} loop

	
	end do ! over inu={-nbfrq+1,nbfrq-1}
	end do ! over m={-ktime+1,ktime}
	end do ! over l={-ktime+1,ktime} 	
	close(mytmp)
	close(mytmp1)
	
	open(mytmp, file='analyze.twop.phi.ph.c.m.frequency.dat', form='formatted', status='unknown')
	do i=-ktime+1,ktime
	do j=-ktime+1,ktime
	  kx = rDk*pi*(real(i,dp)-half)
	  ky = rDk*pi*(real(j,dp)-half)
	    write(mytmp,'(6f16.8)') kx, ky, phi_ph_c(i,j,0,0,0), phi_ph_m(i,j,0,0,0)
	end do ! over n={-ktime+1,ktime} loop
	end do ! over m={-ktime+1,ktime} loop
	close(mytmp)
	 
	open(mytmp, file='analyze.twop.phi.pp.s.t.frequency.dat', form='formatted', status='unknown')
	do i=-ktime+1,ktime
	do j=-ktime+1,ktime
	  kx = rDk*pi*(real(i,dp)-half)
	  ky = rDk*pi*(real(j,dp)-half)
	    write(mytmp,'(6f16.8)') kx, ky, phi_pp_s(i,j,-1,1,1), phi_pp_t(i,j,-1,1,1)
	end do ! over n={-ktime+1,ktime} loop
	end do ! over m={-ktime+1,ktime} loop
	close(mytmp)
	
	phi_ph_c(:,:,:,:,:) = phi_ph_c(:,:,:,:,:)/(real(nffrq*nffrq,dp))
	phi_ph_m(:,:,:,:,:) = phi_ph_m(:,:,:,:,:)/(real(nffrq*nffrq,dp))
	phi_pp_s(:,:,:,:,:) = phi_pp_s(:,:,:,:,:)/(real(nffrq*nffrq,dp))
	phi_pp_t(:,:,:,:,:) = phi_pp_t(:,:,:,:,:)/(real(nffrq*nffrq,dp))
   
	open(mytmp, file='analyze.twop.phi.ph.c.m.dat', form='formatted', status='unknown')
	do i=-ktime+1,ktime
	do j=-ktime+1,ktime
	  kx = rDk*pi*(real(i,dp)-half)
	  ky = rDk*pi*(real(j,dp)-half)
	  write(mytmp,'(6f16.8)') kx, ky, phi_ph_c(i,j,0,1,1), phi_ph_m(i,j,0,1,1)
	end do ! over n={-ktime+1,ktime} loop
	end do ! over m={-ktime+1,ktime} loop
	 close(mytmp)

	 open(mytmp, file='analyze.twop.phi.pp.s.t.dat', form='formatted', status='unknown')
	do i=-ktime+1,ktime
	do j=-ktime+1,ktime
	   kx = rDk*pi*(real(i,dp)-half)
	   ky = rDk*pi*(real(j,dp)-half)
	  write(mytmp,'(6f16.8)') kx, ky, phi_pp_s(i,j,0,1,1), phi_pp_t(i,j,0,1,1)
	end do ! over n={-ktime+1,ktime} loop
	end do ! over m={-ktime+1,ktime} loop
	close(mytmp)
	
! deallocate memory
        deallocate(ph_magnet)
        deallocate(ph_charge)
        deallocate(ph_singlet)
	deallocate(ph_triplet)
        write(*,*) "finish analyze_calculate_vertexladder!!!"

	return
	end subroutine analyze_calculate_vertexladder


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     subroutine analyze_calculate_Gamma_1st_order()
	use constants
	use control
	use context

	implicit none

! local variables
! loop index over momentum, frequency, orbital spaces
      integer :: i, j, k, q, l, m, ip, lp, p1, p2, p3, q1, &
                 iqx,iqy,ikx1,iky1,ikx2,iky2,ikx3,iky3,iqx1,iqy1,iqx2,iqy2,iqx3,iqy3,inu
! function add_k_1D
      integer :: add_k_1D
      real(dp) :: kx, ky, rDk
!    Fix q=0 and nu=0, and
!    Gamma_s(k,k',w,w') = Fc_pp_s(k,k',w,w')+(1/2)(Phid1 + Phid2)-(3/2)(Phim1 + Phim2),
!    where
!    Phi(d/m)1(q=0,nu=0)_p,p' = Phi(d/m)1(w'-w,k'-k)_-w',w
!    Phi(d/m)2(q=0,nu=0)_p,p' = Phi(d/m)2(w'+w,k'+k)_-w',-w
 
     do i=1,4*ktime**2         ! k = (kx,ky)
     ikx1=(mod(i-1,2*ktime) + 1) - ktime ! kx: (-ktime+1:ktime)
     iky1=((i-1)/(2*ktime) + 1) - ktime  ! ky: (-ktime+1:ktime)
     do j=1,4*ktime**2           ! k' = (kx',ky')
     ikx2=(mod(j-1,2*ktime) + 1) - ktime ! kx': (-ktime+1:ktime)
     iky2=((j-1)/(2*ktime) + 1) - ktime  ! ky': (-ktime+1:ktime)
     
!     do k=1,4*ktime**2		! q = (qx,qy)
!	 k=4*ktime**2		! this is for q=(pi,pi)
!	 k=2*ktime**2 - ktime	! this is for q=(0,0)
!      ikx3 = (mod(k-1,2*ktime) + 1) - ktime ! qx : (-ktime+1:ktime)
!      iky3 = ((k-1)/(2*ktime) + 1) - ktime  ! qy : (-ktime+1:ktime)

       iqx1=add_k_1D(ikx2 - ikx1) ! kx'-kx
       iqy1=add_k_1D(iky2 - iky1) ! ky'-ky
       iqx2=add_k_1D(ikx2 + ikx1) ! kx'+kx
       iqy2=add_k_1D(iky2 + iky1) ! ky'+kx
       
!      q1= (k-1)*(2*nbfrq-1) + 0 + nbfrq  ! q1 = (q,v=0) (1,2,...,4*ktime**2*(2*nbfrq-1))
       
       do inu=-nbfrq+1,nbfrq-1 ! nu = (-7,-6,...,-1,0,1,...,6,7)
         
          do l=-nffrq/2+1,nffrq/2   ! w = (-3,-2,...,3,4)
          p1= (i-1)*nffrq + l + nffrq/2 ! p1 = (k,w) = (1,2,...,4*ktime**2*nffrq)
          do m=-nffrq/2+1,nffrq/2   ! w'= (-3,-2,...,3,4)
          p2= (j-1)*nffrq + m + nffrq/2 ! p2 = (k',w')=(1,2,...,4*ktime**2*nffrq)

!!	construct the particle-hole magnetic irreducible vertex
          if(m-l.eq.inu) then         ! Phi(d/m)1(w'-w,k'-k)_-w',w ! no restriction in momentum
          gamat_m(p1,p2)= lambda_m(l,m) + 0.5d0*phi_ph_c(iqx1,iqy1,inu,l,l) - 0.5d0*phi_ph_m(iqx1,iqy1,inu,l,l)
          gamat_c(p1,p2)= lambda_c(l,m) + 0.5d0*phi_ph_c(iqx1,iqy1,inu,l,l) + 1.5d0*phi_ph_m(iqx1,iqy1,inu,l,l)
          else if(m+l-1.eq.inu) then  ! Phi(d/m)2(w'+w,k'+k)_-w',-w ! Note that the index of w+w' is (l+m-1)
          gamat_m(p1,p2)= lambda_m(l,m) + 0.5d0*phi_pp_s(iqx2,iqy2,inu,-l+1,-l+1) - 0.5d0*phi_pp_t(iqx2,iqy2,inu,-l+1,-l+1) 
          gamat_c(p1,p2)= lambda_c(l,m) - 0.5d0*phi_pp_s(iqx2,iqy2,inu,-l+1,-l+1) - 1.5d0*phi_pp_t(iqx2,iqy2,inu,-l+1,-l+1) 
          end if 
          
!!	construct the triplet pairing irreducible vertex
!!          if(m-l.eq.inu) then         ! Phi(d/m)1(w'-w,k'-k)_-w',w ! no restriction in momentum
!!          gamat_s(p1,p2)= lambdaF(l,m) - 0.5d0*phi_ph_c(iqx1,iqy1,inu,-m+1,l) - 0.5d0*phi_ph_m(iqx1,iqy1,inu,-m+1,l) 
!!          else if(m+l-1.eq.inu) then  ! Phi(d/m)2(w'+w,k'+k)_-w',-w ! Note that the index of w+w' is (l+m-1)
!!          gamat_s(p1,p2)= lambdaF(l,m) + 0.5d0*phi_ph_c(iqx2,iqy2,inu,-m+1,-l+1) + 0.5d0*phi_ph_m(iqx2,iqy2,inu,-m+1,-l+1)
!!          end if

          end do ! over m={-nffrq/2+1,nffrq/2} loop
          end do ! over l={-nffrq/2+1,nffrq/2} loop

       end do ! over inu={-nbfrq+1,nbfrq-1} loop

     end do ! over j={1, 4*ntime**2} loop
     end do ! over i={1, 4*ntime**2} loop

        write(*,*) "finish analyze_calculate 1st order: Gamma_m, Gamma_c !!!"
	
!!	write(*,*) "check the momentum dependence of gamat_m"
!!     rDk = 1.0d0/real(ktime,dp)
!!     open(mytmp, file='analyze.twop.gamat_m.kp.dat', form='formatted', status='unknown')
!!     do i=1,4*ktime**2         ! k = (kx,ky)
!!       i = 1
!!       ikx1=(mod(i-1,2*ktime) + 1) - ktime ! kx: (-ktime+1:ktime)
!!       iky1=((i-1)/(2*ktime) + 1) - ktime  ! ky: (-ktime+1:ktime)
!!     do j=1,4*ktime**2           ! k' = (kx',ky')
!!       j = ktime
!!       ikx2=(mod(j-1,2*ktime) + 1) - ktime ! kx': (-ktime+1:ktime)
!!       iky2=((j-1)/(2*ktime) + 1) - ktime  ! ky': (-ktime+1:ktime)
	 
!!          do l=-nffrq/2+1,nffrq/2   ! w = (-3,-2,...,3,4)
!!	  l=1
!!          p1= (i-1)*nffrq + l + nffrq/2 ! p1 = (k,w) = (1,2,...,4*ktime**2*nffrq)
!!          do m=-nffrq/2+1,nffrq/2   ! w'= (-3,-2,...,3,4)
!!	  m=1
!!          p2= (j-1)*nffrq + m + nffrq/2 ! p2 = (k',w')=(1,2,...,4*ktime**2*nffrq)
            
!!          kx = rDk*pi*(real(ikx2,dp)-half)
!!	  ky = rDk*pi*(real(iky2,dp)-half)
!!	  write(mytmp,'(4f16.8)') kx, ky, gamat_s(p1,p2)

!!          end do ! over m={-nffrq/2+1,nffrq/2} loop
!!          end do ! over l={-nffrq/2+1,nffrq/2} loop

!!     end do ! over j={1, 4*ntime**2} loop
!!     end do ! over i={1, 4*ntime**2} loop
     
!!     close(mytmp)
     return
     end subroutine analyze_calculate_Gamma_1st_order

     
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subroutine analyze_calculate_vertexladder_2nd_order
	use constants
	use control
	use context
	implicit none
	
!
! phi_ph_k_kp(k,w,k',w',q,v) = sum_{(k1,w1),(k2,w2)} &
!	{ gamma_ph_X(k,w,k1,w1)*chi_ph_X(k1,w1,k2,w2,q,v)*gamma_ph_X(k2,w2,k',w') }
! where chi_ph_X(k,w,k',w',q,v) = (1/chi0_dm(k,w,q,v) - gamma_ph_X(k,w,k',w'))^-1
     
! local variables
! loop index over momentum, frequency, orbital spaces
      integer :: i, j, k, q, l, m, ip, lp, p1,p1m,p2,p2m,q1,q2, &
      &          ikx1,iky1,ikx2,iky2,iqx,iqy,iqx1,iqy1,iqx2,iqy2,iqx3,iqy3,inu
      integer :: i3,j4,l3,m4,p3,p4,ikx3,iky3,ikx4,iky4
      
! function add_k_1D
      integer :: add_k_1D
! status flag
		integer :: istat
		real(dp) :: kx, ky, rDk
     
! dummy two-particle, intermediate, ph_magnet
		complex(dp), allocatable :: ph_magnet(:,:), phi_ph_m_k_kp(:,:), phi_ph_m_k_kp_1(:,:)
! dummy two-particle, intermediate, ph_charge
		complex(dp), allocatable :: ph_charge(:,:), phi_ph_c_k_kp(:,:), phi_ph_c_k_kp_1(:,:)
! vertex ladder, phi_ph_m_k_kp_q
		complex(dp), allocatable :: phi_ph_m_k_kp_q(:,:,:)
! vertex ladder, phi_ph_c_k_kp_q
		complex(dp), allocatable :: phi_ph_c_k_kp_q(:,:,:)
     		
! allocate memory
	allocate(ph_magnet(4*ktime*ktime*nffrq,4*ktime*ktime*nffrq), stat=istat)
	allocate(phi_ph_m_k_kp(4*ktime*ktime*nffrq,4*ktime*ktime*nffrq), stat=istat)
	allocate(phi_ph_m_k_kp_1(4*ktime*ktime*nffrq,4*ktime*ktime*nffrq), stat=istat)
	allocate(phi_ph_m_k_kp_q(4*ktime*ktime*nffrq,4*ktime*ktime*nffrq,4*ktime*ktime*(2*nbfrq-1)), stat=istat)
	
	allocate(ph_charge(4*ktime*ktime*nffrq,4*ktime*ktime*nffrq), stat=istat)
	allocate(phi_ph_c_k_kp(4*ktime*ktime*nffrq,4*ktime*ktime*nffrq), stat=istat)
	allocate(phi_ph_c_k_kp_1(4*ktime*ktime*nffrq,4*ktime*ktime*nffrq), stat=istat)
	allocate(phi_ph_c_k_kp_q(4*ktime*ktime*nffrq,4*ktime*ktime*nffrq,4*ktime*ktime*(2*nbfrq-1)), stat=istat)
	
	if ( istat /= 0 ) then
	    call analyze_print_error('analyze_calculate_vertexladder_2nd_order','can not allocate enough memory')
	endif

	ph_magnet(:,:) = czero
	phi_ph_m_k_kp(:,:) = czero
	phi_ph_m_k_kp_1(:,:) = czero
	phi_ph_m_k_kp_q(:,:,:) = czero
	
	ph_charge(:,:) = czero
	phi_ph_c_k_kp(:,:) = czero
	phi_ph_c_k_kp_1(:,:) = czero
	phi_ph_c_k_kp_q(:,:,:) = czero
	
!! prepare the momentum-dependence of ladder phi_ph_c_k_kp(k,w,k',w',q,v)
!! prepare the momentum-dependence of ladder phi_ph_m_k_kp(k,w,k',w',q,v)

	rDk = 1.0d0/real(ktime,dp)
!!	open(mytmp, file='analyze.twop.chi.ph.m.q.dat', form='formatted', status='unknown')
!!	open(mytmp1, file='analyze.twop.phi.ph.m.q.dat', form='formatted', status='unknown')
	
	do q=1,4*ktime**2		 ! q  = (qx,qy)
	iqx=(mod(q-1,2*ktime) + 1) - ktime ! qx : (-ktime+1:ktime)
	iqy=((q-1)/(2*ktime) + 1) - ktime  ! qy : (-ktime+1:ktime)
	
	  kx = rDk*pi*(real(iqx,dp)-half)
	  ky = rDk*pi*(real(iqy,dp)-half)

	do inu=-nbfrq+1,nbfrq-1 ! nu = (-7,-6,...,-1,0,1,...,6,7)
!!	  inu=0
	  q1= (q-1)*(2*nbfrq-1) + inu + nbfrq  ! q1 = (q,v) = (1,2,...,4*ktime**2*(2*nbfrq-1))
	  
!  Step 1: Prepare 1/chi_ph_k_kp(k,w,k',w',q,nu) = 1/G(k+q,w+nu)G(k,w) - gamma_ph_X(k,w,k'w',q=0,nu=0)

	    do i=1,4*ktime**2         ! k = (kx,ky)
		ikx1=(mod(i-1,2*ktime) + 1) - ktime ! kx: (-ktime+1:ktime)
		iky1=((i-1)/(2*ktime) + 1) - ktime  ! ky: (-ktime+1:ktime)
				
		iqx1=add_k_1D(ikx1 + iqx) ! kx1 + qx
		iqy1=add_k_1D(iky1 + iqy) ! ky1 + qy
		
	    do j=1,4*ktime**2           ! k' = (kx',ky')
		ikx2=(mod(j-1,2*ktime) + 1) - ktime ! kx': (-ktime+1:ktime)
		iky2=((j-1)/(2*ktime) + 1) - ktime  ! ky': (-ktime+1:ktime)
		
	      do l=-nffrq/2+1,nffrq/2   ! w = (-3,-2,...,3,4)
		  p1= (i-1)*nffrq + l + nffrq/2 ! p1 = (k,w) = (1,2,...,4*ktime**2*nffrq)
	      do m=-nffrq/2+1,nffrq/2   ! w'= (-3,-2,...,3,4)
		  p2= (j-1)*nffrq + m + nffrq/2 ! p2 = (k',w')=(1,2,...,4*ktime**2*nffrq)

		  ph_magnet(p1,p2) = -gamat_m(p1,p2)
		  ph_charge(p1,p2) = -gamat_c(p1,p2)
		  
		if ((i.eq.j).and.(l.eq.m)) then
		    ph_magnet(p1,p2) = -cone/(grnf_lattice(iqx1,iqy1,l+inu)*grnf_lattice(ikx1,iky1,l)) + ph_magnet(p1,p2)
		    ph_charge(p1,p2) = -cone/(grnf_lattice(iqx1,iqy1,l+inu)*grnf_lattice(ikx1,iky1,l)) + ph_charge(p1,p2)
		end if
		
	      end do ! over l={-nffrq/2+1, nffrq/2} loop
	      end do ! over m={-nffrq/2+1, nffrq/2} loop

	     end do ! over j={1,4*ktime**2} loop
	     end do ! over i={1,4*ktime**2} loop
	      
	    call analyze_zmat_inv(4*ktime*ktime*nffrq,ph_magnet)
	    call analyze_zmat_inv(4*ktime*ktime*nffrq,ph_charge)
	   
!!	    write(6,*) "obtian the chi_ph_m_k_kp(k,w,k',w',q,v) for one (q,v)"
!!	    write(mytmp,'(6f16.8)') kx, ky, ph_magnet(1,1)
	
!  Step 2: sum_{k1,w1,k2,w2} gamma_ph_X{k,w,k1,w1} * chi_ph_X(inu)_{qx,qy,k1,w1,k2,w2} * gamma_ph_X(k2,w2,k,w')

	    call analyze_zmat_prd(4*ktime**2*nffrq,gamat_m,ph_magnet,phi_ph_m_k_kp_1)
	    call analyze_zmat_prd(4*ktime**2*nffrq,phi_ph_m_k_kp_1,gamat_m,phi_ph_m_k_kp)
	    
	    call analyze_zmat_prd(4*ktime**2*nffrq,gamat_c,ph_charge,phi_ph_c_k_kp_1)
	    call analyze_zmat_prd(4*ktime**2*nffrq,phi_ph_c_k_kp_1,gamat_c,phi_ph_c_k_kp)
	    
	    phi_ph_m_k_kp_q(:,:,q1) = phi_ph_m_k_kp(:,:)
	    phi_ph_c_k_kp_q(:,:,q1) = phi_ph_c_k_kp(:,:)
	    
	end do ! over inu={-nbfrq+1,nbfrq-1}
	end do ! over q={1,4*ktime**2} loop
	
	! deallocate memory
        deallocate(ph_magnet)
        deallocate(phi_ph_m_k_kp)
        deallocate(phi_ph_m_k_kp_1)
        deallocate(ph_charge)
        deallocate(phi_ph_c_k_kp)
        deallocate(phi_ph_c_k_kp_1)
                
        write(*,*) "finish analyze_calculate_vertexladder_2nd_order!!!"

!    Fix q=0 and nu=0, and
!    Gamma_s(k,k',w,w') = Lambda_pp_s(w,w',nu)+(1/2)(Phid1 + Phid2)-(3/2)(Phim1 + Phim2),
!    where
!    Phi(d/m)1(q=0,nu=0)_p,p' = Phi(d/m)1(w'-w,k'-k)_-k',-w',k,w
!    Phi(d/m)2(q=0,nu=0)_p,p' = Phi(d/m)2(w'+w,k'+k)_-k',-w',-k,-w
 
     do i=1,4*ktime**2         ! k = (kx,ky)
     ikx1=(mod(i-1,2*ktime) + 1) - ktime ! kx: (-ktime+1:ktime)
     iky1=((i-1)/(2*ktime) + 1) - ktime  ! ky: (-ktime+1:ktime)
     do j=1,4*ktime**2           ! k' = (kx',ky')
     ikx2=(mod(j-1,2*ktime) + 1) - ktime ! kx': (-ktime+1:ktime)
     iky2=((j-1)/(2*ktime) + 1) - ktime  ! ky': (-ktime+1:ktime)
!!     do k=1,4*ktime**2		! q = (qx,qy)
!!	 k=4*ktime**2		! this is for q=(pi,pi)
!!	 k=2*ktime**2 - ktime	! this is for q=(0,0)
!!      ikx3 = (mod(k-1,2*ktime) + 1) - ktime ! qx : (-ktime+1:ktime)
!!      iky3 = ((k-1)/(2*ktime) + 1) - ktime  ! qy : (-ktime+1:ktime)

       iqx1=add_k_1D(ikx2 - ikx1) ! kx'-kx
       iqy1=add_k_1D(iky2 - iky1) ! ky'-ky
       iqx2=add_k_1D(ikx2 + ikx1) ! kx'+kx
       iqy2=add_k_1D(iky2 + iky1) ! ky'+kx
       
!!      q1= (k-1)*(2*nbfrq-1) + 0 + nbfrq  ! q1 = (q,v=0) (1,2,...,4*ktime**2*(2*nbfrq-1))
       
!!       do inu=-nbfrq+1,nbfrq-1 ! nu = (-7,-6,...,-1,0,1,...,6,7)
         
          do l=-nffrq/2+1,nffrq/2   ! w = (-3,-2,...,3,4)
          p1 = (i-1)*nffrq + l + nffrq/2 ! p1 = (k,w) = (1,2,...,4*ktime**2*nffrq)
          p1m = (4*ktime**2+1-i-1)*nffrq -l+1 + nffrq/2 ! -p1 = (-k,-w)
          do m=-nffrq/2+1,nffrq/2   ! w'= (-3,-2,...,3,4)
          p2 = (j-1)*nffrq + m + nffrq/2 ! p2 = (k',w')=(1,2,...,4*ktime**2*nffrq)
          p2m = (4*ktime**2+1-j-1)*nffrq -m+1 + nffrq/2 ! -p2 = (-k',-w')

!!	construct the particle-particle irreducible vertex
!!          if(m-l.eq.inu) then         ! Phi(d/m)1(w'-w,k'-k)_-k',-w',k,w ! no restriction in momentum
	    if(((m-l).le.(nbfrq-1)).and.((m-l).ge.(-nbfrq+1))) then
	    q1 = (((iqx1+ktime) + (iqy1+ktime-1)*2*ktime)-1)*(2*nbfrq-1) + (m-l) + nbfrq	! q1 = p2-p1 = (k'-k, w'-w)
!!	singlet
!!	    gamat_s(p1,p2)= lambda_m(l,m) - 0.5d0*phi_ph_c_k_kp_q(p2m,p1,q1) + 1.5d0*phi_ph_m_k_kp_q(p2m,p1,q1)
!!	triplet
	    gamat_s(p1,p2)= lambda_m(l,m) - 0.5d0*phi_ph_c_k_kp_q(p2m,p1,q1) - 0.5d0*phi_ph_m_k_kp_q(p2m,p1,q1)
	    
	    else if (((m+l-1).le.(nbfrq-1)).and.((m+l-1).ge.(-nbfrq+1))) then
	    q2 = (((iqx2+ktime) + (iqy2+ktime-1)*2*ktime)-1)*(2*nbfrq-1) + (m+l-1) + nbfrq  ! q2 = p2+p1 = (k'+k, w'+w)
!!	singlet	    
!!          gamat_s(p1,p2)= lambda_m(l,m)- 0.5d0*phi_ph_c_k_kp_q(p2m,p1m,q2) + 1.5d0*phi_ph_m_k_kp_q(p2m,p1m,q2)
!!	triplet
	  gamat_s(p1,p2)= lambda_m(l,m) + 0.5d0*phi_ph_c_k_kp_q(p2m,p1m,q2) + 0.5d0*phi_ph_m_k_kp_q(p2m,p1m,q2)
          
          end if
!!          else if(m+l-1.eq.inu) then  ! Phi(d/m)2(w'+w,k'+k)_-w',-w ! Note that the index of w+w' is (l+m-1)
!!          gamat_s(p1,p2,q1)= lambda_m(l,m) - 0.5d0*phi_ph_c(iqx2,iqy2,inu,-m+1,-l+1) + 1.5d0*phi_ph_m(iqx2,iqy2,inu,-m+1,-l+1)
!!          end if 

          end do ! over m={-nffrq/2+1,nffrq/2} loop
          end do ! over l={-nffrq/2+1,nffrq/2} loop

!!       end do ! over inu={-nbfrq+1,nbfrq-1} loop

     end do ! over j={1, 4*ntime**2} loop
     end do ! over i={1, 4*ntime**2} loop
    
    ! deallocate memory
        deallocate(phi_ph_m_k_kp_q)
        deallocate(phi_ph_c_k_kp_q)
        
        write(*,*) "finish analyze_calculate 2nd order: Gamma_s !!!"
        
	return
	end subroutine analyze_calculate_vertexladder_2nd_order
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function add_k_1D(ik)
      use control, only : ktime
      implicit none
      integer :: ik,tmpk
      integer :: add_k_1D

      add_k_1D=ik
      tmpk = ik
       if(tmpk.le.-ktime) then
         add_k_1D = ik + 2*ktime
       else if(tmpk.gt.ktime) then
         add_k_1D = ik - 2*ktime
       end if

      return
      end function


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     subroutine analyze_pairing_matrix()
	use constants
	use control
	use context

	implicit none

     integer :: p1,p2,jw,j,k,ikx2,iky2,ikx3,iky3,iqx,iqy
     real(dp) :: rDk,kx,ky
! status flag
     integer :: istat
     integer :: add_k_1D

! eigenvalues and eigenvectors
     external  zgeev
     intrinsic min,int
     integer :: i,ieig,il,ind,lwork
     integer :: neigs = 5
     integer, allocatable :: ilead(:)
     real(dp), allocatable :: rwork(:)
     real(dp), allocatable :: Sval(:)
     complex(dp) :: c1
     complex(dp), allocatable :: Eval(:)
     complex(dp), allocatable :: Umat(:,:)
!     complex(dp), allocatable :: Gamat(:,:)
     complex(dp), allocatable :: cwork(:)

!     allocate (Gamat(nt_2p,nt_2p),stat=istat)
!     Gamat = czero

!    Redefine matrix gamat_s(p1,p2) as pairing matrix:
!    gamat_s(p1,p2) = sum_p3 gamat_s(p1,p3)*chi0(p3,p3)

     do p1=1,nffrq*4*ktime**2
        do j=1,4*ktime**2           ! k' = (kx',ky')
        ikx2=(mod(j-1,2*ktime) + 1) - ktime ! kx': (-ktime+1:ktime)
        iky2=((j-1)/(2*ktime) + 1) - ktime  ! ky': (-ktime+1:ktime)

!     	do k=1,4*ktime**2		! q = (qx,qy)
!!	 k=4*ktime**2			! this is for q = (pi,pi)
!!	 k=2*ktime**2 - ktime	! this is for q=(0,0)
	 
!!	 ikx3 = (mod(k-1,2*ktime) + 1) - ktime ! qx : (-ktime+1:ktime)
!!	 iky3 = ((k-1)/(2*ktime) + 1) - ktime  ! qy : (-ktime+1:ktime)
	 
!!	 iqx=add_k_1D(ikx2 + ikx3) ! kx'+qx
!!	 iqy=add_k_1D(iky2 + iky3) ! ky'+qy
	 
          do jw=-nffrq/2+1,nffrq/2   ! w' = (-3,-2,...,3,4)

          p2= (j-1)*nffrq + jw + nffrq/2 ! p2 = (k',w')=(1,2,...,4*ktime**2*nffrq)

!!	this is the calculation for magnetic LEVs
!!          Gamat_m(p1,p2) =- Gamat_sm(p1,p2)*grnf_lattice(iqx,iqy,jw)*grnf_lattice(ikx2,iky2,jw)/(beta*real((2*ktime)**2,dp))
	  Gamat_s(p1,p2) = -Gamat_s(p1,p2)*grnf_lattice(ikx2,iky2,jw)*grnf_lattice(-ikx2+1,-iky2+1,-jw+1)/(beta*real((2*ktime)**2,dp))

!!          (grnf(jw,2) - grnf_lattice(ikx2,iky2,jw))*(grnf_lattice(-ikx2+1,-iky2+1,-jw+1)-grnf(-jw+1,2)) 
!!          (grnf_lattice(ikx2,iky2,jw) - grnf(jw,2))*(grnf_lattice(-ikx2+1,-iky2+1,-jw+1)-grnf(-jw+1,2)) 

!        write(6001,*) p1,p2,Gamat_s(p1,p2)
          end do ! jw

       end do ! over j={1,4*ktime**2} loop
     end do ! over p1={1,nffrq*4*ktime**2} loop


        write(*,*) "Obtained pairing matrix!!!"

!    Calculate eigenvalues of the pairing matrix:


! eigenvelue and eigenvector
         allocate (ilead(nt_2p),stat=istat)
         allocate (rwork(2*nt_2p),stat=istat)
         allocate (Sval(nt_2p),stat=istat)
         allocate (Eval(nt_2p),stat=istat)
         allocate (Umat(nt_2p,nt_2p),stat=istat)
         allocate (cwork(3*nt_2p),stat=istat)

! Initialize matrices
         ilead = 0
         rwork = zero
         Sval = zero
         Eval = czero
         Umat = czero
         cwork = czero
         rDk = 1.0d0/real(ktime,dp)


! check the status
         if ( istat /= 0 ) then
             write(*,*) "istat=",istat
             call analyze_print_error('analyze_pairing_matrix','can not allocate enough memory')
         endif
        write(*,*) "Diagonizing the pairing matrix....."
! Calculate the eigenvalues and vectors of Gamat_s: 

          lwork = 3*nt_2p


          call zgeev('N','V',nt_2p,Gamat_m,nt_2p,Eval,Umat,nt_2p, &
                     Umat,nt_2p,cwork,lwork,rwork,istat)


          if( istat /= 0 ) then
            write(6,*) 'zgeev error: istat = ',istat
            stop
          endif
! Sort the eigenvalues of Gamat_s:   
          do i=1,nt_2p ! the smallest Sval is the lead

            Sval(i)=(one-real(Eval(i)))**2+(aimag(Eval(i)))**2
          end do

          do ieig=1,nt_2p
            ilead(ieig)=1
            do i=2,nt_2p
               if (Sval(i).lt.Sval(ilead(ieig))) ilead(ieig)=i
            end do
            Sval(ilead(ieig))=1.0E8
          end do

          do i=1,nt_2p ! now the largest Sval is the leading eigenvalues
            if(real(Eval(i)).lt.one) then
              Sval(i)=one-sqrt((one-real(Eval(i)))**2 +(aimag(Eval(i)))**2)
            else
              Sval(i)=one+sqrt((one-real(Eval(i)))**2 +(aimag(Eval(i)))**2)
            end if            
          end do

! Output the leading eigenvalues and eigenvectors:
	open(mytmp, file='analyze.eig.dat', form='formatted', status='unknown')
         do ieig=1,neigs
            il=ilead(ieig)
!            eigs(itype,isus,ieig)=Sval(il)
            write(mytmp,*) ' '
            write(mytmp,"('# |eigenvalue|',i2,' =',f8.4)") ieig,Sval(il)
            write(mytmp,"('# eigenvalue',i2,' =',2f8.4)")  ieig,Eval(il)
            write(mytmp,"('#  Kx      Ky   eigvecR(K,n=0)')") 
            do ikx2=-ktime+1,ktime
            do iky2=-ktime+1,ktime
		kx = rDk*pi*(real(ikx2,dp)-half)
		ky = rDk*pi*(real(iky2,dp)-half)
!              ind=ick + (nwn_2p)*Nc  ! n=0 => nwn_2p+1
!   Take values only at w = w' = pi*T:
              ind = (ikx2+ktime-1)*(2*ktime)*nffrq + (iky2+ktime-1)*nffrq + nffrq/2 +1
              c1=Umat(ind,il)
              write(mytmp,"(4(f16.8,2x))") kx,ky,c1
           
            end do ! over iky2={-ktime+1,ktime} loop             
            end do ! over ikx2={-ktime+1,ktime} loop


           end do ! over ieig={1,neigs} loop
        close(mytmp)

! Deallocate matrices
	 if ( allocated(ilead))    deallocate(ilead)
	 if ( allocated(rwork))    deallocate(rwork)
	 if ( allocated(Sval))    deallocate(Sval)
	 if ( allocated(Eval))    deallocate(Eval)
	 if ( allocated(Umat))    deallocate(Umat)
!	 if ( allocated(Gamat))    deallocate(Gamat)
	 if ( allocated(cwork))    deallocate(cwork)

     return
     end subroutine analyze_pairing_matrix
