!!!-----------------------------------------------------------------------
!!! project : azalea
!!! program : dt_setup_param
!!!           dt_setup_model <<<---
!!!           dt_input_mesh_
!!!           dt_input_dmft_
!!!           dt_input_latt_
!!!           dt_input_dual_
!!!           dt_input_vert_ <<<---
!!!           dt_alloc_array
!!!           dt_reset_array
!!!           dt_final_array <<<---
!!! source  : dt_stream.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 09/16/2009 by li huang (created)
!!!           01/02/2018 by li huang (last modified)
!!! purpose :
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!
!! @sub dt_setup_param
!!
!!
!!
  subroutine dt_setup_param()
     use constants, only : dp

     use control

     implicit none

! only for debug
     nffrq = 16
     nbfrq = 7

     part  = 1.0_dp 
     beta  = 1.0_dp

     return
  end subroutine dt_setup_param

!!
!! @sub dt_setup_model
!!
!!
!!
  subroutine dt_setup_model()
     implicit none

     call dt_input_mesh_()
     call dt_input_dmft_()
     call dt_input_latt_()
     call dt_input_dual_()
     call dt_input_vert_()

     return
  end subroutine dt_setup_model

!!
!! @sub dt_input_mesh_
!!
!!
!!
  subroutine dt_input_mesh_()
     use constants, only : dp
     use constants, only : one, two, pi

     use control, only : nffrq, nbfrq
     use control, only : nkpts, nkp_x, nkp_y, nkp_z
     use control, only : beta
     use control, only : part
     use context, only : kx, ky, kz, ek
     use context, only : fmesh, bmesh

     implicit none

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k

     do i=1,nkp_x
         kx(i) = (two * pi) * float(i - 1)/ float(nkp_x)
     enddo ! over i={1,nkp_x} loop

     do i=1,nkp_y
         ky(i) = (two * pi) * float(i - 1)/ float(nkp_y)
     enddo ! over i={1,nkp_y} loop

     do i=1,nkp_z
         kz(i) = (two * pi) * float(i - 1)/ float(nkp_z)
     enddo ! over i={1,nkp_z} loop

! build a 2d lattice
     k = 0
     do i=1,nkp_x
         do j=1,nkp_y
             k = k + 1
             ek(k) = -two * part * ( cos( kx(i) ) + cos( ky(j) ) )
         enddo ! over j={1,nkp_y} loop
     enddo ! over i={1,nkp_x} loop     
     call s_assert(k == nkpts)

     do i=1,nffrq
         fmesh(i) = (two * i - one - nffrq) * pi / beta
     enddo ! over i={1,nffrq} loop

     do i=1,nbfrq
         bmesh(i) = (two * i - one - nbfrq) * pi / beta
     enddo ! over i={1,nbfrq} loop

     return
  end subroutine dt_input_mesh_

!!
!! @sub dt_input_dmft_
!!
!!
!!
  subroutine dt_input_dmft_()
     use constants, only : dp
     use constants, only : mytmp

     use control, only : nffrq, nbfrq
     use context, only : dmft_g, dmft_h
     use context, only : vert_d, vert_m

     implicit none

! local variables
! loop index
     integer  :: i
     integer  :: if1, if2
     real(dp) :: r1, r2
     real(dp) :: c1, c2
     real(dp) :: d1, d2
     real(dp) :: v1, v2

! read in impurity green's function
     open(mytmp, file = 'df.dmft_g.in', form = 'formatted', status = 'unknown')
     do i=1,nffrq
         read(mytmp,*) r1, r2, c1, c2
         dmft_g(i,1) = dcmplx(c1, c2)
         dmft_g(i,2) = dcmplx(c1, c2)
     enddo ! over i={1,nffrq} loop
     close(mytmp)

! read in hybridization function
     open(mytmp, file = 'df.dmft_h.in', form = 'formatted', status = 'unknown')
     do i=1,nffrq
         read(mytmp,*) r1, r2, c1, c2
         dmft_h(i,1) = dcmplx(c1, c2)
         dmft_h(i,2) = dcmplx(c1, c2)
     enddo ! over i={1,nffrq} loop
     close(mytmp)

! read in vertex function, density channel
     open(mytmp, file = 'df.vert_d.in', form = 'formatted', status = 'unknown')
     do i=1,nbfrq
         do if1=1,nffrq
             do if2=1,nffrq
                 read(mytmp,*) r1, r2, c1, c2, d1, d2, v1, v2
                 vert_d(if2,if1,i) = dcmplx(v1, v2)
             enddo ! over if2={1,nffrq} loop
             read(mytmp,*) ! skip one line
         enddo ! over if1={1,nffrq} loop
     enddo ! over i={1,nbfrq} loop
     close(mytmp)

! read in vertex function, magentic channel
     open(mytmp, file = 'df.vert_m.in', form = 'formatted', status = 'unknown')
     do i=1,nbfrq
         do if1=1,nffrq
             do if2=1,nffrq
                 read(mytmp,*) r1, r2, c1, c2, d1, d2, v1, v2
                 vert_m(if2,if1,i) = dcmplx(v1, v2)
             enddo
             read(mytmp,*) ! skip one line
         enddo
     enddo
     close(mytmp)

     return
  end subroutine dt_input_dmft_

!!
!! @sub dt_input_latt_
!!
!!
!!
  subroutine dt_input_latt_()
     use constants, only : dp
     use constants, only : one

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts
     use context, only : ek
     use context, only : dmft_g, dmft_h
     use context, only : latt_g

     implicit none

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k

     do k=1,nkpts
         do j=1,norbs
             do i=1,nffrq
                 latt_g(i,j,k) = one / ( one / dmft_g(i,j) + dmft_h(i,j) - ek(k) ) 
             enddo ! over i={1,nffrq} loop
         enddo ! over j={1,norbs} loop
     enddo ! over k={1,nkpts} loop

     return
  end subroutine dt_input_latt_

!!
!! @sub dt_input_dual_
!!
!!
!!
  subroutine dt_input_dual_()
     use constants, only : dp
     use constants, only : czero

     use control, only : norbs
     use control, only : nffrq
     use control, only : nkpts
     use context, only : dmft_g
     use context, only : latt_g
     use context, only : dual_g, dual_s, dual_b

     implicit none

! local variables
     integer :: i
     integer :: j
     integer :: k

     do k=1,nkpts
         do j=1,norbs
             do i=1,nffrq
                 dual_b(i,j,k) = latt_g(i,j,k) - dmft_g(i,j) 
                 dual_g(i,j,k) = dual_b(i,j,k)
                 dual_s(i,j,k) = czero
             enddo ! over i={1,nffrq} loop
         enddo ! over j={1,norbs} loop
     enddo ! over k={1,nkpts} loop

     return
  end subroutine dt_input_dual_

!!
!! @sub dt_input_vert_
!!
!!
!!
  subroutine dt_input_vert_()
     implicit none

     return
  end subroutine dt_input_vert_

!!
!! @sub dt_alloc_array
!!
!!
!!
  subroutine dt_alloc_array()
     use context ! ALL

     implicit none

! allocate memory for context module
     call cat_alloc_mesh()
     call cat_alloc_dmft()
     call cat_alloc_dual()
     call cat_alloc_latt()
     call cat_alloc_vert()
     
     return
  end subroutine dt_alloc_array

!!
!! @sub dt_reset_array
!!
!!
!!
  subroutine dt_reset_array()
     implicit none

     return
  end subroutine dt_reset_array

!!
!! @sub dt_final_array
!!
!!
!!
  subroutine dt_final_array()
     use context ! ALL

     implicit none

! deallocate memory for context module
     call cat_free_mesh()
     call cat_free_dmft()
     call cat_free_dual()
     call cat_free_latt()
     call cat_free_vert()

     return
  end subroutine dt_final_array
