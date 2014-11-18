!!!-----------------------------------------------------------------------
!!! project : hibiscus/stoch
!!! program : sac_make_grid
!!!           sac_make_fphi
!!!           sac_make_const
!!!           sac_make_gauss
!!!           sac_make_alpha
!!!           sac_make_rgamm
!!!           sac_make_delta
!!!           sac_make_ppleg
!!!           sac_warp_image
!!!           sac_make_normal
!!!           sac_make_hamil0
!!!           sac_make_hamil1
!!!           sac_make_kernel
!!! source  : sac_toolbox.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 10/01/2008 by li huang
!!!           12/13/2011 by li huang
!!! purpose : to provide utility functions and subroutines for stochastic
!!!           analytic continuation code
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> sac_make_grid: build wgrid and xgrid, very dense grids
  subroutine sac_make_grid(wgrid, xgrid)
     use constants, only : dp, zero, one, two, half

     use control, only : nwmax, ngrid
     use control, only : sigma, wstep

     implicit none

! external arguments
! very dense grid on [ -wstep*nwmax, +wstep*nwmax ]
     real(dp), intent(out) :: wgrid(ngrid)

! very dense grid on [0, 1]
     real(dp), intent(out) :: xgrid(ngrid)

! local variables
! loop index
     integer  :: i

! prefactor for gaussian type model
     real(dp) :: f

! default model on wgrid
     real(dp) :: model(ngrid)

! build wgrid
     call s_linspace_d( -wstep*nwmax, +wstep*nwmax, ngrid, wgrid )

! build default model according to sigma
! case 1: gaussian model
     if ( sigma > zero ) then
         f = half / ( sigma * sigma )
         do i=1,ngrid
             model(i) = exp( - f * wgrid(i) * wgrid(i) )
         enddo ! over i={1,ngrid} loop
! case 2: flat model
     else
         model = one
     endif ! back if ( sigma > zero ) block

! normalization default model
     call sac_make_normal( ngrid, one, model )

! build xgrid
     call s_cumsum_d( ngrid, model, xgrid )

     return
  end subroutine sac_make_grid

!!>>> sac_make_fphi: calculate \phi(\omega), it is belong to [0, 1]
! please refer to equation (16)
  subroutine sac_make_fphi(model, fphi)
     use constants, only : dp

     use control, only : nwmax
     use control, only : wstep

     implicit none

! external arguments
! default model
     real(dp), intent(in) :: model(-nwmax:nwmax)

! \phi(\omega)
     real(dp), intent(out) :: fphi(-nwmax:nwmax)

     call s_cumsum_d( 2*nwmax+1, model, fphi )
     fphi = fphi * wstep

     return
  end subroutine sac_make_fphi

!!>>> sac_make_const: build default model (flat model)
  subroutine sac_make_const(model)
     use constants, only : dp, one

     use control, only : nwmax
     use control, only : wstep

     implicit none

! external arguments
! default model
     real(dp), intent(out) :: model(-nwmax:nwmax)

! build model
     model = one

! enforce it to be normalized
     call sac_make_normal( 2*nwmax+1, wstep, model )

     return
  end subroutine sac_make_const

!!>>> sac_make_gauss: build default model: gaussian
  subroutine sac_make_gauss(model)
     use constants, only : dp, half

     use control, only : nwmax
     use control, only : sigma, wstep
     use context, only : wmesh

     implicit none

! external arguments
! default model
     real(dp), intent(out) :: model(-nwmax:nwmax)

! local variables
! loop index
     integer  :: i

! prefactor
     real(dp) :: f

! build model
     f = half / ( sigma * sigma )
     do i=-nwmax,nwmax
         model(i) = exp( -f * wmesh(i) * wmesh(i) )
     enddo ! over i={-nwmax,nwmax} loop

! enforce it to be normalized
     call sac_make_normal( 2*nwmax+1, wstep, model )

     return
  end subroutine sac_make_gauss

!!>>> sac_make_alpha: build alpha parameters list
  subroutine sac_make_alpha(alpha)
     use constants, only : dp

     use control, only : nalph
     use control, only : ainit, ratio

     implicit none

! external arguments
! alpha parameters list
     real(dp), intent(out) :: alpha(nalph)

! local variables
! loop index
     integer :: i

! alpha_{p+1} / alpha_p = ratio = two (default)
     do i=1,nalph
         alpha(i) = ainit * ( ratio**( i - 1 ) )
     enddo ! over i={1,nalph} loop

     return
  end subroutine sac_make_alpha

!!>>> sac_make_rgamm: build rgamm and igamm, randomly
  subroutine sac_make_rgamm(igamm, rgamm)
     use constants, only : dp, one
     use spring, only : spring_sfmt_stream

     use control, only : ngrid, ngamm, nalph

     implicit none

! external arguments
! configurations of a_{\gamma}
     integer, intent(out)  :: igamm(nalph,ngamm)

! configurations of r_{\gamma}
     real(dp), intent(out) :: rgamm(nalph,ngamm)

! local variables
! loop index
     integer :: i
     integer :: j

! loop over alpha
     do j=1,nalph

! generate rgamm and igamm using random numbers
         do i=1,ngamm
             rgamm(j,i) = spring_sfmt_stream()
             igamm(j,i) = ceiling( spring_sfmt_stream() * ngrid )
         enddo ! over i={1,ngamm} loop

! enforce \sum_{\gamma} r_{\gamma} = 1, it is necessary
         call sac_make_normal( ngamm, one, rgamm(j,:) )

     enddo ! over j={1,nalph} loop

     return
  end subroutine sac_make_rgamm

!!>>> sac_make_delta: build the delta function in advance
  subroutine sac_make_delta(xgrid, F_phi, delta)
     use constants, only : dp

     use control, only : nwmax, ngrid
     use control, only : eta1, eta2

     implicit none

! external arguments
! very dense grid on [0,1]
     real(dp), intent(in)  :: xgrid(ngrid)

! \phi(\omega) function
     real(dp), intent(in)  :: F_phi(-nwmax:nwmax)

! precalculated delta function
     real(dp), intent(out) :: delta(-nwmax:nwmax,ngrid)

! local variables
! loop index
     integer  :: i
     integer  :: j

! real(dp) dummy variables
     real(dp) :: s

     do j=-nwmax,nwmax
         do i=1,ngrid
             s = F_phi(j) - xgrid(i)
             delta(j,i) = eta1 / ( s * s + eta2 )
         enddo ! over i={1,ngrid} loop
     enddo ! over j={-nwmax,nwmax} loop

     return
  end subroutine sac_make_delta

!!>>> sac_warp_image: convert the image function from legendre polynomial
!!>>> representation to normal representation
  subroutine sac_warp_image(image_l, image_t)
     use constants, only : dp, zero, one, two

     use control, only : nwmax, nalph
     use control, only : lemax, legrd
     use context, only : ppleg, wmesh

     implicit none

! external arguments
! image function in legendre polynomial representation
     real(dp), intent(in)  :: image_l(nalph,lemax)

! image function in normal representation
     real(dp), intent(out) :: image_t(nalph,-nwmax:nwmax)

! local variables
! loop index
     integer  :: i
     integer  :: j

! loop index over legendre polynomial
     integer  :: fleg

! index for x in [0,2]
     integer  :: curr

! interval for x
     real(dp) :: step

! dummy variables
     real(dp) :: raux

     image_t = zero

     step = real(legrd - 1) / two
     do i=1,nalph
         do j=-nwmax,nwmax
             raux = wmesh(j) / maxval(wmesh) + one
             curr = nint(raux * step) + 1
             do fleg=1,lemax
                 raux = sqrt(two * fleg - 1)
                 image_t(i,j) = image_t(i,j) + raux * image_l(i,fleg) * ppleg(curr,fleg)
             enddo ! over fleg={1,lemax} loop
         enddo ! over j={-nwmax,nwmax} loop
     enddo ! over i={1,nalph} loop

     return
  end subroutine sac_warp_image

!!>>> sac_make_normal: normalize a given function
  subroutine sac_make_normal(ndim, weight, vector)
     use constants, only : dp

     implicit none

! external arguments
! size of given function
     integer, intent(in)     :: ndim

! weight of the given function
     real(dp), intent(in)    :: weight

! the function
     real(dp), intent(inout) :: vector(ndim)

! local variables
! norm of given function
     real(dp) :: norm

     norm = sum(vector) * weight
     vector = vector / norm

     return
  end subroutine sac_make_normal

!>>> calculate hamiltonian, h_C(\tau)
! please refer to equation (39)
  subroutine sac_make_hamil0(rg, ig, hc)
     use constants
     use control
     use context

     implicit none

! external arguments
! a_{\gamma}
     integer ,intent(in)  :: ig(ngamm)

! r_{\gamma}
     real(dp),intent(in)  :: rg(ngamm)

! h_C(\tau)
     real(dp),intent(out) :: hc(ntime)

! local variables
! loop index
     integer :: i
     integer :: j

! for a given configuration, we calculate h_C
     hc = zero
     do i=1,ntime
         do j=1,ngamm
             hc(i) = hc(i) + rg(j) * fkern( i, ig(j) )
         enddo ! over j={1,ngamm} loop
         hc(i) = hc(i) / G_dev(i) - G_tau(i)
     enddo ! over i={1,ntime} loop

     return
  end subroutine sac_make_hamil0

!>>> calculate hamiltonian: H_C
! please refer to equation (35)
  subroutine sac_make_hamil1(hc, hh)
     use constants
     use control
     use context

     implicit none

! external arguments
! h_C term
     real(dp), intent(in ) :: hc(ntime)

! H_C term
     real(dp), intent(out) :: hh

! local variables
! loop index
     integer :: i

! H_C = \int_0^{\beta} ( h_C(\tau) )^2
     hh = zero
     do i=1,ntime
         hh = hh + hc(i) * hc(i)
     enddo ! over i={1,ntime} loop
     hh = hh * ( tmesh(2) - tmesh(1) )

     return
  end subroutine sac_make_hamil1

!>>> calculate kernel function
! please refer to equation (7)
  subroutine sac_make_kernel(kern)
     use constants, only : dp, one, zero
     use control, only : beta, ntime, ngrid
     use context, only : tmesh, wgrid

     implicit none

! external arguments
! kernel function
     real(dp), intent(out) :: kern(ntime,ngrid)

! local variables
! loop index
     integer  :: i
     integer  :: j

! imaginary time point
     real(dp) :: tau

! frequency point
     real(dp) :: omega

     do j=1,ngrid
         omega = wgrid(j)
         do i=1,ntime
             tau = tmesh(i)
             if ( omega >= zero )then
                 kern(i,j) = exp(        - tau   * omega ) / ( one + exp( - beta * omega ) )
             else
                 kern(i,j) = exp( ( beta - tau ) * omega ) / ( one + exp(   beta * omega ) )
             endif
         enddo ! over i={1,ntime} loop
     enddo ! over j={1,ngrid} loop

     return
  end subroutine sac_make_kernel
