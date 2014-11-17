!!!-----------------------------------------------------------------------
!!! project : hibiscus
!!! program : entropy_make_smooth
!!!           entropy_make_normal
!!!           entropy_make_model
!!!           entropy_make_fnorm
!!!           entropy_make_srule
!!!           entropy_make_sterm
!!!           entropy_make_chihc
!!!           entropy_make_trace
!!!           entropy_make_akern
!!!           entropy_make_ckern
!!!           entropy_make_fkern
!!! source  : entropy_toolbox.f90
!!! type    : subroutines
!!! author  : li huang (email:huangli712@gmail.com)
!!! history : 10/01/2008 by li huang
!!!           01/26/2011 by li huang
!!!           11/17/2014 by li huang
!!! purpose : to provide utility functions and subroutines for the classic
!!!           maximum entropy method code
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> entropy_make_smooth: to smooth the image function. its principle is
!!>>> very simple. the value of every points in the curve is equal to the
!!>>> avarage value of 2[naver]-near-neighbors points.
  subroutine entropy_make_smooth(naver, image)
     use constants, only : dp, zero

     use control, only : nwmax

     implicit none

! external arguments
! number of neighbor points
     integer, intent(in) :: naver

! specctrum function
     real(dp), intent(inout) :: image(-nwmax:nwmax)

! local variables
! loop index
     integer  :: i
     integer  :: j

! array index, used to bracket the smoothing zone
     integer  :: i1
     integer  :: i2

! smoothed spectrum function
     real(dp) :: image_s(-nwmax:nwmax)

! initialize it
     image_s = zero

     do i=-nwmax,nwmax
! bracket the smoothing zone
         i1 = i - naver
         if ( i1 < -nwmax ) then
             i1 = -nwmax
         endif ! back if ( i1 < -nwmax ) block

         i2 = i + naver
         if ( i2 > +nwmax ) then
             i2 = nwmax
         endif ! back if ( i2 > +nwmax ) block

! smooth it in the smoothing zone
         do j=i1,i2
            image_s(i) = image_s(i) + image(j)
         enddo ! over j={i1,i2} loop
         image_s(i) = image_s(i) / real( i2 - i1 + 1 )
     enddo ! over i={-nwmax,nwmax} loop

! copy image1 to image
     image = image_s

     return
  end subroutine entropy_make_smooth

!!>>> entropy_make_normal: used to perform normalization on image function
  subroutine entropy_make_normal(npara, fnorm, image)
     use constants, only : dp

     use control, only : nwmax

     implicit none

! external arguments
! normalized parameter
     real(dp), intent(in) :: npara

! normalized function
     real(dp), intent(in) :: fnorm(-nwmax:nwmax)

! spectrum function
     real(dp), intent(inout) :: image(-nwmax:nwmax)

! local variables
! normalized factor
     real(dp) :: f

     f = npara / dot_product(fnorm, image)
     image = f * image

     return
  end subroutine entropy_make_normal

!!>>> entropy_make_model: to build the default model, here for simplicity,
!!>>> we only implement the flat model---model(\omega) = constant, and
!!>>> gaussian model
  subroutine entropy_make_model(wmesh, model)
     use constants, only : dp, one, two, pi

     use control, only : nwmax, ntype
     use control, only : sigma

     implicit none

! external arguments
! real frequency grid
     real(dp), intent(in)  :: wmesh(-nwmax:nwmax)

! reference model function
     real(dp), intent(out) :: model(-nwmax:nwmax)

! local variables
! loop index
     integer :: i

! gaussian model
! please refer to equation (5.9) in the reference
     if ( ntype == 0 ) then
         do i=-nwmax,nwmax
             model(i) = one / ( sigma * sqrt(pi) ) * exp( -( wmesh(i) / sigma )**two )
         enddo ! over i={-nwmax,nwmax} loop
! flat model
     else
         do i=-nwmax,nwmax
             model(i) = one
         enddo ! over i={-nwmax,nwmax} loop
     endif ! back if ( ntype == 0 ) block

     return
  end subroutine entropy_make_model

!!>>> entropy_make_fnorm: to calculate f-norms, which is necessary in the
!!>>> calculation of sum rules for spectrum function
  subroutine entropy_make_fnorm(wmesh, fnorm)
     use constants, only : dp

     use control, only : nwmax
     use control, only : wstep

     implicit none

! external arguments
! real frequency grid
     real(dp), intent(in)  :: wmesh(-nwmax:nwmax)

! normalization function
     real(dp), intent(out) :: fnorm(-nwmax:nwmax,3)

! local variables
! loop index
     integer :: i

     do i=-nwmax,nwmax
         fnorm(i,1) = wstep
         fnorm(i,2) = wstep * wmesh(i)
         fnorm(i,3) = wstep * wmesh(i) * wmesh(i)
     enddo ! over i={-nwmax,nwmax} loop

     return
  end subroutine entropy_make_fnorm

!>>> to calculate the following integrations
!    s0 = \int A dw
!    s1 = \int A w dw
!    s2 = \int A w^{2} dw
! here A means the spectrum function, and w means frequency grid
  subroutine entropy_make_srule(fnorm, image, srule)
     use constants
     use control

     implicit none

! external arguments
! sum rules values
     real(dp), intent(out) :: srule(3)

! normalization function
     real(dp), intent(in)  :: fnorm(-nwmax:nwmax,3)

! spectrum function
     real(dp), intent(in)  :: image(-nwmax:nwmax)

! local variables
! loop index
     integer :: i

     srule = zero
     do i=-nwmax,nwmax
         srule(1) = srule(1) + image(i) * fnorm(i,1)
         srule(2) = srule(2) + image(i) * fnorm(i,2)
         srule(3) = srule(3) + image(i) * fnorm(i,3)
     enddo ! over i={-nwmax,nwmax} loop

     return
  end subroutine entropy_make_srule

!>>> to calculate the entropy term S
!    S = \int dw (A(\omage)-m(\omega)-A(\omage)*ln [A(\omega)/m(\omega)])
  subroutine entropy_make_sterm(sterm, image, model)
     use constants
     use control

     implicit none

! external arguments
! entropy term
     real(dp), intent(out) :: sterm

! spectrum function
     real(dp), intent(in)  :: image(-nwmax:nwmax)

! reference model
     real(dp), intent(in)  :: model(-nwmax:nwmax)

! local variables
! loop index
     integer :: i

! please refer to equations (3.18) and (3.17) in the reference
     sterm = zero
     do i=-nwmax,nwmax
         if ( image(i) > epss .and. model(i) > epss ) then
             sterm = sterm + ( image(i) - model(i) - image(i) * log( image(i) / model(i) ) ) * wstep
         endif
     enddo ! over i={-nwmax,nwmax} loop

     return
  end subroutine entropy_make_sterm

!>>> to calculate \chi^{2}
!    \chi^{2} = \sum_{l=1}^{L} (\frac{ G_{l}-\sum_{j} K_{lj}A_{j} }{ \sigma_{l} })^{2}
  subroutine entropy_make_chihc(chi2, akern, G_qmc, G_dev)
     use constants
     use control

     implicit none

! external arguments
! denote as \chi^{2}
     real(dp), intent(out) :: chi2

! \sum_{j} K_{lj} . A_{j} term
     real(dp), intent(in)  :: akern(ntime)

! imaginary time green's function
     real(dp), intent(in)  :: G_qmc(ntime)

! error bar data
     real(dp), intent(in)  :: G_dev(ntime)

! local variables
! loop index
     integer :: i

! please refer to equation (5.4) in the reference
     chi2 = zero
     do i=1,ntime
         chi2 = chi2 + G_dev(i) * ( G_qmc(i) - akern(i) )**2
     enddo ! over i={1,ntime} loop

     return
  end subroutine entropy_make_chihc

!>>> to calculate the trace:
!    Tr \Lambda [ \Lambda + \alpha *I ]^{-1}
! it is used to solve the classic maximum entropy method equation
  subroutine entropy_make_trace(trace, alpha, image, ckern)
     use constants
     use control

     implicit none

! external arguments
! see the code, trace value
     real(dp), intent(out) :: trace

! alpha factor
     real(dp), intent(in)  :: alpha

! see entropy_make_ckern() subroutine
! \frac { \partial^{2} L } { \partial A_i \partial A_j } term
     real(dp), intent(in)  :: ckern(2*nwmax+1,2*nwmax+1)

! spectrum function
     real(dp), intent(in)  :: image(-nwmax:nwmax)

! local variables
! loop index
     integer  :: i
     integer  :: j

! index for frequency grid
     integer  :: iw
     integer  :: jw

! \Lambda matrix
     real(dp) :: lamb1(2*nwmax+1,2*nwmax+1)

! inversion of (\Lambda + \alpha I) matrix
     real(dp) :: lamb2(2*nwmax+1,2*nwmax+1)

! calculate lamb1 and lamb2
! lamb1 = \Lambda
! lamb2 = \Lambda + \alpha * I
! \Lambda = \sqrt{A_{i}} ckern_{ij} \sqrt{A_{j}}
! please refer to equation (4.8) in the reference
     do j=1,2*nwmax+1
         do i=1,2*nwmax+1
             iw = i - nwmax - 1
             jw = j - nwmax - 1
             lamb1(i,j) = sqrt( image(iw) * image(jw) ) * ckern(i,j) * wstep
             lamb2(i,j) = lamb1(i,j)
         enddo ! over i={1,2*nwmax+1} loop
! alpha only add up to diagonal element
         lamb2(j,j) = lamb1(j,j) + alpha
     enddo ! over j={1,2*nwmax+1} loop

! calculate lamb2^{-1}
     call s_inv_d( 2*nwmax+1, lamb2(1:2*nwmax+1,1:2*nwmax+1) )

! calculate trace = Tr \Lambda [ \Lambda + \alpha * I ]^{-1}
! please refer to equation (4.28) in the reference
     trace = zero
     do i=1,2*nwmax+1
         do j=1,2*nwmax+1
             trace = trace + lamb1(i,j) * lamb2(j,i)
         enddo ! over j={1,2*nwmax+1} loop
     enddo ! over i={1,2*nwmax+1} loop

     return
  end subroutine entropy_make_trace

!>>> to calculate akern, an important immediate variable
!    akern_{l} = \sum_{j} K_{lj} A_{j}
  subroutine entropy_make_akern(akern, image, fkern)
     use constants
     use control

     implicit none

! external arguments
! akern = A . K
     real(dp), intent(out) :: akern(ntime)

! spectrum function
     real(dp), intent(in)  :: image(-nwmax:nwmax)

! kernel function
     real(dp), intent(in)  :: fkern(-nwmax:nwmax,ntime)

! local variables
! loop index for time slice
     integer :: i

! loop index for frequency grid
     integer :: j

! please refer to equation (5.4) in the reference
     akern = zero
     do i=1,ntime
         do j=-nwmax,nwmax
             akern(i) = akern(i) + image(j) * fkern(j,i)
         enddo ! over j={-nwmax,nwmax} loop
     enddo ! over i={1,ntime} loop

     return
  end subroutine entropy_make_akern

!>>> to calculate the following quantity:
!    ckern = \frac{ \partial ^{2} L }{ \partial A_{i} \partial A_{j} }
!          = [ K^{T} . C^{-1} K ]_{ij}
!          = \sum_{kl} K_{ki} [ C^{-1} ]_{kl} K_{lj}
! here C means covariance matrix, and it is a diagonal matrix if the
! measurements at different values of \tau are uncorrelated. K is kernel
! and L means likelihood function. A means spectrum function.
  subroutine entropy_make_ckern(G_dev, fkern, ckern)
     use constants
     use control

     implicit none

! external arguments
! error bar data
     real(dp), intent(in)  :: G_dev(ntime)

! kernel function
     real(dp), intent(in)  :: fkern(-nwmax:nwmax,ntime)

! contain the results
     real(dp), intent(out) :: ckern(2*nwmax+1,2*nwmax+1)

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: k

! index for frequency grid
     integer :: iw
     integer :: jw

! please refer to equation (4.9) in the reference
     ckern = zero
     do i=1,2*nwmax+1
         do j=1,2*nwmax+1
             iw = i - nwmax - 1
             jw = j - nwmax - 1
             do k=1,ntime
                 ckern(j,i) = ckern(j,i) + fkern(iw,k) * G_dev(k) * fkern(jw,k) / ( wstep**2 )
             enddo ! over k={1,ntime} loop
         enddo ! over j={1,2*nwmax+1} loop
     enddo ! over i={1,2*nwmax+1} loop

     return
  end subroutine entropy_make_ckern

!>>> to calculate fermion kernel function
!    fkern = \frac{ \exp{-\tau\omega} }{ 1.0 + \exp{-\beta\omega} }
  subroutine entropy_make_fkern(tmesh, wmesh, fkern)
     use constants
     use control

     implicit none

! external arguments
! imaginary time slice
     real(dp), intent(in)  :: tmesh(ntime)

! real frequency grid
     real(dp), intent(in)  :: wmesh(-nwmax:nwmax)

! kernel function
     real(dp), intent(out) :: fkern(-nwmax:nwmax,ntime)

! local variables
! loop index for time slice
     integer :: i

! loop index for frequency grid
     integer :: j

! please refer to equation (2.11) in the reference
     do i=1,ntime
         do j=-nwmax,nwmax
             if ( wmesh(j) >= zero ) then
                 fkern(j,i) = wstep * exp(        - tmesh(i)   * wmesh(j) ) / ( one + exp( -beta * wmesh(j) ) )
             else
                 fkern(j,i) = wstep * exp( ( beta - tmesh(i) ) * wmesh(j) ) / ( one + exp(  beta * wmesh(j) ) )
             endif
         enddo ! over j={-nwmax,nwmax} loop
     enddo ! over i={1,ntime} loop

     return
  end subroutine entropy_make_fkern
