!-------------------------------------------------------------------------
! project : daisy
! program : hfqmc_fourier_t2w
!           hfqmc_fourier_w2t
!           hfqmc_fourier_nderiv
!           hfqmc_fourier_forward
!           hfqmc_fourier_backward
! source  : hfqmc_fourier.f90
! type    : subroutine
! author  : li huang (email:huangli712@yahoo.com.cn)
! history : 12/28/2005 by li huang
!           05/15/2007 by li huang
!           10/28/2008 by li huang
!           12/20/2008 by li huang
!           01/04/2009 by li huang
!           04/18/2009 by li huang
!           08/10/2009 by li huang
!           08/23/2009 by li huang
!           12/24/2009 by li huang
!           02/26/2010 by li huang
!           03/26/2010 by li huang
! purpose : perform forward fourier transformation (\tau to \omega) and
!           backward fourier transformation (\omega to \tau))
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> wrapper subroutines for fourier transformation from imaginary time
! space to matsubara frequency space
  subroutine hfqmc_fourier_t2w(grnt, grnw)
     use constants
     use control

     implicit none

! external arguments
! green's function in imaginary time space
     real(dp), intent(in) :: grnt(ntime,norbs)

! green's function in matsubara frequency space
     complex(dp), intent(out) :: grnw(mfreq,norbs)

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k

! first order derivate of G(0), used by fourier subroutine
     real(dp) :: deriv1

! first order derivate of G(\beta), used by fourier subroutine
     real(dp) :: deriv2

! dummmy arrays
     real(dp) :: raux(ntime)
     complex(dp) :: caux(mfreq)

     QMC_FOURIER_FORWARD: do i=1,norbs

! copy tau data at first, grnt ---> raux
         raux = zero
         do j=1,ntime
             raux(j) = grnt(j,i)
         enddo ! over j={1,ntime} loop

! and then fourier it, raux ---> caux
         caux = czero
         call hfqmc_fourier_nderiv(raux, deriv1, deriv2)
         call hfqmc_fourier_forward(raux, caux, deriv1, deriv2)

! copy omega data at last, caux ---> grnw
         do k=1,mfreq
             grnw(k,i) = caux(k)
         enddo ! over k={1,mfreq} loop

     enddo QMC_FOURIER_FORWARD ! over i={1,norbs} loop

     return
  end subroutine hfqmc_fourier_t2w

!>>> wrapper subroutines for fourier transformation from matsubara frequency
! space to imaginary time space
  subroutine hfqmc_fourier_w2t(grnw, grnt)
     use constants
     use control

     implicit none

! external arguments
! green's function in imaginary time space
     real(dp), intent(out) :: grnt(ntime,norbs)

! green's function in matsubara frequency space
     complex(dp), intent(in) :: grnw(mfreq,norbs)

! local variables
! loop index
     integer  :: i
     integer  :: j
     integer  :: k

! dummy arrays
     real(dp) :: raux(ntime)
     complex(dp) :: caux(mfreq)

     QMC_FOURIER_BACKWARD: do i=1,norbs

! copy omega data at first, grnw ---> caux
         caux = czero
         do j=1,mfreq
             caux(j) = grnw(j,i)
         enddo ! over j={1,mfreq} loop

! and then invert fourier it, caux ---> raux
         raux = zero
         call hfqmc_fourier_backward(caux, raux)

! copy tau data at last, raux ---> grnt
         do k=1,ntime
             grnt(k,i) = raux(k)
         enddo ! over k={1,ntime} loop

     enddo QMC_FOURIER_BACKWARD ! over i={1,norbs} loop

     return
  end subroutine hfqmc_fourier_w2t

!>>> to calculate the first-order derivate of G(0) and G(\beta)
  subroutine hfqmc_fourier_nderiv(gt, d1, d2)
     use constants
     use control

     implicit none

! external arguments
! green's function data
     real(dp), intent(in)  :: gt(ntime)

! first order derivated of G(0)
     real(dp), intent(out) :: d1

! first order derivated of G(\beta)
     real(dp), intent(out) :: d2

! local variables
! real(dp) dummy variable, means \delta \tau
     real(dp) :: deltau

! calculate deltau
     deltau = beta / real(ntime)

! calculate first-order derivate of G(0)
     d1 = ( -25.0_dp * gt(1) + 48.0_dp * gt(2) - 36.0_dp * gt(3)  &
         + 16.0_dp * gt(4) - 3.0_dp * gt(5) ) / ( 12.0_dp * deltau )

! calculate first-order derivate of G(\beta)
     d2 = ( 25.0_dp * ( -1.0_dp - gt(1) ) - 48.0_dp * gt(ntime)   &
         + 36.0_dp * gt(ntime-1) - 16.0_dp * gt(ntime-2)          &
         + 3.0_dp * gt(ntime-3) ) / ( 12.0_dp * deltau )

     return
  end subroutine hfqmc_fourier_nderiv

!>>> fourier transformation, from imaginary time to matsubara freqency
  subroutine hfqmc_fourier_forward(taudat, omegadat, deriv1, deriv2)
     use constants
     use control

     implicit none

! external arguments
! first order derivated of G(0)
     real(dp), intent(in) :: deriv1

! first order derivated of G(\beta)
     real(dp), intent(in) :: deriv2

! original data in imaginary-time
     real(dp), intent(in) :: taudat(ntime)

! fouriered data in matsubara frequency
     complex(dp), intent(out) :: omegadat(mfreq)

! local variables
! loop index
     integer  :: i
     integer  :: j

! be equal to ntime
     integer  :: l

! frequency
     real(dp) :: omega

! \delta \tau = \bete/(time slice)
     real(dp) :: delta

! auxliary variable for solving xm vector
     real(dp) :: p

! complex(dp) dummy variables
     complex(dp) :: ex
     complex(dp) :: explus
     complex(dp) :: cdummy

! spline interpolation coefficience alpha_{i}
     real(dp) :: a(ntime)

! spline interpolation coefficience beta_{i}
     real(dp) :: b(ntime)

! spline interpolation coefficience gamma_{i}
     real(dp) :: c(ntime)

! spline interpolation coefficience delta_{i}
     real(dp) :: d(ntime)

! auxliary vector for solving linear equations to obtain xm vector
     real(dp) :: u(ntime+1)

! auxliary vector for solving linear equations to obtain xm vector
     real(dp) :: q(ntime+1)

! second derivatives at knots x_{i} of the desired spline function
     real(dp) :: xm(ntime+1)

! a copy of original data, in addition, the G(l+1) is supplemented
     real(dp) :: tcopy(ntime+1)

! evaluate l, a integer dummy variables
     l = ntime

! evaluate \delta \tau
     delta = beta / real(l)

! build tcopy vectors
     do i=1,l
         tcopy(i) = taudat(i)
     enddo ! over i={1,l} loop
     tcopy(l+1) = - one - taudat(1)

! spline interpolation: the spline is given by
!    G(tau) = a(i) + b(i)*(tau-tau_i) + c(i)*(tau-tau_i)^2 + d(i)*(tau-tau_i)^3
! the following formules are taken directly from the book written
! by Stoer and Bulirsch, p. 102
     if ( deriv1 > .99E30 ) then
         q(1) = zero
         u(1) = zero
     else
         q(1) = -half
         u(1) = (3.0_dp / delta) * ( ( tcopy(2) - tcopy(1) ) / delta - deriv1 )
     endif

     do i=2,l
         p = q(i-1) / two + two
         q(i) = - one / two / p
         u(i) = 3.0_dp / delta**2 * ( tcopy(i+1) + tcopy(i-1) - two * tcopy(i) )
         u(i) = ( u(i) - u(i-1) / two ) / p
     enddo ! over i={2,l} loop

     if ( deriv2 > .99E30 ) then
         q(l+1) = zero
         u(l+1) = zero
     else
         q(l+1) = half
         u(l+1) = ( 3.0_dp / delta ) * ( deriv2 - ( tcopy(l+1) - tcopy(l) ) / delta )
     endif

     xm(l+1) = ( u(l+1) - q(l+1) * u(l) ) / ( q(l+1) * q(l) + 1 )
     do i=l,1,-1
         xm(i) = q(i) * xm(i+1) + u(i)
     enddo ! over i={l,1} loop

! the following formulas are taken directly from the book written
! by Stoer and Bulirsch, p. 98
     do i=1,l
         a(i) = tcopy(i)
         c(i) = xm(i) / two
         b(i) = ( tcopy(i+1) - tcopy(i) ) / delta - ( two * xm(i) + xm(i+1) ) * delta / 6.0_dp
         d(i) = ( xm(i+1) - xm(i) ) / ( 6.0_dp * delta )
     enddo ! over i={1,l} loop

! the spline multiplied by the exponential can now be exlicitely
! integrated. the following formules were obtained using mathematica
     do i=1,mfreq
         omega = ( two * real(i-1) + one ) * pi / beta
         omegadat(i) = zero
         do j=1,l
             cdummy = czi * omega * delta * (j-1)
             ex     = exp(cdummy)

             cdummy = czi * omega * delta * j
             explus = exp(cdummy)

             cdummy = 6.0_dp * d(j) / (omega**4)
             cdummy = cdummy - two * czi * c(j) / (omega**3)
             cdummy = cdummy - b(j) / (omega**2)
             cdummy = cdummy + czi * a(j) / omega

             omegadat(i) = omegadat(i) + ex * cdummy

             cdummy = ( zero - 6.0_dp * d(j) ) / (omega**4)
             cdummy = cdummy + czi * ( two * c(j) + 6.0_dp * delta * d(j) ) / (omega**3)
             cdummy = cdummy + ( b(j) + two * delta * c(j) + 3.0_dp * (delta**2) * d(j) ) / (omega**2)
             cdummy = cdummy + czi * ( -a(j) - delta * b(j) - (delta**2) * c(j) - (delta**3) * d(j) ) / omega

             omegadat(i) = omegadat(i) + explus * cdummy
         enddo ! over j={1,l} loop
     enddo ! over i={1,mfreq} loop

     return
  end subroutine hfqmc_fourier_forward

!>>> invert fourier green's or weiss's function from matsubara frequency
! representation to imaginary time representation
  subroutine hfqmc_fourier_backward(omegadat, taudat)
     use constants
     use control

     implicit none

! external arguments
! invfouriered data in imaginary time
     real(dp), intent(out) :: taudat(ntime)

! original data in matsubara frequency
     complex(dp), intent(in) :: omegadat(mfreq)

! local variables
! loop index
     integer  :: i
     integer  :: j

! $\tau$ variable
     real(dp) :: tau

! frequency point
     real(dp) :: omega

! real(dp) dummy variable
     real(dp) :: dummy

     do i=1,ntime
         dummy = zero
         tau = real(i-1) * beta / real(ntime)
         do j=1,mfreq
             omega = ( two * real(j-1) + one ) * pi / beta
             dummy = dummy + cos(omega * tau) *    real( omegadat(j) )
             dummy = dummy + sin(omega * tau) * ( aimag( omegadat(j) ) + one / omega )
         enddo ! over j={1,mfreq} loop
         taudat(i) = two * dummy / beta - half
     enddo ! over i={1,ntime} loop

     return
  end subroutine hfqmc_fourier_backward
