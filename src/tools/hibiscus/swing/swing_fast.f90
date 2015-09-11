!!!-----------------------------------------------------------------------
!!! project : hibiscus/swing
!!! program : fchi
!!!           matsum
!!!           kramskron
!!! source  : swing_fast.f90
!!! type    : subroutines
!!! author  : li huang (email:lihuang.dmft@gmail.com)
!!! history : 01/10/2011 by li huang (created)
!!!           08/17/2015 by li huang (last modified)
!!! purpose : implement several time-consuming subroutines to accelerate
!!!           the hibiscus/swing code
!!! status  : unstable
!!! comment :
!!!-----------------------------------------------------------------------

!!>>> fchi: calculate \chi partly
  subroutine fchi( chi2, chi4, nrm, gweigh, vary, gwfix, fixed, &
                   sqmc, ifunr, ifuni, ders, expand, expand_sig,&
                   nvary, nfix, nal, nom, nds )
     implicit none

! local parameters
     integer, parameter  :: dp = kind(1.0d0)
     real(dp), parameter :: zero = 0.0_dp
     real(dp), parameter :: half = 0.5_dp

! external arguments
     integer, intent(in)   :: nvary, nfix, nal, nom, nds
     real(dp), intent(out) :: chi2, chi4, nrm
     real(dp), intent(in)  :: gweigh(nvary)
     real(dp), intent(in)  :: vary(nvary)
     real(dp), intent(in)  :: gwfix(nfix)
     real(dp), intent(in)  :: fixed(nfix)
     real(dp), intent(in)  :: sqmc(nom,4)
     real(dp), intent(in)  :: ifunr(nal,nom)
     real(dp), intent(in)  :: ifuni(nal,nom)
     real(dp), intent(in)  :: ders(nal,nds)
     real(dp), intent(in)  :: expand(nds)
     real(dp), intent(in)  :: expand_sig(nds)

!f2py integer intent(hide), depend(gweigh) :: nvary=shape(gweigh,0)
!f2py integer intent(hide), depend(gwfix)  :: nfix=shape(gwfix,0)
!f2py integer intent(hide), depend(sqmc)   :: nom=shape(sqmc,0)
!f2py integer intent(hide), depend(ifunr)  :: nal=shape(ifunr,0)
!f2py integer intent(hide), depend(expand) :: nds=shape(expand,0)

! local varibales
     integer  :: i
     integer  :: j

     real(dp) :: sig0
     real(dp) :: gi
     real(dp) :: gr
     real(dp) :: tders(nds)

     sig0 = sqrt( half * ( sqmc(1,3)**2 + sqmc(1,4) ) )

! calculate chi2
     chi2 = zero
     do j=1,nom
         gi = zero
         gr = zero
         do i=1,nvary
             gi = gi + ifuni(vary(i),j) * gweigh(i)
             gr = gr + ifunr(vary(i),j) * gweigh(i)
         enddo ! over i={1,nvary} loop
         do i=1,nfix
             gi = gi + ifuni(fixed(i),j) * gwfix(i)
             gr = gr + ifunr(fixed(i),j) * gwfix(i)
         enddo ! over i={1,nfix} loop
         chi2 = chi2 + ( ( sqmc(j,1) - gr ) * sig0 / sqmc(j,3) )**2
         chi2 = chi2 + ( ( sqmc(j,2) - gi ) * sig0 / sqmc(j,4) )**2
     enddo ! over j={1,nom} loop

! calculate chi4
     tders = zero
     do i=1,nvary
         tders = tders + ders(vary(i),:) * gweigh(i)
     enddo ! over i={1,nvary} loop
     do i=1,nfix
         tders = tders + ders(fixed(i),:) * gwfix(i)
     enddo ! over i={1,nfix} loop
     chi4 = zero
     do i=1,nds
         chi4 = chi4 + ( ( tders(i) - expand(i) ) / expand_sig(i) )**2
     enddo ! over i={1,nds} loop

! calculate nrm
     nrm = sum(gweigh) + sum(gwfix)

     return
  end subroutine fchi

!!>>> matsum: which is used to calculate the high frequency self-energy
!!>>> function on matsubara axis
  subroutine matsum(gr, gi, En, iom, x0, dh, wb, nom, nx)
     implicit none

! local parameters
     integer, parameter  :: dp = kind(1.0d0)
     real(dp), parameter :: zero = 0.0_dp

! external arguments
     integer, intent(in)   :: nom, nx
     real(dp), intent(out) :: gr(nom), gi(nom)
     real(dp), intent(in)  :: En
     real(dp), intent(in)  :: iom(nom)
     real(dp), intent(in)  :: x0(nx)
     real(dp), intent(in)  :: dh(nx)
     real(dp), intent(in)  :: wb(nx)

!f2py integer intent(hide), depend(iom) :: nom=shape(iom,0)
!f2py integer intent(hide), depend(x0)  :: nx=shape(x0,0)

! local variables
     integer  :: i
     integer  :: j

     real(dp) :: simps
     real(dp) :: omn
     real(dp) :: w1
     real(dp) :: sumr
     real(dp) :: sumi

     do i=1,nom
         omn  = iom(i)
         sumr = zero
         sumi = zero
         do j=1,nx
             w1 = wb(j) / ( omn**2 + ( x0(j) + En )**2 )
             sumr = sumr + w1 * ( x0(j) + En ) * dh(j)
             sumi = sumi + w1 * omn * dh(j)
         enddo ! over j={1,nx} loop
         gr(i) = -sumr
         gi(i) = -sumi
     enddo ! over i={1,nom} loop

     return
  end subroutine matsum

!!>>> kramskron: calculate the kramas-kronig transformation, please refer
!!>>> to eq.(116)
  subroutine kramskron(Frc, om, F0, wb, x0, dhx, En, nom, nx)
     implicit none

! local parameters
     integer, parameter  :: dp = kind(1.0d0)
     real(dp), parameter :: zero = 0.0_dp
     real(dp), parameter :: pi   = 3.141592653589793238462643383279_dp

! external arguments
     integer, intent(in)  :: nom, nx
     real(dp), intent(in) :: En
     real(dp), intent(in) :: om(nom)
     real(dp), intent(in) :: F0(nom)
     real(dp), intent(in) :: wb(nx)
     real(dp), intent(in) :: x0(nx)
     real(dp), intent(in) :: dhx(nx)
     complex(dp), intent(out) :: frc(nom)

!f2py integer intent(hide), depend(om)  :: nom=shape(om,0)
!f2py integer intent(hide), depend(x0)  :: nx=shape(x0,0)

! local variables
     integer  :: i
     integer  :: j

     real(dp) :: omi
     real(dp) :: wi
     real(dp) :: sums
     real(dp) :: Fre
     real(dp) :: Fri

     do i=1,nom
         omi = om(i)
         wi = F0(i)
         sums = zero
         do j=1,nx
             sums = sums + ( wb(j) - wi ) * dhx(j) / ( omi - x0(j) - En )
         enddo ! over j={1,nx} loop
         Fre = sums - wi * log( abs( ( x0(nx) + En - omi ) / ( omi - x0(1) - En ) ) )
         Fri = -pi * F0(i)
         Frc(i) = dcmplx(Fre, Fri)
     enddo ! over i={1,nom} loop

     return
  end subroutine kramskron
