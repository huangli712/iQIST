!-------------------------------------------------------------------------
! project : jasmine
! program : atomic_print_header
!         : atomic_print_footer
!         : atomic_print_summary
!         : atomic_print_error
!         : atomic_print_exception
! source  : atomic_natural.f90
! type    : subroutines
! author  : yilin wang (email: qhwyl2006@126.com)
! history : 07/09/2014 by yilin wang
! purpose : print information
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

!>>> print header
subroutine atomic_print_header()
    use constants,  only: mystd
    use control,    only: nprocs

    ! string for current date and time
    character (len = 20) :: date_time_string

    ! obtain current date and time
    call atomic_time_builder(date_time_string)

    write(mystd,'(2X,a)') 'JASMINE'
    write(mystd,'(2X,a)') '>>> AN ATOMIC PROGRAM FOR CTQMC'
    
    write(mystd,*)

    write(mystd,'(2X,a)') 'version: 2014.07.08T '//'(built at '//__TIME__//" "//__DATE__//')'
    write(mystd,'(2X,a)') 'develop: by yilin wang @IOP'
    write(mystd,'(2X,a)') 'support: qhwyl2006@126.com'
    write(mystd,'(2X,a)') 'license: GPL2 and later versions'
    write(mystd,*)

    write(mystd,'(2X,a)') 'JASMINE >>> start running at '//date_time_string

    write(mystd,'(2X,a,i4)') 'JASMINE >>> parallelism: No  >>> processors:', 1

    write(mystd,*)

    return
end subroutine atomic_print_header

!>>> print footer
subroutine atomic_print_footer()
    use constants, only: dp, mystd

    implicit none

    ! string for current date and time
    character (len = 20) :: date_time_string

    ! used to record the time usage information
    real(dp) :: tot_time

    ! obtain time usage information
    call cpu_time(tot_time)

    ! obtain current date and time
    call atomic_time_builder(date_time_string)

    write(mystd,'(2X,a,f10.2,a)') 'JASMINE >>> total time spent:', tot_time, 's'
    write(mystd,*)

    write(mystd,'(2X,a)') 'JASMINE >>> I am tired and want to go to bed. Bye!'
    write(mystd,'(2X,a)') 'JASMINE >>> happy ending at '//date_time_string

    return
end subroutine atomic_print_footer

!>>> print summary
subroutine atomic_print_summary()
    use constants, only: mystd
    use control

    implicit none

    write(mystd,'(2X,a)') 'JASMINE >>> parameters list:'

    write(mystd,'(2(4X,a,i10))')   'itask :', itask  , 'ictqmc :', ictqmc
    write(mystd,'(2(4X,a,i10))')   'icf   :', icf  ,   'isoc   :', isoc
    write(mystd,'(1(4X,a,i10))')   'icu   :', icu  

    write(mystd,'(2(4X,a,i10))')   'nband :', nband  , 'nspin  :', nspin
    write(mystd,'(2(4X,a,i10))')   'norbs :', norbs  , 'ncfgs  :', ncfgs

    write(mystd,'(2(4X,a,f10.5))') 'Uc    :', Uc     , 'Uv     :', Uv
    write(mystd,'(2(4X,a,f10.5))') 'Jz    :', Jz     , 'Js     :', Js
    write(mystd,'(1(4X,a,f10.5))') 'Jp    :', Jp     

    write(mystd,'(2(4X,a,f10.5))') 'F0    :', F0   ,   'F2     :', F2
    write(mystd,'(2(4X,a,f10.5))') 'F4    :', F4   ,   'F6     :', F6

    write(mystd,'(1(4X,a,f10.5))') 'lambda:', lambda

    write(mystd,*)

    return
end subroutine atomic_print_summary

!>>> print warning
subroutine atomic_print_exception(sub, msg)
    use constants,  only: mystd

    implicit none

    ! external arguments
    ! subroutine name
    character(len=*), intent(in) :: sub

    ! exception message
    character(len=*), intent(in) :: msg

    ! print error information
    write(mystd,'(2X,4a)') 'runtime exception occurred in ', sub, ': ', msg

    ! CONTINUE/PAUSE THE PROGRAM
    !-------------------------------------------------------------------------
         CONTINUE ! OR PAUSE
    !-------------------------------------------------------------------------

    return
end subroutine atomic_print_exception

!>>> print error
subroutine atomic_print_error(sub, msg)
    use constants, only: mystd

    implicit none

    ! external arguments
    ! subroutine name
    character(len=*), intent(in) :: sub

    ! error message
    character(len=*), intent(in) :: msg

    ! print error information
    write(mystd,'(2X,4a)') 'fatal error occurred in ', sub, ': ', msg

    ! TERMINATE THE PROGRAM
    !-------------------------------------------------------------------------
         STOP
    !-------------------------------------------------------------------------

    return
end subroutine atomic_print_error
