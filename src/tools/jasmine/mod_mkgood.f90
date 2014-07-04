module mod_mkgood
    implicit none

    contains

    subroutine make_good_3band(good)
        implicit none

        integer, intent(out) :: good(6)
        ! for 3band case the order is:
        ! -1, 1, -3, -1, 1, 3
        good(1) = -1
        good(2) =  1
        good(3) = -3
        good(4) = -1
        good(5) =  1
        good(6) =  3

        return
    end subroutine make_good_3band

    subroutine make_good_5band(good)
        implicit none

        integer, intent(out) :: good(10)
        ! for 3band case the order is:
        ! -3, -1, 1, 3, -5, -3, -1, 1, 3, 5
        !good(1)  = -3
        !good(2)  = -1
        !good(3)  =  1
        !good(4)  =  3
        !good(5)  = -5
        !good(6)  = -3
        !good(7)  = -1
        !good(8)  =  1
        !good(9)  =  3
        !good(10) =  5

        good(1)  = -3
        good(2)  = -1
        good(3)  = -5
        good(4)  = -3
        good(5)  = -1
        good(6)  =  3
        good(7)  =  1
        good(8)  =  5
        good(9)  =  3
        good(10) =  1

        return
    end subroutine make_good_5band

    subroutine make_good_7band(good)
        implicit none

        integer, intent(out) :: good(14)
        ! for 3band case the order is:
        ! -5, -3, -1, 1, 3, 5, -7, -5, -3, -1, 1, 3, 5, 7
        good(1)  = -5
        good(2)  = -3
        good(3)  = -1
        good(4)  =  1
        good(5)  =  3
        good(6)  =  5
        good(7)  = -7
        good(8)  = -5
        good(9)  = -3
        good(10) = -1
        good(11) =  1
        good(12) =  3
        good(13) =  5
        good(14) =  7

        return
    end subroutine make_good_7band

end module mod_mkgood
