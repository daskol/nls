!
!   solve.f95
!
!   (c) Daniel Bershatsky <daniel.bershatsky@skolkovotech.ru>, 2015
!

program solve
    use nls

    integer, parameter :: sp = selected_real_kind(6, 37)
    integer, parameter :: dp = selected_real_kind(15, 307)

    call solve_equation()

contains

    subroutine solve_equation()
        implicit none

        integer, parameter :: n = 200, iters = 10000000, order = 5
        complex(sp), dimension(n) :: u
        real(sp), parameter :: dt = 0.000001, dx = 0.001

        call solve_nls(dt, dx, n, order, iters, u)

        print *, real(conjg(u) * u)
    end subroutine solve_equation

end program solve
