!
!   solve.f95
!
!   (c) Daniel Bershatsky <daniel.bershatsky@skolkovotech.ru>, 2015
!

program solve
    use nls

    call test_solve_nls()

contains

    subroutine test_solve_nls()
        implicit none

        integer, parameter :: dp = selected_real_kind(15, 307)
        integer, parameter :: n = 1000, iters = 100000, order = 5
        complex(dp), dimension(n) :: u
        real(dp), parameter :: dt = 0.00001, dx = 0.01

        call solve_nls(dt, dx, n, order, iters, u)

        print *, real(conjg(u) * u)
    end subroutine

end program solve