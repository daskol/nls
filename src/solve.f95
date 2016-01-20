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

        integer, parameter :: n = 1000, iters = 100000, order = 5
        integer :: i
        real(sp), parameter :: dt = 0.001, dx = 0.1
        real(sp), parameter :: R = 0.05, g = 1.0e-3, gamma = 0.566, gamma_R = 10, tilde_g = 0.11
        real(sp), dimension(n) :: pumping
        real(sp), dimension(23) :: coeffs  ! TODO: fix magic number `23`
        complex(sp), dimension(n) :: u0, u

        ! Build pumping profile
        do i = 1, n
            pumping(i) = 3.0 * exp( - ((i - 1) * dx) ** 2.0 / (2.0 * 6.8))
        end do

        coeffs = 0.0

        ! NLS equation coeficients
        coeffs(1) = 1.0  ! \partial_t
        coeffs(1) = 1.0  ! \nabla^2
        coeffs(2) = R / (4.0 * tilde_g)  ! 
        coeffs(3) = gamma * 1.0 / 2.0  ! linear damping
        coeffs(4) = 1.0  ! nonlinearity
        coeffs(5) = 1.0  ! interaction to reservoir

        ! Reservoir equation coefficients
        coeffs(10) = 0.0  ! \parital_t
        coeffs(11) = 2.0 * tilde_g * 1.0 / gamma_R  ! pumping coefficient
        coeffs(12) = 1.0  ! damping
        coeffs(13) = R / (gamma_R * g)  ! interaction term
        coeffs(14) = 0.0  ! diffusive term

        ! Inital solution
        u0 = 0.1

        ! Solve
        call solve_nls(dt, dx, n, order, iters, pumping, coeffs, u0, u)

        print *, real(conjg(u) * u)
    end subroutine solve_equation

end program solve
