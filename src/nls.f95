!
!   nls.f95
!
!   (c) Daniel Bershatsky, 2015
!

module nls
    implicit none

    private

    public :: make_banded_matrix
    public :: make_laplacian_o3, make_laplacian_o5, make_laplacian_o7, make_laplacian
    public :: rgbmv
    public :: hamiltonian
    public :: runge_kutta
    public :: solve_nls

contains

    pure subroutine make_banded_matrix(n, m, row, mat)        
        implicit none

        integer, parameter :: dp = selected_real_kind(15, 307)
        integer, intent(in) :: n, m
        real(dp), intent(in), dimension(m) :: row
        real(dp), intent(out), dimension(m, n) :: mat

        integer :: j, k, l, r

        k = (m - 1) / 2
        l = k + 1
        r = n - k

        ! Fill left triangle part of band
        do j = 1, k
            mat(k + 2 - j + 0:, j) = row(k + 2 - j:)
            mat(:k + 2 - j - 1, j) = 0 ! can be removed
        end do

        ! Fill body
        do j = l, r
            mat(:, j) = row(:)
        end do

        ! Fill right triangle part of band
        do j = 1, k
            mat(:m - j + 0, n - k + j) = row(:m - j)
            mat(m - j + 1:, n - k + j) = 0 ! can be removed
        end do
    end subroutine make_banded_matrix

    pure subroutine make_laplacian_o3(n, h, op)
        implicit none

        integer, parameter :: dp = selected_real_kind(15, 307)
        integer, intent(in) :: n
        integer, parameter :: m = 3
        real(dp), intent(in) :: h
        real(dp), intent(out), dimension(m, n) :: op

        op = h
    end subroutine make_laplacian_o3

    subroutine make_laplacian_o5(n, h, op)
        implicit none

        integer, parameter :: dp = selected_real_kind(15, 307)
        integer, intent(in) :: n
        integer, parameter :: m = 5
        real(dp), intent(in) :: h
        real(dp), intent(out), dimension(m, n) :: op

        integer :: i, j
        real(dp), dimension(m) :: row1d, row2d
        real(dp), dimension(m, n) :: L1, L2
        real(dp) :: dx1, dx2

        dx1 = 12 * h ** 1
        dx2 = 24 * h ** 2

        row1d = (/ -1, 8, 0, -8, 1/) / dx1
        row2d = (/ -1, 16, -30, 16, -1 /) / (24 * h ** 2)

        call make_banded_matrix(n, m, row1d, L1)
        call make_banded_matrix(n, m, row2d, L2)

        L1(3, 2) = L1(3, 2) + 1.0 / dx1
        L2(3, 2) = L2(3, 2) - 1.0 / dx2
        L2(3, 1) = 2 * L2(3, 1)
        L2(2, 2) = 4 * L2(2, 2)
        L2(1, 3) = 4 * L2(1, 3)

        ! Fill in origin
        do j = 1, m - 2
            L1(3 - j + 1, j) = 0.0!L2(3 - j + 1, j)
        end do
        do j = 1, m - 1
            L1(4 - j + 1, j) = L1(4 - j + 1, j) / (1 * h)
        end do

        ! Fill far away from origin and infinity
        do i = 1, n - m + 1
            do j = 1, m
                L1(m - j + 1, i + j - 1) = L1(m - j + 1, i + j - 1) / ((1 + i) * h)
            end do
        end do

        ! Fill tail
        do i = n - m + 2, n
            do j = 1, n - i + 1
                L1(m - j + 1, i + j - 1) = L1(m - j + 1, i + j - 1) / ((i + 1)* h)
            end do
        end do

        op = L1 + L2
    end subroutine make_laplacian_o5

    pure subroutine make_laplacian_o7(n, h, op)
        implicit none

        integer, parameter :: dp = selected_real_kind(15, 307)
        integer, intent(in) :: n
        integer, parameter :: m = 7
        real(dp), intent(in) :: h
        real(dp), intent(out), dimension(m, n) :: op

        op = h
    end subroutine make_laplacian_o7

    subroutine make_laplacian(n, m, h, op)
        implicit none

        integer, parameter :: dp = selected_real_kind(15, 307)
        integer, intent(in) :: n
        integer, intent(in) :: m ! order
        real(dp), intent(in) :: h
        real(dp), intent(out), dimension(m, n) :: op

        if (m == 3) then
            call make_laplacian_o3(n, h, op)
        else if (m == 5) then
            call make_laplacian_o5(n, h, op)
        else if (m == 7) then
            call make_laplacian_o7(n, h, op)
        else
        end if
    end subroutine make_laplacian

    subroutine rgbmv(x, u, sign, op, klu, n) ! reduced gbmv()
        implicit none

        integer, parameter :: dp = selected_real_kind(15, 307)
        integer, intent(in) :: n, klu
        real(dp), intent(in), dimension(n) :: x
        real(dp), intent(out), dimension(n) :: u
        real(dp), intent(in), dimension(2 * klu + 1, n) :: op
        real(dp), intent(in) :: sign

        call dgbmv('N', n, n, klu, klu, sign, op, 2 * klu + 1, x, 1, 1.0, u, 1)
    end subroutine rgbmv

    subroutine pumping(x_exp, u_sqr, p, n)
        implicit none

        integer, parameter :: dp = selected_real_kind(15, 307)
        integer, intent(in) :: n
        real(dp), intent(in), dimension(n) :: x_exp
        real(dp), intent(in), dimension(n) :: u_sqr
        real(dp), intent(out), dimension(n) :: p ! actually n(r)

        real(dp), parameter :: d = 5.0

        p = x_exp / (1.0 + d * u_sqr)
    end subroutine

    subroutine hamiltonian(x, x_exp, u, v, op, klu, n)
        implicit none

        integer, parameter :: dp = selected_real_kind(15, 307)
        integer, intent(in) :: n, klu
        complex(dp), intent(in), dimension(n) :: u
        complex(dp), intent(out), dimension(n) :: v
        real(dp), intent(in), dimension(n) :: x, x_exp
        real(dp), intent(in), dimension(2 * klu + 1, n) :: op

        complex(dp), parameter :: i = (0.0, 1.0)
        real(dp), parameter :: sign = 1.0, a = 1.13636363636, b = 0.283
        real(dp), dimension(n) :: p, u_sqr
        real(dp), dimension(n) :: v_real, v_imag, u_real, u_imag

        u_real = real(u)
        u_imag = aimag(u)
        u_sqr = real(conjg(u) * u)
        v_real = (a * p - b) * u_real + (u_sqr + p) * u_imag
        v_imag = (a * p - b) * u_imag - (u_sqr + p) * u_real

        call pumping(x_exp, u_sqr, p, n)
        call rgbmv(u_imag, v_real, -sign, op, klu, n)
        call rgbmv(u_real, v_imag, +sign, op, klu, n)

        v = cmplx(v_real, v_imag, dp)
    end subroutine

    subroutine runge_kutta(dt, dx, t0, u0, op, n, order, iters, u)
        implicit none

        integer, parameter :: dp = selected_real_kind(15, 307)
        integer, intent(in) :: n, order, iters
        real(dp), intent(in) :: dt, dx, t0
        real(dp), intent(in), dimension(order, n) :: op
        complex(dp), intent(in), dimension(n) :: u0
        complex(dp), intent(out), dimension(n) :: u

        integer :: i, klu
        complex(dp) :: t
        complex(dp), dimension(n) :: k1, k2, k3, k4
        real(dp), dimension(n) :: x, x_exp
        real(dp), parameter :: mu = 6.84931506849, p_0 = 1.87

        do i = 1, n
            x(i) = (i - 1) * dx
            x_exp(i) = p_0 * exp(- x(i) ** 2 / (2 * mu))
        end do

        klu = (order - 1) / 2
        u = u0
        t = t0

        do i = 1, iters
            call hamiltonian(x, x_exp, u + 0. * dt / 2, k1, op, 2, n)
            call hamiltonian(x, x_exp, u + k1 * dt / 2, k2, op, 2, n)
            call hamiltonian(x, x_exp, u + k2 * dt / 2, k3, op, 2, n)
            call hamiltonian(x, x_exp, u + k3 * dt / 1, k4, op, 2, n)

            u = u + (k1 + 2 * k2 + 2 * k3 + k4) * dt / 6
            t = t + dt
        end do
    end subroutine runge_kutta

    subroutine solve_nls(dt, dx, n, order, iters, u)
        implicit none

        integer, parameter :: dp = selected_real_kind(15, 307)
        real(dp), intent(in) :: dt, dx
        integer, intent(in) :: n, order, iters
        complex(dp), intent(out), dimension(n) :: u

        complex(dp), dimension(n) :: u0
        real(dp), parameter :: t0 = 0.0
        real(dp), dimension(order, n) :: op

        u0 = 0.01

        call make_laplacian(n, order, dx, op)
        call runge_kutta(dt, dx, t0, u0, op, n, order, iters, u)
    end subroutine solve_nls

end module nls
