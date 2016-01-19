!
!   nls.f95
!
!   (c) Daniel Bershatsky, 2015
!

module nls
    implicit none

    private

    integer, parameter :: sp = selected_real_kind(6, 37)
    integer, parameter :: dp = selected_real_kind(15, 307)

    public :: make_banded_matrix
    public :: make_laplacian_o3, make_laplacian_o5, make_laplacian_o7, make_laplacian
    public :: rgbmv
    public :: revervoir
    public :: hamiltonian
    public :: runge_kutta
    public :: solve_nls

contains

    pure subroutine make_banded_matrix(n, m, row, mat)        
        implicit none

        integer, parameter :: sp = selected_real_kind(6, 37)
        integer, intent(in) :: n, m
        real(sp), intent(in), dimension(m) :: row
        real(sp), intent(out), dimension(m, n) :: mat

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

        integer, parameter :: sp = selected_real_kind(6, 37)
        integer, intent(in) :: n
        integer, parameter :: m = 3
        real(sp), intent(in) :: h
        real(sp), intent(out), dimension(m, n) :: op

        op = h
    end subroutine make_laplacian_o3

    subroutine make_laplacian_o5(n, h, op)
        implicit none

        integer, parameter :: sp = selected_real_kind(6, 37)
        integer, intent(in) :: n
        integer, parameter :: m = 5
        real(sp), intent(in) :: h
        real(sp), intent(out), dimension(m, n) :: op

        integer :: i, j
        real(sp), dimension(m) :: row1d, row2d
        real(sp), dimension(m, n) :: L1, L2
        real(sp) :: dx1, dx2

        dx1 = 12 * h ** 1
        dx2 = 12 * h ** 2

        row1d = (/ -1, 8, 0, -8, 1/) / dx1
        row2d = (/ -1, 16, -30, 16, -1 /) / dx2

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

        integer, parameter :: sp = selected_real_kind(6, 37)
        integer, intent(in) :: n
        integer, parameter :: m = 7
        real(sp), intent(in) :: h
        real(sp), intent(out), dimension(m, n) :: op

        op = h
    end subroutine make_laplacian_o7

    subroutine make_laplacian(n, m, h, op)
        implicit none

        integer, parameter :: sp = selected_real_kind(6, 37)
        integer, intent(in) :: n
        integer, intent(in) :: m ! order
        real(sp), intent(in) :: h
        real(sp), intent(out), dimension(m, n) :: op

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

        integer, parameter :: sp = selected_real_kind(6, 37)
        integer, intent(in) :: n, klu
        real(sp), intent(in), dimension(n) :: x
        real(sp), intent(out), dimension(n) :: u
        real(sp), intent(in), dimension(2 * klu + 1, n) :: op
        real(sp), intent(in) :: sign

        call sgbmv('N', n, n, klu, klu, sign, op, 2 * klu + 1, x, 1, 1.0, u, 1)
    end subroutine rgbmv

    subroutine revervoir(pumping, coeffs, u_sqr, r, n)
        implicit none

        integer, parameter :: sp = selected_real_kind(6, 37)
        integer, intent(in) :: n
        real(sp), intent(in), dimension(n) :: pumping
        real(sp), intent(in), dimension(23) :: coeffs
        real(sp), intent(in), dimension(n) :: u_sqr
        real(sp), intent(out), dimension(n) :: r ! actually n(r)

        r = coeffs(12) * pumping  / (coeffs(13) + coeffs(14) * u_sqr)
    end subroutine

    subroutine hamiltonian(pumping, coeffs, u, v, op, klu, n)
        implicit none

        integer, parameter :: sp = selected_real_kind(6, 37)
        integer, intent(in) :: n, klu
        complex(sp), intent(in), dimension(n) :: u
        complex(sp), intent(out), dimension(n) :: v
        real(sp), intent(in), dimension(n) :: pumping
        real(sp), intent(in), dimension(23) :: coeffs
        real(sp), intent(in), dimension(2 * klu + 1, n) :: op

        complex(sp), parameter :: i = (0.0, 1.0)
        real(sp), parameter :: sign = 1.0
        real(sp), dimension(n) :: r, u_sqr
        real(sp), dimension(n) :: v_real, v_imag, u_real, u_imag

        u_real = real(u)
        u_imag = aimag(u)
        u_sqr = real(conjg(u) * u)

        call revervoir(pumping, coeffs, u_sqr, r, n)

        v_real = (coeffs(3) * r - coeffs(4)) * u_real + (coeffs(5) * u_sqr + coeffs(6) * r) * u_imag
        v_imag = (coeffs(3) * r - coeffs(4)) * u_imag - (coeffs(5) * u_sqr + coeffs(6) * r) * u_real

        call rgbmv(u_imag, v_real, -sign, op, klu, n)
        call rgbmv(u_real, v_imag, +sign, op, klu, n)

        v = cmplx(v_real, v_imag, sp)
    end subroutine

    subroutine runge_kutta(dt, t0, u0, op, n, order, iters, u, pumping, coeffs)
        implicit none

        integer, parameter :: sp = selected_real_kind(6, 37)
        integer, intent(in) :: n, order, iters
        real(sp), intent(in) :: dt, t0
        real(sp), intent(in), dimension(order, n) :: op
        complex(sp), intent(in), dimension(n) :: u0
        complex(sp), intent(out), dimension(n) :: u
        real(sp), intent(in), dimension(n) :: pumping
        real(sp), intent(in), dimension(23) :: coeffs

        integer :: i, klu
        complex(sp) :: t
        complex(sp), dimension(n) :: k1, k2, k3, k4

        klu = (order - 1) / 2
        u = u0
        t = t0

        do i = 1, iters
            call hamiltonian(pumping, coeffs, u + 0. * dt / 2, k1, op, 2, n)
            call hamiltonian(pumping, coeffs, u + k1 * dt / 2, k2, op, 2, n)
            call hamiltonian(pumping, coeffs, u + k2 * dt / 2, k3, op, 2, n)
            call hamiltonian(pumping, coeffs, u + k3 * dt / 1, k4, op, 2, n)

            u = u + (k1 + 2 * k2 + 2 * k3 + k4) * dt / 6
            t = t + dt
        end do
    end subroutine runge_kutta

    !
    !  \brief Solve Non-Linear Schrodinger equation in axial symmentric geometry
    !  \author Daniel Bershatsky
    !
    !   \param[in] dt
    !   \verbatim
    !       Time step.
    !   \endverbatim
    !
    !   \param[in] dx
    !   \param[in] n
    !   \param[in] order
    !   \param[in] iters
    !   \param[out] u
    !   \param[in] pumping
    !   \param[in] coeffs
    !   \verbatim
    !       Contains coefficients of NLS equation with reservoir. Each coefficient coresponds to each term of NLS
    !       equation or exciton-polariton revervoir.
    !       Coefficient 1 corresponds to i \partial_t term.
    !       Coefficient 2 corresponds to - \nabla^2 term.
    !       Coefficient 3 corresponds to i n \psi term.
    !       Coefficient 4 corresponds to - i \psi term.
    !       Coefficient 5 corresponds to |\psi|^3 term.
    !       Coefficient 6 corresponds to n \psi term.
    !       Coefficient 11 corresponds to \partial_t term of reservoire equation.
    !       Coefficient 12 corresponds to P term of reservoire equation.
    !       Coefficient 13 corresponds to - n term of reservoire equation.
    !       Coefficient 14 corresponds to - n |\psi|^2 term of reservoire equation.
    !       Coefficient 15 corresponds to \nabla^2 term of reservoire equation.
    !       The coresspondence is induced by Wouters&Carusotto, 2007.
    !   \endverbatim
    !
    subroutine solve_nls(dt, dx, n, order, iters, u, pumping, coeffs)
        implicit none

        integer, parameter :: sp = selected_real_kind(6, 37)
        real(sp), intent(in) :: dt, dx
        integer, intent(in) :: n, order, iters
        complex(sp), intent(out), dimension(n) :: u
        real(sp), intent(in), dimension(n) :: pumping
        real(sp), intent(in), dimension(23) :: coeffs

        complex(sp), dimension(n) :: u0
        real(sp), parameter :: t0 = 0.0
        real(sp), dimension(order, n) :: op

        u0 = 0.1

        call make_laplacian(n, order, dx, op)
        call runge_kutta(dt, t0, u0, op, n, order, iters, u, pumping, coeffs)
    end subroutine solve_nls

end module nls
