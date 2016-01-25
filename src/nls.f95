!   nls.f95
!   Computational core of the project. Define routines that solve Non-Linear Schrodinger equation.
!   (c) Daniel Bershatsky, 2015-2016
!   See LISENCE for details.
module nls
    implicit none

    private

    integer, parameter :: sp = selected_real_kind(6, 37)
    integer, parameter :: dp = selected_real_kind(15, 307)

    public :: make_banded_matrix
    public :: clear_first_row_of_derivative
    public :: divide_derivative_on_radius
    public :: make_laplacian_o3, make_laplacian_o5, make_laplacian_o7, make_laplacian, make_laplacian_2d
    public :: rgbmv, rbbmv
    public :: revervoir, revervoir_2d
    public :: hamiltonian, hamiltonian_2d
    public :: runge_kutta, runge_kutta_2d
    public :: solve_nls_1d, solve_nls_2d, solve_nls

contains

    !   \brief Make banded matrix from a row. Banded matrix representation corresponds to BLAS-2 documentation. Memory
    !   usage O(nm), time complexity O(nm).
    !
    !   \param[in] n
    !   \verbatim
    !       Size of square banded matrix `mat`.
    !   \endverbatim
    !
    !   \param[in] m
    !   \verbatim
    !       Length of `row`.
    !   \endverbatim
    !
    !   \param[in] row
    !   \verbatim
    !       Row that composes banded matrix.
    !   \endverbatim
    !
    !   \param[out] mat
    !   \verbatim
    !       Banded matrix as a Fortran array of shape (m, n).
    !   \endverbatim
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

    pure subroutine clear_first_row_of_derivative(n, m, L1)
        implicit none

        integer, parameter :: sp = selected_real_kind(6, 37)
        integer, intent(in) :: n
        integer, intent(in) :: m
        real(sp), intent(out), dimension(m, n) :: L1

        integer :: i

        ! Here (m + 1) / 2 is middle row of banded materix L1
        do i = 1, (m + 1) / 2
            L1(i, (m + 1) / 2 - i + 1) = 0.0
        end do
    end subroutine clear_first_row_of_derivative

    pure subroutine divide_derivative_on_radius(n, m, h, L1)
        implicit none

        integer, parameter :: sp = selected_real_kind(6, 37)
        integer, intent(in) :: n
        integer, intent(in) :: m
        real(sp), intent(in) :: h
        real(sp), intent(out), dimension(m, n) :: L1

        integer :: i, j, row

        ! Row k = i + j - k0: k = 0,1,2,...,n - 1, where k0: m = 2 k0 - 3. Cause is sum (i + j) is invarian of row.
        do i = 1, m
            do j = 1, n
                row = i + j - (m + 3) / 2

                if (row > 0) then
                    L1(i, j) = L1(i, j) / (row * h)
                end if
            end do
        end do
    end subroutine divide_derivative_on_radius

    pure subroutine make_laplacian_o3(n, h, op)
        implicit none

        integer, parameter :: sp = selected_real_kind(6, 37)
        integer, intent(in) :: n
        integer, parameter :: m = 3
        real(sp), intent(in) :: h
        real(sp), intent(out), dimension(m, n) :: op

        real(sp), dimension(m) :: row1d, row2d
        real(sp), dimension(m, n) :: L1, L2
        real(sp) :: dx1, dx2

        dx1 = 2 * h ** 1
        dx2 = 1 * h ** 2

        row1d = (/ 1, 0, -1 /) / dx1
        row2d = (/ 1, -2, 1 /) / dx2

        call make_banded_matrix(n, m, row1d, L1)
        call make_banded_matrix(n, m, row2d, L2)

        L2(2, 1) = 2 * L2(2, 1)  ! node #0, row #0
        L2(1, 2) = 4 * L2(1, 2)  ! node #0, row #0

        call clear_first_row_of_derivative(n, m, L1)
        call divide_derivative_on_radius(n, m, h, L1)

        op = L1 + L2
    end subroutine make_laplacian_o3

    subroutine make_laplacian_o5(n, h, op)
        implicit none

        integer, parameter :: sp = selected_real_kind(6, 37)
        integer, intent(in) :: n
        integer, parameter :: m = 5
        real(sp), intent(in) :: h
        real(sp), intent(out), dimension(m, n) :: op

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

        call clear_first_row_of_derivative(n, m, L1)
        call divide_derivative_on_radius(n, m, h, L1)

        op = L1 + L2
    end subroutine make_laplacian_o5

    pure subroutine make_laplacian_o7(n, h, op)
        implicit none

        integer, parameter :: sp = selected_real_kind(6, 37)
        integer, intent(in) :: n
        integer, parameter :: m = 7
        real(sp), intent(in) :: h
        real(sp), intent(out), dimension(m, n) :: op

        real(sp), dimension(m) :: row1d, row2d
        real(sp), dimension(m, n) :: L1, L2
        real(sp) :: dx1, dx2

        dx1 = 60 * h ** 1
        dx2 = 180 * h ** 2

        row1d = (/ 1, -9, 45, 0, -45, 9, -1/) / dx1
        row2d = (/ 2, -27, 270, -490, 270, -27, 2 /) / dx2

        call make_banded_matrix(n, m, row1d, L1)
        call make_banded_matrix(n, m, row2d, L2)

        ! The second row
        L1(5, 1) = L1(5, 1)
        L1(4, 2) = L1(4, 2) + 9.0 / dx1
        L1(3, 3) = L1(3, 3) - 1.0 / dx1
        L1(2, 4) = L1(2, 4)
        L1(1, 5) = L1(1, 5)

        ! The third row
        L1(6, 1) = L1(6, 1)
        L1(5, 2) = L1(5, 2) - 1.0 / dx1
        L1(4, 3) = L1(4, 3)
        L1(3, 4) = L1(3, 4)
        L1(2, 5) = L1(2, 5)
        L1(1, 6) = L1(1, 6)

        ! The first row
        L2(4, 1) = 2 * L2(4, 1)
        L2(3, 2) = 4 * L2(3, 2)
        L2(2, 3) = 4 * L2(2, 3)
        L2(1, 4) = 4 * L2(1, 4)

        ! The second row
        L2(5, 1) = L2(5, 1)
        L2(4, 2) = L2(4, 2) - 27.0 / dx2
        L2(3, 3) = L2(3, 3) + 2.0 / dx2
        L2(2, 4) = L2(2, 4)
        L2(1, 5) = L2(1, 5)

        ! The third row
        L2(6, 1) = L2(6, 1)
        L2(5, 2) = L2(5, 2) + 2.0 / dx2
        L2(4, 3) = L2(4, 3)
        L2(3, 4) = L2(3, 4)
        L2(2, 5) = L2(2, 5)
        L2(1, 6) = L2(1, 6)

        call clear_first_row_of_derivative(n, m, L1)
        call divide_derivative_on_radius(n, m, h, L1)

        op = L1 + L2
    end subroutine make_laplacian_o7

    !   \brief Make laplacian matrix that approximates laplacian operator in axial symmetric case on a given grid and
    !   given approximation. TODO: implement approximation of order 3 and 7.
    !
    !   \param[in] n
    !   \verbatim
    !       Size of square banded matrix `op`.
    !   \endverbatim
    !
    !   \param[in] m
    !   \verbatim
    !       Order of approximation.
    !   \endverbatim
    !
    !   \param[out] op
    !   \verbatim
    !       Banded matrix of size `n` that approximates laplacian operator in order `m`.
    !   \endverbatim
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

    ! TODO: implement high order approximations.
    subroutine make_laplacian_2d(n, m, h, blocks, orders)
        implicit none

        integer, intent(in) :: n, m
        integer, intent(out), dimension(m) :: orders
        real(sp), intent(in) :: h
        real(sp), intent(out), dimension(2 * m - 1, n) :: blocks

        real(sp) :: dx2
        real(sp), dimension(1) :: left, right
        real(sp), dimension(3) :: middle

        dx2 = 1.0 / h ** 2

        left = (/ 1 /) / dx2
        middle = (/ 1, -4, 1 /) / dx2
        right = (/ 1 /) / dx2

        call make_banded_matrix(n, 1, left, blocks(1:1, :))
        call make_banded_matrix(n, 3, middle, blocks(2:4, :))
        call make_banded_matrix(n, 1, right, blocks(5:5, :))

        orders = (/ 0, 1, 0 /)
    end subroutine make_laplacian_2d

    !   \brief Perform one of the matrix-vector operations   y := alpha*A*x + beta*y, or y := alpha*A'*x + beta*y in
    !   case of block band matrix.
    !   TODO: support high-order approximations.
    subroutine rbbmv(x, y, sign, blocks, ms, m, n)
        implicit none

        integer, intent(in) :: m, n
        integer, intent(in), dimension(m) :: ms
        real(sp), intent(in) :: sign
        real(sp), intent(in), dimension(n * n) :: x
        real(sp), intent(out), dimension(n * n) :: y
        real(sp), intent(in), dimension(2 * m - 1, n) :: blocks

        integer :: i
        integer, dimension(m) :: orders

        orders = (/ 0, 1, 0 /)

        y = 0.0

        i = 1
        call rgbmv(x((i - 1) * n + 1:(i + 0) * n), y((i - 1) * n + 1:i * n), sign, blocks(2:4, :), ms(2), n)
        call rgbmv(x((i + 0) * n + 1:(i + 1) * n), y((i - 1) * n + 1:i * n), sign, blocks(5:5, :), ms(3), n)

        do i = 2, n - 1
            call rgbmv(x((i - 2) * n + 1:(i - 1) * n), y((i - 1) * n + 1:i * n), sign, blocks(1:1, :), ms(1), n)
            call rgbmv(x((i - 1) * n + 1:(i + 0) * n), y((i - 1) * n + 1:i * n), sign, blocks(2:4, :), ms(2), n)
            call rgbmv(x((i - 0) * n + 1:(i + 1) * n), y((i - 1) * n + 1:i * n), sign, blocks(5:5, :), ms(3), n)
        end do

        i = n
        call rgbmv(x((i - 2) * n + 1:(i - 1) * n), y((i - 1) * n + 1:i * n), sign, blocks(5:5, :), ms(1), n)
        call rgbmv(x((i - 1) * n + 1:(i - 0) * n), y((i - 1) * n + 1:i * n), sign, blocks(2:4, :), ms(2), n)
    end subroutine rbbmv

    !   \brief Wrapper of BLAS-2 function `sgbmv` which is matvec implementation for banded matrix.
    subroutine rgbmv(x, u, sign, op, klu, n) ! reduced gbmv()
        implicit none

        integer, parameter :: sp = selected_real_kind(6, 37)
        integer, intent(in) :: n, klu
        real(sp), intent(in), dimension(n) :: x
        real(sp), intent(inout), dimension(n) :: u
        real(sp), intent(in), dimension(2 * klu + 1, n) :: op
        real(sp), intent(in) :: sign

        call sgbmv('N', n, n, klu, klu, sign, op, 2 * klu + 1, x, 1, 1.0, u, 1)
    end subroutine rgbmv

    !   \brief Calculate density of revervoir particles with given pumping profile and condensate density.
    !
    !   \param[in] pumping
    !   \verbatim
    !       Pumping profile on spacial grid.
    !   \endverbatim
    !
    !   \param[in] coeffs
    !   \verbatim
    !       Contains coefficients of NLS equation with reservoir. Each coefficient coresponds to each term of NLS
    !       equation or exciton-polariton revervoir. See `solve_nls`.
    !   \endverbatim
    !
    !   \param[in] u_sqr
    !   \verbatim
    !       Squared absolute value of psi-function.
    !   \endverbatim
    !
    !   \param[out] r
    !   \verbatim
    !       Density of revervoir particles.
    !   \endverbatim
    !
    !   \param[in] n
    !   \verbatim
    !       Number of nodes of spacial grid.
    !   \endverbatim
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

    !   \brief Calculate action of hamiltonian operator on function defined on a grid with a given pumping and
    !   condensate density.
    !
    !   \param[in] pumping
    !   \verbatim
    !       Pumping profile on spacial grid.
    !   \endverbatim
    !
    !   \param[in] coeffs
    !   \verbatim
    !       Contains coefficients of NLS equation with reservoir. Each coefficient coresponds to each term of NLS
    !       equation or exciton-polariton revervoir. See `solve_nls`.
    !   \endverbatim
    !
    !   \param[in] u
    !   \verbatim
    !       Psi-function on which hamiltonian operator acts.
    !   \endverbatim
    !
    !   \param[out] v
    !   \verbatim
    !       Result of application of hamiltonian operator on the right hand side.
    !   \endverbatim
    !
    !   \param[out] op
    !   \verbatim
    !       Banded laplacian representation of order `klu` and size `n`. See `make_laplacian`.
    !   \endverbatim
    !
    !   \param[in] klu
    !   \verbatim
    !       Order of laplacian approximation.
    !   \endverbatim
    !
    !   \param[in] n
    !   \verbatim
    !       Number of nodes of spacial grid.
    !   \endverbatim
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

    !   \brief Yet another Runge-Kutta of order 4 implementation with matrix-valeud right hand side. The criterion of
    !   break is steady state solution which is not passed explicitly.
    !
    !   \param[in] dt
    !   \verbatim
    !       Time step.
    !   \endverbatim
    !
    !   \param[in] t0
    !   \verbatim
    !       Initial time.
    !   \endverbatim
    !
    !   \param[in] u0
    !   \verbatim
    !       Initial solution.
    !   \endverbatim
    !
    !   \param[in] op
    !   \verbatim
    !       Banded laplacian representation of order `order` and size `n`. See `make_laplacian`.
    !   \endverbatim
    !
    !   \param[in] n
    !   \verbatim
    !       Number of grid nodes.
    !   \endverbatim
    !
    !   \param[in] order
    !   \verbatim
    !       Order of laplacian approximation.
    !   \endverbatim
    !
    !   \param[in] iters
    !   \verbatim
    !       Number of iterations.
    !   \endverbatim
    !
    !   \param[out] u
    !   \verbatim
    !       Psi-function that solves NLS equation with reservoir.
    !   \endverbatim
    !
    !   \param[in] pumping
    !   \verbatim
    !       Pumping profile on spacial grid.
    !   \endverbatim
    !
    !   \param[in] coeffs
    !   \verbatim
    !       Contains coefficients of NLS equation with reservoir. Each coefficient coresponds to each term of NLS
    !       equation or exciton-polariton revervoir. See `solve_nls`.
    !   \endverbatim
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
            call hamiltonian(pumping, coeffs, u + 0. * dt / 2, k1, op, klu, n)
            call hamiltonian(pumping, coeffs, u + k1 * dt / 2, k2, op, klu, n)
            call hamiltonian(pumping, coeffs, u + k2 * dt / 2, k3, op, klu, n)
            call hamiltonian(pumping, coeffs, u + k3 * dt / 1, k4, op, klu, n)

            u = u + (k1 + 2 * k2 + 2 * k3 + k4) * dt / 6
            t = t + dt
        end do
    end subroutine runge_kutta

    !   \brief Solve Non-Linear Schrodinger equation in axial symmentric geometry. It is based on Runge-Kutta method of
    !   order 4 with laplacian approximation order that passed as argument. Initial solution is uniform everywhere on
    !   the grid.
    !
    !   \param[in] dt
    !   \verbatim
    !       Time step.
    !   \endverbatim
    !
    !   \param[in] dx
    !   \verbatim
    !       Spartial step.
    !   \endverbatim
    !
    !   \param[in] n
    !   \verbatim
    !       Number of grid nodes.
    !   \endverbatim
    !
    !   \param[in] order
    !   \verbatim
    !       Order of laplacian approximation.
    !   \endverbatim
    !
    !   \param[in] iters
    !   \verbatim
    !       Number of iterations.
    !   \endverbatim
    !
    !   \param[in] pumping
    !   \verbatim
    !       Pumping profile on spacial grid.
    !   \endverbatim
    !
    !   \param[in] coeffs
    !   \verbatim
    !       Contains coefficients of NLS equation with reservoir. Each coefficient coresponds to each term of NLS
    !       equation or exciton-polariton revervoir.
    !       Coefficient 1 corresponds to `i \partial_t` term.
    !       Coefficient 2 corresponds to `-\nabla^2` term.
    !       Coefficient 3 corresponds to `i n \psi` term.
    !       Coefficient 4 corresponds to `-i \psi` term.
    !       Coefficient 5 corresponds to `|\psi|^2 \psi` term.
    !       Coefficient 6 corresponds to `n \psi` term.
    !       Coefficient 11 corresponds to `\partial_t` term of reservoire equation.
    !       Coefficient 12 corresponds to `P term of reservoire equation.
    !       Coefficient 13 corresponds to `-n` term of reservoire equation.
    !       Coefficient 14 corresponds to `-n |\psi|^2` term of reservoire equation.
    !       Coefficient 15 corresponds to `\nabla^2` term of reservoire equation.
    !       The coresspondence is induced sequensially with equations in Wouters&Carusotto, 2007.
    !   \endverbatim
    !
    !   \param[in] u_0
    !   \verbatim
    !       Initial psi-function that solves NLS equation with reservoir.
    !   \endverbatim
    !
    !   \param[out] u
    !   \verbatim
    !       Psi-function that solves NLS equation with reservoir.
    !   \endverbatim
    subroutine solve_nls(dt, dx, n, order, iters, pumping, coeffs, u0, u)
        implicit none

        integer, parameter :: sp = selected_real_kind(6, 37)
        real(sp), intent(in) :: dt, dx
        integer, intent(in) :: n, order, iters
        real(sp), intent(in), dimension(n) :: pumping
        real(sp), intent(in), dimension(23) :: coeffs
        complex(sp), intent(in), dimension(n) :: u0
        complex(sp), intent(out), dimension(n) :: u

        real(sp), parameter :: t0 = 0.0
        real(sp), dimension(order, n) :: op

        call make_laplacian(n, order, dx, op)
        call runge_kutta(dt, t0, u0, op, n, order, iters, u, pumping, coeffs)
    end subroutine solve_nls

    subroutine solve_nls_1d(dt, dx, n, order, iters, pumping, coeffs, u0, u)
        implicit none

        integer, parameter :: sp = selected_real_kind(6, 37)
        real(sp), intent(in) :: dt, dx
        integer, intent(in) :: n, order, iters
        real(sp), intent(in), dimension(n, n) :: pumping
        complex(sp), intent(in), dimension(n, n) :: u0
        complex(sp), intent(out), dimension(n, n) :: u
        real(sp), intent(in), dimension(23) :: coeffs

        call solve_nls(dt, dx, n, order, iters, pumping, coeffs, u0, u)
    end subroutine solve_nls_1d

    subroutine revervoir_2d(pumping, coeffs, u_sqr, r, n)
        implicit none

        integer, intent(in) :: n
        real(sp), intent(in), dimension(n * n) :: pumping
        real(sp), intent(in), dimension(23) :: coeffs
        real(sp), intent(in), dimension(n * n) :: u_sqr
        real(sp), intent(out), dimension(n * n) :: r ! actually n(r)

        call revervoir(pumping, coeffs, u_sqr, r, n * n)
    end subroutine revervoir_2d

    subroutine hamiltonian_2d(pumping, coeffs, u, v, blocks, orders, order, n)
        implicit none

        integer, intent(in) :: n, order
        integer, intent(in), dimension(order) :: orders
        complex(sp), intent(in), dimension(n * n) :: u
        complex(sp), intent(out), dimension(n * n) :: v
        real(sp), intent(in), dimension(n * n) :: pumping
        real(sp), intent(in), dimension(23) :: coeffs
        real(sp), intent(in), dimension(2 * order - 1, n) :: blocks

        complex(sp), parameter :: i = (0.0, 1.0)
        real(sp), parameter :: sign = 1.0
        real(sp), dimension(n) :: r, u_sqr
        real(sp), dimension(n) :: v_real, v_imag, u_real, u_imag

        u_real = real(u)
        u_imag = aimag(u)
        u_sqr = real(conjg(u) * u)

        call revervoir_2d(pumping, coeffs, u_sqr, r, n)

        v_real = (coeffs(3) * r - coeffs(4)) * u_real + (coeffs(5) * u_sqr + coeffs(6) * r) * u_imag
        v_imag = (coeffs(3) * r - coeffs(4)) * u_imag - (coeffs(5) * u_sqr + coeffs(6) * r) * u_real

        call rbbmv(u_imag, v_real, -sign, blocks, orders, order, n)
        call rbbmv(u_real, v_imag, +sign, blocks, orders, order, n)

        v = cmplx(v_real, v_imag, sp)
    end subroutine hamiltonian_2d

    ! Code is almost the same as `runge_kutta`.
    subroutine runge_kutta_2d(dt, t0, u0, n, blocks, orders, order, iters, u, pumping, coeffs)
        implicit none

        integer, intent(in) :: n, order, iters
        integer, intent(in), dimension(order) :: orders
        real(sp), intent(in) :: dt, t0
        real(sp), intent(in), dimension(23) :: coeffs
        real(sp), intent(in), dimension((order + 1) / 2, n) :: blocks
        real(sp), intent(in), dimension(n * n) :: pumping
        complex(sp), intent(in), dimension(n * n) :: u0
        complex(sp), intent(out), dimension(n * n) :: u

        integer :: i
        real(sp) :: t
        complex(sp), dimension(n) :: k1, k2, k3, k4

        u = u0
        t = t0

        do i = 1, iters
            call hamiltonian_2d(pumping, coeffs, u + 0. * dt / 2, k1, blocks, orders, order, n)
            call hamiltonian_2d(pumping, coeffs, u + k1 * dt / 2, k2, blocks, orders, order, n)
            call hamiltonian_2d(pumping, coeffs, u + k2 * dt / 2, k3, blocks, orders, order, n)
            call hamiltonian_2d(pumping, coeffs, u + k3 * dt / 1, k4, blocks, orders, order, n)

            u = u + (k1 + 2 * k2 + 2 * k3 + k4) * dt / 6
            t = t + dt
        end do
    end subroutine runge_kutta_2d

    subroutine solve_nls_2d(dt, dx, n, order, iters, pumping, coeffs, u0, u)
        implicit none

        integer, intent(in) :: n, order, iters
        real(sp), intent(in) :: dt, dx
        real(sp), intent(in), dimension(n, n) :: pumping
        real(sp), intent(in), dimension(23) :: coeffs
        complex(sp), intent(in), dimension(n * n) :: u0
        complex(sp), intent(out), dimension(n * n) :: u

        integer, dimension(order) :: orders
        real(sp), parameter :: t0 = 0.0
        real(sp), dimension((order - 1) / 2) :: blocks

        call make_laplacian_2d(n, 3, dx, blocks, orders)
        call runge_kutta_2d(dt, t0, u0, n, blocks, orders, order, iters, u, pumping, coeffs)
    end subroutine solve_nls_2d

end module nls