!
!   test/nls.f95
!
!   (c) Daniel Bershatsky, 2015
!

program test_nls
    use nls

    call test()

contains

    subroutine test()
        implicit none

        logical, dimension(10) :: success

        success = .true.

        call test_make_banded_matrix(success(1))
        call test_make_laplacian_o3(success(2))
        call test_make_laplacian_o5(success(3))
        call test_make_laplacian_o7(success(4))
        call test_make_laplacian(success(5))
        call test_rgbmv(success(6))
        call test_pumping(success(7))
        call test_hamiltonian(success(8))
        call test_runge_kutta(success(9))
        call test_solve_nls(success(10))

        print *
        print *, 'Success: ', success
    end subroutine test

    subroutine test_make_banded_matrix(success)
        implicit none

        logical, intent(inout) :: success
        integer, parameter :: dp = selected_real_kind(15, 307)
        integer, parameter :: m = 5, n = 7
        real(dp), parameter :: tolerance = 1.0e-6
        real(dp), dimension(m) :: row
        real(dp), dimension(m, n) :: ans, mat
        logical, dimension(m, n) :: equlity
        integer :: i, j

        print *, 'Testing make_banded_matrix():'

        row = (/ 1, 2, 3, 4, 5 /)
        ans = reshape((/ &
            0, 0, 3, 4, 5, &
            0, 2, 3, 4, 5, &
            1, 2, 3, 4, 5, &
            1, 2, 3, 4, 5, &
            1, 2, 3, 4, 5, &
            1, 2, 3, 4, 0, &
            1, 2, 3, 0, 0  &
        /), shape(ans))

        call make_banded_matrix(n, m, row, mat)

        equlity = abs(mat - ans) <= tolerance

        print *, 'Row: ', row

        do i = 1, m
            print *, equlity(i, :)
            do j = 1, n
                success = success .and. equlity(i, j)
            end do
        end do

        print *, 'Success: ', success
        print *
    end subroutine test_make_banded_matrix

    subroutine test_make_laplacian_o3(success)
        implicit none

        logical, intent(inout) :: success

        success = success .and. .true.

        print *, 'Testing make_laplacian_o3():', ' TODO: test body'
        print *, 'Success: ', success
        print *
    end subroutine test_make_laplacian_o3

    subroutine test_make_laplacian_o5(success)
        implicit none

        logical, intent(inout) :: success

        success = success .and. .true.

        print *, 'Testing make_laplacian_o5():'
        print *, 'Success: ', success
        print *
    end subroutine test_make_laplacian_o5

    subroutine test_make_laplacian_o7(success)
        implicit none

        logical, intent(inout) :: success

        success = success .and. .true.

        print *, 'Testing make_laplacian_o7():', ' TODO: test body'
        print *, 'Success: ', success
        print *
    end subroutine test_make_laplacian_o7

    subroutine test_make_laplacian(success)
        implicit none

        logical, intent(inout) :: success

        success = success .and. .true.

        print *, 'Testing test_make_laplacian():', ' TODO: test body'
        print *, 'Success: ', success
        print *
    end subroutine test_make_laplacian

    subroutine test_rgbmv(success)
        implicit none

        logical, intent(inout) :: success
        integer, parameter :: dp = selected_real_kind(15, 307)
        integer, parameter :: m = 5, n = 7, klu = 2
        real(dp), parameter :: sign = 1.0, tolerance = 1.0e-6
        real(dp), dimension(m) :: row
        real(dp), dimension(m, n) :: op
        real(dp), dimension(n, n) :: x, y, z
        logical, dimension(n) :: equlity
        integer :: i, j

        print *, 'Testing test_rgbmv()'

        row = (/ 1, 2, 3, 4, 5 /)
        x = reshape((/            &
            1, 0, 0, 0, 0, 0, 0, &
            0, 1, 0, 0, 0, 0, 0, &
            0, 0, 1, 0, 0, 0, 0, &
            0, 0, 0, 1, 0, 0, 0, &
            0, 0, 0, 0, 1, 0, 0, &
            0, 0, 0, 0, 0, 1, 0, &
            0, 0, 0, 0, 0, 0, 1  &
        /), shape(x))
        z = transpose(reshape((/   &
            3, 4, 5, 0, 0, 0, 0, &
            2, 3, 4, 5, 0, 0, 0, &
            1, 2, 3, 4, 5, 0, 0, &
            0, 1, 2, 3, 4, 5, 0, &
            0, 0, 1, 2, 3, 4, 5, &
            0, 0, 0, 1, 2, 3, 4, &
            0, 0, 0, 0, 1, 2, 3  &
        /), shape(z)))

        call make_banded_matrix(n, m, row, op)

        do i = 1, n
            call rgbmv(x(i, :), y(i, :), sign, op, klu, n)
        end do

        do i = 1, n
            equlity = abs(y(i, :) - z(i, :)) <= tolerance
            print *, 'Row #', i, equlity
            do j = 1, n
                success = success .and. equlity(j)
            end do
        end do

        print *, 'Success: ', success
        print *
    end subroutine test_rgbmv

    subroutine test_pumping(success)
        implicit none

        logical, intent(inout) :: success

        success = success .and. .true.

        print *, 'Testing test_pumping():', ' TODO: test body'
        print *, 'Success: ', success
        print *
    end subroutine test_pumping

    subroutine test_hamiltonian(success)
        implicit none

        logical, intent(inout) :: success

        success = success .and. .true.

        print *, 'Testing test_hamiltonian():', ' TODO: test body'
        print *, 'Success: ', success
        print *
    end subroutine test_hamiltonian

    subroutine test_runge_kutta(success)
        implicit none

        logical, intent(inout) :: success

        success = success .and. .true.

        print *, 'Testing test_runge_kutta():', ' TODO: test body'
        print *, 'Success: ', success
        print *
    end subroutine test_runge_kutta

    subroutine test_solve_nls(success)
        implicit none

        logical, intent(inout) :: success
        integer, parameter :: dp = selected_real_kind(15, 307)
        integer, parameter :: n = 7, iters = 1, order = 5
        complex(dp), dimension(n) :: u
        real(dp), parameter :: dt = 0.00001, dx = 0.01

        print *, 'Testing solve_nls():'

        call solve_nls(dt, dx, n, order, iters, u)

        success = success .and. .true.

        print *, 'Success: ', success
        print *
    end subroutine

end program test_nls