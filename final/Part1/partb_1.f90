program gauss_iterative_methods 

    use LinAl, only: dp, prettyprint, from_file, gauss_jacobi, gauss_seidel, generate_test_matrix

    real (dp), dimension(10, 10) :: A 
    real (dp), dimension(10) :: b, x 
    real (dp) :: val 
    real (dp), dimension(5) :: vals 
    real (dp), allocatable, dimension(:) :: errors 
    integer :: i, max_iter, m 
    real (dp) :: tolerance 
    character(len=4096) :: filename

    tolerance = 10.**(-5)
    print *, 'Tolerance is ', tolerance 
    
    max_iter = 500 
    m = 10
    allocate(errors(max_iter))

    vals = (/2., 5., 10., 100., 1000./)

    do i=1, m
        b(i) = i
    end do 

    do i=1, 5
        val = vals(i)
        write(filename, "(A5,I2)") "GJ_err", i
    
        call generate_test_matrix(A, 10, val)
        call gauss_jacobi(A, m, x, b, max_iter, tolerance, 'errors/' // trim(filename))

        print *, 'Solution reached with Gauss-Jacobi for D=', val, 'is'
        print *, x 
    end do 

    ! Print a couple newlines to make this more readable
    print *, ''
    print *, ''

    do i=1, 5
        val = vals(i)
        write(filename, "(A5,I2)") "GS_err", i

        call generate_test_matrix(A, 10, val)
        call gauss_seidel(A, m, x, b, max_iter, tolerance, 'errors/' // trim(filename))

        print *, 'Solution reached with Gauss-Jacobi for D=', val, 'is'
        print *, x 
    end do

    A = 1.
    do i=1, m
        A(i, i) = i
    end do

    call gauss_jacobi(A, m, x, b, max_iter, tolerance, 'GJ_err_i_matrix.dat')
    call gauss_seidel(A, m, x, b, max_iter, tolerance, 'GS_err_i_matrix.dat')
    
end program gauss_iterative_methods