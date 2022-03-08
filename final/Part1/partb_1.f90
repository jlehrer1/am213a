program gauss_iterative_methods 

    use LinAl, only: dp, prettyprint, from_file, gauss_jacobi, gauss_seidel, generate_test_matrix

    real (dp), dimension(10, 10) :: A 
    real (dp), dimension(10) :: b, x 
    real (dp) :: val 
    real (dp), dimension(5) :: vals 
    real (dp), allocatable, dimension(:) :: errors 
    integer :: i, max_iter, m 

    max_iter = 1000 
    m = 10
    allocate(errors(max_iter))

    vals = (/2., 5., 10., 100., 1000./)

    do i=1, 10
        b(i) = i
    end do 

    ! do i=1, 5
    !     val = vals(i)
    !     call generate_test_matrix(A, 10, val)
    !     call gauss_jacobi(A, m, x, b, errors, max_iter)

    !     print *, 'Solution reached with Gauss-Jacobi for D=', val, 'is'
    !     print *, x 
    ! end do 

    do i=1, 5
        val = vals(i)
        call generate_test_matrix(A, 10, val)
        print * ,'max_iter is', max_iter
        call gauss_seidel(A, m, x, b, errors, max_iter)

        print *, 'Solution reached with Gauss-Jacobi for D=', val, 'is'
        print *, x 
    end do 
    
end program gauss_iterative_methods