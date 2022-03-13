program conjugate_gradient_solver

    use LinAl, only: dp, prettyprint, from_file, conjugate_gradient, generate_test_matrix, printvec

    real (dp), dimension(10, 10) :: A 
    real (dp), dimension(10) :: b, x 
    real (dp) :: val 
    real (dp), dimension(5) :: vals 
    real (dp), allocatable, dimension(:) :: errors 
    integer :: i, max_iter, m 
    real (dp) :: tolerance

    tolerance = 10.**(-5)
    max_iter = 10
    m = 10
    allocate(errors(max_iter))

    vals = (/2., 5., 10., 100., 1000./)

    do i=1, 10
        b(i) = real(i)
    end do

    do i=1, 5
        val = vals(i)
        call generate_test_matrix(A, m, val)
        call conjugate_gradient(A, m, x, b, max_iter, tolerance)

        print *, 'Solution reached with Conjugate Gradient for D=', val, 'is'
        call printvec(x, m)
    end do 
    
end program conjugate_gradient_solver