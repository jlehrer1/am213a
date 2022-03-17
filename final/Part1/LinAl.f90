module LinAl
  implicit none
  integer, save :: msize, nsize
  integer, parameter :: dp = kind(0.d0)
contains

subroutine from_file(matrix, m, n)
  ! Reads a matrix from a file 
  ! First line of file must contain two integers m & n
  ! definining the number of rows and columns in the matrix, respectively.

  ! Parameters:
  ! filename: Path to file to read matrix from 
  ! matrix: allocatable output matrix to write into
  integer, intent(in) :: m, n 
  real (dp), intent(out), dimension(m, n) :: matrix 

  integer :: i, j

  open(10, file='dog_bw_data.dat')
  do i=1, m
    read(10, *) (matrix(i, j), j=1,n)
  end do 

end subroutine from_file

subroutine gauss_jacobi(A, m, x, b, max_iter, tolerance, filename) 
  integer, intent(in) :: m, max_iter 
  character(len=*), intent(in) :: filename
  real (dp), intent(in), dimension(m, m) :: A
  real (dp), intent(in), dimension(m) :: b 
  real (dp), intent(out), dimension(:) :: x 
  real (dp), intent(in) :: tolerance 

  real (dp), dimension(m) :: y
  real (dp) :: sum, r
  integer :: i, j, k 
  real (dp), dimension(max_iter) :: errors 

  errors = 0. 
  x = 1. 
  do k=1, max_iter 
    do i=1, m 
      sum = 0 
      do j=1, m 
        if (i .ne. j) then 
          sum = sum + A(i, j)*x(j)
        end if 
      end do 
      y(i) = 1./A(i,i)*(b(i) - sum)
    end do 

    r = norm2(matmul(A, x) - b)
    errors(k) = r
    if (r .le. tolerance) then 
      print *, 'Jacobi algorithm converged on iteration', i 
      exit
    end if
    x = y
  end do

  open(1, file=trim(filename))
  do i=1, max_iter
    write(1, *) (errors(i))
  end do 

  if (k .eq. max_iter + 1) then 
    print *, 'Convergence not reached with', max_iter, 'iterations'
  end if 

end subroutine gauss_jacobi

subroutine gauss_seidel(A, m, x, b, max_iter, tolerance, filename)
  integer, intent(in) :: m, max_iter 
  character(len=*), intent(in) :: filename
  real (dp), intent(in), dimension(m, m) :: A 
  real (dp), intent(in), dimension(m) :: b 
  real (dp), intent(out), dimension(:) :: x 
  real (dp), intent(in) :: tolerance 

  real (dp) :: sum1, r 
  integer :: i, j, k 
  real (dp), dimension(max_iter) :: errors 

  x = 1. 

  do k=1, max_iter 
    do i=1, m 
      sum1 = 0 
      do j=1, m
        if (j .ne. i) then 
          sum1 = sum1 + A(i, j)*x(j)
        end if 
      end do 
      x(i) = (b(i) - sum1) / A(i, i)
    end do 

    r = norm2(matmul(A, x) - b)
    errors(k) = r
    if (r < tolerance) then 
      print *, 'Convergence reached on iteration', k 
      exit 
    end if 
  end do 

  open(1, file=trim(filename))
  do i=1, k
    write(1, *) (errors(i))
  end do 

  if (k .eq. max_iter + 1) then 
    print *, 'Convergence not reached with', max_iter, 'iterations'
  end if 

end subroutine gauss_seidel

subroutine conjugate_gradient(A, m, x, b, max_iter, tolerance)
  integer, intent(in) :: m, max_iter 
  real (dp), intent(in) :: tolerance 
  real (dp), intent(in), dimension(m, m) :: A 
  real (dp), intent(in), dimension(m) :: b 
  real (dp), intent(out), dimension(:) :: x 
  
  real (dp), dimension(m) :: r, p 
  real (dp) :: alpha, denom, beta, err 
  real (dp), dimension(max_iter) :: errors 
  integer :: k 

  call frobenius_norm(A - transpose(A), m, m, err)
  if (err .ge. tolerance) then 
    print *, 'Error: A is not symmetric.'
    return 
  end if 

  x = 1.
  r = b - matmul(A, x)

  if (norm2(r) < tolerance) then
    print *, 'Convergence reached on iteration 0'
    return
  end if

  p = r
  k = 0
  do k=1, max_iter
    alpha = dot_product(r, r)/dot_product(p, matmul(A, p))
    ! print *, 'alpha is', alpha
    x = x + alpha*p

    ! Calculate denominator of beta before r is updated 
    denom = dot_product(r, r)
    ! Now continue and update r 
    r = r - alpha*matmul(A, p)
    ! print *, 'r is ', r 
    errors(k) = norm2(r)
    if (norm2(r) < tolerance) then 
      print *, 'Convergence reached on iteration', k 
      return
    end if

    beta = dot_product(r, r) / denom 
    p = r + beta*p 
  end do 
  print *, 'Warning: convergence not reached with', k, 'iterations'

end subroutine conjugate_gradient

subroutine preconditioned_conjugate_gradient(A, m, x, b, max_iter, tolerance)
  integer, intent(in) :: m, max_iter 
  real (dp), intent(in) :: tolerance 
  real (dp), intent(in), dimension(m, m) :: A
  real (dp), intent(in), dimension(m) :: b 
  real (dp), intent(out), dimension(:) :: x 
  
  real (dp), dimension(m) :: r, p 
  real (dp) :: alpha, denom, beta, err 
  real (dp), dimension(max_iter) :: errors 
  real (dp), dimension(m, m) :: precond, precond_inv
  integer :: k 

  call frobenius_norm(A - transpose(A), m, m, err)
  if (err .ge. tolerance) then 
    print *, 'Error: A is not symmetric.'
    return 
  end if 

  x = 1.
  r = b - matmul(A, x)
  
  if (norm2(r) < tolerance) then
    print *, 'Convergence reached on iteration 0'
    return
  end if

  p = r
  k = 0
  do k=1, max_iter
    alpha = dot_product(r, r)/dot_product(p, matmul(A, p))
    ! print *, 'alpha is', alpha
    x = x + alpha*p

    ! Calculate denominator of beta before r is updated 
    denom = dot_product(r, r)
    ! Now continue and update r 
    r = r - alpha*matmul(A, p)
    ! print *, 'r is ', r 
    errors(k) = norm2(r)
    if (norm2(r) < tolerance) then 
      print *, 'Convergence reached on iteration', k 
      return
    end if

    beta = dot_product(r, r) / denom 
    p = r + beta*p 
  end do 
  print *, 'Warning: convergence not reached with', k, 'iterations'

end subroutine preconditioned_conjugate_gradient

subroutine generate_test_matrix(A, m, val)
  real (dp), intent(in) :: val 
  integer, intent(in) :: m 
  real (dp), intent(inout), dimension(m, m) :: A 

  integer :: i 

  A = 1. 
  do i=1,m 
    A(i, i) = val 
  end do 
end subroutine generate_test_matrix

subroutine prettyprint(A, m, n)
  ! Prints a 2D array in a human-readable way 

  ! Parameters:
  ! A: Matrix to print 
  ! m: Number of rows in matrix 
  ! n: Number of columns in matrix 

  integer, intent(in) :: m, n
  real (dp), intent(in), dimension(:, :) :: A
  integer :: i, j

  do i=1,m
      write(*,"(100g15.5)") ( A(i,j), j=1,n )
  enddo
  print *, ''
end subroutine prettyprint
subroutine printvec(x, m)
  integer, intent(in) :: m 
  real (dp), intent(in), dimension(m) :: x 

  integer :: i 

  do i=1,m 
    write(*,"(100g15.5)") x(i)
  end do 
end subroutine printvec

subroutine frobenius_norm(A, m, n, norm)
  ! Calculates the Frobenius norm of a matrix A 

  ! Parameters:
  ! A: Matrix to calculate norm of 
  ! m: Number of rows in A 
  ! n: Number of columns in A 
  ! norm: Output variable to store Frobenius norm of A to 

  integer, intent(in) :: m, n 
  real (dp), intent(in), dimension(m, n) :: A 
  real (dp), intent(out) :: norm 

  integer :: i, j 

  norm = 0.
  do i=1,m 
    do j=1,n
      norm = norm + abs(A(i, j))
    end do 
  end do 

  norm = sqrt(norm)
end subroutine frobenius_norm

end module LinAl

