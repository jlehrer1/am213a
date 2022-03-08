module LinAl
  implicit none
  integer, save :: msize, nsize
  ! integer, parameter :: dp = SELECTED_REAL_KIND(15) ! to use double precision throughout 
  integer, parameter :: dp = kind(0.d0)
  real (dp), dimension(:, :), allocatable, save :: mat
contains

subroutine from_file(matrix, m, n)
  ! Reads a matrix from a file 
  ! First line of file must contain two integers m & n, 
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

subroutine gauss_jacobi(A, m, x, b, errors, max_iter) 
  integer, intent(in) :: m, max_iter 
  real (dp), intent(in), dimension(m, m) :: A 
  real (dp), intent(in), dimension(m) :: b 
  real (dp), intent(out), dimension(:) :: x 
  real (dp), intent(inout), dimension(max_iter) :: errors 

  real (dp), dimension(m) :: y
  real (dp) :: sum, r
  integer :: i, j, k 

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
    if (r .le. 1e-6) then 
      print *, 'Jacobi algorithm converged on iteration', i 
      return
    end if
    x = y

  end do 
  print *, 'Convergence not reached with', max_iter, 'iterations'
end subroutine gauss_jacobi

subroutine gauss_seidel(A, m, x, b, errors, max_iter)
  integer, intent(in) :: m, max_iter 
  real (dp), intent(in), dimension(m, m) :: A 
  real (dp), intent(in), dimension(m) :: b 
  real (dp), intent(out), dimension(:) :: x 
  real (dp), intent(inout), dimension(max_iter) :: errors 

  real (dp) :: sum1, r 
  real (dp), dimension(m, m) :: L, U, D 
  integer :: i, j, k 

  D = 0. 
  L = 0. 
  U = 0. 

  do i=1, m 
    do j=1, m 
      if (i .le. j) then 
        U(i, j) = A(i, j)
      end if 

      if (i .ge. j) then 
        L(i, j) = A(i, j)
      end if 

      if (i .eq. j) then 
        D(i, j) = A(i, j)
      end if 
    end do 
  end do 

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
    if (r < 1e-6) then 
      print *, 'Convergence reached on iteration', k 
      return 
    end if 
  end do 
  print *, 'Convergence not reached with', max_iter, 'iterations'
end subroutine gauss_seidel

subroutine conjugate_gradient(A, m, x, b, max_iter, eps)
  real (dp), intent(in) :: eps 
  integer, intent(in) :: m, max_iter 
  real (dp), intent(in), dimension(m, m) :: A 
  real (dp), intent(in), dimension(m) :: b 
  real (dp), intent(out), dimension(:) :: x 

  real (dp), dimension(m) :: r, p 
  real (dp) :: alpha, denom, beta  
  integer :: k 

  r = b - matmul(A, x)

  if (norm2(r) < eps) then 
    x = r 
    return 
  end if 

  p = r 
  k = 0 

  do k=1, max_iter
    alpha = dot_product(r, r)/dot_product(p, matmul(A, p))
    x = x + alpha*p 
    ! Calculate denominator of beta before r is updated 
    denom = dot_product(r, r)
    ! Now continue and update r 
    r = r - alpha*matmul(A, p)

    if (norm2(r) < eps) then 
      print *, 'Convergence reached on iteration', k 
      x = r
      return  
    end if 

    beta = dot_product(r, r) / denom 
    p = r + beta*p 
  end do 
end subroutine conjugate_gradient

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

end module LinAl

