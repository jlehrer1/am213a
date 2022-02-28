module LinAl
  implicit none
  integer, save :: msize, nsize
  integer, parameter :: dp = SELECTED_REAL_KIND(15) ! to use double precision throughout 
  real (dp), dimension(:,:), allocatable, save :: mat

contains

subroutine from_file(matrix)
  ! Reads a matrix from a file 
  ! First line of file must contain two integers m & n, 
  ! definining the number of rows and columns in the matrix, respectively.

  ! Parameters:
  ! filename: Path to file to read matrix from 
  ! matrix: allocatable output matrix to write into
  real (dp), intent(out), dimension(21, 2) :: matrix 

  integer :: i, j

  open(10, file='least_squares_data.dat')
  do i=1, 21
    read(10, *) (matrix(i, j), j=1,2)
  end do 

end subroutine from_file 

subroutine make_A_b(A, b)
  real (dp), intent(out), dimension(21, 1) :: A, b
  real (dp), dimension(21, 2) :: matrix 

  integer :: i, j

  open(10, file='least_squares_data.dat')
  do i=1, 21
    read(10, *) (matrix(i, j), j=1,2)
  end do 

  A(:, 1) = matrix(:, 1)
  b(:, 1) = matrix(:, 2)
  
end subroutine make_A_b

subroutine cholesky_factorization(A, flag) 
  real (dp), intent(inout), dimension(:, :) :: A 
  logical, intent(out) :: flag 

  integer :: i, j, k

  do j=1,ubound(A, 1)
    do k=1, j-1 
      A(j, j) = A(j, j) - A(j, k)**2 
    end do 

    A(j, j) = sqrt(A(j, j))
    do i=j+1, ubound(A, 1)
      do k=1, j-1 
        A(i, j) = A(i, j) - A(i, k)*A(j, k)
      end do 
      if (A(j, j) .eq. 0.0) then 
        print *, 'Error, matrix is singular. Cannot continue'
        flag = .true.
        return 
      end if 

      ! Now we know it's safe to divide by A(j, j)
      A(i, j) = A(i, j) / A(j, j)
    end do 
  end do 
end subroutine cholesky_factorization

subroutine cholesky_backsubsitution(A, b, x)
  real (dp), intent(in), dimension(:, :) :: A 
  real (dp), intent(in), dimension(:) ::  b 
  real (dp), intent(inout), dimension(:) :: x 

  integer :: i, j, k 
  real (dp) :: sum 
  real (dp), dimension(ubound(A, 1)) :: y
  y = 0.0

  ! Calculate Ly = b 
  do i=1,ubound(A, 1)
    sum = b(i) 
    do j=1,i-1 
      sum = sum - y(j)*A(i, j)
    end do 
    y(i) = sum / A(i, i)
  end do 

  ! Calculate L^* x = y 
  do i=ubound(A, 1), 1, -1
    if (A(i, i) .eq. 0.0) then 
      print *, 'Error, matrix is sinular'
      return 
    end if 
    do k=i+1, ubound(A, 1)
      y(i) = y(i) - A(k, i)*x(k)
    end do 
    x(i) = y(i)/A(i, i)
  end do
end subroutine cholesky_backsubsitution

subroutine qr_factorization(A, R, Q, m, n)
  ! Implementation of QR Decomposition via householder
  integer, intent(in) :: m, n 
  real (dp), intent(in), dimension(m, n) :: A
  real (dp), intent(inout), dimension(m, m) :: Q
  real (dp), intent(inout), dimension(m, n) :: R

  real (dp), dimension(m) :: v_j
  real (dp), dimension(m, m) :: eye 
  real (dp), dimension(m, m) :: outer_product
  real (dp) :: s_j
  integer :: j, k, i

  ! Setup identity matrix for householder reflections 

  eye = 0.0
  do i=1, m
    eye(i, i) = 1.0
  end do 

  R = A 
  Q = eye
  do j=1, n
    v_j = 0.0
    s_j = 0.0
    
    do k=j, m
      s_j = s_j + R(k, j)**2
    end do 
    s_j = sqrt(s_j)

    if (A(j, j) < 0.0) then
      s_j = -s_j
    end if

    v_j(j) = A(j, j) + s_j
    v_j(j+1:) = A(j+1:, j)
    v_j = v_j / norm2(v_j)

    do i=1,m
      do k=1, m
        outer_product(i, k) = v_j(i)*v_j(k)
      end do
    end do

    Q = matmul(Q, eye - 2*outer_product)
    R = R - 2*matmul(outer_product, R)
  end do
end subroutine qr_factorization

subroutine backsubstitution(U, x, b, m)
  ! Performs backsubstitution for Ux = b

  ! Parameters:
  ! U: Upper triangular matrix 
  ! x: m vector, solutions are stored here
  ! b: m vector, RHS of Ux = b 
  ! m: number of rows in U 

  real (dp), intent(in), dimension(:, :):: U 
  real (dp), intent(in), dimension(:) :: b 
  real (dp), intent(inout), dimension(:) :: x 
  integer, intent(in) :: m 

  integer :: i, k 
  real (dp) :: sum

  if (U(m, m) .eq. 0.0) then
    print *, 'Error: U is singular'
    return 
  end if 

  x(m) = b(m)/U(m, m)
  do i=m-1, 1, -1 
    if (U(i, i) .eq. 0.0) then 
      print *, 'Entry at ', i, 'is zero, matrix is singular'
      exit 
    end if 

    sum = 0.0 
    do k=i+1, m 
      sum = sum + U(i, k)*x(k)
    end do 
    x(i) = (b(i) - sum)/U(i, i)
  end do 
end subroutine backsubstitution

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

subroutine frobenius_norm(A, m, n, norm)
  integer, intent(in) :: m, n 
  real (dp), intent(in), dimension(:, :) :: A 
  real (dp), intent(out) :: norm 

  integer :: i, j 

  norm = 0

  do i=1,m 
    do j=1,n
      norm = norm + abs(A(i, j))
    end do 
  end do 

  norm = sqrt(norm)
end subroutine frobenius_norm

end module LinAl
