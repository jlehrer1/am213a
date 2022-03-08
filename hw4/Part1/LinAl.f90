module LinAl
  implicit none
  integer, save :: msize, nsize
  integer, parameter :: dp = SELECTED_REAL_KIND(6) ! to use double precision throughout 
  real (dp), dimension(:, :), allocatable, save :: mat
contains

subroutine qr_factorization(A, R, Q, m, n)
  ! Implementation of QR Decomposition via Householder reflections 
  
  ! Parameters:
  ! A: Input matrix to factorize 
  ! R: Matrix to write upper triangular part of QR factorization to 
  ! Q: Matrix to write orthogonal part Q of QR factorization to 
  ! m: Number of rows in A 
  ! n: Number of columns in A

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

    if (norm2(v_j) .ne. 0.0) then 
      v_j = v_j / norm2(v_j)
    end if

    do i=1,m
      do k=1, m
        outer_product(i, k) = v_j(i)*v_j(k)
      end do
    end do

    Q = matmul(Q, eye - 2*outer_product)
    R = R - 2*matmul(outer_product, R)
  end do
end subroutine qr_factorization

subroutine hessenberg(A, m)
  ! Converts a matrix into it's Hessenberg form via modified Householder reflections 
  ! Parameters:
  ! A: Input matrix, modified to Hessenberg form 
  ! m: Number of rows/columns of A

  real (dp), intent(inout), dimension(:, :) :: A 
  integer, intent(in) :: m

  integer :: i, j, k
  real (dp), dimension(m) :: v_j
  real (dp), dimension(m, m) :: outer_product
  real (dp) :: s_j
  ! Define first standard basis vector for algorithm

  do j=1, m - 2
    v_j = 0.0
    s_j = 0.0
    
    do i=j+1, m
      s_j = s_j + A(i, j)**2
    end do 
    s_j = sqrt(s_j)

    ! s_j = sign(A_{j+1, j})*s_j
    if (A(j+1, j) < 0.0) then
      s_j = -s_j
    end if

    v_j(j+1) = A(j+1, j) + s_j
    v_j(j+2:) = A(j+2:, j)
    v_j = v_j / norm2(v_j)

    do i=1,m
      do k=1, m
        outer_product(i, k) = v_j(i)*v_j(k)
      end do
    end do

    A = A - 2*matmul(outer_product, A)
    A = A - 2*matmul(A, outer_product)
  end do
end subroutine hessenberg

subroutine QR_without_shift(A, m, max_iter)
  ! Reveals the eigenvalues of A via the QR algorithm without shift. 

  ! Parameters:
  ! A: Input matrix. A is transformed to an upper triangular matrix with eigenvalues along the diagonal 
  ! m: Number of rows/columns of A 
  ! max_iter: Max number of iterations that the QR algorithm is allowed to run 

  real (dp), intent(inout), dimension(:, :) :: A 
  integer, intent(in) :: m, max_iter 
  
  real (dp), dimension(m, m) :: Q
  real (dp), dimension(m, m) :: R
  real (dp) :: eps = 1e-6
  real (dp), dimension(m) :: diag, prev_diag 
  integer :: i, k

  do i=1, m 
    prev_diag(i) = A(i, i)
  end do 

  k = 0 
  do while (k .le. max_iter)
    call qr_factorization(A, R, Q, m, m)
    A = matmul(R, Q)

    do i=1, m 
      diag(i) = A(i, i)
    end do 
  
    if (norm2(diag - prev_diag) .le. eps) then 
      print *, 'Convergence criterion reached, exiting...'
      print *, 'Finished in', k, 'iterations'
      return 
    end if 
    prev_diag = diag 
    k = k + 1
  end do 

  print *, 'QR with shift: Max iterations reached without convergence criterion'
end subroutine QR_without_shift

subroutine QR_with_shift(A, m, max_iter)
  ! Reveals the eigenvalues of A via the QR algorithm without shift 

  ! Parameters:
  ! A: Input matrix, transformed to an upper triangular matrix with the eigenvalues of A along the diagonal 
  ! m: Number of rows/columns of A 
  ! max_iter: Maximum number of iterations that the algorithm is allowed to run 
  real (dp), intent(inout), dimension(:, :) :: A 
  integer, intent(in) :: m, max_iter 
  
  real (dp), dimension(m, m) :: Q, eye 
  real (dp), dimension(m, m) :: R 
  integer :: i, k
  real (dp) :: eps = 1e-5
  real (dp) :: mu 
  real (dp), dimension(m) :: diag, prev_diag 

  eye = 0.0
  do i=1,m 
    eye(i, i) = 1.
    diag(i) = A(i, i)
  end do 

  k = 0
  do while (k .le. max_iter)
    ! Run alg.
    mu = A(m, m)
    call qr_factorization(A - mu*eye, R, Q, m, m)

    A = matmul(R, Q) + mu*eye 
    ! Calculate error 
    do i=1, m 
      diag(i) = A(i, i)
    end do

    print *,'error is', norm2(diag - prev_diag)
    if (norm2(diag - prev_diag) .le. eps) then 
      print *, 'Convergence criterion reached, exiting...'
      print *, 'Finished in', k, 'iterations'
      return 
    end if 

    prev_diag = diag 
    k = k + 1 
  end do 

  print *, 'QR with shift: Max iterations reached without convergence criterion'
end subroutine QR_with_shift

subroutine inverse_iteration(A, x, m, mu, max_iter)
  ! Calculates the eigenvector associated with a given eigenvalue 
  
  ! Parameters:
  ! A: Input matrix to calculate eigenvector for 
  ! x: Output vector to write estimation of eigenvector to 
  ! m: Number of rows/columns of A 
  ! mu: Estimation of eigenvalue of A 
  ! max_iter: Maximum number of iterations the inverse iteration algorithm is allowed to run for 
  integer, intent(in) :: m, max_iter
  real (dp), intent(in), dimension(m, m) :: A 
  real (dp), intent(out), dimension(m) :: x
  real, intent(in) :: mu

  integer :: k
  real (dp) :: eps = 1e-6
  real (dp), dimension(m) :: r, y
  integer (dp), dimension(m) :: s
  real (dp), dimension(m, m) :: B, eye 
  logical :: flag = .false.

  ! initialize I_m 
  eye = 0. 
  do k=1,m 
    eye(k, k) = 1. 
  end do 

  ! Initialize B 
  B = A - mu*eye

  x = 0.
  x(1) = 1.
  k = 0

  do while (k .le. max_iter)
    call gaussian_elimination(B, x, flag, m, m)
    call single_backsolve(B, y, x, m)

    y = y / norm2(y)

    if (norm2(y - x) .le. eps) then 
      print *, 'Convergence reached on iteration', k 
      return
    end if

    x = y
    k = k + 1
  end do
  print *, 'Convergence not reached'
end subroutine inverse_iteration


subroutine gaussian_elimination(A, b, flag, m, n)
  ! Gaussian elimination operation on A

  ! Parameters:
  ! A: Matrix corresponding to coefficients of linear equations 
  ! B: Matrix of RHS b's, [b_1 ... b_n]
  ! flag: Boolean, .true. if matrix is singular, .false. otherwise
  ! m: Number of rows & columns in A
  ! n: Number of columns in B

  real (dp), intent(inout), dimension(:, :) :: A
  real (dp), intent(inout), dimension(:) :: b
  logical, intent(out) :: flag 
  integer, intent(in) :: m, n

  real (dp) :: p, divisor
  real (dp), allocatable, dimension(:) :: temp, r ! keeps vector when we're doing row swaps 
  integer :: i, j, K

  allocate(temp(n))
  allocate(r(n))

  temp = 0.0
  r = 0.0 
  flag = .false. !begin by assuming matrix is not singular

  do j=1, m - 1 
    ! print *, "On step", j, "of Gaussian elimination, the matrix looks like this"
    ! call prettyprint(A, m, m)
    ! find maximum pivot 
    call find_pivot(A, j, K, p, m)
    
    ! if matrix is singular exit 
    if (A(j, j) .eq. 0.0) then 
      print *, 'ERROR: Matrix is singular'
      flag = .true.
      return 
    end if 

    if (K .ne. j) then !swap the rows 
      temp = A(K, :)
      A(K, :) = A(j, :)
      A(j, :) = temp 
    end if 

    ! print *, "On step", j, "of Gaussian elimination, after partial pivoting the matrix looks like this"
    ! call prettyprint(A, m, m)

    do i=j+1, m 
      divisor =  A(i, j)/A(j, j) !since this gets overwritten
      A(i, :) = A(i, :) - A(j, :)*divisor
      b(i) = b(i) - b(j)*divisor
    end do 
  end do
  ! if we made it here then the matrix is nonsingular
  flag = .false.
end subroutine gaussian_elimination

subroutine find_pivot(A, j, K, p, m)
  ! Find the index K and pivot P such that 
  ! p = max_{k=j...,msize} |a_{kj}|
  ! This will be used in gaussian elimination with partial pivoting 

  ! Parameters:
  ! A: Matrix to find max pivot of 
  ! j: current row index in Gaussian elimination 
  ! K: output index for pivot row 
  ! p: output value for maximum pivot 
  ! m: Number of rows of A 

  real (dp), intent(in), dimension(:, :) :: A 
  integer, intent(in) :: j, m
  
  integer, intent(out) :: K 
  real (dp), intent(out) :: p

  integer :: i 
  p = 0.0 
  
  do i=j,m
    if (abs(A(i, j)) > p) then 
      p = abs(A(i, j))
      K = i 
    end if
  end do

end subroutine find_pivot 

subroutine single_backsolve(U, x, b, m)
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
end subroutine single_backsolve

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

