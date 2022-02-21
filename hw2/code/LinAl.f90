module LinAl
  implicit none
  integer, save :: msize, nsize
  integer, parameter :: dp = SELECTED_REAL_KIND(15) ! to use double precision throughout 
  real (dp), dimension(:,:), allocatable, save :: mat

contains

! ********************************************************
subroutine readMat(filename)
  character(len=*) :: filename

  integer :: i,j

  ! Reads a file containing the matrix A 
  ! Sample file:
  !
  ! 4 4 
  ! 2.0 1.0 1.0 0.0
  ! 4.0 3.0 3.0 1.0
  ! 8.0 7.0 9.0 5.0
  ! 6.0 7.0 9.0 8.0
  !
  ! Note that the first 2 numbers in the first line are the matrix dimensions, i.e., 4x4,
  ! then the next msize lines are the matrix entries. This matrix is found in Eq. 2.18 of the lecture note.
  ! Note that entries must be separated by a tab.

  open(10,file=filename)

  ! Read the matrix dimensions
  read(10,*) i,j

  ! Read matrix
  do i=1,msize
      read(10,*) ( mat(i,j), j=1,nsize )
  enddo

  close(10)
  
end subroutine readMat

subroutine from_file(filename, matrix)
  ! Reads a matrix from a file 
  ! First line of file must contain two integers m & n, 
  ! definining the number of rows and columns in the matrix, respectively.

  ! Parameters:
  ! filename: Path to file to read matrix from 
  ! matrix: allocatable output matrix to write into

  character(len=*) :: filename 
  real (dp), intent(out), allocatable, dimension(:, :) :: matrix 
  integer :: i, j, m, n

  open(10,file=filename)
  read(10,*) m, n
  close(10)

  open(10, file=filename)
  read(10,*) i, j

  allocate(matrix(m, n))

  do i=1,m
    read(10, *) (matrix(i, j), j=1,n)
  end do

end subroutine from_file 

subroutine matrix_trace(A, m, trace)
  ! Calculates the trace of a matrix A 

  ! Parameters:
  ! A: Matrix to calculate norm of 
  ! m: Number of rows in A 
  ! trace: output to store trace value in 

  integer, intent(in) :: m
  real (dp), intent(in), dimension(m, m) :: A
  real (dp), intent(out) :: trace 
  integer :: i 

  trace = 0.0
  do i=1,m
      trace = trace + A(i, i)
  end do

end subroutine matrix_trace

subroutine two_norm(vec, m, norm)
  ! Calculates the 2 norm of a vector 

  ! Parameters:
  ! vec: Vector to calculate norm of 
  ! m: Length of vector 
  ! norm: Output to store 2 norm value in 
  integer, intent(in) :: m
  real (dp), intent(in), dimension(m) :: vec 
  real (dp), intent(out) :: norm 


  integer :: i
  norm = 0.0

  do i=1,m
    norm = norm + vec(i) ** 2
  end do
  norm = sqrt(norm)

end subroutine two_norm

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

subroutine gaussian_elimination(A, B, flag, m, n)
  ! Gaussian elimination operation on A

  ! Parameters:
  ! A: Matrix corresponding to coefficients of linear equations 
  ! B: Matrix of RHS b's, [b_1 ... b_n]
  ! flag: Boolean, .true. if matrix is singular, .false. otherwise
  ! m: Number of rows & columns in A
  ! n: Number of columns in B

  real (dp), intent(inout), dimension(:, :) :: A, B 
  logical, intent(out) :: flag 
  integer, intent(in) :: m, n

  real (dp) :: p
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
      r = A(i, j)*B(j, :)/A(j, j) !since this gets overwritten
      A(i, :) = A(i, :) - A(i, j)*A(j, :)/A(j, j)
      B(i, :) = B(i, :) - r
      print *, 'B(i, :) is ', B(i, :)
    end do 
  end do
  ! if we made it here then the matrix is nonsingular
  flag = .false.
end subroutine gaussian_elimination

subroutine lu_decomp(A, m, flag, s)
  ! LU Decomposition on A

  ! Parameters:
  ! A: Matrix of coefficients corresponding to linear equations 
  ! m: Number of rows & columns of A 
  ! flag: Boolean, indicated if A is singular 
  ! s: Vector of length m, equivalent to the permutation matrix P
  integer, intent(in) :: m 
  real (dp), intent(inout), dimension(m, m) :: A 

  logical, intent(inout) :: flag 
  integer, intent(inout), dimension(m) :: s 

  integer :: i, j, K, k2, permtemp ! permutation temp variable for swaps 
  real (dp) :: p 
  real (dp), dimension(m) :: temp ! vector for swaps 

  temp = 0.0
  s = 0.0 

  ! Initialize permutation vector 
  do j=1,m 
    s(j) = j
  end do 

  ! Perform LU loop
  do j=1, m 

    ! Find index k and pivot p
    call find_pivot(A, j, K, p, m)

    ! Swap rows K and j or A, swap entries K and j of s
    if (K .ne. j) then 
      temp = A(K, :)
      A(K, :) = A(j, :)
      A(j,: ) = temp 

      permtemp = s(K)
      s(K) = s(j)
      s(j) = permtemp
    end if 
  
    if (A(j,j) .eq. 0.0) then 
      print *, 'Error: Pivot is zero'
      flag = .true.
      return 
    end if 

    do i=j+1, m 
      A(i, j) = A(i, j) / A(j, j)

      do k2=j+1,m 
        A(i, k2) = A(i, k2) - A(i, j)*A(j, k2)
      end do 
    end do
  end do
  flag = .false.
end subroutine lu_decomp

subroutine LU_backsolve(A, m, b, s, x)
  ! Backsubsitution for LU 

  ! Parameters:
  ! A: LU decomposition, stored in the format of 2.55 from the text 
  ! m: Number of rows & columns of A 
  ! b: RHS of (LU)x=b, m vector 
  ! s: Permutation vector from LU decomposition 
  ! x: Output vector, where the solutions are stored 

  ! In this case A is LU where the diagonal is l_{ii}
  integer, intent(in) :: m 
  real (dp), intent(inout), dimension(m, m) :: A 
  real (dp), intent(inout), dimension(m) :: b
  real (dp), intent(inout), dimension(m) :: x 

  integer, intent(inout), dimension(m) :: s 

  real (dp), dimension(m) :: y 
  integer :: i, j
  
  x = 0.0 
  y = 0.0 

  do i=1, m
    y(i) = b(s(i))
  end do 

  ! forward substitution for y=L^{-1}Pb 
  do j=1, m-1
    do i=j+1, m 
      y(i) = y(i) - y(j)*A(i, j)
    end do
  end do 

  call single_backsolve(A, x, y, m)
end subroutine LU_backsolve 

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

subroutine backsolve(U, B, X, m, n)
  ! Performs backsubstitution UX=B on Gaussian eliminated matrix U

  ! Parameters:
  ! U: Upper triangular matrix 
  ! B: Matrix of RHS vectors 
  ! X: Matrix where solutions are written to 
  ! m: Number of rows of U 
  ! n: Number of columns of B 

  real (dp), intent(in), dimension(:, :) :: U, B
  real (dp), intent(inout), dimension(:, :) :: X
  integer, intent(in) :: m, n
  integer :: i 

  real (dp), dimension(m) :: x_temp, b_temp  
  X = 0.0 

  do i=1, n
    x_temp = X(:, i)
    b_temp = B(:, i)

    call single_backsolve(U, x_temp, b_temp, m)

    X(:, i) = x_temp
  end do 
end subroutine backsolve

end module LinAl
