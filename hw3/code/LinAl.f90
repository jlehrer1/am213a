module LinAl
  implicit none
  integer, save :: msize, nsize
  integer, parameter :: dp = SELECTED_REAL_KIND(15) ! to use double precision throughout 
  real (dp), dimension(:,:), allocatable, save :: mat

contains

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

subroutine qr_factorization(A, m, n)
  ! Implementation of QR Decomposition via householder
  real (dp), intent(inout), dimension(:, :) :: A 
  integer, intent(in) :: m, n 

  real (dp), dimension(m) :: v_j 
  integer :: j, k 

  do j=1,n 

  end do 

end subroutine qr_factorization
end module LinAl
