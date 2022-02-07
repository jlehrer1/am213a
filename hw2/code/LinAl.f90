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
  integer, intent(in) :: m
  real (dp), intent(in), dimension(m, m) :: A
  real (dp), intent(out) :: trace 
  integer :: i 

  do i=1,m
      trace = trace + A(i, i)
  end do

end subroutine matrix_trace

subroutine n_norm(vec, m, norm)
    integer, intent(in) :: m
    real (dp), intent(in), dimension(m) :: vec 
    real (dp), intent(out) :: norm 

    integer :: i
    norm = 0.0

    do i=1,m
      norm = norm + vec(i) ** 2
    end do
    norm = sqrt(norm)

end subroutine n_norm

subroutine prettyprint(A, m, n)
  integer, intent(in) :: m, n
  real (dp), intent(in), dimension(:, :) :: A
  integer :: i, j

  do i=1,m
      write(*,"(100g15.5)") ( A(i,j), j=1,n )
  enddo

  ! write(*,"100g15.5") ( matrix(i,j), j=1,n )
  ! do i=1,m 
  !   print *, A(i, :)
  ! end do 

  print *, ''
end subroutine prettyprint

subroutine find_pivot(A, j, K, p, m)
  ! Find the index K and pivot P such that 
  ! p = max_{k=j...,msize} |a_{kj}|

  ! This will be used in gaussian elimination with partial pivoting 
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
  real (dp), intent(inout), dimension(:, :) :: A, B 
  logical, intent(out) :: flag 
  integer, intent(in) :: m, n

  real (dp) :: p, maxval
  real (dp), allocatable, dimension(:) :: temp, r ! keeps vector when we're doing row swaps 
  integer :: i, j, K

  allocate(temp(n))
  allocate(r(n))

  temp = 0.0
  r = 0.0 
  flag = .false. !begin by assuming matrix is not singular

  do j=1, m - 1 
    ! find maximum pivot 
    call find_pivot(A, j, K, p, m)
    
    ! if matrix is singular exit 
    if (A(j, j) .eq. 0.0) then 
      print *, 'ERROR: Matrix is singular'
      flag = .true.
      exit 
    end if 

    if (K .ne. j) then !swap the rows 
      temp = A(K, :)
      A(K, :) = A(j, :)
      A(j, :) = temp 
    end if 

    do i=j+1, m 
      A(i, :) = A(i, :) - A(i,j)*A(j, :)/A(j,j)
      B(i, :) = B(i, :) - A(i, j)*B(j, :)/A(j, j)
    end do 
  end do
  print *, 'hi'
end subroutine gaussian_elimination

subroutine lu_decomp(A, m, flag, s)
  integer, intent(in) :: m 
  real (dp), intent(in), dimension(m, m) :: A 

  logical, intent(out) :: flag 
  integer, intent(out), dimension(m) :: s 

end subroutine lu_decomp

subroutine backsolve(U, B, X)
  real (dp), intent(in), dimension(:, :) :: U, B
  real (dp), intent(out), dimension(:, :) :: X
  real (dp) :: sum 
  integer :: i, k, j

  x = 0.0

  ! Instead of Ux = b, we have UX = B
  do j=1, nsize
    do i=msize-1, 1, -1 
      sum = 0.0 
      do k=i+1, msize 
        sum = sum + U(i, k)*X(j, k)
      end do 
    end do 
    X(i, j) = B(i, j) / U(i, i)
  end do
end subroutine backsolve

end module LinAl
