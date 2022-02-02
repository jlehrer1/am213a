module LinAl
  implicit none
  integer, save :: msize, nsize
  real, dimension(:,:), allocatable, save :: mat

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

subroutine matrix_trace(A, m, trace)
  integer, intent(in) :: m
  real, intent(in), dimension(m, m) :: A
  real, intent(out) :: trace 
  integer :: i 

  do i=1,m
      trace = trace + A(i, i)
  end do
end subroutine matrix_trace

subroutine n_norm(vec, m, norm)
    integer, intent(in) :: m
    real, intent(in), dimension(m) :: vec 
    real, intent(out) :: norm 

    integer :: i

    do i=1,m
      norm = norm + abs(vec(i)) ** 2
    end do
    norm = sqrt(norm)

end subroutine n_norm

subroutine print_matrix(A, m, n):
  integer, intent(in) :: m 
  integer, intent(in) :: n 
  real, intent(in) :: A(m, n)


end subroutine print_matrix

end module LinAl
