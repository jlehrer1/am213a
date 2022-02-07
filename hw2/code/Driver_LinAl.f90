Program Driver_LinAl

  use LinAl, only: mat, msize, nsize, readMat, prettyprint, n_norm, matrix_trace, dp, from_file, gaussian_elimination

  implicit none
  
  character(len=100) :: myFileName, file1, file2 
  integer :: i,j
  real (dp) :: traceA, normval 
  real (dp), allocatable, dimension(:) :: vec
  real (dp), allocatable, dimension(:, :) :: A, B
  logical :: flag 

  myFileName = 'Amat.dat'
  call from_file(myFileName, A)

  myFileName = 'Bmat.dat'
  call from_file(myFileName, B)

  call prettyprint(A, 4, 4)
  call prettyprint(B, 4, 6)

  call matrix_trace(A, 4, traceA)

  print *, 'The trace of A is', traceA

  do i=1, 4
    vec = A(:, i)
    call n_norm(vec, 4, normval)
    print *, 'The norm of column', i, 'is', normval 
  end do 
  
  call gaussian_elimination(A, B, flag, 4, 6)

  print *, 'A after gaussian elimination is'
  call prettyprint(A, 4, 4)

  ! print *, 'B after gaussian elimination is'
  ! call prettyprint(B)

  ! deallocate(mat)

End Program Driver_LinAl
