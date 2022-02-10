Program Driver_LinAl

  use LinAl, only: mat, msize, nsize, readMat, prettyprint, n_norm, matrix_trace, dp, from_file, gaussian_elimination, backsolve

  implicit none
  
  character(len=100) :: myFileName, file1, file2 
  integer :: i,j
  real (dp) :: traceA, normval 
  real (dp), allocatable, dimension(:) :: vec
  real (dp), allocatable, dimension(:, :) :: A, B, X, A_inv, T
  logical :: flag

  myFileName = 'Amat.dat'
  call from_file(myFileName, A)

  myFileName = 'Bmat.dat'
  call from_file(myFileName, B)

  call prettyprint(A, 4, 4)
  call prettyprint(B, 4, 6)

  call matrix_trace(A, 4, traceA)

  ! allocate to be the same size as B
  allocate(X(4, 6))
  allocate(A_inv(4, 4))
  allocate(T(4,4))

  X = 0.0
  T=0.0

  print *, 'The trace of A is', traceA

  do i=1, 4
    vec = A(:, i)
    call n_norm(vec, 4, normval)
    print *, 'The norm of column', i, 'is', normval 
  end do 
  
  call gaussian_elimination(A, B, flag, 4, 6)

  print *, 'A after gaussian elimination is'
  call prettyprint(A, 4, 4)

  print *, 'B after gaussian elimination is'
  call prettyprint(B, 4, 6)

  call backsolve(A, B, X, 4)
  print *, 'Solution matrix is'
  call prettyprint(X, 4, 6)

  ! print *, 'B after gaussian elimination is'
  ! call prettyprint(B)
  T = matmul(A, X)
  print *, 'AX is '
  call prettyprint(T, 4, 4)
  ! deallocate(mat)

End Program Driver_LinAl
