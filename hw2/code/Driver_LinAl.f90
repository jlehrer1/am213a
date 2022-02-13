Program Driver_LinAl

  use LinAl, only: readMat, prettyprint, n_norm, matrix_trace, dp, from_file, gaussian_elimination, backsolve, &
  lu_decomp, LU_backsolve

  implicit none
  
  character(len=100) :: myFileName 
  integer :: i 
  real (dp) :: traceA, normval 
  real (dp), allocatable, dimension(:) :: vec, x_lu

  integer, allocatable, dimension(:) :: s ! permutation vector for LU decomp
  real (dp), allocatable, dimension(:, :) :: A, B, X, T
  logical :: flag

  integer :: m, n 
  m = 4 
  n = 6

  myFileName = 'Amat.dat'
  call from_file(myFileName, A)

  myFileName = 'Bmat.dat'
  call from_file(myFileName, B)

  print *, 'A is '
  call prettyprint(A, m, m)

  print *, 'B is '
  call prettyprint(B, m, n)

  call matrix_trace(A, m, traceA)
  print *, 'The trace of A is', traceA

  ! allocate to be the same size as B
  allocate(X(m, n))
  allocate(T(m, m))
  allocate(vec(m))

  X = 0.0
  T = 0.0
  vec = 0.0 

  do i=1, m
    vec = A(:, i)
    call n_norm(vec, m, normval)
    print *, 'The norm of column', i, 'is', normval 
  end do 
  
  print *, 'A before Gaussian elimination is'
  call prettyprint(A, m, m)

  call gaussian_elimination(A, B, flag, m, n)

  print *, 'A after gaussian elimination is'
  call prettyprint(A, m, m)

  print *, 'B after gaussian elimination is'
  call prettyprint(B, m, n)

  call backsolve(A, B, X, m, n)
  print *, 'Solution matrix is'

  call prettyprint(X, m, n)

  ! print *, 'B after gaussian elimination is'
  ! call prettyprint(B)
  T = matmul(A, X)

  print *, 'AX is '
  call prettyprint(T, m, n)

  print *, 'B is '
  call prettyprint(B, m, n)
  ! deallocate(mat)

  print *, 'Error matrix from Gaussian elimination is '
  call prettyprint(T - B, m, n)

  ! DONE WITH GAUSSIAN ELIM, RESETTING FOR LU
  print *, 'Resetting matrices A & B'
  myFileName = 'Amat.dat'
  call from_file(myFileName, A)

  myFileName = 'Bmat.dat'
  call from_file(myFileName, B)
  ! DONE RESETTING 

  ! LU Decomp.
  ! QUETSTION 4
  allocate(s(m))
  allocate(x_lu(m))
  s = 0
  T = A ! make a copy to keep in LU

  print *, 'Printing matrix A after reset'
  call prettyprint(A, m, m)

  print *, 'Matrix A after LU Decomposition is '
  call lu_decomp(T, m, flag, s)

  print *, 'A after LU decomposition (held as LU) is '
  call prettyprint(T, m, m)

  do i=1, n 
    vec = B(:, i)
    x_lu = B(:, i)
    call LU_backsolve(T, m, vec, s, x_lu)
    print *, 'Solution vector at ', i, 'using LU is '
    print *, x_lu
  end do 

  ! Need to reset A since it is overwritten with LU
  myFileName = 'Amat.dat'
  call from_file(myFileName, A)

  print *, 'Ax is '
  print *, matmul(A, x_lu)

  print *, 'And b is '
  print *, vec

End Program Driver_LinAl
