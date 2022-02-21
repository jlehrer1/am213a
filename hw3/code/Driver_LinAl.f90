Program Driver_LinAl

  use LinAl, only: readMat, prettyprint, two_norm, matrix_trace, dp, from_file, gaussian_elimination, backsolve, &
  lu_decomp, LU_backsolve, single_backsolve

  implicit none
  
  character(len=100) :: myFileName 
  integer :: i 
  real (dp) :: traceA, normval 
  real (dp), allocatable, dimension(:) :: vec, x_lu

  integer, allocatable, dimension(:) :: s ! permutation vector for LU decomp
  real (dp), allocatable, dimension(:, :) :: A, B, X, T
  logical :: flag

  ! For question 6
  real (dp), dimension(4,4) :: plane 
  real (dp), dimension(4,1) :: b_plane
  real (dp), dimension(4,1) :: x_plane 
  integer, dimension(4) :: s_plane ! permutation for s 

  real (dp) :: pi = acos(-1.0_8)
  real (dp) :: e = dexp(1.0_8)

  ! For array dimensions 
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
    call two_norm(vec, m, normval)
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

    print *, 'Error vector is ', matmul(A, x_lu) - vec

    print *, 'Error norm at', i, 'is'
    call two_norm(matmul(A, x_lu) - vec, m, normval)
    print *, normval
  end do 

! End of QUESTION 5.
! QUESTION 6: FITTING A PLANE

plane(1, :) = (/1,2,3,1/)
plane(2, :) = (/-3,2,5,1/)
plane(3, :) = (/pi, e, -sqrt(2.), 1./)
plane(4, :) = (/4, 0,-2, 1/)

b_plane = 0.0
x_plane = 0.0 
s_plane = 0 

m = 4
print *, 'Equations to solve are'
call prettyprint(plane, m, m)

print *, 'Using GE, we have that the matrix becomes the upper triangular' 
call gaussian_elimination(plane, b_plane, flag, m, 1)

print *, b_plane
call prettyprint(plane, m, m)
print *, 'And with backsubstitution we have that our solution is'

call backsolve(A, b_plane, x_plane, 4, 1)

print *, x_plane

End Program Driver_LinAl
