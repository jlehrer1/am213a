Program Driver_LinAl

  use LinAl, only: dp, from_file, cholesky_backsubsitution, cholesky_factorization, qr_factorization, prettyprint

  implicit none
  
  !m=21, n=2
  integer :: m 
  integer :: n 
  real (dp), allocatable, dimension(:, :) :: A, ERR, L, ata
  real (dp), allocatable, dimension(:, :) :: A_T 

  real (dp), allocatable, dimension(:, :) :: Q 
  real (dp), allocatable, dimension(:, :) :: R
  real (dp), allocatable, dimension(:) :: b, x ! just a test 
  logical :: flag 

  m = 21
  n = 2
  allocate(A(m, n), A_T(n, m), Q(m, m), R(m, n), b(m), x(n))
  A = 0.0 
  A_T = 0.0 
  Q = 0.0 
  R = 0.0 
  b=1.0 !test 
  ! For question 6
  call from_file(A)
  call prettyprint(A, m, n)

  call qr_factorization(A, R, Q, m, n)
  print *, 'A-QR is'
  
  ERR = A - matmul(Q, R)
  call prettyprint(ERR, m, n)
  call prettyprint(A, m, n)

  ata = matmul(transpose(A), A)
  a_t = transpose(A)
  b = matmul(a_t, b)

  call prettyprint(ata, n, n)
  call cholesky_factorization(ata, flag)
  call cholesky_backsubsitution(ata, b, x)

  print *, matmul(ata, x)
  ! print *, shape(transpose(A)), shape(A)

End Program Driver_LinAl
