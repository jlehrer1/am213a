Program Driver_LinAl

  use LinAl, only: dp, from_file, cholesky_backsubsitution, cholesky_factorization, qr_factorization, prettyprint

  implicit none
  
  !m=21, n=2
  integer :: m 
  integer :: n 
  real (dp), allocatable, dimension(:, :) :: A, A_ERR, Q_ERR, L, ata, atac
  real (dp), allocatable, dimension(:, :) :: A_T 

  real (dp), allocatable, dimension(:, :) :: Q 
  real (dp), allocatable, dimension(:, :) :: R
  real (dp), allocatable, dimension(:) :: b, x, b_mod ! just a test 
  logical :: flag 

  m = 21
  n = 2
  allocate(A(m, n), A_T(n, m), Q(m, m), R(m, n), x(n), b_mod(n), ata(n, n), atac(n, n), Q_ERR(m, m))

  A = 0.0 
  A_T = 0.0 
  Q = 0.0 
  R = 0.0 
  ata = 0.0 
  atac = 0.0

  ! For question 6
  call from_file(A)
  call prettyprint(A, m, n)

  call qr_factorization(A, R, Q, m, n)
  print *, 'A-QR is'
  
  A_ERR = A - matmul(Q, R)
  Q_ERR = matmul(transpose(Q), Q)
  call prettyprint(A_ERR, m, n)
  call prettyprint(Q_ERR, m, m)
  ata = matmul(transpose(A), A)
  atac = ata 
  a_t = transpose(A)
  b_mod = 1.

  call prettyprint(ata, n, n)

  call cholesky_factorization(ata, flag)
  call cholesky_backsubsitution(ata, b_mod, x)

  print *, matmul(atac, x)
  print *, b_mod
  ! print *, shape(transpose(A)), shape(A)

End Program Driver_LinAl
