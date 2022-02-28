Program question_2

  use LinAl, only: dp, from_file, qr_factorization, prettyprint, frobenius_norm, make_A_b

  implicit none
  
  integer :: m 
  integer :: n 
  real (dp), allocatable, dimension(:, :) :: A, A_ERR, Q_ERR, L, ata, atac, eye
  real (dp), allocatable, dimension(:, :) :: A_T 

  real (dp), allocatable, dimension(:, :) :: Q 
  real (dp), allocatable, dimension(:, :) :: R
  real (dp), allocatable, dimension(:) :: b, x ! just a test 
  real (dp) :: norm
  integer :: i 

  m = 21
  n = 1
  allocate(A(m, n), A_T(n, m), Q(m, m), R(m, n), x(n), b(n), ata(n, n), atac(n, n), Q_ERR(m, m), eye(m, m))

  A = 0.0 
  A_T = 0.0 
  Q = 0.0 
  R = 0.0 
  ata = 0.0 
  atac = 0.0
  eye = 0 

  do i=1,m 
    eye(i, i) = 1. 
  end do 

  ! For question 6
  call make_A_b(A, b)
  call prettyprint(A, m, n)

  call qr_factorization(A, R, Q, m, n)
  
  A_ERR = A - matmul(Q, R)
  Q_ERR = matmul(transpose(Q), Q) - eye 

  print *, 'A-QR is'
  call prettyprint(A_ERR, m, n)

  call frobenius_norm(A_ERR, m, n, norm)

  print *, 'Frobenius norm of A-QR is ', norm

  print *, 'Q^TQ - I is'
  call prettyprint(Q_ERR, m, m)

  call frobenius_norm(Q_ERR, m, m, norm)
  print *, 'Frobenius norm of Q^TQ - I is', norm 

  

End Program question_2
