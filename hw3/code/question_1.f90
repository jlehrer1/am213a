Program Driver_LinAl

  use LinAl, only: dp, make_A_b, cholesky_backsubsitution, cholesky_factorization, qr_factorization, prettyprint

  implicit none
  integer :: m 
  integer :: n 
  real (dp), allocatable, dimension(:, :) :: A, A_ERR, Q_ERR, L, ata, atac, Vand1, Vand2 
  real (dp), allocatable, dimension(:, :) :: A_T 

  real (dp), allocatable, dimension(:, :) :: Q 
  real (dp), allocatable, dimension(:, :) :: R
  real (dp), allocatable, dimension(:) :: b, x, b_mod ! just a test 

  real (dp), dimension(4) :: vand1_x, vand1_b 
  real (dp), dimension(6) :: vand2_x, vand2_b 

  real (dp), dimension(4, 4) :: normal_vand1
  real (dp), dimension(6, 6) :: normal_vand2

  logical :: flag
  integer :: i 

  m = 21
  n = 1
  allocate(A(m, n), A_T(n, m), Q(m, m), R(m, n), x(n), b_mod(n), ata(n, n), atac(n, n), Q_ERR(m, m), b(m))

  allocate(Vand1(m, 4), Vand2(m, 6))
  A = 0.0 
  A_T = 0.0 
  Q = 0.0 
  R = 0.0 
  ata = 0.0 
  atac = 0.0

  ! For question 6
  call make_A_b(A, b)
  print *, 'A is '
  call prettyprint(A, m, n)
  print *, 'b is ', b

  print *, 'Setting up vandermonde for 3rd degree polynomial with intercept'
  Vand1(:, 1) = 1 

  do i=2, 4 
    Vand1(:, i) = A(:, 1)**(i-1)
  end do 

  call prettyprint(Vand1, m, 4)

  print *, 'Setting up Vandermonde for 5th degree polynomial with intercep'

  Vand2(:, 1) = 1
  do i=2, 6 
    Vand2(:, i) = A(:, 1)**(i-1)
  end do 

  call prettyprint(Vand2, m, 6)

  print *, 'Calculating normal equations and performing QR Decomposition for vand1'

  normal_vand1 = matmul(transpose(Vand1), Vand1) 
  vand1_b = matmul(transpose(Vand1), b)

  call cholesky_factorization(normal_vand1, flag)
  call cholesky_backsubsitution(normal_vand1, vand1_b, vand1_x)

  print *, 'Solution vector is', vand1_x
  print *, 'b is ', b 
  print *, 'err is ', matmul(Vand1, vand1_x) - b
  print *, 'error norm is', norm2(matmul(Vand1, vand1_x) - b, m)

  open (unit=10, file='vand1.txt')
  do i=1, 4
    write(10, *) vand1_x(i)
  end do 
  
  print *, 'Calcuating normal equations and performing QR Decomposition for vand2'

  normal_vand2 = matmul(transpose(Vand2), Vand2) 
  vand2_b = matmul(transpose(Vand2), b)

  call cholesky_factorization(normal_vand2, flag)
  call cholesky_backsubsitution(normal_vand2, vand2_b, vand2_x)

  print *, 'Solution vector is ', vand2_x
  print *, 'b is ', b 
  print *, 'err is ', matmul(Vand2, vand2_x) - b
  print *, 'error norm is', norm2(matmul(Vand2, vand2_x) - b, m)

  open (unit=10, file='vand2.txt')
  do i=1, 6
    write(10, *) vand2_x(i)
  end do

End Program Driver_LinAl
