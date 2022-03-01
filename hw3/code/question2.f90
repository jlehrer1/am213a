Program question2

  use LinAl, only: dp, from_file, qr_factorization, prettyprint, frobenius_norm, make_A_b, backsubstitution

  implicit none
  
  integer :: m
  integer :: n
  real (dp), allocatable, dimension(:, :) :: A, A_ERR, Q_ERR, eye, A_ERR_2
  real (dp), allocatable, dimension(:) :: b, b_mod, x ! just a test 
  real (dp) :: norm
  integer :: i

  real (dp), dimension(4) :: vand1_x, vand1_b 
  real (dp), dimension(6) :: vand2_x, vand2_b

  real (dp), dimension(21, 4) :: R_vand1, Vand1
  real (dp), dimension(21, 6) :: R_vand2, Vand2

  real (dp), dimension(21, 21) :: Q_vand1
  real (dp), dimension(21, 21) :: Q_vand2

  m = 21
  n = 4

  allocate(A(m, n), eye(m, m), b(m), Q_ERR(m, m), A_ERR(m, n), A_ERR_2(m, 6))
  eye = 0 

  do i=1,m 
    eye(i, i) = 1. 
  end do 

  ! For question 6
  call make_A_b(A, b)

  print *, 'Setting up vandermonde for 3rd degree polynomial with intercept'
  Vand1(:, 1) = 1 

  do i=2, n 
    Vand1(:, i) = A(:, 1)**(i-1)
  end do 

  call prettyprint(Vand1, m, n)
  call qr_factorization(Vand1, R_vand1, Q_vand1, m, n)
  
  A_ERR = Vand1 - matmul(Q_vand1, R_vand1)

  print *, 'A-QR is'
  call prettyprint(A_ERR, m, n)
  call frobenius_norm(A_ERR, m, n, norm)
  print *, 'Frobenius norm of A-QR is ', norm

  Q_ERR = matmul(transpose(Q_vand1), Q_vand1) - eye
  print *, 'Q^TQ - I is'
  call prettyprint(Q_ERR, m, m)

  call frobenius_norm(Q_ERR, m, m, norm)
  print *, 'Frobenius norm of Q^TQ - I is', norm

  print *, 'Setting up Vandermonde for 5th degree polynomial with intercept'

  Vand2(:, 1) = 1
  do i=2, 6 
    Vand2(:, i) = A(:, 1)**(i-1)
  end do 

  n=6 
  call prettyprint(Vand2, m, n)
  call qr_factorization(Vand2, R_vand2, Q_vand2, m, n)
  
  A_ERR = Vand2 - matmul(Q_vand2, R_vand2)

  print *, 'A-QR is'
  call prettyprint(A_ERR, m, n)
  call frobenius_norm(A_ERR, m, n, norm)
  print *, 'Frobenius norm of A-QR is ', norm

  Q_ERR = matmul(transpose(Q_vand2), Q_vand2) - eye
  print *, 'Q^TQ - I is'
  call prettyprint(Q_ERR, m, m)

  call frobenius_norm(Q_ERR, m, m, norm)
  print *, 'Frobenius norm of Q^TQ - I is', norm

  ! open (unit=10, file='vand1.txt')
  ! do i=1, 4
  !   write(10, *) vand2_x(i)
  ! end do

  ! open (unit=10, file='vand2.txt')
  ! do i=1, 6
  !   write(10, *) vand2_x(i)
  ! end do

End Program question2
