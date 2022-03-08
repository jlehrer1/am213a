Program question1

  use LinAl, only: dp, hessenberg, prettyprint, QR_with_shift, QR_without_shift, inverse_iteration, gaussian_elimination
  
  implicit none
  real (dp), dimension(4, 4) :: A_1, A_3, A_3LU, eye 
  real (dp), dimension(3, 3) :: A_2, A_2C
  real (dp), dimension(4) :: v, b, x
  integer :: m = 4
  integer :: i 
  integer :: max_iter = 10000
  integer, dimension(4) :: s
  logical :: flag = .false. 

  ! Define A_1 
  A_1(1,:) = (/5,4,1,1/)
  A_1(2, :) = (/4,5,1,1/)
  A_1(3, :) = (/1,1,4,2/)
  A_1(4, :) = (/1,1,2,4/)

  ! Define A_2 
  A_2(1, :) = (/3,1,0/)
  A_2(2, :) = (/1,2,1/)
  A_2(3, :) = (/0,1,1/)

  !Define A_3
  A_3(1, :) = (/2,1,3,4/)
  A_3(2, :) = (/1,-3,1,5/)
  A_3(3, :) = (/3,1,6,-2/)
  A_3(4, :) = (/4,5,-2,-1/)

  print *, 'Before converting to Hessenberg form, A_1 is'
  call prettyprint(A_1, m, m)
  call hessenberg(A_1, m)
  print *, 'After converting to Hessenberg form, A_1 is'
  call prettyprint(A_1, m, m)

  A_2C = A_2 
  m = 3
  print *, 'Printing A_2'
  call prettyprint(A_2C, m, m)

  call QR_with_shift(A_2C, m, 10)

  print *, 'Printing result of QR alg. with shift'
  call prettyprint(A_2C, m, m)

  call QR_without_shift(A_2, m, max_iter)
  print *, 'Printing result of QR alg. without shift'

  call prettyprint(A_2, m, m)
  m = 4

  print *, 'Printing eigenvector associated with \lambda=-8.0286'
  call inverse_iteration(A_3, v, m, -8.0286, max_iter)
  print *, 'V is ', v

  print *, 'Printing eigenvector associated with \lambda=7.9329'
  call inverse_iteration(A_3, v, m, 7.9329, max_iter)
  print *, 'V is ', v  

  print *, 'Printing eigenvector associated with \lambda=5.6689'
  call inverse_iteration(A_3, v, m, 5.6689, max_iter)
  print *, 'V is ', v  

  print *, 'Printing eigenvector associated with \lambda=-1.5732'
  call inverse_iteration(A_3, v, m, -1.5732, max_iter)
  print *, 'V is ', v
  
End Program question1
