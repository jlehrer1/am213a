Program question1

  use LinAl, only: dp, prettyprint, from_file, frobenius_norm, printvec

  implicit none
  real (dp), allocatable, dimension(:, :) :: A, U, V_T, recon, full_sigma, error_mat, copy 
  real (dp), allocatable, dimension(:) :: sigma, work, temp
  real (dp) :: dummy(1, 1)

  integer, dimension(9) :: ks
  real (dp), dimension(9) :: errors 
  real (dp) :: err
  integer :: m, n, lwork, info, i, j, l, p, k
  character(len=4096) :: filename

  m = 3355
  n = 5295
  ks = (/20, 40, 80, 160, 320, 640, 1280, 2560, 3355/)

  print *, 'Finding optimal work size'
  lwork = -1
  allocate(A(m, n), sigma(m), U(m, m), V_T(n, n), error_mat(m, n), copy(m, n))

  call from_file(A, m, n)
  copy = A 

  call dgesvd('A', 'A', m, n, A, m, sigma, U, m, V_T, n, dummy, lwork, info)

  print *, 'Calculating SVD for image'
  lwork = nint(dummy(1, 1))
  allocate(work(lwork))
  call dgesvd('A', 'A', m, n, A, m, sigma, U, m, V_T, n, work, lwork, info)
  print *, 'Info is', info 

  print *, 'Sigma vector is'
  call printvec(sigma, m)

  allocate(full_sigma(m, n), recon(m, n), temp(m))
  do i=1, 9 
    full_sigma = 0. 
    recon = 0.

    ! Construct full sigma matrix for multiplications
    do j=1, ks(i)
      full_sigma(j, j) = sigma(j)
    end do
    
    print *, 'Reconstructing image with', ks(i), 'singular values'
    recon = matmul(U, matmul(full_sigma, V_T))

    print *, 'Calculating error'
    error_mat = copy - recon
    call frobenius_norm(error_mat, m, n, err)

    errors(i) = err/real(m*n)

    if (errors(i) .le. 10.**(-3)) then 
      print *, 'Error less than 10^(-3) at with', ks(i), 'singular values'
    end if 

    print *, 'Writing recontructed image to file'

    write(filename, "(A5,I2)") "compressed_image", i
    open(1, file=trim(filename))
    do l=1,m 
      write(1, *) (recon(l, p), p=1,n)
    end do 
  end do

  print *, 'Writing errors to file'
  open(1, file='frobenius_errors.dat')
  do i=1,9 
    write(1, *) errors(i)
  end do

End Program question1
