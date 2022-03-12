Program question1

  use LinAl, only: dp, prettyprint, from_file, frobenius_norm
  
  implicit none
  real (dp), allocatable, dimension(:, :) :: A, U, V_T, recon, full_sigma
  real (kind=dp) :: dummy(1, 1)
  real (kind=dp), allocatable, dimension(:) :: sigma, work, temp
  integer, dimension(9) :: ks
  real (dp), dimension(9) :: errors 
  real (dp) :: err 
  integer :: m, n, lwork, info, i, k, j, l, p
  character(len=4096) :: filename

  m = 3355
  n = 5295
  ks = (/20, 40, 80, 160, 320, 640, 1280, 2560, 3355/)

  print *, 'Finding optimal work size'
  lwork = -1
  allocate(A(m, n), sigma(m), U(m, m), V_T(n, n))
  call from_file(A, m, n)
  call dgesvd('A', 'A', m, n, A, m, sigma, U, m, V_T, n, dummy, lwork, info)

  print *, 'Calculating SVD for image'
  lwork = nint(dummy(1, 1))
  allocate(work(lwork))
  call dgesvd('A', 'A', m, n, A, m, sigma, U, m, V_T, n, work, lwork, info)

  allocate(full_sigma(m, n), recon(m, n), temp(m))
  do i=1,9 
    temp = 0. 
    full_sigma = 0. 
    recon = 0.

    ! Construct full sigma matrix for multiplications
    k = ks(i) 
    temp(1:k) = sigma(1:k)

    do j=1,m
      full_sigma(j, j) = temp(j)
    end do

    print *, 'Reconstructing image with', k, 'singular values'
    recon = matmul(U, matmul(full_sigma, V_T))

    print *, 'Calculating error'
    call frobenius_norm(A-recon, m, n, err)
    print *, 'Error is', err
    errors(i) = err/real(m*n)

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
