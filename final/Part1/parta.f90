Program question1

  use LinAl, only: dp, prettyprint, from_file
  
  implicit none
  real (dp), allocatable, dimension(:, :) :: A, U, V_T, recon, full_sigma
  real (kind=dp) :: dummy(1, 1)
  real (kind=dp), allocatable, dimension(:) :: sigma, work, temp 
  integer, dimension(9) :: ks
  integer :: m, n, lwork, info, i, k, j, l, p
  character(len=4096) :: filename

  m = 3355
  n = 5295
  ks = (/20, 40, 80, 160, 320, 640, 1280, 2560, 3355/)


  do i=1, 9 
    k = ks(i)
    write(filename, "(A5,I2)") "compressed_image", i
    print *, 'filename is'
    print *, trim(filename)

  end do 

  print *, 'Finding optimal work size'
  lwork = -1
  allocate(A(m, n), sigma(m), U(m, m), V_T(n, n))
  call from_file(A, m, n)
  call dgesvd('A', 'A', m, n, A, m, sigma, U, m, V_T, n, dummy, lwork, info)

  print *, 'Calculating SVD for image'
  lwork = nint(dummy(1, 1))
  allocate(work(lwork))
  call dgesvd('A', 'A', m, n, A, m, sigma, U, m, V_T, n, work, lwork, info)

  allocate(full_sigma(m, n), recon(m, n))
  do i=1,9 
    temp = 0. 
    full_sigma = 0. 
    recon = 0.

    ! Construct full sigma matrix for multiplications
    k = ks(i) 
    temp = sigma(0:k)
    
    do j=1,m
      full_sigma(j, j) = temp(j)
    end do

    print *, 'Reconstructing image with', k, 'singular values'
    recon = matmul(U, matmul(full_sigma, V_T))
    
    print *, 'Writing recontructed image to file'
    write(filename, "(A5,I2)") "compressed_image", i
    open(1, file=trim(filename), status='new')
    do l=1,m 
      write(1, *) (recon(l, p), p=1,n)
    end do 
  end do 
  
End Program question1
