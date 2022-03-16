Program question1

  use LinAl, only: dp, prettyprint, from_file
  
  implicit none
  real (kind=dp), dimension(4, 4) :: test, test_U, test_V
  real (kind=dp), dimension(4) :: test_sigma
  real (kind=dp), allocatable, dimension(:) :: work 
  integer :: m, n, lwork, info, i

  m = 4 
  n = 4

  ! Initialize simple matrix for testing 
  test = 0.
  do i=1,m 
    test(i, i) = i
  end do 
  allocate(work(4))
  call dgetrf(m, n, test, m, work, info)

  print *, test 

End Program question1
