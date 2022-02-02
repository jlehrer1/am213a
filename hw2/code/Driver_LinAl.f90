Program Driver_LinAl

  use LinAl, only: mat, msize, nsize, readMat

  implicit none
  
  character(len=100) :: myFileName
  integer :: i,j
  real :: traceA

  
  myFileName = 'Amat.dat'

  open(10,file=myFileName)
  read(10,*) msize,nsize
  close(10)

  allocate(mat(msize,nsize))
  ! Always initialize with zeros
  mat = 0.0

  
  
  call readMat(myFileName)

  do i = 1, msize
     write(*,*) (mat(i,j) , j = 1, nsize )
  end do
  

  deallocate(mat)





End Program Driver_LinAl
