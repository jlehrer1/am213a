module test 
    implicit none
    integer, save :: msize, nsize
    real, dimension(:,:), allocatable, save :: mat
  
contains 

subroutine testroutine(A, m, trace)
    integer, intent(in) :: m
    real, intent(in), dimension(m, m) :: A
    real, intent(out) :: trace 
    integer :: i 

    do i=1,m
        trace = trace + A(i, i)
    end do
end subroutine testroutine

subroutine periodic_boundry(N)
    ! Not technically necessary since the initial conditions are periodic
    ! But makes programming the finite differences easier since we don't have to
    ! Specify which initial condition was used 
    integer, intent(in) :: N 
end subroutine periodic_boundry

! subroutine trace(A, m, trace)
!     integer, intent(in) :: m
!     real, intent(in), dimension(m, m) :: A
!     real, intent(out) :: trace 
!     integer :: i

!     do i=1,m
!         trace = trace + A(i, i)
!     end do
! end subroutine trace 

end module test