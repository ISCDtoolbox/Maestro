subroutine cumul1(a,b,c,n)

  implicit none

  integer :: n,i
  real(8),dimension(n) :: a,b,c

  do i=1,n
    b(i)=b(i)*c(i)+a(i)
    c(i)=c(i)+1.d0
    b(i)=b(i)/c(i)
  enddo

  return

end subroutine
