subroutine cumul(a,b,c,n)
      
  implicit none       
  integer :: k,i,n
  real(8),dimension(n) :: a
  real(8),dimension(2*n) :: b,c
 
  k=0
  do i=1,n
    k=k+2
  
    b(k)=(b(k)+b(k-1)**2)*c(k)+a(i)**2
    b(k-1)=b(k-1)*c(k-1)+a(i)
  
    c(k-1)=c(k-1)+1.d0
    c(k)=c(k)+1.d0
  
    b(k-1)=b(k-1)/c(k-1)
    b(k)=b(k)/c(k)-(b(k-1)**2)
  
  enddo

  return

end subroutine
