subroutine propX(tau)

  use md_variables
  implicit none

  integer :: ind,i,l,k
  real(8) :: tau

  do k=1,nbeadMD 
    do i = 1, N
      ind=indx(i)
        do l=1,ndimMD 
          rpos(l,i,k) = rpos(l,i,k) + tau*pimp(l,i,k)/amas(ind) 
        enddo 
    end do
  end do 

  if(fixcm) then
    do k=1,nbeadMD
      do l=1,ndimMD
        rcm(l,k)=0.d0
        do i=1,n
          rcm(l,k)=rcm(l,k)+amas(indx(i))*rpos(l,i,k)
        enddo
        rcm(l,k)=rcm(l,k)/mtot
!       do i=1,n
!        rpos(l,i)=rpos(l,i)-rcm(l)
!       enddo
      enddo
    enddo
  endif

  return

end subroutine
