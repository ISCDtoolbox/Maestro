subroutine setvel(iff)

  use md_variables

  integer :: i,k,l,iff
  real(8), dimension(:), allocatable :: v2,f
  real(8) :: dummy

! Multiplication by 'renyi_order' added by Miha Srdinsek
  allocate(v2(nbeadMD * renyi_order),f(nbeadMD * renyi_order))

  vcm=0.d0
  velocity=0.d0
  vel=0.d0
  v2=0.d0


! Multiplication by 'renyi_order' added by Miha Srdinsek
  do k=1,(nbeadMD * renyi_order)
    do i=1,n
      if (iff==0) then
        do l=1,ndimMD
          call random_number(dummy)
          vel(l,i,k)=dummy-0.5d0
        enddo
      endif
    enddo
  enddo

! Multiplication by 'renyi_order' added by Miha Srdinsek
  do k=1,(nbeadMD * renyi_order)
    do i=1,n
      do l=1,ndimMD
        vcm(l,k)=vcm(l,k)+(amas(indx(i))*vel(l,i,k))/mtot
      enddo
    enddo

    if(fixcm) then
      do i=1,n
        do l=1,ndimMD
          vel(l,i,k)=vel(l,i,k)-vcm(l,k)
        enddo
      enddo
    endif


    do l=1,ndimMD
      do i=1,n
        v2(k)=v2(k)+amas(indx(i))*vel(l,i,k)**2
      enddo
    enddo
    f(k)=sqrt(gMD*tfakeMD/v2(k))
    do i=1,n
      do l=1,ndimMD
        vel(l,i,k)=vel(l,i,k)*f(k)
      enddo
    enddo

    !!! check for the center of mass velocity
    do l=1,ndimMD
      vcm(l,k)=0.d0
      do i=1,n
        vcm(l,k)=vcm(l,k)+amas(indx(i))*vel(l,i,k)/mtot
      enddo
    enddo

  enddo

  deallocate(v2,f)

  return

end subroutine setvel
