subroutine zeroav(iblock)
      
  use md_variables
      
  integer :: i, iblock
      
  if(iblock == 0) then
! initialization
    allocate(ipot(4))
 
    ipot(1)=1
    ipot(2)=ipot(1)+2
    ipot(3)=ipot(2)+2
    ipot(4)=ipot(3)+2
    ikin=ipot(4)+2
    nav=ikin+2
    if(nav .gt. mav) stop 'zeroav: nav.gt.mav'

    allocate(av(nav))
    allocate(avp(nav))
    allocate(anorm(nav))
    allocate(anormp(nav))
    do i=1,nav
      av(i)=0.d0
      anorm(i)=0.d0
    enddo
  else
    do i=1,nav
      avp(i)=0.d0
      anormp(i)=0.d0
    enddo
  endif

  return

end subroutine
