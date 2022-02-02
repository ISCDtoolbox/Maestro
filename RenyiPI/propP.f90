subroutine propP( tau, enkin)

  use md_variables
  
  real(8) :: tau,scaleV,scaleF,enkin
  real(8) :: dts
  real(8), parameter :: eps=1.d-5
  integer :: i,l,ind,k
        
  if(nh) then
    dts=tau*PS
    scaleV=exp(-dts)
    if(abs(dts).gt.eps) then
      scaleF=(exp(dts)-1.d0)/PS
    else
      scaleF=tau*(1.d0+0.5d0*tau*PS*(1.0d0+tau*PS/3.d0))
    endif

    do k=1,nbeadMD
      do i = 1, N
        ind = indx(i)
        do l=1,ndimMD
          pimp(l,i,k) = scaleV*(pimp(l,i,k) + forceMD(l,i,k)*scaleF)
        enddo
      enddo
    enddo
  else
    do k=1,nbeadMD
      do i = 1, N
        ind=indx(i)
        do l=1,ndimMD
          pimp(l,i,k) = pimp(l,i,k)+tau*forceMD(l,i,k)
        enddo
      enddo
    enddo
  endif

  enkin=0.d0
  vcm=0.d0

  do k=1,nbeadMD
    do i=1,n
      ind=indx(i)
      do l=1,ndimMD
        vel(l,i,k) = pimp(l,i,k)/amas(ind)
        enkin=enkin+amas(ind)*vel(l,i,k)**2
      enddo
    enddo
  enddo
  enkin=0.5d0*enkin/nbeadMD
  
  if(fixcm) then
    do k=1,nbeadMD
      do l=1,ndimMD
        vcm(l,k)=0.d0
        do i=1,n
          vcm(l,k)=vcm(l,k)+amas(indx(i))*vel(l,i,k)
        enddo
        vcm(l,k)=vcm(l,k)/mtot
      enddo

    enddo
  endif

  return

end subroutine
