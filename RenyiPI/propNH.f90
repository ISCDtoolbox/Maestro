SUBROUTINE propNH( tau, enh)

  use md_variables
  real(8) :: tau, enh,pps,sumv2
  integer :: i,l,ind 
       
  sumv2 = 0.d0
       
  do k=1,nbeadMD 
    do i = 1, N
      ind = indx(i)
      do l=1,ndimMD
        sumv2 = sumv2 + amas(ind)*vel(l,i,k)**2
      enddo
    enddo
  enddo 

!  sumv2=sumv2/(gMD*nbeadMD)
!  S = S + 0.5d0*tau*PS
!  pps=(sumv2/tempMD-1.d0)/Q   
!  PS = PS + tau * pps
!  S = S + 0.5d0*tau*PS
!  enh=gMD*tempMD*(0.5d0*Q*PS**2 + S)

  sumv2=sumv2/nbeadMD
  S = S + 0.5d0*tau*PS
  pps=(sumv2-gMD*tempMD)/Q
  PS = PS + tau * pps
  S = S + 0.5d0*tau*PS
  enh=Q*0.5d0*PS**2 +gMD*S*tempMD


!  if (debug) write(98,*) pps,tempMD,sumv2
  RETURN

END
