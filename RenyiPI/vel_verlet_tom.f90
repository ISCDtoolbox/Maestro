subroutine vel_verlet_tom(enk,ep,enh)
  
  use md_variables
  implicit none
  real(8) :: enk,ep,ep_centroid,enh
  
  
  if(nh) call propNH(0.5d0*delt, enh)
  
  call propP(0.5d0*delt,enk)

  call propX(delt)
  call force0(forceMD)

  call propP(0.5d0*delt,enk)
  
  call pot(ep,ep_centroid)
 
  if(nh) call propNH(0.5d0*delt, enh)
  
  
end subroutine vel_verlet_tom
