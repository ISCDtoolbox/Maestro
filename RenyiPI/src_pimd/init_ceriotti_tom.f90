subroutine init_ceriotti_tom
  
  use md_variables
  implicit none
  integer i,j,ii,jj
    
    allocate(cost1(nbeadMD),cost2(nbeadMD))
    allocate(tmes_bead(nbeadMD))
    allocate(mass_ion(ndimMD,n))
    allocate(rpos_old(ndimMD,n,nbeadMD),forceMD_old(ndimMD,n,nbeadMD))
    
    if(sigmacov .ne. 0) then
      allocate(fk(n*ndimMD,n*ndimMD))
      allocate(cov(n*ndimMD*n*ndimMD))
      allocate(alpha_qmc(n*ndimMD*n*ndimMD))
      allocate(alphaqmc_eig(n*ndimMD))
      allocate(gamma_eigen(n*ndimMD))
      cov=0.d0 ! for the first iteration
      gamma_eigen=gammaMD ! for the first iteration
      alpha_qmc=0.d0
      alphaqmc_eig=0.d0
    endif    

    if(yesquantum) then 
      allocate(omega_mode(nbeadMD),friction_mode(nbeadMD)) 
      allocate(cmatrix(nbeadMD,nbeadMD))  
      allocate(ptilde(ndimMD,n,nbeadMD))
      allocate(rtilde_mode(ndimMD,n,nbeadMD)) 
      allocate(fbead(ndimMD,n),rcentroid(ndimMD,n),rcentroidtilde(n,ndimMD))

      do i=1,nbeadMD
        omega_mode(i) = 2.d0*tfakeMD*sin((i-1)*pi/nbeadMD)
        if(i .eq. 1) then 
          friction_mode(i) = gammaMD
        else
          friction_mode(i) = 2.d0*omega_mode(i)
        endif
        if(i .eq. 1 .and. yesglobal) then
          cost1(i) = exp(-friction_mode(i)*delt) ! PILE_G thermostat
          cost2(i) = 1.d0-cost1(i)
        else
          cost1(i) = exp(-friction_mode(i)*delt*0.5d0)
          cost2(i) = sqrt(1.d0-cost1(i)**2)
        endif
        do j=1,nbeadMD
          if(j.eq. 1) then 
            cmatrix(i,j)=1.d0
          elseif (j .ge. 2 .and. j .le. nbeadMD/2) then
            cmatrix(i,j) = sqrt(2.d0)*cos(2*pi*(i-1)*(j-1)/nbeadMD)
          elseif (j .eq. nbeadMD/2+1) then
            cmatrix(i,j) = (-1.d0)**(i-1)
          else 
            cmatrix(i,j) = sqrt(2.d0)*sin(2*pi*(i-1)*(j-1)/nbeadMD)
          endif
        enddo
      enddo
  
      cmatrix=cmatrix*sqrt(1.d0/nbeadMD)

    else ! classical  
      cost1(1) = exp(-gammaMD*delt/2.d0)
      cost2(1) = sqrt(1.d0-cost1(1)**2)
    endif ! yesquantum 

    do ii=1,ndimMD
      do jj=1,n
        mass_ion(ii,jj)=amas(jj)
      enddo
    enddo
   
    deltahtilde=0.d0 

end subroutine init_ceriotti_tom
