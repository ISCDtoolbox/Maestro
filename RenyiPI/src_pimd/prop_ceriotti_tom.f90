subroutine prop_ceriotti_tom(ekin,epot)
   
  use md_variables
     
  implicit none
  integer :: i,l,ind,k,kk,jj,ii
  real(8) :: ekin,epot,energyq,ekinetic

 ! kinetic energy averaged over quantum images
  ekin=0.d0
  do k=1,nbeadMD
    do i=1,n
      ind=indx(i)
      do l=1,ndimMD
        ekin=ekin+amas(ind)*vel(l,i,k)**2
      enddo
    enddo
  enddo
  ekin=0.5d0*ekin/nbeadMD
 
! Instantaneous temperature of each bead
  do k=1,nbeadMD
    tmes_bead(k)=sum(mass_ion(:,:)*vel(:,:,k)**2)/gMD  ! gMD= # degrees of freedom
  enddo

! potential on the configuration (averaged over quantum images)

  call pot(epot,epot_centroid)

  ener_true=epot
  sigma_true=0.d0

  if(yesquantum) then

! Quantum kinetic energy

  ! virial
    cost=0.d0
    do ii=1,nbeadMD
      fbead(:,:)=rcentroid(:,:)-rpos(:,:,ii)
      do jj=1,n
        do kk=1,ndimMD
          cost=cost+0.5d0*fbead(kk,jj)*forceMD(kk,jj,ii) ! Ry units!!
        enddo
      enddo
    enddo
    ekinq=cost/nbeadMD+ndimMD*0.5d0*tempMD*n

  ! primitive
    ekinqp=0.d0
    do k=1,nbeadMD
      if(k.gt.1 .and. k.lt.nbeadMD) then                 
        ekinqp=ekinqp+sum(mass_ion(:,:)*(rpos(:,:,k)-rpos(:,:,k-1))**2)
      elseif(k.eq.1) then
        ekinqp=ekinqp+sum(mass_ion(:,:)*(rpos(:,:,k)-rpos(:,:,nbeadMD))**2)
      else
        ekinqp=ekinqp+sum(mass_ion(:,:)*(rpos(:,:,k)-rpos(:,:,k-1))**2)
      endif
    enddo
    ekinqp=-ekinqp/nbeadMD*tfakeMD**2*0.5d0 + ndimMD*0.5d0*n*tfakeMD

 ! Propagation in the quantum case
    call normal_modes_fw(ptilde,pimp,.true.) ! half step for the friction in normal modes representation
    if (sigmacov .ne. 0.d0) then 
       call normal_modes_bw(ptilde,pimp,.true.)
    else
       call normal_modes_bw(ptilde,pimp,.false.) ! back to the real space
    endif
  
  else 
    call propGamma
  endif ! for yesquantum

  call propP(delt/2.d0,ekinetic) ! half step for the momenta
  if (yesquantum) then 
     call normal_modes_fw(ptilde,pimp,.false.) ! normal modes for momenta
     call normal_modes_fw(rtilde_mode,rpos,.false.) ! normal modes for positions
     call propHarm(delt) ! free ring polymer propagation in normal modes space
     call normal_modes_bw(ptilde,pimp,.false.) ! back to the real space for momenta
     call normal_modes_bw(rtilde_mode,rpos,.false.) ! back to the real space for positions
  else
     call propX(delt)
  endif
      
  call force0(forceMD)                 ! new forces calculation
  
  if(sigmacov .ne. 0.d0) call add_noise ! noisy forces
 
  call propP(delt/2.d0,ekinetic)       ! second half step for the momenta
   
  if (yesquantum) then 
     call normal_modes_fw(ptilde,pimp,.true.) ! second half step for the friction in normal modes representation
     if (sigmacov .ne. 0.d0) then
        call normal_modes_bw(ptilde,pimp,.true.)
     else 
        call normal_modes_bw(ptilde,pimp,.false.) ! back to the real spaces 
     endif
  else
     call propGamma
  endif

! Updating velocities (only momenta are updated during the propagation)
  do k=1,nbeadMD
    do i=1,n  
      ind=indx(i)
      do l=1,ndimMD
         vel(l,i,k)=pimp(l,i,k)/amas(ind)
      enddo
    enddo
  enddo

! Evaluation of deltaHtilde --> conserved quantity 
! (eq. 19 of the paper 'Accurate sampling using Langevin dynamics' - Bussi)
  do k=1,nbeadMD
    do i=1,n
      ind=indx(i)
      do l=1,ndimMD
        deltahtilde=deltahtilde+(rpos(l,i,k)-rpos_old(l,i,k))*(forceMD(l,i,k) + forceMD_old(l,i,k))*0.5d0 &
                    + delt**2/8.d0/amas(ind)*(forceMD(l,i,k)**2 &
                    - forceMD_old(l,i,k)**2)  
      enddo
    enddo 
  enddo 
  deltahtilde=deltahtilde+epot-epot_old

  if(yesquantum) then
  !  tmes already in Ha!!!!
!    true target temperature
    tmes=sum(tmes_bead(:))/nbeadMD**2
!   fake temperature (T_fake = T_physical * nbead)
!   tmes=sum(tmes_bead(:))/nbead
!   ekin=ieskin*tmes/2.d0
    energyq=ener_true+ekinq
    write(unit_dot_sigma,123) energyq,sigma_true,tmes,ekinq,ekinqp 
  else
    tmes=tmes_bead(1)   
    write(unit_dot_sigma,123) ener_true,sigma_true,tmes!,((forceMD_old(l,i,1),l=1,ndimMD),i=1,n)
!     ekin=ieskin*tmes/2.d0  
  endif
  flush(unit_dot_sigma)

123 format(1000000e15.7)
   
  return

end subroutine prop_ceriotti_tom


subroutine normal_modes_fw(coord_mode,coord,first)
  use md_variables
  implicit none
  integer :: i,k,l,kk,ind
  real(8) :: drand1,drand2,xi,alpha2,alpha,kin,sumxi2,xi1,testsign
  real(8), dimension(ndimMD,n,nbeadMD) :: coord_mode,coord
  real(8), dimension(:,:), allocatable :: sov5,sov6
  real(8), dimension(:), allocatable :: pimpion,sov4
  logical :: first
   
  coord_mode=0.d0
       
  if (yesglobal) then
     kin=0.d0
     sumxi2=0.d0
     xi1=0.d0
  endif
      
  if (first .and. sigmacov .ne. 0.d0) then ! Noise on CC forces 

     allocate(sov4(n*ndimMD),pimpion(n*ndimMD),sov6(n*ndimMD,nbeadMD))
     allocate(sov5(n*ndimMD,nbeadMD))
     sov4=0.d0
     sov6=0.d0
! sov4 non ha gli indici bead (così come pimpion), mentre sov5 e sov6 si!!
     do k=1,nbeadMD
       do i=1,n
       ind=indx(i)
         do l=1,ndimMD
! Mass scaling to have mass-independant propagation in the space which diagonalizes gamma matrix
            pimpion(l+(i-1)*ndimMD)=coord(l,i,k)/sqrt(amas(ind))   
         enddo
       enddo
       
       do i=1,n*ndimMD
! Transformation in the basis which diagonalizes the friction
!! perchè fa quest aoperazione dentro un ciclo ????
         call dgemv('T',n*ndimMD,n*ndimMD,1.d0,cov,n*ndimMD,pimpion,1,0.d0,sov4,1) ! Now, momenta are in the sov4 vector
       enddo
! Store the momenta in a bead dependant vector
       sov5(1:n*ndimMD,k)=sov4(1:n*ndimMD)
     enddo 
     
! Second basis transformation to work in normal modes
     do k=1,nbeadMD
       do i=1,n*ndimMD
         do kk=1,nbeadMD
           sov6(i,k)=sov6(i,k)+sov5(i,kk)*cmatrix(kk,k) ! sov6 are the momenta after the double transformation
         enddo
       enddo 
     enddo
     
     do k=1,nbeadMD
       sov4(1:n*ndimMD)=sov6(1:n*ndimMD,k)
       do i=1,n*ndimMD
   !!! gamma_eigen(i) non viene mai ricalcolato e rimane sempre uguale a gamma.....    
         cost1(k)=exp(-(friction_mode(k)+gamma_eigen(i)-gammaMD)*delt/2.d0) 
         cost2(k)=sqrt(1.d0-cost1(k)**2)   
         call random_number(drand1)
         call random_number(drand2)
         xi=dsqrt(-2.d0*dlog(1.d0-drand1))*dcos(2.d0*pi*drand2)
         sov4(i) = cost1(k)*sov4(i) + sqrt(tfakeMD)*cost2(k)*xi
       enddo
       sov5(1:n*ndimMD,k)=sov4(1:n*ndimMD)
       do i=1,n
         do l=1,ndimMD
           coord_mode(l,i,k)=sov5(l+(i-1)*ndimMD,k) ! to store momenta before calling the backward transformation routine
         enddo
       enddo 
     enddo

     deallocate(sov4,sov5,sov6,pimpion)
 
   else ! No noise on CC forces

     do k=1,nbeadMD
       do i=1,n
         ind=indx(i)
         do l=1,ndimMD
           do kk=1,nbeadMD
             coord_mode(l,i,k) = coord_mode(l,i,k) + coord(l,i,kk)*cmatrix(kk,k)
           enddo
           if (first) then
             if (.not. yesglobal .or. k .gt. 1) then ! PILE_L + PILE_G for modes > 0
                call random_number(drand1)
                call random_number(drand2)
                xi=dsqrt(-2.d0*dlog(1.d0-drand1))*dcos(2.d0*pi*drand2)
                coord_mode(l,i,k) = cost1(k)*coord_mode(l,i,k) + sqrt(amas(ind)*tfakeMD)*cost2(k)*xi
             endif
           endif
         enddo
       enddo
     enddo

     if (first) then
       if (yesglobal) then ! PILE_G for mode = 0 => To be checked because I suspect it does not work properly
          do i=1,n
            ind=indx(i)
            do l=1,ndimMD
              kin=kin+coord_mode(l,i,1)**2/amas(ind)
              call random_number(drand1)
              call random_number(drand2)
              xi=dsqrt(-2.d0*dlog(1.d0-drand1))*dcos(2.d0*pi*drand2)
              sumxi2=sumxi2+xi**2
              if (i .eq. 1) xi1=xi1+xi
            enddo 
          enddo
          kin=kin/2.d0 
          alpha2=cost1(1)+cost2(1)*sumxi2*tfakeMD/(2.d0*kin)
          alpha2=alpha2+2.d0*xi1*dsqrt(cost1(1)*cost2(1)*tfakeMD/(2.d0*kin))
          testsign=xi1+dsqrt(2.d0*kin*cost1(1)/(cost2(1)*tfakeMD))
          if (testsign .ge. 0.d0) then 
             alpha=dsqrt(alpha2)
          else
             alpha=-dsqrt(alpha2)
          endif
          coord_mode(l,i,1)=alpha*coord_mode(l,i,1)
        endif
     endif

  endif 

  return

end subroutine normal_modes_fw


subroutine normal_modes_bw(coord_mode,coord,first)
    
  use md_variables
  implicit none

  integer :: i,k,l,kk,ind
  real(8),dimension(ndimMD,n,nbeadMD) :: coord,coord_mode
  real(8),dimension(:,:), allocatable :: sov5,sov6
  real(8),dimension(:), allocatable :: sov4,pimpion
  logical :: first

  coord=0.d0

  if (first .and. sigmacov .ne. 0.d0) then  ! Noisy CC forces
     allocate(sov5(n*ndimMD,nbeadMD),sov6(n*ndimMD,nbeadMD))
     allocate(sov4(n*ndimMD),pimpion(n*ndimMD))
     sov6=0.d0
     pimpion=0.d0
     do k=1,nbeadMD
       do i=1,n
         do l=1,ndimMD
           sov5(l+(i-1)*ndimMD,k)=coord_mode(l,i,k)
         enddo
       enddo
     enddo
! Backward transformation to leave from normal modes representation
     do k=1,nbeadMD
       do i=1,n*ndimMD
         do kk=1,nbeadMD
           sov6(i,k)=sov6(i,k)+cmatrix(k,kk)*sov5(i,kk) ! sov6 are now the momenta in the space which diagonalizes gamma (1 transformation remaining)
         enddo
       enddo
       sov4(1:n*ndimMD)=sov6(1:n*ndimMD,k)
! backward tranformation with respect to gamma (coordinate space)
       call dgemv('N',n*ndimMD,n*ndimMD,1.d0,cov,n*ndimMD,sov4,1,0.d0,pimpion,1)
! Mass-scaling to recover the physical momenta
       do i=1,n
         ind=indx(i)
         do l=1,ndimMD
           coord(l,i,k)=pimpion(l+(i-1)*ndimMD)*sqrt(amas(ind))
         enddo
       enddo
     enddo
     
     deallocate(sov4,sov5,sov6,pimpion)
         
  else 
 
     do k=1,nbeadMD
       do i=1,n
         do l=1,ndimMD
           do kk=1,nbeadMD
             coord(l,i,k) = coord(l,i,k) + cmatrix(k,kk)*coord_mode(l,i,kk)
           enddo
         enddo
       enddo
     enddo

  endif

  return

end subroutine normal_modes_bw   


subroutine propHarm(timestep)
      
  use md_variables
  implicit none

  integer :: i,k,l,ind
  real(8) :: cosit,sinut,timestep
  real(8),dimension(ndimMD,n,nbeadMD) :: rtilde0,ptilde0

  ptilde0=ptilde
  rtilde0=rtilde_mode
      
  do k=1,nbeadMD
    cosit=cos(omega_mode(k)*timestep)
    sinut=sin(omega_mode(k)*timestep)
    do i=1,n
      ind=indx(i)
      do l=1,ndimMD
        ptilde(l,i,k) = cosit*ptilde0(l,i,k) - amas(ind)*omega_mode(k)*sinut*rtilde0(l,i,k)
        if (k .eq. 1) then 
           rtilde_mode(l,i,k) = cosit*rtilde0(l,i,k) + 1.d0/amas(ind)*timestep*ptilde0(l,i,k) 
        else
           rtilde_mode(l,i,k) = cosit*rtilde0(l,i,k) + 1.d0/(amas(ind)*omega_mode(k))*sinut*ptilde0(l,i,k)
        endif      
      enddo 
    enddo
  enddo
    
  return

end subroutine propHarm


subroutine propGamma
     
  use md_variables
  implicit none
   
  integer :: i,l,ind
  real(8) :: drand1,drand2,xi,correct_noise,cost0
  real(8), dimension(:), allocatable :: pimpion,sov4
      
  if (sigmacov .ne. 0.d0) then  ! Noisy forces
     allocate(sov4(n*ndimMD),pimpion(n*ndimMD))
     sov4=0.d0
     do i=1,n
       ind=indx(i)
       do l=1,ndimMD
! Mass scaling to have mass-independant propagation in the space which diagonalizes gamma matrix       
         pimpion(l+(i-1)*ndimMD)=pimp(l,i,1)/sqrt(amas(ind))  
       enddo 
     enddo

! Transformation in the basis which diagonalizes the friction
     call dgemv('T',n*ndimMD,n*ndimMD,1.d0,cov,n*ndimMD,pimpion,1,0.d0,sov4,1) ! Now, momenta are in the sov4 vector

     do i=1,n*ndimMD
        !correct_noise=gamma/gamma_eigen(i)
       cost0=exp(-alphaqmc_eig(i)*delt/(8.d0*tempMD))
       cost1(1)=exp(-(gamma_eigen(i)-alphaqmc_eig(i)/(2.d0*tempMD))*delt/2.d0)
       cost2(1)=sqrt(1.d0-cost1(1)**2)
       call random_number(drand1)
       call random_number(drand2)
       xi=dsqrt(-2.d0*dlog(1.d0-drand1))*dcos(2.d0*pi*drand2)
       !xi=xi*sqrt(correct_noise)
       sov4(i) = cost0*sov4(i)
       sov4(i) = cost1(1)*sov4(i) + sqrt(tempMD)*cost2(1)*xi
       sov4(i) = cost0*sov4(i)
     enddo

! Back to the original basis 
     call dgemv('N',n*ndimMD,n*ndimMD,1.d0,cov,n*ndimMD,sov4,1,0.d0,pimpion,1)

     do i=1,n
       ind=indx(i)
       do l=1,ndimMD
         pimp(l,i,1)=pimpion(l+(i-1)*ndimMD)*sqrt(amas(ind))
       enddo
     enddo

     deallocate(sov4,pimpion)

  else ! No noise on forces
      
     do i=1,n
       ind=indx(i)
       do l=1,ndimMD
         call random_number(drand1)
         call random_number(drand2)
         xi=dsqrt(-2.d0*dlog(1.d0-drand1))*dcos(2.d0*pi*drand2)
         pimp(l,i,1) = cost1(1)*pimp(l,i,1) + sqrt(amas(ind)*tempMD)*cost2(1)*xi
       enddo
     enddo
 
  endif 

  return
 
end subroutine propGamma


subroutine add_noise
    
  use md_variables
  implicit none
      
  integer :: i,k,l,kk,ii,ll
  real(8) :: delta_noise
  real(8), dimension (:), allocatable :: alphaqmc_inter,alphaqmc_diag
  logical :: yeswrite_loc

! Force fluctuations
  fk=0.d0
! Covariance matrix
  cov=0.d0
    
  allocate(alphaqmc_inter(n*ndimMD*n*ndimMD),alphaqmc_diag(n*ndimMD*n*ndimMD))

  alphaqmc_inter=0.d0
  alphaqmc_diag=0.d0
 
! stochastic model based on dynamical matrix assumption
! fill fk (stochastic displacement diagonal in the normal harmonic coordinates)
  do k=1,nbeadMD
    do i=1,n
      do l=1,ndimMD
        forcedyn(l+(i-1)*ndimMD,k)=forceMD(l,i,k)
      enddo
    enddo
    if(k.eq.1) yeswrite_loc=.true.

    call dynmatrix(rpos(:,:,k),yeswrite_loc) !!! perchè gli passo solo un elemento di matrice ??????, ok siccome è un vettore global
                                             !!! è la stessa cosa passargli rpos(1,1,k) o rpos(:,:,k)
! add stochastic noise to forces
    do i=1,n*ndimMD
      delta_noise=0.d0
      do kk=1,n*ndimMD
        delta_noise=delta_noise+fk(kk,i)
      enddo
      forcedyn(i,k)=forcedyn(i,k)+delta_noise
    enddo
    do i=1,n
      do l=1,ndimMD
        forceMD(l,i,k)=forcedyn(l+(i-1)*ndimMD,k)
      enddo
    enddo
! compute covariance (by adding contribution of each bead)
    
    call dgemm('T','N',n*ndimMD,n*ndimMD,n*ndimMD,1.d0,fk,n*ndimMD,fk,n*ndimMD,1.d0,cov,n*ndimMD)
    yeswrite_loc=.false.
         !do i=1,n*ndimMD 
            !do l=1,n*ndimMD
               !cov(l+(i-1)*n*ndimMD)=cov(l+(i-1)*n*ndimMD)+sigmavar*dynmat_force0(l,i)
            !enddo
         !enddo
  enddo
! average covariance over beads
  cov=cov/dble(nbeadMD)
! Scale covariance by the masses 
  do i=1,n*ndimMD
    ii=int((i-1)/ndimMD+1)
    do l=1,n*ndimMD
      ll=int((l-1)/ndimMD+1)
      cov(l+(i-1)*ndimMD)=cov(l+(i-1)*ndimMD)/sqrt(amas(indx(ii))*amas(indx(ll)))
    enddo
  enddo
! Save alpha_qmc before cov is modified  
  alpha_qmc=cov
  
! construct gamma matrix ---> cov/2T 
  if(tempMD.ne.0.d0) then
     cost=delta0*0.5d0/tempMD
  else
     cost=0.d0
  endif
  call dscal(n*ndimMD*n*ndimMD,cost,cov,1)
  
! add diagonal contribution to gamma matrix
  do i=1,n*ndimMD
 ! No bead index since we assume the covariance matrix is the same for all replicas
     cov(n*ndimMD*(i-1)+i)=cov(n*ndimMD*(i-1)+i)+gammaMD
  enddo

! diagonalize gamma matrix
  call dsyev_my('V','L',n*ndimMD,cov,n*ndimMD,gamma_eigen,info)
! Minimum eigenvalue is gamma_0
  do i=1,n*ndimMD
    if (gamma_eigen(i) .lt. gammaMD) gamma_eigen(i)=gammaMD
  enddo
! Diagonalize alpha_qmc matrix => cov diagonalizes gamma so it must also diagonalize alpha_qmc
  call dsyev_my('V','L',n*ndimMD,alpha_qmc,n*ndimMD,alphaqmc_eig,info)
  !call dgemm('T','N',n*ndimMD,n*ndimMD,n*ndimMD,1.d0,cov,n*ndimMD,alpha_qmc,n*ndimMD,0.d0,alphaqmc_inter,n*ndimMD)
  !call dgemm('N','N',n*ndimMD,n*ndimMD,n*ndimMD,1.d0,alphaqmc_inter,n*ndimMD,cov,n*ndimMD,0.d0,alphaqmc_diag,n*ndimMD)

  !call dgemv('T',n*ndimMD,n*ndimMD,1.d0,cov,n*ndimMD               &
  !     ,alpha_qmc,1,0.d0,alphaqmc_inter,1)
  !call dgemv('N',n*ndimMD,n*ndimMD,1.d0,alphaqmc_inter,n*ndimMD               &
  !      ,cov,1,0.d0,alphaqmc_diag,1)
! Mininum eigenvalue of alpha_qmc is 0 since it is a positive matrix
  do i=1,n*ndimMD
     if (alphaqmc_eig(i) .lt. 0.d0) alphaqmc_eig(i)=0.d0
  enddo
! Covariance and alpha_qmc eigenvalues output
  if (verbose) then 
     write(6,*) ' Eigenvalues covariance'
     do i=1,n*ndimMD
       write(6,*) i,gamma_eigen(i)
     enddo
     write(6,*) ' Alpha_qmc diag '
     do i=1,n*ndimMD
       write(6,*) i,alphaqmc_eig(i)
     enddo
  endif

  deallocate(alphaqmc_inter,alphaqmc_diag)

  return
      
end subroutine add_noise
