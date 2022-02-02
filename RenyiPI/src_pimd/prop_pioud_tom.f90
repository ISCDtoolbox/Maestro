subroutine prop_pioud_tom(ekin,epot)

! Hartree units

  use md_variables
  implicit none 

  integer :: i,l,ind,jj,kk,ii,k
  real(8) :: ekin,epot,amassq,energyq
  real(8) :: delta_noise
  logical :: yeswrite_loc

! Force fluctuations
  fk(:,:)=0.d0

! Force calculation
  call force0(forceMD)

! kinetic energy averaged over quantum images
  ekin=0.d0
  do k=1,nbead
     do i=1,n
        ind=indx(i)
        do l=1,ndimMD
           ekin=ekin+amas(ind)*vel(l,i,k)**2
        enddo
     enddo
  enddo
  ekin=0.5d0*ekin/nbead

! potential on the configuration (averaged over quantum images)

  call pot(epot,epot_centroid)
  ener_true=epot
  sigma_true=0.d0

  if(yesquantum) then

! Quantum kinetic energy

! virial
     cost=0.d0
     do ii=1,nbead
        fbead(:,:)=rcentroid(:,:)-rpos(:,:,ii)
        do jj=1,n
           do kk=1,ndimMD
              cost=cost+0.5d0*fbead(kk,jj)*forceMD(kk,jj,ii) 
           enddo
        enddo
     enddo
     ekinq=cost/nbead+ndimMD*0.5d0*temp*nion/nbead

! primitive
     ekinqp=0.d0
     do k=1,nbead
        if(k.gt.1.and.k.lt.nbead) then                 
           ekinqp=ekinqp+sum(mass_ion(:,:)*(rpos(:,:,k)-rpos(:,:,k-1))**2)
        elseif(k.eq.1) then
           ekinqp=ekinqp+sum(mass_ion(:,:)*(rpos(:,:,k)-rpos(:,:,nbead))**2)
        else
           ekinqp=ekinqp+sum(mass_ion(:,:)*(rpos(:,:,k)-rpos(:,:,k-1))**2)
        endif
     enddo
     ekinqp=-ekinqp/nbead*temp**2*0.5d0+ndimMD*0.5d0*nion*temp   

  endif

! upload forces, positions, and velocities in the proper TurboRVB vector!!!!!
! ion positions in a_0
! forces must be converted in Ry
! velocities in Ry* (scaled with the masses)

  cov_pimd=0.d0

  do k=1,nbead

     do i=1,n
        ind=indx(i)
        amassq=sqrt(amas(ind))
        do l=1,ndimMD 
           velocity(3,l+(i-1)*ndimMD,k)=vel(l,i,k)*amassq  !! qui moltiplico v0_i*sqrt(m_i) ===> è come se passo ai momenti nuovi!!
                                                           !! infatti ora il vettore contiene p0_i/sqrt(m_i)
           forcedyn(l+(i-1)*ndimMD,k)=forceMD(l,i,k)
        enddo
     enddo

     if(sigmacov.ne.0.d0) then 
! stochastic model based on dynamical matrix assumption
! fill fk (stochastic displacement diagonal in the normal harmonic coordinates)
       if(k.eq.1) yeswrite_loc=.true.
       
       call dynmatrix(rpos(1,1,k),yeswrite_loc)      
       
! add stochastic noise to forces
       do i=1,ieskin
           delta_noise=0.d0
           do kk=1,ieskin
              delta_noise=delta_noise+fk(kk,i) 
           enddo
           forcedyn(i,k)=forcedyn(i,k)+delta_noise
        enddo
! compute covariance (by adding contribution of each bead)
        call dgemm('T','N',ieskin,ieskin,ieskin,1.d0,fk,ieskin,fk,ieskin,0.d0,cov_pimd(1,k),ieskin)   

        yeswrite_loc=.false.
     endif
  enddo
  
  if(sigmacov.ne.0.d0 .and. yesquantum .and. averaged_cov) then
     cov(:)=0.d0
     do k=1,nbead
        cov(:)=cov(:)+cov_pimd(:,k)
     enddo
     cov(:)=cov(:)/dble(nbead)
     do k=1,nbead
        cov_pimd(:,k)=cov(:)
     enddo
  endif

  call reweight_pioud_tom(psip,ieskin &
                           ,scalpar,iflagerr,rank,temp,friction &
                           ,delta0,delta0q,delta0k,dt,tmes_bead,cov_pimd &
                           ,nion,normcorr,forcedyn,rpos,velocity)

! rescale velocity into atomic units 
  do i=1,n
     ind = indx(i)
     amassq=sqrt(amas(ind))
     do k=1,nbead
        do l=1,ndimMD
           vel(l,i,k)=velocity(3,l+(i-1)*ndimMD,k)/amassq
        enddo
     enddo
  enddo

  
  if(yesquantum) then
! tmes already in Ha!!!!
! true target temperature
     tmes=sum(tmes_bead(:))/nbead**2
! fake temperature (T_fake = T_physical * nbead)
!    tmes=sum(tmes_bead(:))/nbead
!    ekin=ieskin*tmes/2.d0
     energyq=ener_true+ekinq
     write(unit_dot_sigma,123) energyq,sigma_true,tmes,ekinq,ekinqp
  else
     tmes=tmes_bead(1)   
     write(unit_dot_sigma,123) ener_true,sigma_true,tmes
!     ekin=ieskin*tmes/2.d0  
  endif
  
  flush(unit_dot_sigma)

123 format(1000000e15.7)

  return

end subroutine prop_pioud_tom
