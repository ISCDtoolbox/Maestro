subroutine reweight_pioud_tom(psip,ieskin,scalpar,iflagerr,rank &
                              ,temp,friction &
                              ,delta0,delta0q,delta0k,dtr,tmes,cov &
                              ,nion,normcorr,forza,rion,velocity)

 use md_variables, only: kdyn_eig,yessecond &
                         ,kdyn,nbead,ndimMD

! three vectors as fundamental input (bead dependent)
! rion, velion, forza
! forza(ieskin,*),velocity(3,ieskin,*),rion(3,nion,*)
! velocity(1,ieskin,*) and velocity(3,ieskin,*) are "input/output sensitive"
! velocity(1,ieskin,*) "intermediate" positions stored (to be used in the next iteration, updated in output)
! velocity(3,ieskin,*) actual velocities (updated in output)
! velocity(2,ieskin,*) scratch vector
! forza(ieskin,*) ion forces (only input parameter)
! rion(3,nion,*)  ion positions (updated in output)

! cov matrix bead specific!!!

  implicit none 
  integer nion,i,info &
          ,iflagerr,ieskin,ind  &
          ,n2,n3,n4,n5,maxnm,jj,ii   &
          ,indin,ieskin2,kk,kkf,j,ibead,ieskin2_pimd
  real*8 normcorr,psip(*),cost,zeta &
         ,scalpar(*),forza(ieskin,*),velocity(3,ieskin,*) &
         ,temp,pi,dt,friction,tmes(*),dtr &
         ,cov(ieskin*ieskin,*),delta0,delta0q,delta0k &
         ,rion(ndimMD,nion,*),cov_sav(ieskin*ieskin)
  real*8 mnoise(2,2),zetan(2),costh,eta_v,eta_r &
         ,Tn,Gn,Gni(2,2),alphaqmc,alphaall &
         ,dth,Gnh,Gnih
  real*8 drand1,drand2,dnrm2
  integer  errnoise
  real(8), dimension(:,:), allocatable:: sov4
  real(8), dimension(:,:,:), allocatable:: sov5
  real(8), dimension(:,:), allocatable:: velion
  real(8), dimension(:,:,:), allocatable:: vel0
  parameter(pi=3.14159265358979323846D0) 
  integer rank
  logical info_check
                    
  info_check=.false.
                           
  if(iflagerr.ne.0) return 
  
  errnoise=0 
  maxnm=ieskin
  n2=maxnm+1
  n3=maxnm+n2 
  n4=maxnm+n3 
  n5=maxnm+n4 
  ieskin2=ieskin*ieskin 
  ieskin2_pimd=ieskin2*nbead

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! bead independent: GLOBAL variables 

! NOTE : psip(1:ieskin) and psip(n3:n3+ieskin-1) are filled before calling ion_dynamics
! they must not be overwritten !!!!
! psip(1:ieskin) ---> eigenvalues of gamma matrix
! psip(n3:n3+ieskin-1) ---> \sqrt{m_0/m_i} with m_0 mass of the first atom in the geometry

  allocate(vel0(3,ieskin,nbead))
  allocate(velion(3,ieskin))
  allocate(sov4(ieskin,2))
  sov4=0.d0

! time step
  dt=dtr
  if(yessecond) then 
     dth=dt/2.d0
  else
     dth=dt
  endif

! mass scaling factor
  do j=1,ieskin 
     psip(n3+j-1)=dsqrt(scalpar(j))
  enddo

!! gamma deve avere dimensioni di (1/tempo); temp ha dimensioni di energia ==> cov deve avere dimensioni di [energia/tempo]
! construct gamma matrix ---> cov/2T
! covariance computed outside (cov given in input)!!!

  if(temp.ne.0.d0) then 
     cost=delta0*0.5d0/temp !!! cost ha le dimensioni di un'azione
  else 
     cost=0.d0 
  endif
  call dscal(ieskin2_pimd,cost,cov,1) !!! cov ora ha le dimensioni di (gamma=alfa/2T)


  do ibead=1,nbead
 
! forces stored in sov4(1:ieskin,2) 
     call dcopy(ieskin,forza(1,ibead),1,sov4(1,2),1)

! scale the force by the Mass (in realtà per sqrt(mass_ion(1,1)/mass_ion(i,j)) )                               
     do i=1,ieskin 
        sov4(i,2)=sov4(i,2)*psip(n3+i-1)   
     enddo

! velocities stored in velion(1:ieskin,3)
     vel0(:,:,ibead)=velocity(:,:,ibead)
     velion(:,:)=velocity(:,:,ibead)

! each bead has its own covariance matrix

! scale covariance by masses   (in realtà per mass_ion(1,1)/mass_ion(i,j) )                                                                   
     do j=1,ieskin 
        do i=1,ieskin 
           cov(ieskin*(j-1)+i,ibead)=cov(ieskin*(j-1)+i,ibead)*psip(n3+i-1)*psip(n3+j-1) 
        enddo
     enddo

! add diagonal contribution to gamma matrix
     do i=1,ieskin 
        cov(ieskin*(i-1)+i,ibead)=cov(ieskin*(i-1)+i,ibead)+delta0q 
     enddo
    
  
!  do j=1,ieskin
!     write(6,*) j,psip(n3+j-1)
!     do i=1,ieskin
!        write(6,*) i,j,cov(ieskin*(j-1)+i)
!     enddo
!  enddo

! diagonalize gamma matrix
     call dsyev_my('V','L',ieskin,cov(1,ibead),ieskin,psip,info)  !! psip(1:ieskin) contiene gli autovalori di gamma ma 
     ! c'è il fattore di prima psip(n3:4*ieskin)
 
    ! if(ibead.eq.1) then 
    !    write(6,*) ' Eigenvalues covariance'
    !    do i=1,ieskin
    !       write(6,*) i,psip(i)
    !    enddo
!       write(6,*) ' Eigenvectors covariance'
!       do i=1,ieskin
!          write(6,*) i,(cov(ieskin*(i-1)+j),j=1,ieskin)
!       enddo
    ! endif
  
     if(info.ne.0) then 
        write(6,*) ' Error in lapack dsyev dynamic !!! ' 
        errnoise=4 
        info_check=.true.
     endif

!  if(rank.eq.0) write(6,*) ' Ratio dyn =',(psip(ieskin)-friction)/sqrt(temp*157.8873306d0)
     if(ibead.eq.1 .and. rank.eq.0) write(6,*) ' Ratio dyn =',(psip(ieskin)-delta0q)*dt
!  write(6,*) psip(ieskin),friction,temp

     do i=1,ieskin 
        if(psip(i).gt.0.d0) then 
           psip(n5+i-1)=(1.d0-dexp(-psip(i)*dth))/psip(i) 
        else 
           if(rank.eq.0) write(6,*) ' Refused eigenvalue ',i,ibead 
           psip(n5+i-1)=dth 
        endif
     enddo
!         compute the noise correction                                  
     alphaqmc=0.d0
                                                                        
     do i=1,ieskin 

        if(psip(i).gt.0.d0) then 
!         the following is protected for overflow                       
           cost=exp(-dth*psip(i)) 
           alphaqmc=normcorr*2.d0*temp*(psip(i)-delta0q)/delta0
           psip(n4+i-1)=temp*psip(i)**2*(1.d0+cost)/(1.d0-cost)-alphaqmc                          
!          subtracting the noise already present in the forces          
        else 
           psip(n4+i-1)=0.d0 
        endif

        if(psip(n4+i-1).gt.0.d0) then 
           psip(n4+i-1)=dsqrt(psip(n4+i-1)) 
        else 
           psip(n4+i-1)=0.d0 
           if(delta0q.gt.0.d0) then 
              if(rank.eq.0) write(6,*) 'There should be some error in reweight0',i,ibead,psip(n4+i-1)
              errnoise=5 
           endif
        endif

     enddo
                                                                        
! now we are able to update the velocity                        
! first change basis actual velocity and force                  
     call dgemv('T',ieskin,ieskin,1.d0,cov(1,ibead),ieskin,velion(3,1),3,0.d0,sov4,1)                               
     call dgemv('T',ieskin,ieskin,1.d0,cov(1,ibead),ieskin,sov4(1,2),1,0.d0,velion(2,1),3)                        

     if(yessecond) then
! Compute the temperature at half time intervals.
        do i=1,ieskin 
           call random_number(drand1)
           call random_number(drand2)
           zeta=dsqrt(-2.d0*dlog(1.d0-drand1))*dcos(2.d0*pi*drand2) 
           sov4(i,1)=sov4(i,1)*dexp(-psip(i)*dth)                    &
                     +psip(n5+i-1)*(velion(2,i)+psip(n4+i-1)*zeta)                     
        enddo
        tmes(ibead)=dnrm2(ieskin,sov4,1)**2/ieskin 
! Second half time interval.
        do i=1,ieskin 
           call random_number(drand1)
           call random_number(drand2)
           zeta=dsqrt(-2.d0*dlog(1.d0-drand1))*dcos(2.d0*pi*drand2) 
           velion(2,i)=sov4(i,1)*dexp(-psip(i)*dth)                    &
                       +psip(n5+i-1)*(velion(2,i)+psip(n4+i-1)*zeta)                     
        enddo
     else
        do i=1,ieskin
           call random_number(drand1)
           call random_number(drand2)
           zeta=dsqrt(-2.d0*dlog(1.d0-drand1))*dcos(2.d0*pi*drand2) 
           velion(2,i)=sov4(i,1)*dexp(-psip(i)*dth)                    &
                       +psip(n5+i-1)*(velion(2,i)+psip(n4+i-1)*zeta)                     
        enddo
     endif
                                                                        
! go back to the original basis                                 
                                                                        
     call dgemv('N',ieskin,ieskin,1.d0,cov(1,ibead),ieskin,velion(2,1),3,0.d0,velion(3,1),3)                                
                                                                        
     if(.not.yessecond) tmes(ibead)=dnrm2(ieskin,velion(3,1),3)**2/ieskin 

! prepare for the next (harmonic) step

     velocity(3,1:ieskin,ibead)=velion(3,1:ieskin) ! velocities
     velocity(2,1:ieskin,ibead)=0.d0 ! forces 
     do j=1,nion
        do i=1,ndimMD
           ind=i+(j-1)*ndimMD
           velocity(1,ind,ibead)=rion(i,j,ibead)/psip(n3+ind-1) ! positions
        enddo
     enddo

!     write(6,*) 'velocity before'
!     do j=1,nion
!        do i=1,ndimMD
!           ind=i+(j-1)*ndimMD
!           write(6,*) i,j,ibead,velocity(1,ind,ibead)
!           write(6,*) i,j,ibead,velocity(2,ind,ibead)
!        enddo
!     enddo     

  enddo

! harmonic forces step
  allocate(sov5(3,ieskin,nbead))
  sov5=0.d0

! rotation in the kdyn basis
! write the above operations with rank-3 blas!!!!
  call dgemm('N','N',3*ieskin,nbead,nbead,1.d0,velocity,3*ieskin,kdyn,nbead,0.d0,sov5,3*ieskin)

! sov5(1,1:ieskin,ibead)  ! positions
! sov5(2,1:ieskin,ibead)  ! forces
! sov5(3,1:ieskin,ibead)  ! velocities

  do ibead=1,nbead

!!!! START equations of motion integration 

     do i=1,ieskin 

        call random_number(drand1)
        call random_number(drand2)
        zetan(1)=dsqrt(-2.d0*dlog(1.d0-drand1))*dcos(2.d0*pi*drand2) 

        call random_number(drand1)
        call random_number(drand2)
        zetan(2)=dsqrt(-2.d0*dlog(1.d0-drand1))*dcos(2.d0*pi*drand2) 
 
! TEST: freezing the random number generator         
!          zetan(1)=0.55
!          zetan(2)=0.75

! the QMC noisy part has been integrated in the first part
        alphaqmc=0.d0

!!!! START equations of motion integration for harmonic part
! gammaall
        psip(i)=2.d0*delta0k*sqrt(abs(kdyn_eig(ibead)))
        if(psip(i).lt.friction) psip(i)=friction ! cutoff on the lowest eigenvalue. Ceriotti's choice.

        alphaall=2.d0*temp*psip(i)

        call set_turboq(psip(i),kdyn_eig(ibead),dt &
             ,alphaall,alphaqmc,mnoise,Gn,Tn,Gni,Gnh,Gnih)

        call root2mat(mnoise,errnoise)

        if(errnoise.ne.0) write(6,*) ' Error negative definite matrix noise '
                                                                        
        eta_v=mnoise(1,1)*zetan(1)+mnoise(1,2)*zetan(2)                                   
        eta_r=mnoise(2,1)*zetan(1)+mnoise(2,2)*zetan(2)
          
        sov4(i,2)=sov5(3,i,ibead)*Gni(2,1)+Gni(2,2)*sov5(1,i,ibead)+Tn*eta_r
        velion(2,i)=sov5(3,i,ibead)*Gni(1,1)+Gni(1,2)*sov5(1,i,ibead)+Gn*eta_v
       
     enddo

     velocity(2,1:ieskin,ibead)=velion(2,1:ieskin) ! velocities
     velocity(1,1:ieskin,ibead)=sov4(1:ieskin,2) ! coordinates

!     write(6,*) 'velocity'
!     do j=1,nion
!        do i=1,ndimMD
!           ind=i+(j-1)*ndimMD
!           write(6,*) i,j,ibead,velocity(1,ind,ibead)
!           write(6,*) i,j,ibead,velocity(2,ind,ibead)
!        enddo
!     enddo     

  enddo
                                                                        
! Now go back to the original basis                             

! write the above operations with rank-3 blas!!!!
  call dgemm('N','T',3*ieskin,nbead,nbead,1.d0,velocity,3*ieskin,kdyn,nbead,0.d0,sov5,3*ieskin)

! sov5(1,1:ieskin,ibead) coordinates
! sov5(2,1:ieskin,ibead) velocities
     
  do ibead=1,nbead

     call dcopy(ieskin,sov5(2,1,ibead),3,velion(3,1),3)
     call dcopy(ieskin,sov5(1,1,ibead),3,psip(n2),1)

!     write(6,*) 'sov5'
!     do j=1,nion
!        do i=1,ndimMD
!           ind=i+(j-1)*ndimMD
!           write(6,*) i,j,ibead,sov5(1,ind,ibead)
!           write(6,*) i,j,ibead,sov5(2,ind,ibead)
!        enddo
!     enddo     

!       back to the scale of coordinates                                
     do i=1,ieskin 
        psip(n2+i-1)=psip(n2+i-1)*psip(n3+i-1)
        if(psip(n3+i-1).eq.0.d0) then 
           velion(3,i)=0.d0
        endif
     enddo

     if(ibead.eq.1 .and. rank.eq.0) then
        write(6,*) ' Temperature (Ha) = ',tmes(ibead) 
        cost=dnrm2(ieskin,psip(n2),1)
        write(6,*) ' Norm change ions =',cost
        write(6,*) ' Coordinates variation for bead',ibead                              
        do i=1,ieskin 
           write(6,*) i,psip(n2+i-1)   
        enddo
     endif

! update velocities
     velocity(:,:,ibead)=velion(:,:)

     if(info_check) then 
        write(6,*) ' Warning info =',info 
        write(6,*) 'Continuing with previous param. !!!' 
! do not move ions
        call dscalzero(ieskin,0.d0,psip(n2),1) 
! restore previous velocities
        velocity(:,:,ibead)=vel0(:,:,ibead)
     endif

     kkf=0
     do kk=1,nion
        do jj=1,ndimMD
           kkf=kkf+1
           rion(jj,kk,ibead)=rion(jj,kk,ibead)+psip(n2+kkf-1)
        enddo
     enddo
     
  enddo

! over beads

  if(errnoise.ne.0) then                                                           
     write(6,*) ' Error in reweight0  =',errnoise 
     stop                                                                         
  endif

  deallocate(vel0,velion,sov4,sov5)
                                                                              
  return

123 format(1000000e15.7)
     
end subroutine reweight_pioud_tom
