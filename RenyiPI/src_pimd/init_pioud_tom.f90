subroutine init_pioud_tom

! Hartree units

  use md_variables

  implicit none

  ! Variable 'rorder' added by Miha Srdinsek
  integer :: ii,jj,ind,i,rorder

  ieskin=ndimMD*n
  iflagerr=0
  dt=delt
  nbead=nbeadMD
  temp=tempMD
  nion=n
  rorder=renyi_order
  iscramax=6*ieskin ! Perhaps this should be higher number for Renyi - Miha S.
  iscramax=max(3*nbead * rorder,iscramax)
  normcorr=1.d0
  scalecov=1.d0

  if(nbeadMD.eq.1) then
     delta0q=gammaMD
  else
     friction=gammaMD
  endif

! set output control in Turbo Langevin dynamics
  rank=1
  if(verbose) rank=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(6,*) '*************************************'
  write(6,*) 'Turbo Ceriotti integrator parameters'
  write(6,*) 'time step =',dt
  write(6,*) 'friction (gamma) =',gammaMD
  write(6,*) 'delta0 =',delta0
  if(yessecond) write(6,*) 'second order Trotter breakup'
  if(nbeadMD.gt.1) then
     write(6,*) 'Quantum MD: Path integral Langevin!!!!'
     write(6,*) 'Generalized dynamics in the bead eigenmodes'
     write(6,*) 'delta0k, delta0q =',delta0k,delta0q
     if(averaged_cov.and.sigmacov.ne.0.d0) &
          write(6,*) 'Warning: averaging the covariance over the beads!!!'
  endif

  if (.not. restart_pimd) then
    write(unit_dot_out,*) '*************************************'
    write(unit_dot_out,*) 'Turbo Ceriotti integrator parameters'
    write(unit_dot_out,*) 'time step =',dt
    write(unit_dot_out,*) 'friction (gamma) =',gammaMD
    write(unit_dot_out,*) 'delta0 =',delta0
    if(yessecond) write(unit_dot_out,*) 'second order Trotter beakup'
    if(nbeadMD.gt.1) then
      write(unit_dot_out,*) 'Quantum MD: Path integral Langevin!!!!'
      write(unit_dot_out,*) 'Generalized dynamics in the bead eigenmodes'
      write(unit_dot_out,*) 'delta0k, delta0q =',delta0k,delta0q
      if(averaged_cov.and.sigmacov.ne.0.d0) &
          write(unit_dot_out,*) 'Warning: averaging the covariance over the beads!!!'
    endif
  end if
  allocate(psip(iscramax))
  allocate(scalpar(ieskin))
  ! Multiplication by '* order' added by Miha Srdinsek
  allocate(cov_pimd(ieskin*ieskin,nbead * rorder))
  allocate(fk(n*ndimMD,n*ndimMD))

  if(delta0 .lt. dt*normcorr) then
     delta0=dt*normcorr
     write(6,*) ' Warning  delta_0 >= dt !!! Changed to ',delta0
     if (.not. restart_pimd) write(unit_dot_out,*) ' Warning  delta_0 >= dt !!! Changed to ',delta0
  endif
  write(6,*) '*************************************'
  if (.not. restart_pimd) write(unit_dot_out,*) '*************************************'

  if(yesquantum) then
     temp=temp*nbead ! Miha: Here I don't modify
     allocate(tmes_bead(nbead)) ! Miha: This is irrelavant
     tmes_bead=0.d0
     allocate(fbead(ndimMD,nion),rcentroid(ndimMD,nion),rcentroidtilde(nion,ndimMD))
     allocate(mass_ion(ndimMD,nion))
     write(6,*) ' Mass particles (unit m_e) '
     do ii=1,ndimMD
        do jj=1,nion
           mass_ion(ii,jj)=amas(jj)
           if(ii.eq.1) write(6,*) jj,mass_ion(ii,jj)
        enddo
     enddo
     ! '* rorder' was added by Miha Srdinsek
     allocate(kdyn(nbead * rorder, nbead * rorder), kdyn_eig(nbead * rorder))
     kdyn=0.d0
     ! the following loop was modified by Miha Srdinsek
     do ii=1,(nbead * rorder)
       if (mod(ii, nbead).eq.1) then
         kdyn(ii,ii + 1)=-1.d0
         kdyn(ii + 1,ii)=-1.d0
         kdyn(ii, ii)=1.d0 + renyi_lam_s + renyi_lam_j
       elseif (mod(ii, nbead).eq.0) then
         kdyn(ii,(ii/nbead-1)*nbead + 1)=-1.d0 * renyi_lam_s
         kdyn((ii/nbead-1)*nbead + 1,ii)=-1.d0 * renyi_lam_s

         kdyn(ii,mod(ii,nbead*rorder)+1)=-1.d0 * renyi_lam_j
         kdyn(mod(ii,nbead*rorder)+1,ii)=-1.d0 * renyi_lam_j
         kdyn(ii, ii)=1.d0 + renyi_lam_s + renyi_lam_j
       else
         kdyn(ii,ii + 1)=-1.d0
         kdyn(ii + 1,ii)=-1.d0
         kdyn(ii, ii)=2.d0
       endif
     enddo
!     kdyn=kdyn*temp**2
     kdyn=kdyn*temp**2*mass_ion(1,1)  !! attenzione: qui moltiplico per mass_ion(1,1), ma in realtà scrivo
                                      !! la matrice K come se fosse già divisa per la massa, per essere consistente con
                                      !! le unità (come P che è diviso per sqrt(m) )
     do ii=1,ndimMD
        do jj=1,nion
           ind=ii+(jj-1)*ndimMD
!           scalpar(ind)=1.d0/mass_ion(ii,jj)
           scalpar(ind)=mass_ion(1,1)/mass_ion(ii,jj)
        enddo
     enddo

! Multiplication by '* rorder' added by Miha Srdinsek
     call dsyev('V','L',nbead * rorder,kdyn,nbead * rorder,kdyn_eig,psip,3*nbead * rorder,info)

     write(6,*) ' Eigenvalues elastic quantum term '
! Multiplication by '* rorder' added by Miha Srdinsek
     do ii=1,(nbead*rorder)
        write(6,*) ii,kdyn_eig(ii)
     enddo

     if(averaged_cov .and. sigmacov.ne.0.d0) allocate(cov(ieskin*ieskin))


  elseif(nbead.eq.1) then ! classical molecular dynamic case

     write(6,*) 'classical dynamics not implemented yet!!'
     stop

   !  allocate(tmes_bead(1))
   !  tmes_bead=0.d0
   !  allocate(rcentroid(ndimMD,nion),rcentroidtilde(nion,ndimMD))
   !  allocate(mass_ion(ndimMD,nion))
   !  write(6,*) ' Mass particles (unit m_e) '
   !  do jj=1,nion
   !     do ii=1,ndimMD
   !        mass_ion(ii,jj)=amas(jj)
   !        if(ii.eq.1) write(6,*) jj,mass_ion(ii,jj)
   !     enddo
   !  enddo
     ! define scalpar as inverse mass in Ry
   !  do ii=1,ndimMD
   !     do jj=1,nion
   !        ind=ii+(jj-1)*ndimMD
           !   scalpar(ind)=1.d0/mass_ion(ii,jj)
   !        scalpar(ind)=dt
   !     enddo
   !  enddo

  endif

end subroutine init_pioud_tom
