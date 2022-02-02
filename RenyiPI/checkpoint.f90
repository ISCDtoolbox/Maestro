subroutine checkpoint(ttk)
        
  use md_variables, only : nbeadMD,nunitcells,unit_dot_xyz,&
                             forceMD,vel,natoms,ion_name,rcentroid,&
                             ndimMD,unit_dot_positions,&
                             unit_dot_velocities,unit_dot_forces,&
                             unit_dot_localtemp,rpos,&
                             unit_dot_positions_cen,&
                             unit_dot_forces_cen,&
                             unit_dot_velocities_cen
  
  implicit none
!    *******************************************************************
!    ** subroutine to write out the configurations and trajectories   **
!    *******************************************************************

  integer :: i,l,k
  real(8) :: ttk
  real(8), allocatable :: fcentroid(:,:),vcentroid(:,:)
!    ********************************************************************
  allocate(fcentroid(ndimMD,natoms),vcentroid(ndimMD,natoms))

  if(nbeadMD.gt.1) then

     !if(nunitcells.gt.1) then
       fcentroid=0.0; vcentroid=0.0; rcentroid=0.0
       do k=1,nbeadMD
         fcentroid(:,:)=fcentroid(:,:)+forceMD(:,:,k)
         vcentroid(:,:)=vcentroid(:,:)+vel(:,:,k)
         rcentroid(:,:)=rcentroid(:,:)+rpos(:,:,k)
       end do
       fcentroid=fcentroid/nbeadMD; vcentroid=vcentroid/nbeadMD; rcentroid=rcentroid/nbeadMD
     !end if

!!! -------------------write the xyz file for visualizing trajectories------     
     write(unit_dot_xyz,* ) natoms
     write(unit_dot_xyz,* ) 
     do i=1, natoms
        write(unit_dot_xyz,'(a2,3g15.5)' ) ion_name(i),(rcentroid(l,i), l=1,ndimMD)
     enddo
     flush(unit_dot_xyz)
!!!--------------------writes velocity,force and position files for postprocessing----------
!!!--------------------if the number of unit cells > 1 (crystal system)---------------------
!!!--------------------writes only the centroid coordinatomses---------------------------------
     !if(nunitcells.gt.1) then
         write(unit_dot_positions_cen,'(400e15.5)') ((rcentroid(l,i),l=1,ndimMD),i=1,natoms)
         flush(unit_dot_positions_cen)
         write(unit_dot_velocities_cen,'(400e15.5)') ((vcentroid(l,i),l=1,ndimMD),i=1,natoms)
         flush(unit_dot_velocities_cen)
         write(unit_dot_forces_cen,'(400e15.5)') ((fcentroid(l,i),l=1,ndimMD),i=1,natoms)
         flush(unit_dot_forces_cen)
       !  write(unit_dot_localtemp,'(e15.5)') ttk
       !  flush(unit_dot_localtemp)
     !end if
     
     do k=1,nbeadMD
           write(unit_dot_positions,'(400e15.5)') ((rpos(l,i,k),l=1,ndimMD),i=1,natoms)
           flush(unit_dot_positions)
           write(unit_dot_velocities,'(400e15.5)') ((vel(l,i,k),l=1,ndimMD),i=1,natoms)
           flush(unit_dot_velocities)
           write(unit_dot_forces,'(400e15.5)') ((forceMD(l,i,k),l=1,ndimMD),i=1,natoms)
           flush(unit_dot_forces)
           if(k.eq.1) write(unit_dot_localtemp,'(e15.5)') ttk
           if(k.eq.1) flush(unit_dot_localtemp)
     enddo
     
  else !!! nbeadMD=1 => classical particles
     
     write (unit_dot_xyz,*) natoms
     write (unit_dot_xyz,*) 
     do i=1, natoms
        write (unit_dot_xyz,'(a2,3g15.5)' ) ion_name(i),(rpos(l,i,1), l=1,ndimMD)
     enddo
     flush(unit_dot_xyz)
     
     write(unit_dot_positions,'(400e15.5)') ((rpos(l,i,1),l=1,ndimMD),i=1,natoms)
     flush(unit_dot_positions)
     write(unit_dot_velocities,'(400e15.5)') ((vel(l,i,1),l=1,ndimMD),i=1,natoms)
     flush(unit_dot_velocities)
     write(unit_dot_forces,'(400e15.5)') ((forceMD(l,i,1),l=1,ndimMD),i=1,natoms)
     flush(unit_dot_forces)
     write(unit_dot_localtemp,'(e15.5)') ttk
     flush(unit_dot_localtemp)

  endif

  deallocate(fcentroid,vcentroid) 

  return
end subroutine checkpoint
