subroutine input_tom

!!! *********************************************************************!!!
!!! initialize all the global variables - not the specific (to irun) ones!!!
!!! *********************************************************************!!!

  use md_variables
  implicit none   
  integer :: i,j,l,k,ii,jj,ind,ntest
  integer :: myfind_free_unit
  logical :: ifex
  double precision :: dxvar
  character(len=20) :: dummy1,varname1,format_string1
   

! 'system' is the only namelist the QE does not need
  namelist /system/ potential,kspring,c0_k,c1_k,c2_k,    &
                    nspecies,output_dir,verbose

!!! 'dynamics' should be equal to QE interface                     
  namelist /dynamics/ restart_pimd,nbeadMD,nunitcells,nblocks,nstep_block  &
                     ,iprint,irun   &
                     ,delt,tempMD,gammaMD,sigmacov     &
                     ,delta_force,delta_harm     &
                     ,delta0,delta0k,delta0q,yessecond         &
                     ,yesturboq,yesglobal

!!! 'cell_and_atoms' should be equal to QE interface                     
  namelist /cell_and_atoms/ natoms,ndimMD,ncellsx,ncellsy,ncellsz &
                            ,ibrav,avec,numqsymmpoints !!.....amas(in fondo),qsymmpoints(in fondo)

  open(223,file='output_run_dir',status='old',form='formatted')
  read(223,*) output_dir
  lonl=index(output_dir,' ')-1
  close(223)
  
  open(222,file=output_dir(1:lonl)//'/prefix.input',status='old',form='formatted')

  write(*,'(''            *** PIL Program ***                 '' )')
  write(*,'('' Classical/Quantum (N,V,T) molecular dynamics       '')')
  write(*,'('' with Path Integral Langevin approach               '')')

  read(222,*) filen
  lfnl=index(filen,' ')-1
  write(*,*) 'name of the run: ',filen(1:lfnl)
  close(222)
! does this file exist?
  inquire(file=output_dir(1:lonl)//'/'//filen(1:lfnl)//'.in',exist=ifex)
  if(ifex) then
     open(1,file=output_dir(1:lonl)//'/'//filen(1:lfnl)//'.in',status='old',form='formatted')
  else 
     stop ' input file does not exist '
  endif
   

!    ***   System parameters    ***

! Default values defined
  n=1 ! one atom (however this does not belong to 'system' block)
  nspecies=0 ! only 1 type of particle in the system
  kspring=0.d0
  c0_k=0.d0
  c1_k=0.d0
  c2_k=0.d0
  avec(1,1)=100.d0; avec(2,1)=0.d0;   avec(3,1)=0.d0
  avec(1,2)=0.d0;   avec(2,2)=100.d0; avec(3,2)=0.d0
  avec(1,3)=0.d0;   avec(2,3)=0.d0;   avec(3,3)=100.d0
  ibrav=1
  numqsymmpoints=1
  verbose=.false.
  natoms=1
  ndimMD=3
  
! Namelist 'system' call
  read(1,nml=system) 

! Controlling 'system' parameters and stop if they are not correct
  if(nspecies .eq. 0) then
    write(6,*) 'number of species "nspecies" not defined'
    stop
  endif
  if(trim(potential) .eq. 'zundel') then
    write(6,*) 'we are working on Zundel ion'
    write(6,*) 'potential energy parametrized from CC'
  elseif(trim(potential) .eq. 'harmonic_chain1p_1d') then
    write(*,*) 'we are working with a 1d chain coupled by the harmonic constant', kspring
  elseif(trim(potential) .eq. 'harmonic_chain2p_1d') then
    write(*,*) 'we are working with a 2particles- 1d chain coupled by the harmonic constant', kspring
  elseif(trim(potential) .eq. 'harmonic_1d') then 
    write(6,*) 'we are working with the harmonic potential'
    write(6,*) 'harmonic constant =',kspring
  elseif(trim(potential) .eq. 'coupled_harmonic_1d') then 
    write(6,*) 'we are working with two coupled harmonic oscillators in 1d'
    write(6,*) 'harmonic constant =',kspring
  elseif(trim(potential) .eq. 'coupled_harmonic_2d') then
    write(6,*) 'we are working with two coupled harmonic oscillators'
    write(6,*) 'harmonic constants =',kspring*c1_k,kspring*c2_k
  elseif(trim(potential) .eq. 'coupled_sdw_ha_x2y2_2d') then
    write(6,*) 'we are working with a symmetc double well coupled with an harmonic one'
  elseif(trim(potential) .eq. 'coupled_sdw_qu_xy2_2d') then
    write(6,*) 'we are working with a symmetric double well couple dwith a quartic oscillator'
  elseif(trim(potential) .eq. 'coupled_sdw_qu_xy3_2d') then
    write(6,*) 'we are working with a symmetric double well couple dwith a quartic oscillator'      
  elseif(trim(potential) .eq. 'symmetric_double_well_1d') then 
    write(6,*) 'we are working with the double well potential'
  elseif(trim(potential) .eq. 'morse_1d') then
    write(6,*) 'we are working with the Morse potential'
  elseif(trim(potential) .eq. 'quartic_1d') then
    write(6,*) 'we are working with the quartic potential'
  elseif(trim(potential) .eq. 'quartic_plus_harmonic_1d') then
    write(6,*) 'we are working with a general potential'
  elseif(trim(potential) .eq. 'asymmetric_double_well_1d') then 
    write(6,*) 'we are working with the asymmetric double well potential'
  else
    write(6,*) 'error in the choice of the potential!!!'
    stop
  endif


!    ***   Dynamic parameters    ***

! Default values defined
  restart_pimd=.false.
  nh = .false.
  sigmacov = 0.d0
  fixcm = .false.
  nbeadMD=1
  averaged_cov=.false.
  dynmat_cutoff=0.d0
  delta0 = 0.d0
  delta0k = 1.d0
  delta0q = 0.d0
  epstion = 1.d0
  yessecond = .true.
  yesturboq = .true.
  yesglobal = .false. ! PILE_L is the default thermostat for irun=3
  delta_force = 1.d-4       
  delta_harm = 5.d-3
  nunitcells=1
  nblocks=10
  nstep_block=100
  iprint=2
  irun=3
  delt=1.d0
  tempMD=100.d0
  gammaMD=0.00146d0
  ncellsx=1
  ncellsy=1
  ncellsz=1
  
  
  read(1,nml=dynamics)
  read(1,nml=cell_and_atoms)
  
  avecsp(:,1)=avec(:,1)*dble(ncellsx)
  avecsp(:,2)=avec(:,2)*dble(ncellsy)
  avecsp(:,3)=avec(:,3)*dble(ncellsz)
  
  n=natoms
  write(6,*) 'total number of particles in the system =',n

  allocate(amas(n))
  allocate(ion_name(n))
  
  read(1,*)
  do i=1,n/nunitcells
    read(1,*) amas(i),ion_name(i)
  enddo
  
  
  ii=0
  do i=n/nunitcells+1,n
    ii=ii+1  
    amas(i)=amas(ii)
    ion_name(i)=ion_name(ii)
    if(ii.eq.n/nunitcells)ii=0
  end do

  if(irun .eq. 0) then   
    nbeadMD=1
    PS=0.5 ! set initial value for NH
    !!!!!!!!!!!!!!!!!!!!!!!!!! q now is set in the input file
  endif 
  
  if(nbeadMD .eq. 1) then  
    yesquantum=.false.
    yesturboq=.false.   
  else
    yesquantum=.true.
  endif

  if(fixcm) verbose=.true.

 
  
!    ** write input data **
  if(fixcm) then 
    gMD=ndimMD*(n-1)      !   # degrees of freedom
  else 
    gMD=ndimMD*n
  endif
  tempMD=tempMD/kbm1
  tfakeMD=tempMD*nbeadMD
  sigmavar=sigmacov*sigmacov
       
       
  write(*,'('' number of particles       '',i10   )') n
  do i=1,n
    write(*,'('' masses of particles        '',f10.4 )') amas(i)
  enddo
  write(*,'('' number of blocks          '',i10   )') nblocks
  write(*,'('' number of cycles          '',i10   )') nstep
  write(*,'('' output frequency          '',i10   )') iprint
  write(*,'('' initial temperature       '',f10.8 )') tempMD
  write(*,'('' Nose-Hoover               '',l )')     nh
  write(*,'('' Maximum single part. disp.'',f10.6 )') delt
  if(irun == 0) then
    write(*,*) 'Class. MD with Nosé-Hoover thermostat'
    write(unit_dot_out,*) 'Class. MD with Nosé-Hoover thermostat'
  elseif(irun == 3) then
    write(*,*) 'Langevin dyn. with classical/quantum protons (Ceriotti integrator)'
    write(unit_dot_out,*) 'Langevin dyn. with classical/quantum protons (Ceriotti integrator)'
  elseif(irun == 4) then
    write(*,*) 'Langevin dyn. with classical/quantum protons (PIOUD integrator)'
    write(unit_dot_out,*) 'Langevin dyn. with classical/quantum protons (PIOUD integrator)'
  else
    write(*,*) ' type of run unknown: irun = ',irun
    write(unit_dot_out,*) ' type of run unknown: irun = ',irun
    stop
  endif
  write(*,*)

! Allocation of variables
  allocate(indx(n))
  allocate(el(ndimMD,nbeadMD))
  allocate(rpos(ndimMD,n,nbeadMD))
  allocate(rpos_init(ndimMD,n,nbeadMD))
  allocate(forceMD(ndimMD,n,nbeadMD))
  allocate(forcedyn(ndimMD*n,nbeadMD))
  allocate(vel(ndimMD,n,nbeadMD))
  allocate(velocity(3,ndimMD*n,nbeadMD))
  allocate(pimp(ndimMD,n,nbeadMD))
  allocate(rtilde(n,ndimMD))
  allocate(vcm(ndimMD,nbeadMD))
  allocate(rcm(ndimMD,nbeadMD))
  
  
  if(sigmacov .ne. 0.d0) then 
    allocate(dynmat(n*ndimMD,n*ndimMD))
    allocate(dynmat_eig(n*ndimMD),dynmatforce_eig(n*ndimMD))
  endif 


! set indexes of the different atoms for interactions
  do i=1,n
    indx(i)=i 
  enddo
        

! set the total mass of the zundel ion
  mtot=0.d0
  do i=1,n
    mtot=mtot+amas(indx(i))
  enddo 

!    ** read initial configuration **
  if (potential .eq. 'zundel') then
    write(*,*)' restart file not found'
    write(*,*)' starting from input file configuration'
    if(.not.restart_pimd) call readc2
    if(.not.restart_pimd) call setvel(0)
    S=0.d0
    PS=0.d0
    ifl=0
  else if(potential .eq. 'harmonic_chain1p_1d') then
    if(n.lt.10) then          
        format_string1 = "(I1)"
    else if (n.gt.9 .and. n.lt.100) then     
        format_string1 = "(I2)"
    else
        format_string1 = "(I3)"  
    endif
    write (varname1,format_string1) n

    dxvar=avecsp(1,1)/dble(n)
    if(.not.restart_pimd) then 
      open(89,file='input_files_dir/configinit_chain1p'//trim(varname1)//'.xyz',status='replace')
      write(89,*) n
      write(89,*) avecsp(1,1),avecsp(2,2),avecsp(3,3)
      do k=1,nbeadMD
        do i=1,n
          do l=1,ndimMD
            if(l.eq.1) then
               rpos(l,i,k) = (i-1)*dxvar+dxvar*0.5
            else
               rpos(l,i,k) = 0.d0
            end if
          enddo
          if(k.eq.1) write(89,*) ion_name(i),rpos(:,i,1)
        enddo
      enddo
      close(89)

      rpos_init=rpos
      call setvel(0)
    
      S=0.d0
      PS=0.d0
    end if
  else if(potential .eq. 'harmonic_chain2p_1d') then
    if(n.lt.10) then
        format_string1 = "(I1)"
    else if (n.gt.9 .and. n.lt.100) then
        format_string1 = "(I2)"
    else
        format_string1 = "(I3)"
    endif
    write (varname1,format_string1) n

    dxvar=avecsp(1,1)/dble(n)
    if(.not.restart_pimd)  then
      open(89,file='input_files_dir/configinit_chain2p'//trim(varname1)//'.xyz',status='replace')
      write(89,*) n
      write(89,*) avecsp(1,1),avecsp(2,2),avecsp(3,3)
      do k=1,nbeadMD
        do i=1,n
          do l=1,ndimMD
            if(l.eq.1) then
               rpos(l,i,k) = (i-1)*dxvar+dxvar*0.5
            else
               rpos(l,i,k) = 0.d0
            end if
          enddo
          if(k.eq.1 .and. mod(i,2).ne.0) write(89,*) ion_name(i),rpos(:,i,1)
          if(k.eq.1 .and. mod(i,2).eq.0) write(89,*) ion_name(i),rpos(:,i,1)
        enddo
      enddo
      close(89)

      rpos_init=rpos
      call setvel(0)
      S=0.d0
      PS=0.d0
    end if
  else
    if (.not.restart_pimd) then
      rpos=0.d0
      rpos_init=rpos
      call setvel(0)
      S=0.d0
      PS=0.d0
    end if
  endif


  if(.not.restart_pimd) then
     
     unit_dot_out  = myfind_free_unit()
     open(unit_dot_out,file=output_dir(1:lonl)//'/pimd.out',form='formatted')
     unit_dot_ep   = myfind_free_unit()
     open(unit_dot_ep,file=output_dir(1:lonl)//'/pimd.ep',form='formatted')
     unit_dot_ep2  = myfind_free_unit()
     open(unit_dot_ep2,file=output_dir(1:lonl)//'/pimd.ep2',form='formatted')
     unit_dot_epsr = myfind_free_unit()
     open(unit_dot_epsr,file=output_dir(1:lonl)//'/pimd.epsr',form='formatted')
     unit_dot_eplr = myfind_free_unit()
     open(unit_dot_eplr,file=output_dir(1:lonl)//'/pimd.eplr',form='formatted')
     unit_dot_ek   = myfind_free_unit()
     open(unit_dot_ek,file=output_dir(1:lonl)//'/pimd.ek',form='formatted')
  
     write(unit_dot_out,'('' number of blocks          '',i10   )') nblocks
     write(unit_dot_out,'('' number of per block       '',i10   )') nstep_block
     write(unit_dot_out,'('' output frequency          '',i10   )') iprint
     write(unit_dot_out,'('' temperature               '',f10.8 )') tempMD
     write(unit_dot_out,'('' Nose-Hoover               '',l )') nh
     write(unit_dot_out,'('' Maximum single part. disp.'',f10.6 )') delt

     write(unit_dot_out,*)
     
     if(irun == 0) then
       write(unit_dot_out,*) 'Classical MD with Nosé-Hoover thermostat'
     elseif(irun == 3) then
       write(unit_dot_out,*) 'Classical/Quantum Langevin dyn. with Ceriotti integrator'
     elseif(irun == 4) then
       write(unit_dot_out,*) 'Quantum Langevin dyn. with PIOUD integrator'
     else
       write(unit_dot_out,*) ' type of run unknown: irun = ',irun
       stop
     endif
     write(unit_dot_out,*)  
     flush(unit_dot_out)
  
     unit_dot_xyz  = myfind_free_unit()
     open(unit_dot_xyz,file=output_dir(1:lonl)//'/pimd.xyz',form='formatted')
     unit_dot_positions   = myfind_free_unit()
     open(unit_dot_positions,file=output_dir(1:lonl)//'/positions.dat',form='formatted')
     unit_dot_positions_cen   = myfind_free_unit()
     open(unit_dot_positions_cen,file=output_dir(1:lonl)//'/positions_cen.dat',form='formatted')
     unit_dot_velocities  = myfind_free_unit()
     open(unit_dot_velocities,file=output_dir(1:lonl)//'/velocities.dat',form='formatted')
     unit_dot_velocities_cen  = myfind_free_unit()
     open(unit_dot_velocities_cen,file=output_dir(1:lonl)//'/velocities_cen.dat',form='formatted')
     unit_dot_forces = myfind_free_unit()
     open(unit_dot_forces,file=output_dir(1:lonl)//'/forces.dat',form='formatted')
     unit_dot_forces_cen = myfind_free_unit()
     open(unit_dot_forces_cen,file=output_dir(1:lonl)//'/forces_cen.dat',form='formatted')
     unit_dot_localtemp = myfind_free_unit()
     open(unit_dot_localtemp,file=output_dir(1:lonl)//'/local_temp.dat',form='formatted')
     unit_dot_sigma   = myfind_free_unit()
     open(unit_dot_sigma,file=output_dir(1:lonl)//'/sigma.dat',form='formatted')
     rewind(unit_dot_xyz)
     rewind(unit_dot_positions)
     rewind(unit_dot_positions_cen)
     rewind(unit_dot_velocities)
     rewind(unit_dot_velocities_cen)
     rewind(unit_dot_forces)
     rewind(unit_dot_forces_cen)
     rewind(unit_dot_localtemp)
     rewind(unit_dot_sigma)
  
  else
     
     unit_dot_out  = myfind_free_unit()
     open(unit_dot_out,file=output_dir(1:lonl)//'/pimd.out',position='APPEND',form='formatted')
     unit_dot_ep  = myfind_free_unit()
     open(unit_dot_ep,file=output_dir(1:lonl)//'/pimd.ep',position='APPEND',form='formatted')
     unit_dot_ep2  = myfind_free_unit()
     open(unit_dot_ep2,file=output_dir(1:lonl)//'/pimd.ep2',position='APPEND',form='formatted')
     unit_dot_epsr  = myfind_free_unit()
     open(unit_dot_epsr,file=output_dir(1:lonl)//'/pimd.epsr',position='APPEND',form='formatted')
     unit_dot_eplr  = myfind_free_unit()
     open(unit_dot_eplr,file=output_dir(1:lonl)//'/pimd.eplr',position='APPEND',form='formatted')
     unit_dot_ek  = myfind_free_unit()
     open(unit_dot_ek,file=output_dir(1:lonl)//'/pimd.ek',position='APPEND',form='formatted')
     unit_dot_xyz  = myfind_free_unit()
     open(unit_dot_xyz,file=output_dir(1:lonl)//'/pimd.xyz',position='APPEND',form='formatted')
     
     unit_dot_positions  = myfind_free_unit()
     open(unit_dot_positions,file=output_dir(1:lonl)//'/positions.dat',form='formatted') !!! not append because I need to read last pos
     unit_dot_velocities  = myfind_free_unit()
     open(unit_dot_velocities,file=output_dir(1:lonl)//'/velocities.dat',form='formatted') !!! same as pos

     unit_dot_positions_cen  = myfind_free_unit()
     open(unit_dot_positions_cen,file=output_dir(1:lonl)//'/positions_cen.dat',position='APPEND',form='formatted')      
     unit_dot_velocities_cen  = myfind_free_unit()
     open(unit_dot_velocities_cen,file=output_dir(1:lonl)//'/velocities_cen.dat',position='APPEND',form='formatted') 
     unit_dot_forces  = myfind_free_unit()
     open(unit_dot_forces,file=output_dir(1:lonl)//'/forces.dat',position='APPEND',form='formatted')
     unit_dot_forces_cen  = myfind_free_unit()
     open(unit_dot_forces_cen,file=output_dir(1:lonl)//'/forces_cen.dat',position='APPEND',form='formatted')
     unit_dot_localtemp  = myfind_free_unit()
     open(unit_dot_localtemp,file=output_dir(1:lonl)//'/local_temp.dat',position='APPEND',form='formatted')
     unit_dot_sigma  = myfind_free_unit()
     open(unit_dot_sigma,file=output_dir(1:lonl)//'/sigma.dat',position='APPEND',form='formatted')
      
  endif
  
  if (restart_pimd) call pimd_restart_traj()

  do k=1,nbeadMD
    do i=1,n
      ind=indx(i)
      do l=1,ndimMD
        pimp(l,i,k)=amas(ind)*vel(l,i,k)
      enddo
    enddo
  enddo

  return

end subroutine input_tom



FUNCTION myfind_free_unit()
    !--------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER :: myfind_free_unit
    INTEGER :: iunit
    LOGICAL :: opnd
    !
    !
    myfind_free_unit = -1
    unit_loop: DO iunit = 999, 1, -1
       !
       INQUIRE( UNIT = iunit, OPENED = opnd )
       !
       IF ( .NOT. opnd ) THEN
          !
          myfind_free_unit = iunit
          !
          RETURN
          !
       END IF
       !
    END DO unit_loop
    !
    RETURN
    !
END FUNCTION myfind_free_unit



subroutine pimd_restart_traj()
 
 ! use path_input_parameters_module, only : pos
  use md_variables, only : rpos,rpos_init,vel,nbeadMD,natoms,ndimMD,output_dir,lonl
  use md_variables, only : unit_dot_positions, unit_dot_velocities
  
  implicit none
  integer cc,k,iat,i,iflagerr,ngen
  real(8), allocatable :: vec_tmp(:)
  
  allocate(vec_tmp(natoms*ndimMD))

  !count the number of total available snapshots
  ngen=0
  rewind(unit_dot_positions)
  rewind(unit_dot_velocities)
  do while(.true.)
    read(unit_dot_positions,*,iostat=iflagerr)
    if(iflagerr/=0) exit
    ngen=ngen+1
  enddo
  rewind(unit_dot_positions)
  
  do i=1,ngen-nbeadMD
    READ(unit_dot_positions,*)
    READ(unit_dot_velocities,*)
  end do
  
  rpos=0.0
  do k=1,nbeadMD
    vec_tmp=0.0
    read(unit_dot_positions,*) vec_tmp(:)
    cc=0
    DO iat=1,natoms
      DO i=1,ndimMD
        cc=cc+1
        rpos(i,iat,k)=vec_tmp(cc)
      END DO
    END DO
  end do
  rpos_init=rpos

  vel=0.0
  do k=1,nbeadMD
    vec_tmp=0.0
    read(unit_dot_velocities,*) vec_tmp(:)
    cc=0
    DO iat=1,natoms
      DO i=1,ndimMD
        cc=cc+1
        vel(i,iat,k)=vec_tmp(cc)
      END DO
    END DO
  end do
  
  close(unit_dot_positions)
  close(unit_dot_velocities)
  open(unit_dot_positions,file=output_dir(1:lonl)//'/positions.dat',position='APPEND',form='formatted')
  open(unit_dot_velocities,file=output_dir(1:lonl)//'/velocities.dat',position='APPEND',form='formatted')
  

  deallocate (vec_tmp)

  
  return
end subroutine pimd_restart_traj
