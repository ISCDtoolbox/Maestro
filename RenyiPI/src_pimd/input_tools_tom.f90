subroutine input_tools_tom

!!! *********************************************************************!!!
!!! initialize all the global variables for the tools                    !!!
!!! *********************************************************************!!!

  use md_variables
  implicit none   
  integer :: i,i1,j,l,k,ii,jj,ind,ntest,natucell
  logical :: ifex
  character(len=20) :: dummy1

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
  read(223,'(a)') output_directory
  lonl=index(output_directory,' ')-1
  close(223)

  write(*,*) 'name of the directory run: ',output_directory(1:lonl)
  write(*,*)

  open(223,file='type_calc',status='old',form='formatted')
  read(223,'(a)') pot_type
  lonll=index(pot_type,' ')-1
  close(223)

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
  
!    ***   Dynamic parameters    ***

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
  

!------------------------------------------------------------------------------
! definition des variables collectives and reading of the MD running parameters
!------------------------------------------------------------------------------
  if (pot_type(1:9).eq.'ab_initio') then

     open(1,file=output_directory(1:lonl)//'/pimd.dat',status='old',form='formatted')
     read(1,nml=dynamics)
     close(1)

     open(1,file=output_directory(1:lonl)//'/cell.dat',status='old',form='formatted')
     read(1,nml=cell_and_atoms)

  else

     open(1,file=output_directory(1:lonl)//'/'//pot_type(1:lonll)//'.in',status='old',form='formatted')
     read(1,nml=system)
     read(1,nml=dynamics)
     read(1,nml=cell_and_atoms)

  end if
  
  avecsp(:,1)=avec(:,1)*dble(ncellsx)
  avecsp(:,2)=avec(:,2)*dble(ncellsy)
  avecsp(:,3)=avec(:,3)*dble(ncellsz)
  
  n=natoms
  write(6,*) 'total number of particles in the system =',n
  
  natucell = n / nunitcells


  allocate(amas(natucell))
  allocate(ion_name(natucell))
  
  amas=0.0d0
  read(1,*)
  nspecies=0
  do i1=1,natucell
    read(1,*) amas(i1),ion_name(i1)
    if (i1.eq.1) then
      nspecies=nspecies+1
    else
      if(ion_name(i1).ne.(ion_name(i1-1))) nspecies=nspecies+1
    end if
  enddo 
  
  close(1)
  
  allocate(indx(n))
  do i=1,n
    indx(i)=i 
  enddo      

  return

end subroutine input_tools_tom
