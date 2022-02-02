!    **************************************************************************
!    ** molecular dynamics simulation program in the constant-nvt ensemble.  **
!    **************************************************************************

program mdnvt_tom
  
  use md_variables
   
  implicit none 

  real(8) :: tin,tstart,tend
  real(8) :: epot, ekin, en, h
  real(8) :: tt, ttk,vs
  integer :: i, l, iblockMD,k
  integer :: ii,jj,kk,ind
  
  integer :: nseed, clock
  integer, allocatable :: seed(:)
  
  call random_seed(size=nseed)
  allocate(seed(nseed))
  call system_clock(count=clock)
  do i=1,nseed
    seed(i) = clock + 17*i
  end do
  call random_seed(put=seed)

! ****************************************************************

! ** read input data **
  call scnd(tin)
  call input_tom   

! ** set pointers for the averages and zerod accumulators **

  call zeroav(0)         ! set to zero the cumulators
  
!  call fire_k

! ** irun==0 --->   classical MD velocity Verlet algorithm with Nosé-Hoover thermostat !!
! ** irun==3 --->   classical and quantum Langevin MD with Ceriotti integrator (DOI: https://doi.org/10.1063/1.3489925)
! ** irun==4 --->   classical and quantum Langevin MD with TurboCeriotti
!                   the quantum harmonic part is done WITHOUT Trotter breakup 
!                   (exact integration of both harmonic dynamics AND Langevin thermostatting) ---> PIOUD 

      
   if(irun .eq. 3) then
      call init_ceriotti_tom
         
   elseif(irun .eq. 4) then
      call init_pioud_tom

   endif


!!!! **** compute the initial potential & forces... **** !!!!
   
   call force0(forceMD)
   call pot(epot,epot_centroid)
  
   !write(*,*) forceMD, epot
   !STOP

   if (irun .eq. 3) forceMD_old=forceMD
   write(*,*) epot,epot_centroid,'++++++++++++++++++++++++++++++++++++++++++++++++++++'

!!!! **** ...and kinetic energy averaged over quantum images **** !!!!
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

  en=ekin+epot
  h=en

  ttk=2.d0*ekin/gMD*kbm1 ! Temperature in Kelvin
  tt=2.d0*ekin/gMD    ! Temperature in Hartree

  ttk=ttk/nbeadMD
  tt=tt/nbeadMD

  if(yesquantum) then
! initial EXPECTED value of Quantum kinetic energy calculated in two ways

! virial
     cost=0.d0
     do ii=1,nbeadMD
        fbead(:,:)=rcentroid(:,:)-rpos(:,:,ii)
        do jj=1,n
           do kk=1,ndimMD
              cost=cost+0.5d0*fbead(kk,jj)*forceMD(kk,jj,ii) ! Ha units!!
           enddo
        enddo
     enddo
     ekinq=cost/nbeadMD+ndimMD*0.5d0*tempMD*n/nbeadMD

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
     ekinqp=-ekinqp/nbeadMD*tempMD**2*0.5d0+ndimMD*0.5d0*n*tempMD ! Ha units!!!
        
  endif


  call scnd(tstart)
  write (*,*) ' initialization time: ',tstart-tin
  write (unit_dot_out,*) ' initialization time: ',tstart-tin
      
  vs = epot/n
  write(*,'('' initial v =  '', g20.10 )' ) vs
  write(*,*) 'quantities after initialization'
  if(nbeadMD.eq.1) then
     write(*,'(//''1)block   2)nmove   3)H   4)H_Nosé   5)E_pot'' &
         ''   6)E_kin   7)Temp(Ha)   8)Temp(K)''//)')
     write(*,'(a3,i2,a6,i2,a6,g15.7,a6,g15.7,a6,g15.7,a6,g15.7,a6,g15.7,a6,g15.7)') &
              '1)',0,'2)',0,'3)',h,'4)',en,'5)',epot,'6)',ekin,'7)',tt,'8)',ttk
  else
     write(*,'(//''1)block   2)nmove   3)H   4)H_Nosé   5)E_pot'' &
           ''   6)E_kin   7)Temp(Ha)   8)Temp(K)'' &
           ''   9)quantum_kin_virial   10)quantum_kin_primitive  ''//)')
     write(*,'(a3,i2,a6,i2,a6,g15.7,a6,g15.7,a6,g15.7,a6,g15.7,a6,g15.7,a6,g15.7,a6,g15.7,a6,g15.7)') & 
               '1)',0,'2)',0,'3)',h,'4)',en,'5)',epot,'6)',ekin,'7)',tt,&
                             '8)',ttk,'9)',ekinq,'10)',ekinqp
  endif
  
  write(*,*)
  write(*,*)'***********************************************************'
  write(*,*)'**************** start of dynamics ************************'
  
  if(nbeadMD.eq.1) then
    write(*,'(//''1)block   2)nmove   3)H   4)H_Nosé   5)E_pot'' &
         ''   6)E_kin   7)Temp(Ha)   8)Temp(K)''//)')
  else
    write(*,'(//''1)block   2)nmove   3)H   4)H_Nosé   5)E_pot'' &
           ''   6)E_kin   7)Temp(Ha)   8)Temp(K)'' &
           ''   9)quantum_kin_virial   10)quantum_kin_primitive  ''//)')
  endif
  
  if(.not.restart_pimd) then
    write(unit_dot_out,'('' initial potential =  '', f10.4 )' ) epot
    write(unit_dot_out,'(//'' start of dynamics ''//)')
  
    if(nbeadMD.eq.1) then
      write(unit_dot_out,'(//''1)block   2)nmove   3)H   4)H_Nosé   5)E_pot'' &
         ''   6)E_kin   7)Temp(Ha)   8)Temp(K)''//)')
    else
      write(unit_dot_out,'(//''1)block   2)nmove   3)H   4)H_Nosé   5)E_pot'' &
           ''   6)E_kin   7)Temp(Ha)   8)Temp(K)'' &
           ''   9)quantum_kin_virial   10)quantum_kin_primitive  ''//)')
    endif
  end if

!    *******************************************************************
!    **               loops over blocks                               **
!    *******************************************************************
  do iblockMD = 1, nblocks
     call zeroav(1)                    ! zeroed the block cumulators
     call blockMD_tom(iblockMD,epot)
 
! ** block analysis  
     call cumul(avp(ipot(1)),av(ipot(1)),anorm(ipot(1)),1)
     call cumul(avp(ipot(2)),av(ipot(2)),anorm(ipot(2)),1)
     call cumul(avp(ipot(3)),av(ipot(3)),anorm(ipot(3)),1)
     call cumul(avp(ipot(4)),av(ipot(4)),anorm(ipot(4)),1)
     call cumul(avp(ikin),av(ikin),anorm(ikin),1)
     call sprint(unit_dot_ep,ipot(1),1,'epot')
     flush(unit_dot_ep)
     call sprint(unit_dot_ep2,ipot(2),1,'ept2')
     flush(unit_dot_ep2)
     call sprint(unit_dot_epsr,ipot(3),1,'ept2')
     flush(unit_dot_epsr)
     call sprint(unit_dot_eplr,ipot(4),1,'ept2')
     flush(unit_dot_eplr)
     call sprint(unit_dot_ek,ikin,1,'ekin')
     flush(unit_dot_ek)

!    *******************************************************************
!    ** ends the loop over cycles                                     **
!    *******************************************************************
  enddo              !    over blocks


  call scnd(tend)

  write(*,'(//'' end of markov chain    ''//)')
  write(*,'(//'' total CPU time (sec)   '',g20.10,//)') tend-tstart
  write(unit_dot_out,'(//'' end of markov chain    ''//)')
  write(unit_dot_out,'(//'' total CPU time (sec)   '',g20.10,//)') tend-tstart

!    ** Deallocation of dynamical variables allocated in input.f90 and zeroav.f90**
  deallocate(indx)
  deallocate(amas)
  deallocate(rpos)
  deallocate(rpos_init)
  deallocate(forceMD)
  deallocate(vel)
  deallocate(pimp)
  deallocate(rtilde)
  deallocate(av)
  deallocate(avp)
  deallocate(anorm)
  deallocate(anormp)
  deallocate(ipot)
  deallocate(el)
  deallocate(vcm)
  deallocate(rcm)
  deallocate(ion_name)
  deallocate(seed)
 
  if (sigmacov .ne. 0.d0) then
     deallocate(dynmat,dynmat_eig,dynmatforce_eig)
  endif 
  
  
  call dealloc  !! deallocate all the other vectors specific to different iruns
  
  
!    ** calculate and write out running averages **

  write(*,'(/'' end of simulation '')')
  write(unit_dot_out,'(/'' end of simulation '')')

  close(unit_dot_ek)
  close(unit_dot_ep)
  close(unit_dot_ep2)
  close(unit_dot_eplr)
  close(unit_dot_epsr)
  close(unit_dot_forces)
  close(unit_dot_forces_cen)
  close(unit_dot_localtemp)
  close(unit_dot_out)
  close(unit_dot_positions)
  close(unit_dot_positions_cen)
  close(unit_dot_sigma)
  close(unit_dot_velocities)
  close(unit_dot_velocities_cen)
  close(unit_dot_xyz)

 ! call rmaget(0)

  stop
  
end program mdnvt_tom
