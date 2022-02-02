subroutine blockMD_tom(iblockMD,ep)
      
  use md_variables
  implicit none

  integer :: step,ii,iblockMD
  real(8) :: tt,enh,h,en,vs,ep,elr,enk,ttk,ep_centroid
!  real(8) :: eq, eqp

  enh=0.d0
  do step = 1,nstep_block                ! do one block

     if (irun == 0) then 
        call vel_verlet_tom(enk,ep,enh)
     elseif (irun == 3) then
        call prop_ceriotti_tom(enk,ep)  
     elseif (irun == 4) then
          call prop_pioud_tom(enk,ep)
     endif

     en = ep + enk  

     if (irun .eq. 0) then  
        h = en + enh
     !   write(*,*) step,'the conserved quantity is',h
     !   write(*,*) 'vcm',vcm
     !   write(*,*) 'rcm',rcm
     !   write(*,*) '----------------------------------------'
     elseif (irun .eq. 3) then 
        h = deltahtilde
     else 
        h = en
     endif
     

!  ** calculate potential energy per dumbell **
     vs=ep/n !! questa energia in teoria è l'energia potenziale media per atomo (atomo vero si intende)
      
!  ** accumulate step averages **
     call cumul1(vs,avp(ipot(1)),anormp(ipot(1)),1)
     call cumul1(vs**2,avp(ipot(2)),anormp(ipot(2)),1)
     call cumul1(ep,avp(ipot(3)),anormp(ipot(3)),1)
     call cumul1(elr,avp(ipot(4)),anormp(ipot(4)),1)
     call cumul1(enk/n,avp(ikin),anormp(ikin),1)
 
!  ** perform periodic operations  **
     if ( mod ( step, iprint ) .eq. 0 ) then 
           
        ttk=2.d0*enk/(gMD*nbeadMD)*kbm1 ! Temperature in Kelvin
        tt=2.d0*enk/(gMD*nbeadMD)  ! Temperature in Hartree

        call checkpoint(ttk)
        
     !   if (mod(step,100) .eq. 0) call checkpoint_tdep

        if(nbeadMD.gt.1) then
!  ** write out runtime information in quantum case**
          write(*,'(a3,i4,a6,i6,a6,g15.7,a6,g15.7,a6,g15.7,a6,g15.7,a6,g15.7,a6,g15.7,a6,g15.7,a6,g15.7)') & 
                    '1)',iblockMD,'2)',step,'3)',h,'4)',enh,'5)',ep,'6)',enk,'7)',tt,&
                    '8)',ttk,'9)',ekinq,'10)',ekinqp
          write(unit_dot_out,'(2i8,8g20.8)') iblockMD,step,h,enh,ep,&
                enk,tt,ttk,ekinq,ekinqp
        else
!  ** write out runtime information in classical case**
          write(*,'(a3,i4,a6,i6,a6,g15.7,a6,g15.7,a6,g15.7,a6,g15.7,a6,g15.7,a6,g15.7)') &
                   '1)',iblockMD,'2)',step,'3)',h,'4)',enh,'5)',ep,'6)',enk,'7)',tt,'8)',ttk
          write(unit_dot_out,'(2i8,6g20.8)') iblockMD,step,h,enh,ep,&
                enk,tt,ttk
        endif

     endif

     if(irun.eq.3) then
! htilde is computed only for Ceriotti integrator
! Saving potential energy, forces and positions to evaluate deltaHtilde
        rpos_old=rpos
        forceMD_old=forceMD
        epot_old=ep
     endif 
     
  enddo            ! over steps
  
  return

end subroutine blockMD_tom
