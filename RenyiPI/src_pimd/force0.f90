subroutine force0(force)
  
  use md_variables
  implicit none
  integer :: i,l,k
  real(8) :: pot_en_minus,pot_en_plus
  real(8) :: delta
  real(8) :: force(ndimMD,n,nbeadMD)
  
  force=0.d0
  
  if(trim(potential) .eq. 'zundel') then

! *******************************************************************
! *******************ZUNDEL POTENTIAL********************************
!    calcpot(V,xx)
!    xx(7,3) is cartesian array containing O O H H H H H coor, in bohr
!    returned V is potential in hartree, taking C2-sym minimum as zero potential reference

    call prepot() 
    delta=delta_force
    do k=1,nbeadMD
       do i=1,n 
          do l=1,ndimMD
             rtilde(i,l) = rpos(l,i,k)
          enddo
       enddo
       do i=1,n 
          do l=1,ndimMD
             rtilde(i,l) = rpos(l,i,k) + delta
             call calcpot(pot_en_plus,rtilde)
             rtilde(i,l) = rpos(l,i,k) - delta
             call calcpot(pot_en_minus,rtilde)
             force(l,i,k) = (pot_en_minus-pot_en_plus)/2.d0/delta
             rtilde(i,l) = rpos(l,i,k)
          enddo
       enddo
    enddo

  elseif(trim(potential) .eq. 'harmonic_chain1p_1d') then
    
    do k=1,nbeadMD
      do i=1,n
        do l=1,ndimMD
          if(l.eq.1) then
            if(i.ne.1 .and. i.ne.n)  force(l,i,k) = -kspring*(2*rpos(l,i,k)-rpos(l,i+1,k)-rpos(l,i-1,k))
            if(i.eq.1) force(l,i,k) = -kspring*(2*rpos(l,i,k)-rpos(l,i+1,k)-(rpos(l,n,k)-avecsp(1,1)))
            if(i.eq.n) force(l,i,k) = -kspring*(2*rpos(l,i,k)-rpos(l,i-1,k)-(rpos(l,1,k)+avecsp(1,1)))
          else
            force(l,i,k) = -0.7d0*kspring*rpos(l,i,k)
          end if
        enddo
      enddo
    enddo

  elseif(trim(potential) .eq. 'harmonic_chain2p_1d') then

    do k=1,nbeadMD
      do i=1,n
        do l=1,ndimMD
          if(l.eq.1) then
            if(i.ne.1 .and. i.ne.n)  force(l,i,k) = -kspring*(2*rpos(l,i,k)-rpos(l,i+1,k)-rpos(l,i-1,k))
            if(i.eq.1) force(l,i,k) = -kspring*(2*rpos(l,i,k)-rpos(l,i+1,k)-(rpos(l,n,k)-avecsp(1,1)))
            if(i.eq.n) force(l,i,k) = -kspring*(2*rpos(l,i,k)-rpos(l,i-1,k)-(rpos(l,1,k)+avecsp(1,1)))
        !!  else
        !!    force(l,i,k) = -0.7d0*kspring*rpos(l,i,k)
          end if
        enddo
      enddo
    enddo


  elseif(trim(potential) .eq. 'harmonic_1d') then
    
    do k=1,nbeadMD
      do i=1,n 
        do l=1,ndimMD
        if(l.eq.1) force(l,i,k) = -kspring*rpos(l,i,k)
        if(l.ne.1) force(l,i,k) = -0.7d0*kspring*rpos(l,i,k)
        enddo
      enddo
    enddo 
  
  elseif(trim(potential) .eq. 'coupled_harmonic_1d') then
        
    do k=1,nbeadMD
      do i=1,n 
        do l=1,ndimMD
           if (i.eq.1) then
          !!   force(l,i,k) = -c0_k*kspring*rpos(l,i,k) - c1_k*sqrt(c0_k*kspring**2)*rpos(l,1,k)
             if(l.eq. 1) then 
               force(l,i,k) = -kspring*rpos(l,i,k) + 2*c0_k*c1_k*(rpos(l,1,k)-rpos(l,n,k))*&
                                                     exp(-c0_k*(rpos(l,n,k)-rpos(l,i,k))**2)
             else
               force(l,i,k) = -kspring*rpos(l,i,k)
             end if
           else
          !!   force(l,i,k) = -kspring*rpos(l,i,k) - c1_k*sqrt(c0_k*kspring**2)*rpos(l,i+1,k)
             if(l.eq. 1) then
               force(l,i,k) = -kspring*rpos(l,i,k) + 2*c0_k*c1_k*(rpos(l,n,k)-rpos(l,1,k))*&
                                                     exp(-c0_k*(rpos(l,1,k)-rpos(l,i,k))**2)
             else
               force(l,i,k) = -kspring*rpos(l,i,k)
             end if
              
            end if
        enddo
      enddo
    enddo 

  elseif(trim(potential) .eq. 'coupled_harmonic_2d') then

    do k=1,nbeadMD
      do i=1,n
        do l=1,ndimMD
          
           if(i.eq.1) then
             if(l.eq. 1) then
               force(l,i,k) = -c1_k*0.8*kspring*rpos(1,1,k) - (c0_k*0.8*kspring)*rpos(2,2,k)  !!! here c0_k is a scale factor with respect
                                                                                              !!! of kspring
             else         
               force(l,i,k) = -c1_k*kspring*rpos(l,i,k)
             end if
           else
             if(l.eq. 2) then
               force(l,i,k) = -c2_k*0.8*kspring*rpos(2,2,k) - (c0_k*0.8*kspring)*rpos(1,1,k)
             else
               force(l,i,k) = -c2_k*kspring*rpos(l,i,k)
             end if
           end if
        
        enddo
      enddo
    enddo

  elseif(trim(potential) .eq. 'coupled_sdw_ha_x2y2_2d') then

    do k=1,nbeadMD
      do i=1,n
        do l=1,ndimMD

           if(i.eq.1) then
             if(l.eq. 1) then
               force(l,i,k) = -(kspring*(4*(c2_k**2)*rpos(1,1,k)**3 + 4*c1_k*c2_k*rpos(1,1,k))) &
                              - 2.d0*c0_k*rpos(1,1,k)*rpos(2,2,k)**2
             else
               force(l,i,k) = -kspring*rpos(l,i,k)
             end if
           else
             if(l.eq. 2) then
               force(l,i,k) = -(0.89114835d0)*kspring*rpos(2,2,k) - 2.d0*c0_k*(rpos(1,1,k)**2)*rpos(2,2,k)
             else
               force(l,i,k) = -kspring*rpos(l,i,k)  !!! the factor 0.89114835 is to distinguish the 2 ksprings
             end if
           end if

        enddo
      enddo
    enddo
  
  elseif(trim(potential) .eq. 'coupled_sdw_qu_xy2_2d') then

    do k=1,nbeadMD
      do i=1,n
        do l=1,ndimMD

           if(i.eq.1) then
             if(l.eq. 1) then
               force(l,i,k) = -(kspring*(4*(c2_k**2)*rpos(1,1,k)**3 + 4*c1_k*c2_k*rpos(1,1,k))) &
                              - c0_k*rpos(2,2,k)**2
             else
               force(l,i,k) = -kspring*rpos(l,i,k)
             end if
           else
             if(l.eq. 2) then
               force(l,i,k) = -16.d0*kspring*rpos(2,2,k)**3 - 2.d0*c0_k*rpos(1,1,k)*rpos(2,2,k)  !!! quartic potential with 4*kspring
             else
               force(l,i,k) = -kspring*rpos(l,i,k)
             end if
           end if

        enddo
      enddo
    enddo

  elseif(trim(potential) .eq. 'coupled_sdw_qu_xy3_2d') then

    do k=1,nbeadMD
      do i=1,n
        do l=1,ndimMD

           if(i.eq.1) then
             if(l.eq. 1) then
               force(l,i,k) = -(kspring*(4*(c2_k**2)*rpos(1,1,k)**3 + 4*c1_k*c2_k*rpos(1,1,k))) &
                              - c0_k*rpos(2,2,k)**3
             else
               force(l,i,k) = -kspring*rpos(l,i,k)
             end if
           else
             if(l.eq. 2) then
               force(l,i,k) = -16.d0*kspring*rpos(2,2,k)**3 - 3.d0*c0_k*rpos(1,1,k)*rpos(2,2,k)**2  !!! quartic potential with 4*kspring
             else
               force(l,i,k) = -kspring*rpos(l,i,k)
             end if
           end if

        enddo
      enddo
    enddo
    


  elseif(trim(potential) .eq. 'symmetric_double_well_1d') then
    
    do k=1,nbeadMD
      do i=1,n 
        do l=1,ndimMD
           if (l.eq.1) then
             force(l,i,k) = -(kspring*(4*(c2_k**2)*rpos(l,i,k)**3 + 4*c0_k*c2_k*rpos(l,i,k)))
           else
             force(l,i,k) = -kspring*rpos(l,i,k)
           endif  
        enddo
      enddo
    enddo

  elseif(trim(potential) .eq. 'morse_1d') then

    do k=1,nbeadMD
      do i=1,n
        do l=1,ndimMD
           if (l.eq.1) then
             force(l,i,k) =  2*(kspring*0.5/c0_k**2)*c0_k*exp(-2.0*c0_k*rpos(l,i,k))&
                             -2*(kspring*0.5/c0_k**2)*c0_k*exp(-1.0*c0_k*rpos(l,i,k))
           else
             force(l,i,k) = -kspring*rpos(l,i,k)
           endif
        enddo
      enddo
    enddo

  elseif(trim(potential) .eq. 'quartic_1d') then

    do k=1,nbeadMD
      do i=1,n
        do l=1,ndimMD
           if (l.eq.1) then
             force(l,i,k) = -4.d0*c2_k*kspring*rpos(l,i,k)**3 
           else
             force(l,i,k) = -kspring*rpos(l,i,k)
           endif
        enddo
      enddo
    enddo
    
  elseif(trim(potential) .eq. 'quartic_plus_harmonic_1d') then

    do k=1,nbeadMD
      do i=1,n
        do l=1,ndimMD
           if (l.eq.1) then
             force(l,i,k) = -(kspring*(4*(c2_k**2)*rpos(l,i,k)**3 + 4*c0_k*c2_k*rpos(l,i,k)))
           else
             force(l,i,k) = -kspring*rpos(l,i,k)
           endif
        enddo
      enddo
    enddo

  elseif(trim(potential) .eq. 'asymmetric_double_well_1d') then
    
    do k=1,nbeadMD
      do i=1,n 
        do l=1,ndimMD
           if (l.eq.1) then
             force(l,i,k) = -kspring*(4*c2_k*rpos(l,i,k)**3 + 3*c1_k*rpos(l,i,k)**2 + 2*c0_k*rpos(l,i,k)) 
           else
             force(l,i,k) = -kspring*rpos(l,i,k)
           endif  
        enddo
      enddo
    enddo

  endif
  
  return
  
end subroutine force0
