subroutine pot(ep,ep_centroid)
  
  use md_variables
  implicit none
  integer :: i,l,k
  real(8) :: ep,pot_bead,ep_centroid,pot_en
  
! ****************** ep_centroid ***********************

    ep_centroid=0.d0
  
    if(nbeadMD.gt.1) then
      rcentroidtilde=0.d0 
      do k=1,nbeadMD
        do i=1,n 
          do l=1,ndimMD
             rcentroidtilde(i,l) = rcentroidtilde(i,l) + rpos(l,i,k)
          enddo
        enddo
      enddo
      rcentroidtilde(:,:)=rcentroidtilde(:,:)/nbeadMD
 
      do i=1,n
        do l=1,ndimMD
          rcentroid(l,i)=rcentroidtilde(i,l)
        enddo 
      enddo 
      
      if(trim(potential) .eq. 'zundel') then
        call prepot()
        call calcpot(pot_en,rcentroidtilde)
        ep_centroid = pot_en
     
      elseif(trim(potential) .eq. 'harmonic_chain1p_1d' .or. trim(potential) .eq. 'harmonic_chain2p_1d') then
        do i=1,n
          do l=1,ndimMD
            if(l.eq.1) then
              if(i.eq.n) then
                 ep_centroid=ep_centroid+0.5*kspring*(rcentroid(l,i)-(rcentroid(l,1)+avecsp(1,1)))**2
              else
                 ep_centroid=ep_centroid+0.5*kspring*(rcentroid(l,i)-rcentroid(l,i+1))**2
              end if
            end if
        !!    if(l.ne.1) ep_centroid=ep_centroid+0.5*0.7d0*kspring*rcentroid(l,i)**2
          enddo
        enddo

      elseif(trim(potential) .eq. 'harmonic_1d') then
        do i=1,n
          do l=1,ndimMD
            if(l.eq.1) ep_centroid=ep_centroid+0.5*kspring*rcentroid(l,i)**2
            if(l.ne.1) ep_centroid=ep_centroid+0.5*0.7d0*kspring*rcentroid(l,i)**2 
          enddo
        enddo
     
      elseif(trim(potential) .eq. 'coupled_harmonic_1d') then
        do i=1,n
          do l=1,ndimMD
            if(i.eq.1) then          
              if(l.eq.1) then
                !!ep_centroid=ep_centroid+0.5*kspring*rcentroid(l,i)**2+c1_k*exp(-c0_k*(rcentroid(l,i)-rcentroid(l,n))**2)/n
                ep_centroid=ep_centroid + 0.5 * kspring * rcentroid(l,i)**2 + 0.5 * c1_k * rcentroid(l,i) * rcentroid(l,i + 1)
              else
                ep_centroid=ep_centroid+0.5*kspring*rcentroid(l,i)**2
              endif
            else
              if(l.eq.1) then
                !!ep_centroid=ep_centroid+0.5*kspring*rcentroid(l,i)**2+c1_k*exp(-c0_k*(rcentroid(l,n)-rcentroid(l,1))**2)/n
                ep_centroid=ep_centroid + 0.5 * c0_k * rcentroid(l,i)**2
              else
                ep_centroid=ep_centroid+0.5*kspring*rcentroid(l,i)**2
              endif
            end if
          enddo
        enddo        

      elseif(trim(potential) .eq. 'coupled_harmonic_2d') then
        do i=1,n
          do l=1,ndimMD
            if(i.eq.1) then
              if(l.eq.1) then
                ep_centroid=ep_centroid+0.5*c1_k*0.8*kspring*rcentroid(l,i)**2 &
                             + c0_k*0.8*kspring*rcentroid(1,1)*rcentroid(2,2)/n
              else
                ep_centroid=ep_centroid+0.5*c1_k*kspring*rcentroid(l,i)**2
              endif
            else
              if(l.eq.2) then
                ep_centroid=ep_centroid+0.5*c2_k*0.8*kspring*rcentroid(l,i)**2 &
                             + c0_k*0.8*kspring*rcentroid(1,1)*rcentroid(2,2)/n
              else
                ep_centroid=ep_centroid+0.5*c2_k*kspring*rcentroid(l,i)**2
              endif
            
            
            end if
          enddo
        enddo

      elseif(trim(potential) .eq. 'coupled_sdw_ha_x2y2_2d') then
        do i=1,n
          do l=1,ndimMD
            if(i.eq.1) then
              if(l.eq.1) then
                ep_centroid=ep_centroid + kspring*(c2_k*rcentroid(l,i)**2+c1_k)**2 &
                             + c0_k*(rcentroid(1,1)**2)*rcentroid(2,2)**2/n
              else
                ep_centroid=ep_centroid+0.5*kspring*rcentroid(l,i)**2
              endif
            else
              if(l.eq.2) then
                ep_centroid=ep_centroid+0.5*(0.89114835d0)*kspring*rcentroid(l,i)**2 & 
                            + c0_k*(rcentroid(1,1)**2)*rcentroid(2,2)**2/n
              else
                ep_centroid=ep_centroid+0.5*kspring*rcentroid(l,i)**2
              endif


            end if
          enddo
        enddo

      elseif(trim(potential) .eq. 'coupled_sdw_qu_xy2_2d') then
        do i=1,n
          do l=1,ndimMD
            if(i.eq.1) then
              if(l.eq.1) then
                ep_centroid=ep_centroid + kspring*(c2_k*rcentroid(l,i)**2+c1_k)**2 & 
                             + c0_k*rcentroid(1,1)*rcentroid(2,2)**2/n
              else
                ep_centroid=ep_centroid+0.5*kspring*rcentroid(l,i)**2
              endif
            else
              if(l.eq.2) then
                ep_centroid=ep_centroid+4.0*kspring*rcentroid(l,i)**4 & 
                            + c0_k*rcentroid(1,1)*rcentroid(2,2)**2/n
              else
                ep_centroid=ep_centroid+0.5*kspring*rcentroid(l,i)**2
              endif


            end if
          enddo
        enddo
        if(c0_k.eq.0.1d0) ep_centroid=ep_centroid+0.00035535d0
        if(c0_k.eq.0.2d0) ep_centroid=ep_centroid+0.0016118d0
    
      elseif(trim(potential) .eq. 'coupled_sdw_qu_xy3_2d') then
        do i=1,n
          do l=1,ndimMD
            if(i.eq.1) then
              if(l.eq.1) then
                ep_centroid=ep_centroid + kspring*(c2_k*rcentroid(l,i)**2+c1_k)**2 &
                             + c0_k*rcentroid(1,1)*rcentroid(2,2)**3/n
              else
                ep_centroid=ep_centroid+0.5*kspring*rcentroid(l,i)**2
              endif
            else
              if(l.eq.2) then
                ep_centroid=ep_centroid+4.0*kspring*rcentroid(l,i)**4 &
                            + c0_k*rcentroid(1,1)*rcentroid(2,2)**3/n
              else
                ep_centroid=ep_centroid+0.5*kspring*rcentroid(l,i)**2
              endif


            end if
          enddo
        enddo        
        if(c0_k.eq.0.4d0) ep_centroid=ep_centroid+0.00007049d0

      elseif(trim(potential) .eq. 'symmetric_double_well_1d') then
        do i=1,n
          do l=1,ndimMD
            if(l.eq.1) then
              ep_centroid=ep_centroid + kspring*(c2_k*rcentroid(l,i)**2+c0_k)**2
            else
              ep_centroid=ep_centroid+0.5*kspring*rcentroid(l,i)**2
            endif
          enddo
        enddo 

      elseif(trim(potential) .eq. 'morse_1d') then
        do i=1,n
          do l=1,ndimMD
            if(l.eq.1) then
              ep_centroid=ep_centroid+(kspring*0.5/c0_k**2)*(1.0-exp(-c0_k*rcentroid(l,i)))**2
            else
              ep_centroid=ep_centroid+0.5*kspring*rcentroid(l,i)**2
            endif
          enddo
        enddo
        
      elseif(trim(potential) .eq. 'quartic_1d') then
        do i=1,n
          do l=1,ndimMD
            if(l.eq.1) then
              ep_centroid=ep_centroid + c2_k*kspring*rcentroid(l,i)**4
            else
              ep_centroid=ep_centroid+0.5*kspring*rcentroid(l,i)**2
            endif
          enddo
        enddo

      elseif(trim(potential) .eq. 'quartic_plus_harmonic_1d') then
        do i=1,n
          do l=1,ndimMD
            if(l.eq.1) then
              ep_centroid=ep_centroid + kspring*(c2_k*rcentroid(l,i)**2+c0_k)**2 -kspring*c0_k**2
            else
              ep_centroid=ep_centroid+0.5*kspring*rcentroid(l,i)**2
            endif
          enddo
        enddo        

      elseif(trim(potential) .eq. 'asymmetric_double_well_1d') then
        do i=1,n
          do l=1,ndimMD
            if(l.eq.1) then
              ep_centroid=ep_centroid + kspring*(c2_k*rcentroid(l,i)**4+c1_k*rcentroid(l,i)**3+c0_k*rcentroid(l,i)**2)
            else
              ep_centroid=ep_centroid+0.5*kspring*rcentroid(l,i)**2
            endif
          enddo
        enddo 

      endif
  
    endif

! ****************** ep ***********************  
  
    ep=0.d0
    
    if(trim(potential) .eq. 'zundel') then
      call prepot()
      do k=1,nbeadMD
        do i=1,n
           do l=1,ndimMD
              rtilde(i,l)=rpos(l,i,k)
           enddo
        enddo    
        call calcpot(pot_bead,rtilde)
        ep=ep+pot_bead
      enddo
      ep=ep/nbeadMD
    
    elseif(trim(potential) .eq. 'harmonic_chain1p_1d' .or. trim(potential) .eq. 'harmonic_chain2p_1d') then

      do k=1,nbeadMD
        do i=1,n
          do l=1,ndimMD
            if(l.eq.1) then
              if(i.eq.n) then
                 ep=ep+0.5*kspring*(rpos(l,n,k)-(rpos(l,1,k)+avecsp(1,1)))**2
              else
                 ep=ep+0.5*kspring*(rpos(l,i,k)-rpos(l,i+1,k))**2
              end if
            end if  
           !! if(l.ne.1) ep=ep+0.5*0.7d0*kspring*rpos(l,i,k)**2
          enddo
        enddo
      enddo
      ep=ep/nbeadMD

    elseif(trim(potential) .eq. 'harmonic_1d') then
      
      do k=1,nbeadMD
        do i=1,n
          do l=1,ndimMD
            if(l.eq.1) ep=ep+0.5*kspring*rpos(l,i,k)**2
            if(l.ne.1) ep=ep+0.5*0.7d0*kspring*rpos(l,i,k)**2
          enddo
        enddo
      enddo
      ep=ep/nbeadMD 

    elseif(trim(potential) .eq. 'coupled_harmonic_1d') then
       
        do k=1,nbeadmD
          do i=1,n
            do l=1,ndimMD           
             if (i.eq.1) then

               if(l.eq.1) then
                 ep = ep +0.5*kspring*rpos(l,i,k)**2+c1_k*exp(-c0_k*(rpos(l,i,k)-rpos(l,n,k))**2)/n
               else
                 ep = ep + 0.5*kspring*rpos(l,i,k)**2
               endif 
             else
               if(l.eq.1) then
                 ep = ep +0.5*kspring*rpos(l,i,k)**2+c1_k*exp(-c0_k*(rpos(l,i,k)-rpos(l,1,k))**2)/n
               else
                 ep = ep + 0.5*kspring*rpos(l,i,k)**2
               endif

             end if
            enddo
          enddo  
        enddo
        ep=ep/nbeadMD

    elseif(trim(potential) .eq. 'coupled_harmonic_2d') then

        do k=1,nbeadmD
          do i=1,n
            do l=1,ndimMD
             if (i.eq.1) then

               if(l.eq.1) then
                 ep = ep + 0.5*c1_k*0.8*kspring*rpos(l,i,k)**2 + (c0_k*0.8*kspring)*rpos(1,1,k)*rpos(2,2,k)/n
               else
                 ep = ep + 0.5*c1_k*kspring*rpos(l,i,k)**2
               endif
             else
               if(l.eq.2) then
                 ep = ep + 0.5*c2_k*0.8*kspring*rpos(l,i,k)**2 + (c0_k*0.8*kspring)*rpos(1,1,k)*rpos(2,2,k)/n
               else
                 ep = ep + 0.5*c2_k*kspring*rpos(l,i,k)**2
               endif


             end if
            enddo
          enddo
        enddo
        ep=ep/nbeadMD

    elseif(trim(potential) .eq. 'coupled_sdw_ha_x2y2_2d') then

        do k=1,nbeadmD
          do i=1,n
            do l=1,ndimMD
             if (i.eq.1) then

               if(l.eq.1) then
                 ep = ep + kspring*(c2_k*rpos(l,i,k)**2+c1_k)**2 + c0_k*(rpos(1,1,k)**2)*rpos(2,2,k)**2/n
               else
                 ep = ep + 0.5*kspring*rpos(l,i,k)**2
               endif
             else
               if(l.eq.2) then
                 ep = ep + 0.5*(0.89114835d0)*kspring*rpos(l,i,k)**2+c0_k*(rpos(1,1,k)**2)*rpos(2,2,k)**2/n
               else
                 ep = ep + 0.5*kspring*rpos(l,i,k)**2
               endif


             end if
            enddo
          enddo
        enddo
        ep=ep/nbeadMD

    elseif(trim(potential) .eq. 'coupled_sdw_qu_xy2_2d') then

        do k=1,nbeadmD
          do i=1,n
            do l=1,ndimMD
             if (i.eq.1) then

               if(l.eq.1) then
                 ep = ep + kspring*(c2_k*rpos(l,i,k)**2+c1_k)**2 + c0_k*rpos(1,1,k)*rpos(2,2,k)**2/n
               else
                 ep = ep + 0.5*kspring*rpos(l,i,k)**2
               endif
             else
               if(l.eq.2) then
                 ep = ep + 4.d0*kspring*rpos(l,i,k)**4 + c0_k*rpos(1,1,k)*rpos(2,2,k)**2/n
               else
                 ep = ep + 0.5*kspring*rpos(l,i,k)**2
               endif


             end if
            enddo
          enddo
        enddo
        ep=ep/nbeadMD
        if(c0_k.eq.0.1d0) ep=ep+0.00035535d0
        if(c0_k.eq.0.2d0) ep=ep+0.0016118d0

    elseif(trim(potential) .eq. 'coupled_sdw_qu_xy3_2d') then

        do k=1,nbeadmD
          do i=1,n
            do l=1,ndimMD
             if (i.eq.1) then

               if(l.eq.1) then
                 ep = ep + kspring*(c2_k*rpos(l,i,k)**2+c1_k)**2 + c0_k*rpos(1,1,k)*rpos(2,2,k)**3/n
               else
                 ep = ep + 0.5*kspring*rpos(l,i,k)**2
               endif
             else
               if(l.eq.2) then
                 ep = ep + 4.d0*kspring*rpos(l,i,k)**4 + c0_k*rpos(1,1,k)*rpos(2,2,k)**3/n
               else
                 ep = ep + 0.5*kspring*rpos(l,i,k)**2
               endif


             end if
            enddo
          enddo
        enddo
        ep=ep/nbeadMD
        if(c0_k.eq.0.4d0) ep=ep+0.00007049d0



    elseif(trim(potential) .eq. 'symmetric_double_well_1d') then
      
      do k=1,nbeadMD
        do i=1,n
          do l=1,ndimMD
              if(l.eq.1) then
                ep=ep + kspring*(c2_k*rpos(l,i,k)**2+c0_k)**2
              else
               ep=ep+0.5*kspring*rpos(l,i,k)**2
              endif 
         enddo
        enddo
      enddo
      ep=ep/nbeadMD 

    elseif(trim(potential) .eq. 'morse_1d') then

      do k=1,nbeadMD
        do i=1,n
          do l=1,ndimMD
              if(l.eq.1) then
                ep=ep+(kspring*0.5/c0_k**2)*(1.0-exp(-c0_k*rpos(l,i,k)))**2
              else
               ep=ep+0.5*kspring*rpos(l,i,k)**2
              endif
         enddo
        enddo
      enddo
      ep=ep/nbeadMD

    elseif(trim(potential) .eq. 'quartic_1d') then

      do k=1,nbeadMD
        do i=1,n
          do l=1,ndimMD
              if(l.eq.1) then
                ep=ep+c2_k*kspring*rpos(l,i,k)**4
              else
               ep=ep+0.5*kspring*rpos(l,i,k)**2
              endif
         enddo
        enddo
      enddo
      ep=ep/nbeadMD

    elseif(trim(potential) .eq. 'quartic_plus_harmonic_1d') then

      do k=1,nbeadMD
        do i=1,n
          do l=1,ndimMD
              if(l.eq.1) then
                ep=ep + kspring*(c2_k*rpos(l,i,k)**2+c0_k)**2 -kspring*c0_k**2
              else
               ep=ep+0.5*kspring*rpos(l,i,k)**2
              endif
         enddo
        enddo
      enddo
      ep=ep/nbeadMD

    elseif(trim(potential) .eq. 'asymmetric_double_well_1d') then
      
      do k=1,nbeadMD
        do i=1,n
          do l=1,ndimMD
              if(l.eq.1) then
                ep=ep + kspring*(c2_k*rpos(l,i,k)**4 + c1_k*rpos(l,i,k)**3 + c0_k*rpos(l,i,k)**2)
              else
               ep=ep+0.5*kspring*rpos(l,i,k)**2
              endif 
         enddo
        enddo
      enddo
      ep=ep/nbeadMD 
            
     endif
  
  return
  
end subroutine pot
