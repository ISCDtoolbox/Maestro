  double precision function func(energy,eig,width,typ,weight)
    implicit none     
    integer typ
    double precision energy,eig
    double precision width,weight
    double precision,parameter :: pig = 3.141592654d0
    if(typ.eq.1) then
      func=weight*2.d0/Sqrt(pig)/width*exp(-((energy-eig)/width)**2.0)   
    else if(typ.eq.2) then
      func=weight*2.d0/pig*width/((energy-eig)**2.0+width**2.0)
    endif    
  end function

  subroutine force_fire(force,r_fire)
  
    use md_variables
    implicit none
    integer :: i,l
    double precision :: pot_en_minus,pot_en_plus
    double precision :: delta
    double precision :: force(ndimMD,n),r_fire(n,ndimMD)
    double precision, allocatable :: r_fire_pp(:,:)
  
! *******************************************************************
!    calcpot(V,xx)
!    xx(7,3) is cartesian array containing O O H H H H H coor, in bohr
!    returned V is potential in hartree, taking C2-sym minimum as zero potential reference
    call prepot() 
    delta=delta_force
    force=0.d0
    allocate(r_fire_pp(n,ndimMD))
    r_fire_pp(:,:)=r_fire(:,:)
    do i=1,n 
      do l=1,ndimMD
        r_fire_pp(i,l) = r_fire(i,l) + delta
        call calcpot(pot_en_plus,r_fire_pp)
        r_fire_pp(i,l) = r_fire(i,l) - delta
        call calcpot(pot_en_minus,r_fire_pp)
        force(l,i) = (pot_en_minus-pot_en_plus)/2.d0/delta
        r_fire_pp(i,l) = r_fire(i,l)
      enddo
    enddo
    deallocate(r_fire_pp)
    return
  end subroutine force_fire
    
  program fire_k
    use md_variables
    implicit none
    integer i,iter,nitmax,n_min,nstep_fire
    double precision v1_fire,v2_fire,dt_fire,dt_md_fire,dt_max_fire
    double precision P_fire,ene_fire,ene_t_fire,af
    double precision, parameter :: af_start=0.1d0,f_inc=1.1d0,f_dec=0.5d0,f_a=0.99d0
    double precision, allocatable :: v_fire(:,:),enex_fire(:),r_fire_old(:,:),r_fire(:,:),xpr(:)
    double precision, allocatable :: dv_fire(:,:),F_fire(:,:),foverm(:,:) 

    call input_tools_tom
    
    open(2,file='output_dir_phfd/'//filen(1:lfnl)//'_phonon_eigv.out',form='formatted')
    open(3,file='output_dir_phfd/'//filen(1:lfnl)//'_phonon_dos.out',form='formatted')

    write(*,'('' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    '')')
    write(*,'(''            *** PHONONS with FINITE DIFFERENCES ***                 '' )')
    write(*,'('' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    '')')
    
    write(2,'('' PHONONS OF '//filen(1:lfnl)//' ION with FINITE DIFFERENCES AT 0 K       '')') 
    write(2,*)

    nitmax=100000
    dt_md_fire=0.02d0
    n_min=5
   !  Inizia i cicli di ottimizzazione
    dt_fire=dt_md_fire
    dt_max_fire=dt_md_fire*10.d0
    af=af_start
    nstep_fire=0
    allocate(dv_fire(3,n),F_fire(3,n),r_fire_old(n,3),r_fire(n,3),xpr(6),foverm(3,n))
    
    open(unit=39,file='configinit',status='old')
    write(*,*) 'reading the initial configuration file for minimization'
    do i=1,n  
      !call random_number(xpr)
      read(39,*) r_fire(i,1),r_fire(i,2),r_fire(i,3)
      !r_fire(i,1)=r_fire(i,1)!+xpr(1)*2
      !r_fire(i,2)=r_fire(i,2)!+xpr(2)*2
      !r_fire(i,3)=r_fire(i,3)!-xpr(3)*2
    enddo
    close(39)
    
    call force_fire(dv_fire,r_fire)
    call prepot()
    call calcpot(v1_fire,r_fire)

    allocate(v_fire(3,n),enex_fire(nitmax))
    v_fire=0.d0
    do i=1,n
      call random_number(xpr)
      v_fire(1,i)=sqrt(-2.d0*tempMD*2.d0*log(xpr(1)))*sin(2*pi*xpr(2))
      v_fire(2,i)=sqrt(-2.d0*tempMD*2.d0*log(xpr(3)))*sin(2*pi*xpr(4))
      v_fire(3,i)=sqrt(-2.d0*tempMD*2.d0*log(xpr(5)))*sin(2*pi*xpr(6))
    end do

    enex_fire=0.d0
    iter=1
    do while(iter .le. nitmax)
    !   Step MD  ----------
    !   Prende il gradiente precedentemente calcolato il gradiente di V
    !  F_fire=dv_fire
      r_fire_old=r_fire
      
    !   Aggiorna la v a dt/2
      do i=1,n
        v_fire(1,i)=v_fire(1,i)+(1/amas(i))*F_fire(1,i)*dt_fire*0.5d0
        v_fire(2,i)=v_fire(2,i)+(1/amas(i))*F_fire(2,i)*dt_fire*0.5d0
        v_fire(3,i)=v_fire(3,i)+(1/amas(i))*F_fire(3,i)*dt_fire*0.5d0
      end do
 
      !   Aggiorna le posizioni
      r_fire(:,1)=r_fire(:,1)+v_fire(1,:)*dt_fire
      r_fire(:,2)=r_fire(:,2)+v_fire(2,:)*dt_fire
      r_fire(:,3)=r_fire(:,3)+v_fire(3,:)*dt_fire
      
        !   Ricalcola accel.
      
      call force_fire(F_fire,r_fire)
      call prepot()
      v2_fire=0.d0
      call calcpot(v2_fire,r_fire)
      
      enex_fire(iter)=v2_fire

      do i=1,n
        v_fire(1,i)=v_fire(1,i)+(1/amas(i))*F_fire(1,i)*dt_fire*0.5d0
        v_fire(2,i)=v_fire(2,i)+(1/amas(i))*F_fire(2,i)*dt_fire*0.5d0
        v_fire(3,i)=v_fire(3,i)+(1/amas(i))*F_fire(3,i)*dt_fire*0.5d0
      end do

!      write(*,*) v2_fire+0.5*sum(amas(:)*v_fire(1,:)**2)+&
!                         0.5*sum(amas(:)*v_fire(2,:)**2)+&
!                         0.5*sum(amas(:)*v_fire(3,:)**2)

        !   Fine step MD  --------------
        !   Calcola P=F.v
      foverm(1,:)=F_fire(1,:)/amas(:)
      foverm(2,:)=F_fire(2,:)/amas(:)
      foverm(3,:)=F_fire(3,:)/amas(:)
      P_fire=sum(F_fire(1,:)*v_fire(1,:))+sum(F_fire(2,:)*v_fire(2,:))+sum(F_fire(3,:)*v_fire(3,:))
      v_fire(1,:)=(1.d0-af)*v_fire(1,:)+af*sqrt(sum(v_fire(:,:)**2.d0))*F_fire(1,:)/sqrt(sum(F_fire(:,:)**2.d0))
      v_fire(2,:)=(1.d0-af)*v_fire(2,:)+af*sqrt(sum(v_fire(:,:)**2.d0))*F_fire(2,:)/sqrt(sum(F_fire(:,:)**2.d0))
      v_fire(3,:)=(1.d0-af)*v_fire(3,:)+af*sqrt(sum(v_fire(:,:)**2.d0))*F_fire(3,:)/sqrt(sum(F_fire(:,:)**2.d0))
!!!        if ((P .gt. 0.d0) .and. (tvv1-tvv2 .ge. 0.d0)) then
      if ((P_fire .gt. 0.d0) .and. (v1_fire-v2_fire .ge. 0.d0)) then
        if (nstep_fire .gt. n_min) then
          dt_fire=min(dt_fire*f_inc,dt_max_fire)
          af=af*f_a
        end if
        nstep_fire=nstep_fire+1
        dv_fire=F_fire
        v1_fire=v2_fire 
      else   
        dt_fire=dt_fire*f_dec
        r_fire=r_fire_old
        v_fire=0.d0
        af=af_start
      end if
      if (iter .gt. n_min) then
        if ((dt_fire .lt. 1.d-8) .or. (abs(enex_fire(iter)-enex_fire(iter-1)) .lt. 1.d-12)) exit
      end if
      write(*,*) iter,'energy minimum of zundel ion with CC potential',enex_fire(iter)
      iter=iter+1
    end do
    
    open(81,file='output_dir_phfd/'//filen(1:lfnl)//'_opt.xyz',form='formatted')
    write(81,*) n
    write(81,*)
    do i=1,n
      write(81,*) ion_name(i),r_fire(i,:)
    end do
    close(81)
    
    call phonon_fd(r_fire)
            
    deallocate(v_fire,enex_fire,dv_fire,F_fire,r_fire_old,r_fire,xpr,foverm)
    deallocate(amas,ion_name,indx)
  end program
  

  subroutine phonon_fd(rpos_loc)
  
    use md_variables 

    implicit none 
    integer :: i,j,k,l,ll,indi,indj,lworkd,elim_tr_rot
    double precision :: ep,ep_plus,ep_minus,deltasq
    double precision :: ep_plxply,ep_plxmy,ep_mxply,ep_mxmy
    double precision :: delta
    double precision :: rpos_loc(n,ndimMD)
    double precision, allocatable :: dynmat_force(:,:),work(:)
  
    integer, parameter :: num_m = 1000
  
    double precision  t_dos,temp_eig
    double precision, allocatable :: dos(:,:), int_dos(:,:)
    double precision, allocatable :: weight(:)
    double precision  energy,e_l,min_eig,max_eig
    double precision  width,e_sum,sum_weight,func,co_tom
    integer typ

    delta=delta_harm  !!! attention: very important parameter
    deltasq=delta**2
    elim_tr_rot=6

    allocate(dynmat_force(n*ndimMD,n*ndimMD))
    allocate(dynmat(n*ndimMD,n*ndimMD))
    allocate(dynmat_eig(n*ndimMD),dynmatforce_eig(n*ndimMD))
    allocate(rtilde(n,ndimMD))
!    allocate(indx(n))
  
    do i=1,n 
       do l=1,ndimMD
          rtilde(i,l) = rpos_loc(i,l)
       enddo
    enddo
  
    call prepot()
    call calcpot(ep,rtilde) 
        
! Computing dynamical matrix
! lower triangle

    dynmat=0.d0
    dynmat_force=0.d0
    
    do j=1,n
       do ll=1,ndimMD

         indj=(j-1)*ndimMD+ll
 
! diagonal part (equal particles, equal coordinates)

         rtilde(j,ll) = rpos_loc(j,ll) + delta
         !call prepot()
         call calcpot(ep_plus,rtilde)
 
         rtilde(j,ll) = rpos_loc(j,ll) - delta
         !call prepot()
         call calcpot(ep_minus,rtilde)
        
         rtilde(j,ll) = rpos_loc(j,ll)
         dynmat_force(indj,indj) = (ep_plus+ep_minus-2.d0*ep)/deltasq
         dynmat(indj,indj) =  dynmat_force(indj,indj)/amas(indx(j))

! off-diag part (equal particles)

         do l=ll+1,ndimMD
           indi=(j-1)*ndimMD+l 
           
           rtilde(j,l) = rpos_loc(j,l) + delta
           rtilde(j,ll) = rpos_loc(j,ll) + delta
           !call prepot()
           call calcpot(ep_plxply,rtilde)
                        
           rtilde(j,l) = rpos_loc(j,l) - delta
           !call prepot()
           call calcpot(ep_mxply,rtilde)
                 
           rtilde(j,ll) = rpos_loc(j,ll) - delta
           !call prepot()
           call calcpot(ep_mxmy,rtilde)
                 
           rtilde(j,l) = rpos_loc(j,l) + delta                 
           !call prepot()
           call calcpot(ep_plxmy,rtilde)
                 
                 
           dynmat_force(indi,indj) = (ep_plxply+ep_mxmy-ep_plxmy-ep_mxply)/4.d0/deltasq
           dynmat(indi,indj) = dynmat_force(indi,indj)/sqrt(amas(indx(j))*amas(indx(j))) 

           rtilde(j,l) = rpos_loc(j,l)
           rtilde(j,ll) = rpos_loc(j,ll)
         enddo

! off-diag part (unequal particles, unequal coordinates)
         do i=j+1,n
           do l=1,ndimMD
              indi=(i-1)*ndimMD+l 

              rtilde(i,l) = rpos_loc(i,l) + delta
              rtilde(j,ll) = rpos_loc(j,ll) + delta
              !call prepot()
              call calcpot(ep_plxply,rtilde)
                        
              rtilde(i,l) = rpos_loc(i,l) - delta
              !call prepot()
              call calcpot(ep_mxply,rtilde)
                 
              rtilde(j,ll) = rpos_loc(j,ll) - delta
              !call prepot()
              call calcpot(ep_mxmy,rtilde)
                 
              rtilde(i,l) = rpos_loc(i,l) + delta                 
              !call prepot()
              call calcpot(ep_plxmy,rtilde)
          
               
              dynmat_force(indi,indj) = (ep_plxply+ep_mxmy-ep_plxmy-ep_mxply)/4.d0/deltasq
              dynmat(indi,indj) = dynmat_force(indi,indj)/sqrt(amas(indx(i))*amas(indx(j)))    

              rtilde(i,l) = rpos_loc(i,l)
              rtilde(j,ll) = rpos_loc(j,ll)

           enddo
         enddo

       enddo
    enddo

! Write out the dynamical matrix
!  open(12,file='zz_dinmat')
  do indi=1,n*ndimMD
     do indj=1,n*ndimMD        
        write(*,*) indi,indj,dynmat(indi,indj)   
     enddo
   enddo

! Diagonalization of the dynamical matrix
    lworkd = 3*n*ndimMD
  
    if(allocated(work)) deallocate(work)
    allocate(work(lworkd))
  
    call dsyev('V','L', n*ndimMD, dynmat, n*ndimMD, dynmat_eig, work, lworkd, info)  

    if(info.ne.0) then
      write(6,*) 'some problem in dsyev from dynmatrix routine: info =',info
      stop
    endif

    dynmat_eig(:)=abs(dynmat_eig(:))   ! The matrix is positive definite
!  write(6,*) 'Dynamical matrix eigenvalues are'
!  do i=1,n*ndimMD
!     write(6,*) i,dynmat_eig(i)
!  enddo

    write(6,*) 'Harmonic frequencies in cm^-1'
    do i=1,n*ndimMD
       write(6,*) i,sqrt(dynmat_eig(i))*ha2cm1
       write(2,*) i,0.d0,sqrt(dynmat_eig(i))*ha2cm1
    enddo
    close(2)


    allocate(weight(n*ndimMD))
    allocate(dos(2,num_m+1),int_dos(2,num_m+1))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!  read instruction for broadening algorithm  !    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    typ=1
    width=1.d0
    weight=1.d0
  !  write(*,*),'broadening function =',1
  !  write(*,'(''energy broadening'',f15.8)')width     
    min_eig=sqrt(dynmat_eig(1))*ha2cm1
    max_eig=sqrt(dynmat_eig(n*ndimMD))*ha2cm1
    co_tom=10.d0 !! cm^-1
    e_l = (max_eig+co_tom) - (min_eig-co_tom)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!           construct DOS                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    e_sum=0.0
    do i=0,num_m  
      t_dos=0.0
      energy=(min_eig-co_tom)+(i*1.d0/num_m)*e_l
      do j=1+elim_tr_rot,n*ndimMD
        t_dos=t_dos+func(energy,sqrt(dynmat_eig(j))*ha2cm1,width,typ,weight(j))
      end do
      dos(1,i+1)=energy
      dos(2,i+1)=t_dos
      e_sum=e_sum+t_dos
      int_dos(1,i+1)=energy
      int_dos(2,i+1)=e_sum*1.2*e_l/num_m
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
!     check integral of DOS                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     output of DOS                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
    do i=1,num_m+1
      write(3,'(2f16.9)')dos(1,i),dos(2,i)*100.d0
    end do
    close(3)

    deallocate(work)
    deallocate(dos,int_dos,weight)
    deallocate(dynmat_force,dynmat)
    deallocate(dynmat_eig,dynmatforce_eig) 
    deallocate(rtilde)

    return         
  end subroutine phonon_fd
