program readgr
  
  use estimator
  use md_variables, only : filen, nbeadMD, ndimMD, n, pi, lfnl, yesquantum, amas
  implicit none

  integer i, j, b, ig, skip, flag_pair_force,&
       flag_select, done, icut, flag_mol, &
       n_atom1, n_atom2, nskip, repetition, ibinit, iend, &
       lbin, iflagerr, nbin, ngen, nion, ieskinr, d

  real(8) tmp(5), MaxL, Delta, norm, rrr, den, rr, &
       atom1, atom2, r_cut, fant, molfrac,diff,diffs, &
       cellscale(3), omega
  
  logical :: read_cov, iespbc, ifex

  real(8), dimension(:), allocatable :: grr_one, forces, forces_dev, covmat, atomic_weight
  real(8), dimension(:,:), allocatable :: my_iond, rion
  real(8), dimension(:,:,:), allocatable :: iond_cart
  type(avg_scalar), dimension(:), allocatable :: pair_force, pair_force_var
  type(avg_array) :: grr,grrJN
  type(esti_array) :: grrbin
  
  
  call input_tools_tom     

  nion=n
  ieskinr=ndimMD*n

  allocate(rion(ndimMD,nion),my_iond(nion,nion))
  allocate(iond_cart(ndimMD,nion,nion))
  allocate(atomic_weight(nion))
  
  iespbc=.false. ! period world non implemented
  repetition=0
  cellscale=0.d0
  atomic_weight=amas
  
! get nion, atomic_weight
! allocate rion
! allocate my_iond(nion,nion),iond_cart(3,nion,nion)
  
  open(unit=22,file='positions.dat',form='formatted',status='old',iostat=iflagerr)
  if(iflagerr.ne.0) then
    write(*,*) "Opening position.dat error!"
    stop
  endif

  if(yesquantum) then  ! quantum case
     nskip=5
  else
     nskip=3
  endif

  open(unit=25,file='pair_gr.out',form='formatted',status='unknown')

  !count the number of total available snapshots
  ngen=0
  do while(.true.)
    read(22,*,iostat=iflagerr)
    if(iflagerr/=0) exit
    ngen=ngen+1
  enddo
  rewind(22)

  write(*,*)
  write(*,*)
  write(*,*) ' Number of generation = ',ngen

  ! all default options
  flag_select=0
  flag_mol=0
  iflagerr=0

  MaxL=10.0
  omega=1.d0

  write(6,*) ' Max dist =',MaxL
  write(*,*)

  
  write(*,*) ' Histogram points: (if < 0 set rcut molecular frac;'
  write(*,*) '                    if > 0 normal) '
  read(5,*) nbin
  write(*,*) 'bin lenght, ibinit, iend (if <= 0 ngen),  skip (if < 0 select atom): '
  read(5,*) lbin, ibinit,iend,skip

  
  if(skip .lt. 0) then
    write(*,*) 'Compute g(r) only between selected atoms type'
    write(*,*) ' Type 1, Type 2 (if Type<0: position in input file geometry, if Type>0: Z number)'
    read(5,*) atom1,atom2
    skip = -skip

    if(atom1.lt.0 .or. atom2.lt.0) then

       flag_select = 2
       atom1=abs(atom1)
       atom2=abs(atom2)
       n_atom1=1
       n_atom2=1
       write(*,*) 'Compute g(r) only between selected atoms'
       write(*,*) 'atom 1 position',atom1
       write(*,*) 'atom 2 position',atom2
       if(atom1.eq.atom2) then
          write(*,*) 'warning: I cannot compute distance between the same atom'
          stop
       endif

    else

       write(*,*) 'Compute g(r) only between selected species'
       write(*,*) 'species 1',atom1
       write(*,*) 'species 2',atom2

       flag_select = 1
       n_atom1=0
       n_atom2=0
       do i=1,nion
          if(atomic_weight(i) .eq. atom1) n_atom1=n_atom1+1
          if(atomic_weight(i) .eq. atom2) n_atom2=n_atom2+1
       enddo
       write(6,*) 'The number of ion type 1,2 : ', n_atom1, n_atom2
    endif

  else
    n_atom1=nion
    n_atom2=nion
  endif

  r_cut=0.d0

  if(nbin .lt. 0) then
    write(*,*) 'Compute molecular fraction'
    write(*,*) 'Set rcut. (r < rcut -> r is a molecule)'
    read(5,*) r_cut
    nbin = -nbin
    flag_mol = 1
  endif

  write(*,*) ' Compute the pair force projection? ==1 yes, \=1 no'
  read(5,*) flag_pair_force
  if(flag_pair_force==1) then
    if(ieskinr.ne.nion*3) then
      write(6,*) ' ERROR! All the forces components should be computed!'
      stop
    else
      open(unit=23,file='forces.dat',form='formatted',status='old',iostat=iflagerr)
      if(iflagerr.ne.0) then
        write(*,*) ' ERROR! No forces.dat!'
        stop
      else ! the only right case
         open(unit=18,file='covmat.dat',status='old',form='formatted',iostat=iflagerr)
        if(iflagerr.eq.0) then
          read_cov=.true.
          write(*,*) " Read covmat.dat!"
        else
          read_cov=.false.
          write(*,*) " No read covmat.dat!"
        endif
        iflagerr=0
        open(unit=26,file='pair_force.dat',form='formatted',status='unknown')
      endif
    endif
  else
    flag_pair_force=0
  endif

  if(iend>0.and.iend<ngen) ngen=iend

  if(mod(ngen-ibinit+1,lbin).ne.0) then
   ibinit=ibinit+mod(ngen-ibinit+1,lbin)
   write(6,*) ' Warning changing ibinit to match bin length &
   &  of last iterations, new ibinit',ibinit
  endif


  Delta=MaxL/dble(nbin)

  icut = int(r_cut/Delta)

  !all the gr allocation here
  allocate(grr_one(nbin))
  call alloc(grr,nbin)
  if(flag_pair_force==1) then
    allocate(forces(ieskinr),forces_dev(ieskinr),pair_force(nbin),pair_force_var(nbin))
    do i=1,nbin
      call reset(pair_force(i))
      call reset(pair_force_var(i))
    enddo
    if(read_cov) allocate(covmat(ieskinr*ieskinr))
  endif
  if(.not.allocated(forces)) allocate(forces(1))
  if(.not.allocated(forces_dev)) allocate(forces_dev(1))
  if(.not.allocated(pair_force)) allocate(pair_force(1))
  if(.not.allocated(pair_force_var)) allocate(pair_force_var(1))
  if(.not.allocated(covmat)) allocate(covmat(1))

  !skip the initial ibinit-1 records
  do i=1,ibinit-1
    read(22,*,iostat=iflagerr)
    if(iflagerr/=0) exit
    if(flag_pair_force==1) then
      read(23,*,iostat=iflagerr) 
      if(iflagerr/=0) exit
      if(read_cov) then
        read(18,*,iostat=iflagerr) 
        if(iflagerr/=0) exit
      endif
    endif
  enddo

  done=0

  !main loop for all records
  do ig=ibinit,ngen,skip

    call read_one_rec
    if(iflagerr/=0) exit
    if(mod(ig,100).eq.0) write(*,*) ' Generation Number = ',ig
    done=done+1

    call eval_iond(my_iond,rion,nion,iond_cart)

    grr_one(:)=0.d0
    !distances between ions
    do i=1,nion
      do j=i+1,nion ! avoid self-pair counting for pbc
        if (flag_select .eq. 0) then
          call upgrr(i,j,nbin,grr_one,nion,delta,iond_cart,forces,forces_dev,read_cov,covmat,cellscale,&
                     repetition,flag_pair_force,pair_force,pair_force_var)
        elseif(flag_select.eq.1) then
          if( ((atomic_weight(i) .eq. atom1) .and. (atomic_weight(j) .eq. atom2)) &
            & .or. (atomic_weight(i) .eq. atom2) .and. (atomic_weight(j) .eq. atom1)) then
            call upgrr(i,j,nbin,grr_one,nion,delta,iond_cart,forces,forces_dev,read_cov,covmat,cellscale,&
                       repetition,flag_pair_force,pair_force,pair_force_var)
          endif
        else
          if(i.eq.atom1.and.j.eq.atom2) then
             call upgrr(i,j,nbin,grr_one,nion,delta,iond_cart,forces,forces_dev,read_cov,covmat,cellscale,&
                        repetition,flag_pair_force,pair_force,pair_force_var)
          endif
        endif
      enddo
    enddo
    call push(grr, grr_one)
  enddo !main loop ends

  call calc(grr)

  !apply the norm and output the grr
  if(flag_select==1 .and. atom1/=atom2 ) then
    norm=3.d0*omega/(4.d0*PI*dble(n_atom1)*dble(n_atom2))/2.d0
  elseif(flag_select==2 .and. atom1/=atom2 ) then
    norm=3.d0*omega/(4.d0*PI*dble(n_atom1)*dble(n_atom2))/2.d0
  else
    norm=3.d0*omega/(4.d0*PI*dble(n_atom1)*dble(n_atom2-1))
  endif

  write(*,*) done, 'generations done, norm = ', norm

  write(*,*) ' Begin Jacknife!'
  rewind(22)
  if(flag_pair_force==1) then
    rewind(23)
    if(read_cov) rewind(18)
  endif

  !skip the initial ibinit-1 records
  do i=1,ibinit-1
    read(22,*,iostat=iflagerr)
    if(iflagerr/=0) exit
    if(flag_pair_force==1) then
      read(23,*,iostat=iflagerr) 
      if(iflagerr/=0) exit
      if(read_cov) then
        read(18,*,iostat=iflagerr) 
        if(iflagerr/=0) exit
      endif
    endif
  enddo


  call alloc(grrJN,nbin)
  call alloc(grrbin,nbin)


  done=0

  do ig=ibinit,ngen,skip

    call read_one_rec
    if(iflagerr/=0) exit
    write(*,*) ' Generation Number = ',ig
    done=done+1

    call eval_iond(my_iond,rion,nion,iond_cart)

    grr_one(:)=0.d0
    !distances between ions
    do i=1,nion
      do j=i+1,nion ! include self-pair counting for pbc
        if (flag_select .eq. 0) then
          call upgrr(i,j,nbin,grr_one,nion,delta,iond_cart,forces,forces_dev,read_cov,covmat,cellscale,&
                     repetition,0,pair_force,pair_force_var)
        elseif(flag_select.eq.1) then
          if( ((atomic_weight(i) .eq. atom1) .and. (atomic_weight(j) .eq. atom2)) &
            & .or. (atomic_weight(i) .eq. atom2) .and. (atomic_weight(j) .eq. atom1)) then
            call upgrr(i,j,nbin,grr_one,nion,delta,iond_cart,forces,forces_dev,read_cov,covmat,cellscale,&
                       repetition,0,pair_force,pair_force_var)
          endif
       else
          if(i.eq.atom1.and.j.eq.atom2) then
             call upgrr(i,j,nbin,grr_one,nion,delta,iond_cart,forces,forces_dev,read_cov,covmat,cellscale,&
                        repetition,flag_pair_force,pair_force,pair_force_var)
          endif
        endif
      enddo
    enddo
    call push(grrJN, grr_one)
    if(grrJN%num==lbin) then
      call push(grrbin,(grr%summation-grrJN%summation)/(grr%num-grrJN%num))
      call reset(grrJN)
    endif
  enddo !main loop ends

  call calc(grrbin)
  write(*,*) grrbin%num,' Jacknife bins!'

  write(25,*) '#  Radius, g(r), error, g(r) counts'
!  write(25,*) ' Normalization =',sum(grr%average(1:nbin))
  do b=1,nbin
    rrr=(dble(b)-0.5d0)*Delta
    rr=dble(b-1)*Delta
    write(25,'(F18.10,3E18.8)') rrr, grr%average(b)*norm/((rr+Delta)**3-rr**3), &
      & sqrt(grrbin%vari(b)*max(grrbin%num-1,1))*norm/((rr+Delta)**3-rr**3), grrbin%average(b)
   
  enddo

  ! count molecular fraction
  if(flag_mol .eq. 1) then
    fant =0d0
    do b=1,icut
      fant = fant + grr%average(b)
    enddo
    molfrac = fant/nion/dble(done)
    write(6,*) 'molec frac, rcutreal, icut' , molfrac , icut*Delta , icut, Delta, abs(icut*Delta-r_cut)
  endif

  ! output the pair_force
  if(flag_pair_force==1) then
    write(26,*) '#  Radius, pair force average, error bar, counts'
    do b=1,nbin
      rrr=(dble(b)-0.5d0)*Delta
      call calc(pair_force(b))
      if(pair_force_var(b)%num>1) pair_force_var(b)%summation=pair_force_var(b)%summation/pair_force_var(b)%num
      call calc(pair_force_var(b))
      write(26,'(E11.3,2E11.3,I15)') rrr, pair_force(b)%average , sqrt(pair_force_var(b)%average), pair_force(b)%num
    enddo
  endif

  ! close files and deallocate vectors
  if(flag_pair_force==1) then
    close(23)
    close(26)
  endif
  close(25)
  close(22)

  deallocate(forces,forces_dev,pair_force,pair_force_var,covmat)
  deallocate(grr_one)
  call free(grr)
  call free(grrJN)
  call free(grrbin)

  stop

contains
  subroutine read_one_rec
    implicit none

    iflagerr=0

    !skip skip-1 records
    do i=1,skip-1
      read(22,*,iostat=iflagerr)
      if(iflagerr/=0) return
      if(flag_pair_force==1) then
        read(23,*,iostat=iflagerr) 
        if(iflagerr/=0) return
        if(read_cov) then
          read(18,*,iostat=iflagerr) 
          if(iflagerr/=0) return
        endif
      endif
    enddo

    !read one record
    if(.not.iespbc) then
      read(22,*,iostat=iflagerr) ((rion(d,i),d=1,3),i=1,nion)
      if(iflagerr/=0) return
    else
      read(22,*,iostat=iflagerr) ((rion(d,i),d=1,3),i=1,nion),cellscale(1:3)
      if(iflagerr/=0) return
      omega=cellscale(1)*cellscale(2)*cellscale(3)
    endif
    if(flag_pair_force==1) then
      read(23,*,iostat=iflagerr) tmp(1:nskip),(forces(i),forces_dev(i),i=1,ieskinr)
      if(iflagerr/=0) return
      if(read_cov) then
        read(18,*,IOSTAT=iflagerr) ((covmat(i+ieskinr*(j-1)),i=j,ieskinr),j=1,ieskinr)
        if(iflagerr/=0) return
        do i=2,ieskinr
          do j=1,i-1
            covmat(j+(i-1)*ieskinr)=covmat(i+(j-1)*ieskinr)
          enddo
        enddo
      endif
    endif
    return
  end subroutine
END

subroutine upgrr(i,j,nbin,grr_one,nion,delta,iond_cart,forces,forces_dev,read_cov,covmat,cellscale,&
                 rep,flag_pair_force,pair_force,pair_force_var)
  use estimator
  ! this subroutine makes image repetition in 3d if the system is pbc.
  ! the i==j case is not counted in both open and pbc system.
  ! if flag_pair_force==1, compute the pair force projection
  ! The two level loops to call this subroutine should be i=1,nion j=i+1,nion

  implicit none
  integer nbin,nion,ix,i,j
  integer rep
  ! rep=0,1,2,3,... (rep*2+1) times of images in each direction, 0 is open system, >0 is pbc
  integer flag_pair_force
  logical :: read_cov
  type(avg_scalar) :: pair_force(*), pair_force_var(*)
  integer rep_x,rep_y,rep_z
  real*8 delta,r_imag,grr_one(*),iond_cart(3,nion,nion),cellscale(3)
  real*8 forces(3,*),forces_dev(3,*),one_pair_force(3)
  real*8 aux(3,nion),aux_temp(3,nion),aux2(3,3)
  real*8 r_ij(3),r_ij_unit(3),covmat(3,nion,3,*)
  real*8 proj_force, proj_force_var
  real*8 dis_x2,dis_x2y2
  real(8), external :: ddot

  if(rep<0) then
    write(6,*) 'repetition should be positive!'
  else
    ! both open and pbc case

    ! compute the pair force
    if(flag_pair_force==1) then
      one_pair_force=forces(:,i)-forces(:,j)
    endif

    do rep_x=-rep,rep
      r_ij(1)=iond_cart(1,i,j)+cellscale(1)*rep_x
      dis_x2=r_ij(1)**2
      do rep_y=-rep,rep
        r_ij(2)=iond_cart(2,i,j)+cellscale(2)*rep_y
        dis_x2y2=dis_x2+r_ij(2)**2
        do rep_z=-rep,rep
          r_ij(3)=iond_cart(3,i,j)+cellscale(3)*rep_z
          r_imag=dsqrt(dis_x2y2+r_ij(3)**2)
          ix=int(r_imag/Delta)+1
          if(ix.le.nbin) then
            grr_one(ix)=grr_one(ix)+2.d0 !double counting (i,j) (j,i) pair
            if(flag_pair_force==1) then
              !calculate the projected force
              r_ij_unit(:)=r_ij(:)/r_imag ! r_ij_unit is unit vector now
              proj_force=sum(one_pair_force(1:3)*r_ij_unit(1:3))
              if(read_cov) then

                aux2(:,:)=covmat(:,i,:,i)
                aux2(:,:)=aux2(:,:)+covmat(:,j,:,j)
                aux2(:,:)=aux2(:,:)-2*covmat(:,i,:,j)
                call dgemv('N',3,3,1.d0,aux2,3,r_ij_unit,1,0.d0,aux_temp,1)
                proj_force_var=ddot(3,r_ij_unit,1,aux_temp,1)
              else
                proj_force_var=sum((r_ij_unit(:)*forces_dev(:,i))**2)+sum((r_ij_unit(:)*forces_dev(:,j))**2)
              endif
              call push(pair_force(ix),proj_force)
              call push(pair_force_var(ix),proj_force_var)
            endif
          endif
        enddo !rep_z
      enddo !rep_y
    enddo !rep_x

  endif

  return
end


subroutine eval_iond(iond,rion,nion,iond_cart)
!      use Cell
  implicit none
  integer nion,i,j,d,counter,inds,kkk
  real(8) iond(nion,nion),rion(3,*),iond_cart(3,nion,nion)

  !open system case
  do i=1,nion
    do j=i+1,nion
    !distances between ions
      iond(i,j)=dsqrt(sum((rion(:,i)-rion(:,j))**2))
      iond_cart(:,i,j)=rion(:,i)-rion(:,j)
    enddo
  enddo
! endif
       ! put the correspondent (j,i) from (i,j)
  do i=1,nion
    iond(i,i)=0.d0
    iond_cart(:,i,i)=0.d0
    !distances between ions
    do j=i+1,nion
      iond(j,i)=iond(i,j)
      iond_cart(:,j,i)=-iond_cart(:,i,j)
    enddo
  enddo
  return
END subroutine
